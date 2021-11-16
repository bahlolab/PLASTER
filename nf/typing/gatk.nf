
workflow gatk {
    take:
        par // ?ploidy, ?qd, ?targ (am, vcf, tbi)
        bams // am, sm, ?ps, nr, bam, bai
        amps // am, reg
        ref // fasta, fai, dict

    main:
        par = [ploidy:2, qd: '2.0', targ: false] + par
        bams = bams.map { it.size() == 6 ? it : it.take(2) + [null] + it.takeRight(3) }

        (amps |
            combine(bams, by:0))
            .with { haplotype_caller(it, ref, par.ploidy) }

        (haplotype_caller.out |
            groupTuple(by:[0,1], sort: true) |
            map { it[[0, 1, 4, 5]] })
            .with { genotype_gvcfs(it, ref, par.qd) }

        if (par.targ) {
            targeted =
                genotype_gvcfs.out |
                    join(par.targ, by:0, remainder:true) |
                    branch { yes: it[3] != null; no: true }

            get_targ_sites(targeted.yes, ref)

            ((get_targ_sites.out | // am, sites, tbi
                combine(amps, by: 0) | //am, sites, tbi, reg
                combine(bams, by: 0)) // am, sites, tbi, reg,  sm, ps, nr, bam, bai
                .with { call_targ_sites(it, ref, par.ploidy) } |
                groupTuple(by: 0, sort: true) |
                map { it[[0, 3, 4]] } |
                combine(get_targ_sites.out, by:0) |
                combine(genotype_gvcfs.out, by:0))
                .with { merge_targ_sites(it, ref) }
            out = targeted.no |
                map {it.take(3) } |
                mix(merge_targ_sites.out)
        } else {
            out = genotype_gvcfs.out
        }
    emit:
        out//am, vcf, tbi
}

process haplotype_caller {
    label 'M2_NR'
    publishDir "progress/haplotype_caller", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am${ps ? ':' + ps : ''}" }

    input:
        tuple val(am), val(reg), val(sm), val(ps), val(nr), path(bam), path(bai)
        tuple path(ref), path(fai), path(dict)
        val(ploidy)

    output:
        tuple val(am), val(reg), val(sm), val(ps), path(gvcf), path("${gvcf}.tbi")

    script:
        gvcf = "SM-${sm}.AM-${am}" + (ps ? ".PS-${ps}" : '') + ".gvcf.gz"
        """
        gatk HaplotypeCaller \\
            --java-options "-Xmx5G -Djava.io.tmpdir=." \\
            -R $ref \\
            -I $bam \\
            -O $gvcf \\
            -ERC GVCF \\
            --intervals $reg \\
            --sample-ploidy $ploidy \\
            --pcr-indel-model AGGRESSIVE
        """
}

process genotype_gvcfs {
    label 'M2'
    publishDir "progress/genotype_gvcfs", mode: "$params.intermediate_pub_mode"
    tag { am }

    input:
        tuple val(am), val(reg), path(gvcfs), path(tbis)
        tuple path(ref), path(fai), path(dict)
        val(qd)

    output:
        tuple val(am), path(vcf), path("${vcf}.tbi")

    script:
        vcf = "${am}.gatk.vcf.gz"
        """
        gatk CombineGVCFs \\
            --java-options "-Xmx4G -Djava.io.tmpdir=." \\
            -R $ref \\
            -L $reg \\
            -O combined.g.vcf.gz \\
            --disable-sequence-dictionary-validation \\
            --variant ${gvcfs.sort{ it.name }.join(' --variant ')}

        gatk GenotypeGVCFs \\
            --java-options "-Xmx4G -Djava.io.tmpdir=." \\
            -R $ref \\
            -V combined.g.vcf.gz\\
            -O tmp.vcf.gz \\
            -G StandardAnnotation \\
            -G AS_StandardAnnotation \\
            --use-new-qual-calculator

        bcftools filter tmp.vcf.gz -Ou -e 'QD<$qd' -s lowQD |
            bcftools view -f PASS -Ou |
            bcftools norm -f $ref -Ou |
            bcftools view -i 'GT="alt"' -Oz -o $vcf

        bcftools index -t $vcf
        """
}

process get_targ_sites {
    label 'S'
    publishDir "progress/get_targ_sites", mode: "$params.intermediate_pub_mode"
    tag { "$am" }

    input:
        tuple val(am), path(vcf), path(tbi), path(sites), path(sites_tbi)
        tuple path(ref), path(fai), path(dict)

    output:
        tuple val(am), path(out), path("${out}.tbi")

    script:
        out = "${am}.absent.vcf.gz"
        """
        bcftools view -G -Ou $vcf |
            bcftools norm -m-any -f $ref -Oz -o calls.norm.vcf.gz
        bcftools index -t calls.norm.vcf.gz
        bcftools isec $sites calls.norm.vcf.gz -C -w 1 -Ou |
            bcftools norm -m+any -f $ref -Oz -o $out
        bcftools index -t $out
        """
}

process call_targ_sites {
    label 'M2_NR'
    publishDir "progress/call_targ_sites", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am:$ps" }

    input:
        tuple val(am), path(sites), path(tbi), val(reg), val(sm), val(ps), val(nr), path(bam), path(bai)
        tuple path(ref), path(fai), path(dict)
        val(ploidy)

    output:
        tuple val(am), val(sm), val(ps), path(out), path("${out}.tbi")

    script:
        out = "SM-${sm}.AM-${am}" + (ps ? ".PS-${ps}" : '') + ".targ.vcf.gz"
        """
        gatk HaplotypeCaller \\
            --java-options "-Xmx4G -Djava.io.tmpdir=." \\
            -R $ref \\
            -I $bam \\
            -O tmp.vcf.gz \\
            --alleles $sites \\
            --force-call-filtered-alleles \\
            --intervals $reg \\
            --sample-ploidy $ploidy
        bcftools view tmp.vcf.gz -i 'GT="alt"' -Ou |
            bcftools norm -m-any -f $ref -Oz -o $out
        bcftools index -t $out
        """
}

process merge_targ_sites {
    label 'M2'
    publishDir "progress/merge_targ_sites", mode: "$params.intermediate_pub_mode"
    tag { am }

    input:
        tuple val(am), path(vcfs), path(tbis), path(sites), path(stbi), path(disc), path(dtbi)
        tuple path(ref), path(fai), path(dict)

    output:
        tuple val(am), path(out), path("${out}.tbi")

    script:
        out = "${am}.gatk_targeted.vcf.gz"
        """
        bcftools merge -0 -Oz ${vcfs.sort{ it.name }.join(' ')} -o merged.vcf.gz
        bcftools index -t merged.vcf.gz
        bcftools norm -m-any -f $ref $sites -Oz -o sites.norm.vcf.gz
        bcftools index -t sites.norm.vcf.gz
        bcftools isec merged.vcf.gz sites.norm.vcf.gz -p isec -Oz -w 1
        bcftools concat -a isec/0002.vcf.gz $disc -D -Ou |
            bcftools norm -m-any -f $ref -Oz -o $out
        bcftools index -t $out
        """
}
