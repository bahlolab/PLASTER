
workflow gatk {
    take:
        data // am, reg, sm, ps, nr, bam, bai
        ref // fasta, fai, dict
        ploidy
        QD

    main:
        gvcfs = haplotype_caller(data, ref, ploidy) |
            groupTuple(by:[0,1], sort: true)
        genotype_gvcfs(gvcfs, ref, QD)
    emit:
        // am, vcf, tbi
        genotype_gvcfs.out
}

process haplotype_caller {
    label 'M2_NR'
    publishDir "progress/haplotype_caller", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am:$ps" }

    input:
        tuple val(am), val(reg), val(sm), val(ps), val(nr), path(bam), path(bai)
        tuple path(ref), path(fai), path(dict)
        val ploidy

    output:
        tuple val(am), val(reg), path(gvcf), path("${gvcf}.tbi")

    script:
        gvcf = "SM-${sm}.AM-${am}.PS-${ps}.gvcf.gz"
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
        val qd

    output:
        tuple val(am), path(vcf), path("${vcf}.tbi")

    script:
        vcf = "${am}.vcf.gz"
        """
        gatk CombineGVCFs \\
            --java-options "-Xmx4G -Djava.io.tmpdir=." \\
            -R $ref \\
            -L $reg \\
            -O combined.g.vcf.gz \\
            --disable-sequence-dictionary-validation \\
            --variant ${gvcfs.join(' --variant ')}

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
            bcftools norm -m-both -f $ref -Ou |
            bcftools view -i 'GT="alt"' -Oz -o $vcf

        bcftools index -t $vcf
        """
}

