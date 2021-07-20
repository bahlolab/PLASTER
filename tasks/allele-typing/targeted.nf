
workflow targeted {
    take:
        data // am, vcf, tbi, reg, sm, ps, nr, bam, bai
        ref // fasta, fai, dict
    main:
        haplotype_caller(data, ref) |
            groupTuple(by: 0, sort: true) |
            map { it[[0,3,4]] } |
            merge_vcf
    emit:
        merge_vcf.out //am, vcf, tbi
}

process haplotype_caller {
    label 'M2_NR'
    publishDir "progress/haplotype_caller_targ", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am:$ps" }

    input:
        tuple val(am), path(sites), path(tbi), val(reg), val(sm), val(ps), val(nr), path(bam), path(bai)
        tuple path(ref), path(fai), path(dict)

    output:
        tuple val(am), val(sm), val(ps), path(out), path("${out}.tbi")

    script:
        out = "SM-${sm}.AM-${am}.PS-${ps}.vcf.gz"
        """
        bcftools norm -m+any -f $ref $sites -Oz -o multi.vcf.gz
        bcftools index -t multi.vcf.gz
        gatk HaplotypeCaller \\
            --java-options "-Xmx4G -Djava.io.tmpdir=." \\
            -R $ref \\
            -I $bam \\
            -O tmp.vcf.gz \\
            --alleles multi.vcf.gz \\
            --force-call-filtered-alleles \\
            --intervals $reg \\
            --sample-ploidy 1
        bcftools norm -m-any -f $ref tmp.vcf.gz -Oz -o tmp2.vcf.gz
        bcftools index -t tmp2.vcf.gz
        bcftools isec -p tmp -Oz -w 2 $sites tmp2.vcf.gz
        mv tmp/0003.vcf.gz $out
        bcftools index -t $out
        """
}

process merge_vcf {
    label 'M2'
    publishDir "progress/merge_targ", mode: "$params.intermediate_pub_mode"
    tag { am }

    input:
        tuple val(am), path(vcf), path(tbi)

    output:
        tuple val(am), path(out), path("${out}.tbi")

    script:
        out = "${am}.targeted.vcf.gz"
        """
        bcftools merge -Oz ${vcf.join(' ')} -o $out
        bcftools index -t $out
        """
}

