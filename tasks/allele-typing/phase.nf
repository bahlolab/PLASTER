
workflow phase {
    take:
        bams // am, sm, nr, bam, bai
        vcfs // am, vcf, tbi
    main:
        get_snp_pos(vcfs)
        out = get_snp_pos.out
    emit:
        out
}


process get_snp_pos {
    label 'S'
    publishDir "progress/get_snp_pos", mode: params.output_pub_mode
    tag { am }

    input:
    tuple val(am), path(vcf), path(tbi)

    output:
    tuple val(am), path(out)

    script:
        out = "${am}.snp_pos.tsv"
        """
        bcftools view $vcf -i 'GT="alt"' -v snps -M2 -HG | cut -f1,2 > $out
        """
}

