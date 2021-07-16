

process index_bam {
    label 'S_NR'
    publishDir "output/bam", mode: params.output_pub_mode
    tag { "$sm:$am" }

    input:
        tuple val(sm), val(am), val(nr), path(bam)

    output:
        tuple val(sm), val(am), val(nr), path(bam), path("${bam}.bai")

    script:
        """
        samtools index $bam
        """
}