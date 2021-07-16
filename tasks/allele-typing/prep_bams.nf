
workflow prep_bams {
    take:
        // sm, am, nr, bam
        data
    main:
        bams = data |
            groupTuple(by:[0,1], sort: true) |
            map { it[2] = it[2].sum(); it } |
            branch { multi: it[3].size() > 1; single: true }
        bams = bams.multi |
            merge |
            mix(bams.single.map { it[3] = it[3][0]; it }) |
            branch { downsample: it[2] > params.max_reads; pass_through: true }
        bams = bams.downsample |
            downsample |
            mix(bams.pass_through) |
            index
    emit:
        // sm, am, nr, bam, bai
        bams
}

process merge {
    label 'S_NR'
    publishDir "progress/merge_bams", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am" }

    input:
        tuple val(sm), val(am), val(nr), path(bams)

    output:
        tuple val(sm), val(am), val(nr), path(merged)

    script:
        merged = "SM-${sm}.AM-${am}.merged.bam"
        """
        for BAM in *.bam; do samtools index \$BAM; done
        samtools merge $merged *.bam
        """
}

process downsample {
    label 'S_NR'
    publishDir "progress/downsample_bams", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am" }

    input:
        tuple val(sm), val(am), val(nr), path(bam)

    output:
        tuple val(sm), val(am), val(params.max_reads), path(out_bam)

    script:
        out_bam = bam.name.replaceAll('.bam', '.sub.bam')
        """
        samtools index $bam
        downsample_bam.py $bam -u --count $params.max_reads | samtools view -b -o $out_bam
        """
}

process index {
    label 'S_NR'
    publishDir "progress/index_bams", mode: "$params.intermediate_pub_mode"
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