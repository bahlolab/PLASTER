
workflow split_sample_amplicons {
    take:
        data
    main:
        SSA(data)
    emit:
        // sm, am, nr, bam
        bams = SSA.out.bams.flatten().map{ [(it =~ '([^/]+).bam$')[0][1], it] } |
            combine(SSA.out.counts.flatten().map{ [(it =~ '([^/]+).bam.count$')[0][1], it] }, by:0) |
            map { it.drop(1) } |
            map { (it[0] =~ /LB-(.+)\.SM-(.+)\.AM-(.+)\.bam/)[0][2..3] +
                [it[1].toFile().text.trim() as int, it[0]] }
}

process SSA {
    label 'M_NR'
    publishDir "output/bam", mode: params.output_pub_mode, pattern: '*.bam'

    input:
        tuple val(nr), path(bam)

    output:
        path "*LB-*.SM-*.AM-*.bam", emit: bams
        path "*LB-*.SM-*.AM-*.bam.count", emit: counts

    script:
        """
        samtools view -u $bam |
            bam_split_sample_amplicons.py - --lb-tag $params.run_id
        """
}