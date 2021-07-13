

process alignment_stats {
    label 'M_NR'
    publishDir "intermediates/alignment_stats", mode: "$params.intermediate_pub_mode"
    tag { "$rt:$is_bc" }

    input:
        tuple val(rt), val(is_bc), val(nr), path(bam)

    output:
        path stats, emit: stats

    script:
        stats = bam.name.replace('.bam', '.stats.tsv.gz')
        """
        samtools view -u $bam |
            alignment_stats.py - --zmw  --ext read_type:$rt \\
                --tag SM --tag BC --tag AM --tag OL --tag ID --tag PM --tag PO --tag PP |
            gzip > $stats
        """
}