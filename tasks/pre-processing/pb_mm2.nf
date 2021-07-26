
process pb_mm2 {
    label 'XL2_NR'
    publishDir "progress/pb_mm2", mode: "$params.intermediate_pub_mode"
    tag { "$rt:$is_bc" }

    input:
        tuple val(rt), val(is_bc), val(nr), path(bam), path(mmi)

    output:
        tuple val(rt), val(is_bc), val(nr), file(aln), emit: bams

    script:
        aln = bam.name.replace('.bam', '.mm2.bam')
        """
        pbmm2 align $bam $mmi tmp.bam \\
            --num-threads $task.cpus \\
            --preset ${rt ==~ /^CCS/ ? 'CCS' : 'SUBREAD'} \\
            --best-n 1 \\
            --unmapped
        samtools sort tmp.bam -@$task.cpus -m 1G -O BAM -o $aln
        """
}