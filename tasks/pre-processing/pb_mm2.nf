
process pb_mm2 {
    label 'XL2_NR'
    publishDir "intermediates/pb_mm2", mode: "$params.intermediate_pub_mode"
    tag { "$rt:$is_bc" }

    input:
        tuple val(rt), val(is_bc), val(nr), path(bam), val(ref_map), path(ref_files)

    output:
        tuple val(rt), val(is_bc), val(nr), file(aln), emit: bams

    script:
        aln = bam.name.replace('.bam', '.mm2.bam')
        ref_index = rt ==~ /^CCS$/ ? ref_map['ccs.mmi'] : ref_map['subread.mmi']
        """
        pbmm2 align $bam $ref_index tmp.bam \\
            --num-threads $task.cpus \\
            --preset ${rt ==~ /^CCS/ ? 'CCS' : 'SUBREAD'} \\
            --best-n 1 \\
            --unmapped
        samtools sort tmp.bam -@$task.cpus -m 1G -O BAM -o $aln
        """
}