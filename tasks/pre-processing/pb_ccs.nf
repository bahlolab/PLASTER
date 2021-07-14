
process pb_ccs {
    label 'L2'
    publishDir "progress/pb_ccs", mode: "$params.intermediate_pub_mode"
    tag { i }

    input:
        tuple val(i), path(subreads_bam), path(subreads_pbi)

    output:
        tuple val(i), path(ccs_bam), emit: bams
        tuple val(i), path(report), emit: reports

    script:
        ccs_bam = "${params.run_id}_${i}.ccs.bam"
        report = "${params.run_id}_${i}.ccs_report.txt"
        """
        ccs $subreads_bam $ccs_bam \\
            --chunk $i/$params.ccs_n_parallel \\
            --num-threads $task.cpus \\
            --max-length $params.ccs_max_len \\
            --min-length $params.ccs_min_len \\
            --min-rq $params.ccs_min_acc \\
            --min-passes $params.ccs_min_passes \\
            --report-file $report
        """
}