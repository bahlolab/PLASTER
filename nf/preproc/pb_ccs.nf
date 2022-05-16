
workflow pb_ccs {
    take:
        subreads_channel
    main:
        Channel.from((1..params.ccs_n_parallel) as ArrayList) |
            combine(subreads_channel) |
            ccs
        ccs_bam =  ccs.out.bams |
            toSortedList() |
            map { it.collect { it[1] } } |
            merge |
            map { [it[0].toFile().text.trim() as int, it[1]] }
    emit:
        ccs_bam // nr, bam
}

process ccs {
    label 'L2'
    publishDir "progress/pb_ccs", mode: "$params.intermediate_pub_mode"
    tag { i }

    input:
        tuple val(i), path(subreads_bam), path(subreads_pbi)

    output:
        tuple val(i), path(ccs_bam), emit: bams
        tuple val(i), path(report), emit: reports

    script:
        prefix = params.run_id + (params.ccs_n_parallel > 1 ? "_${i}" : "")
        ccs_bam = prefix + ".ccs.bam"
        report = prefix + ".ccs_report.txt"
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

process merge {
    label 'M'
    publishDir "progress/pb_merge", mode: "$params.intermediate_pub_mode"

    input:
    path bams

    output:
    tuple path("${pref}.nr"), path("${pref}.bam")

    script:
    pref = params.run_id + '.ccs_merged'
    """
    pbmerge ${bams.join(' ')} -o ${pref}.bam
    pbindexdump ${pref}.bam.pbi | jq .numReads > ${pref}.nr
    """
}