
workflow pb_merge {
    take:
        data
    main:
        data |
            toSortedList() |
            map { it.collect { it[1] } } |
            PBM
    emit:
        bam = PBM.out.bam
}

process PBM {
    label 'M'
    publishDir "progress/pb_merge", mode: "$params.intermediate_pub_mode"

    input:
        path bams

    output:
        path merged, emit: bam

    script:
        merged = params.run_id + '.ccs_merged.bam'
        """
        pbmerge ${bams.join(' ')} -o $merged --no-pbi
        """
}