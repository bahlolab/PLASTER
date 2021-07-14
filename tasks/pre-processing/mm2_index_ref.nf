

process mm2_index_ref {
    label 'M'
    publishDir "output/ref", mode: params.output_pub_mode

    input:
        path ref

    output:
        path ccs_mmi, emit: ccs_mmi
        path subread_mmi, emit: subread_mmi

    script:
        base = ref.toString().replaceAll('.fa(sta)?$', '')
        ccs_mmi = base + '.ccs.mmi'
        subread_mmi = base + '.subread.mmi'
        """
        pbmm2 index $ref $ccs_mmi --preset CCS --num-threads $task.cpus
        pbmm2 index $ref $subread_mmi --preset SUBREAD --num-threads $task.cpus
        """
}