#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { wget_ref } from './wget_ref'
//include { samtools_index_ref } from './samtools_index_ref'
include { mm2_index_ref } from './mm2_index_ref'

workflow prepare_reference {
    take:
        ref_fasta
    main:
        ref = ref_fasta ==~ '^(ftp|https?)://.+' ?
            wget_ref(ref_fasta) :
            Channel.from(path(ref))
        mm2_index_ref(ref)
    emit:
        mmi = mm2_index_ref.out.ccs_mmi.map {['CCS', it]} |
            mix(mm2_index_ref.out.subread_mmi.map {['SR', it]})
}
