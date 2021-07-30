#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { preproc } from './nf/preproc'
include { typing } from './nf/typing'

workflow {
    if (params.mode == 'preproc') {
        preproc()
    }
    else if (params.mode == 'typing') {
        typing()
    } else {
        println "Error: No stage specified - must specify either 'preproc' or 'typing'"
        System.exit(1)
    }
}
