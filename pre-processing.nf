#!/usr/bin/env nextflow

import static Helpers.*
nextflow.enable.dsl=2

wf_name = "PLASTER: pre-processing"
println "\n------- $wf_name -------"

include { extract_barcode_set } from './tasks/pre-processing/extract_barcode_set'

// default params
params.intermediate_pub_mode = 'symlink'
params.output_pub_mode = 'copy'

// check and load inputs
sample_manifest = file(params.sample_manifest, checkIfExists: true)
barcodes_fasta = file(params.barcodes_fasta, checkIfExists: true)

// main workflow
workflow {
    println 'starting'

    extract_barcode_set(sample_manifest, barcodes_fasta)

    println extract_barcode_set

    println 'finishing'
}
