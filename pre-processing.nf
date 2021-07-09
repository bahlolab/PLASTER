#!/usr/bin/env nextflow

import static Helpers.*
nextflow.enable.dsl=2

wf_name = "PLASTER: pre-processing"
println "\n------- $wf_name -------"

// default params
params.intermediate_pub_mode = 'symlink'
params.output_pub_mode = 'copy'
params.ccs_min_len = 250
params.ccs_max_len = 25000
params.ccs_min_acc = 0.99
params.ccs_min_passes = 3
params.ccs_n_parallel = 10

include { path; checkManiAmps } from './functions'
include { extract_barcode_set } from './tasks/pre-processing/extract_barcode_set'
include { pb_ccs } from './tasks/pre-processing/pb_ccs'

// check and load inputs
subreads_bam = path(params.subreads_bam)
subreads_pbi = path(params.subreads_bam + '.pbi')
sample_manifest = path(params.sample_manifest)
barcodes_fasta = path(params.barcodes_fasta)
amplicons_json = path(params.amplicons_json)
checkManiAmps(sample_manifest, amplicons_json)

// main workflow
workflow {
    extract_barcode_set(sample_manifest, barcodes_fasta)

    Channel.from((1..params.ccs_n_parallel) as ArrayList) |
        map { [it, subreads_bam, subreads_pbi ] } |
        pb_ccs
}
