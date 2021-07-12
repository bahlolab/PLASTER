#!/usr/bin/env nextflow

import static Helpers.*
nextflow.enable.dsl=2

wf_name = "PLASTER: pre-processing"
println "\n------- $wf_name -------\n"

// default params
params.intermediate_pub_mode = 'symlink'
params.output_pub_mode = 'copy'
params.run_id = "plaster-run"
params.ccs_min_len = 250
params.ccs_max_len = 25000
params.ccs_min_acc = 0.99
params.ccs_min_passes = 3
params.ccs_n_parallel = 10

// import functions and processes
include { path; checkManiAmps } from './functions'
include { pb_ccs } from './tasks/pre-processing/pb_ccs'
include { pb_merge } from './tasks/pre-processing/pb_merge'
include { pb_lima } from './tasks/pre-processing/pb_lima'
include { extract_barcode_set } from './tasks/pre-processing/extract_barcode_set'
include { extract_ccs_failed } from './tasks/pre-processing/extract_ccs_failed'

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

    ccs_bam = params.ccs_n_parallel == 1 ?
        pb_ccs.out.bams.map{ it[1] }.first() :
        pb_ccs.out.bams | pb_merge

    lima_in = extract_ccs_failed(subreads_bam, ccs_bam) |
        map { ['SR', it] } |
        mix(ccs_bam.map { ['CCS', it] })

    pb_lima(lima_in, extract_barcode_set.out.fasta)

}
