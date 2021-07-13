#!/usr/bin/env nextflow

//import static Helpers.*
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

// import functions and tasks
include { path; checkManiAmps; refFastaFileMap} from './functions'
include { pb_ccs } from './tasks/pre-processing/pb_ccs'
include { pb_merge } from './tasks/pre-processing/pb_merge'
include { pb_lima } from './tasks/pre-processing/pb_lima'
include { pb_mm2; pb_mm2 as pb_mm2_2 } from './tasks/pre-processing/pb_mm2'
include { merge_lima_smry } from './tasks/pre-processing/merge_lima_smry'
include { extract_barcode_set } from './tasks/pre-processing/extract_barcode_set'
include { extract_ccs_failed } from './tasks/pre-processing/extract_ccs_failed'
include { annotate_samples } from './tasks/pre-processing/annotate_samples'
include { annotate_amplicons } from './tasks/pre-processing/annotate_amplicons'
include { alignment_stats } from './tasks/pre-processing/alignment_stats'
include { split_sample_amplicons } from './tasks/pre-processing/split_sample_amplicons'
include { index_bam } from './tasks/pre-processing/index_bam'

// check and load inputs
subreads_bam = path(params.subreads_bam)
subreads_pbi = path(params.subreads_bam + '.pbi')
sample_manifest = path(params.sample_manifest)
barcodes_fasta = path(params.barcodes_fasta)
amplicons_json = path(params.amplicons_json)
checkManiAmps(sample_manifest, amplicons_json)
path(params.ref_fasta)
ref = refFastaFileMap(params.ref_fasta)

// main workflow
workflow {
    ref_channel = Channel.from([ref]).first()

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
    merge_lima_smry(pb_lima.out.smry)

    pb_lima.out.bams |
        combine (ref_channel) |
        pb_mm2 |
        combine(extract_barcode_set.out.order) |
        map { it + [sample_manifest] } |
        annotate_samples |
        map { it + [amplicons_json, sample_manifest] } |
        annotate_amplicons |
        combine (ref_channel) |
        pb_mm2_2

    pb_mm2_2.out.bams |
        filter { it[0] == 'CCS' & it[1] } |
        map { it.drop(2) } |
        split_sample_amplicons |
        index_bam |
        view

    pb_mm2_2.out.bams |
        alignment_stats |
        toSortedList()

}
