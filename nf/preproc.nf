#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// default params
params.run_id = "plaster-run"
params.ccs_min_len = 250
params.ccs_max_len = 25000
params.ccs_min_acc = 0.99
params.ccs_min_passes = 3
params.ccs_n_parallel = 100

// import functions and tasks
include { path; checkManiAmps } from './functions'
include { prep_ref } from './prep_ref'
include { pb_mm2_index } from './preproc/pb_mm2_index'
include { pb_ccs } from './preproc/pb_ccs'
include { pb_merge } from './preproc/pb_merge'
include { pb_lima } from './preproc/pb_lima'
include { pb_mm2; pb_mm2 as pb_mm2_2 } from './preproc/pb_mm2'
include { extract_barcode_set } from './preproc/extract_barcode_set'
include { extract_ccs_failed } from './preproc/extract_ccs_failed'
include { annotate_samples } from './preproc/annotate_samples'
include { annotate_amplicons } from './preproc/annotate_amplicons'
include { alignment_stats } from './preproc/alignment_stats'
include { split_sample_amplicons } from './preproc/split_sample_amplicons'
include { index_bam } from './preproc/index_bam'
include { pre_processing_report } from './preproc/pre_processing_report'

// main workflow
workflow preproc {

    println "\n------- PLASTER: pre-processing -------\n${params.test ? 'Running test dataset' : ''}"

    // check and load inputs
    subreads_bam = path(params.subreads_bam)
    subreads_pbi = path(params.subreads_bam + '.pbi')
    sample_manifest = path(params.sample_manifest)
    barcodes_fasta = path(params.barcodes_fasta)
    amplicons_json = path(params.amplicons_json)
    checkManiAmps(sample_manifest, amplicons_json)
    rmd = file(workflow.projectDir + '/bin/preproc-report.Rmd')

    // run tasks
    prep_ref(params.ref_fasta, 'mmi')

    extract_barcode_set(sample_manifest, barcodes_fasta)

    Channel.from([[subreads_bam, subreads_pbi]]) |
        pb_ccs |
        map { [subreads_bam, it] } |
        extract_ccs_failed |
        combine(['SR']) |
        mix(pb_ccs.out | combine(['CCS'])) |
        combine(extract_barcode_set.out.fasta) |
        pb_lima

    pb_lima.out.bams |
        combine(prep_ref.out, by:0) |
        pb_mm2 |
        combine(extract_barcode_set.out.order) |
        map { it + [sample_manifest] } |
        annotate_samples |
        map { it + [amplicons_json, sample_manifest] } |
        annotate_amplicons |
        combine(prep_ref.out, by:0) |
        pb_mm2_2 |
        filter { it[0] == 'CCS' & it[1] } |
        map { it.drop(2) } |
        split_sample_amplicons |
        index_bam |
        map { [params.run_id] + it.dropRight(1) } |
        map { it.take(4) + [file(file("./output/bam").toRealPath().toString() + '/'+ it[4].fileName)] } |
        map { it.collect { it.toString() }.join('\t') } |
        collectFile(name: 'sample_amplicon_bam_manifest.tsv', storeDir: './output/', newLine: true,
            seed: ['run_id', 'sample', 'amplicon', 'n_reads', 'bam_file'].join('\t'))

    pb_mm2_2.out.bams |
        alignment_stats |
        toSortedList() |
        map { [it] } |
        combine(pb_lima.out.smry) |
        map { [rmd, amplicons_json, sample_manifest] + it } |
        pre_processing_report
}
