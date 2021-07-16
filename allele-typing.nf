#!/usr/bin/env nextflow

nextflow.enable.dsl=2

wf_name = "PLASTER: allele-typing"
println "\n------- $wf_name -------\n"

// default params
params.intermediate_pub_mode = 'symlink'
params.output_pub_mode = 'copy'
params.run_id = "plaster-run"
params.max_reads = 1000
params.default_ploidy = 2
params.min_reads_phased = 5
params.min_qd_1 = '2.0'
params.min_qd_2 = '20.0'


// import functions and tasks
include { path; parseManifestAT; checkManiAmpsAT } from './functions'
include { prep_ref } from './tasks/prep_ref'
include { prep_bams } from './tasks/allele-typing/prep_bams'

// check and load inputs
manifest = parseManifestAT(params.manifest)
amplicons_json = path(params.amplicons_json)
checkManiAmpsAT(manifest.collect{it.amplicon}.unique(), amplicons_json)

// main workflow
workflow {
    prep_ref(params.ref_fasta, 'fai')

    bams = Channel.fromList(manifest) |
        map { (it.values() as ArrayList)[1..4] } |
        prep_bams |
        view
}
