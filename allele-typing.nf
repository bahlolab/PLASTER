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
params.qd_1 = '2.0'
params.qd_2 = '20.0'


// import functions and tasks
include { path; parseManifestAT; checkManiAmpsAT } from './functions'
include { prep_ref } from './tasks/prep_ref'
include { prep_bams } from './tasks/allele-typing/prep_bams'
include { gatk as gatk_1; gatk as gatk_2 } from './tasks/allele-typing/gatk'

// check and load inputs
manifest = parseManifestAT(params.manifest)
amplicons_json = path(params.amplicons_json)
amplicons = checkManiAmpsAT(manifest.collect{it.amplicon}.unique(), amplicons_json)

// main workflow
workflow {
    prep_ref(params.ref_fasta, 'fai')

    amp_channel = amplicons
        .collect { k, v -> [k, "$v.chrom:$v.start-$v.end"] }
        .with{ Channel.fromList(it) }


    (Channel.fromList(manifest) |
        map { (it.values() as ArrayList)[1..4] } |
        prep_bams |
        combine(amp_channel, by:0) |
        map { it[[0, 5, 1]] + ['NA'] + it[[2, 3, 4]] })
        .with { gatk_1(it, prep_ref.out, 1, params.qd_1) } |
        view
}
