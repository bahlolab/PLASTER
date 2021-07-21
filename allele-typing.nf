#!/usr/bin/env nextflow

nextflow.enable.dsl=2

wf_name = "PLASTER: allele-typing"
println "\n------- $wf_name -------\n"

// default params
params.intermediate_pub_mode = 'symlink'
params.output_pub_mode = 'copy'
params.run_id = "plaster-run"
params.max_reads = 1000
params.min_reads = 20
params.min_reads_phased = 5
params.ploidy = 2
params.qd_1 = '2.0'
params.qd_2 = '20.0'
params.copy_num = null
params.vep_cache_ver = '104'
params.vep_assembly = 'GRCh38'

// import functions and tasks
include { path; parseManifestAT; checkManiAmpsAT; parseCopyNum; readTSV } from './functions'
include { prep_ref } from './tasks/prep_ref'
include { prep_bams } from './tasks/allele-typing/prep_bams'
include { gatk as gatk_1; gatk as gatk_2 } from './tasks/allele-typing/gatk'
include { phase } from './tasks/allele-typing/phase'
include { targeted } from './tasks/allele-typing/targeted'
include { vep } from './tasks/allele-typing/vep'

// check and load inputs
manifest = parseManifestAT(params.manifest)
amplicons_json = path(params.amplicons_json)
amplicons = checkManiAmpsAT(manifest.collect{it.amplicon}.unique(), amplicons_json)
targ = amplicons.any { it.size() == 4 }
copy_num = parseCopyNum(manifest, params.copy_num, params.ploidy)

// main workflow
workflow {
    ref = prep_ref(params.ref_fasta, 'fai')

    amps = Channel.fromList(amplicons)

    bams = Channel.fromList(manifest) |
        map { (it.values() as ArrayList)[1..4] } |
        prep_bams

    gatk_1(bams, amps.map{ it.take(2) }, ref,
           qd: params.qd_1, targ: false)

    phase(bams, gatk_1.out, copy_num, ref)

    gatk_2(phase.out, amps, ref,
           ploidy: 1, qd: params.qd_2, targ: targ) |
        view


//    targ_channel = Channel.fromList(amplicons_targ)
//    (targ_channel |
//        combine(amp_channel, by:0) |
//        combine(phase.out, by: 0))
//        .with { targeted(it, ref) }

//    gatk_2.out |
//        map { ['gatk'] + it} |
//        vep
}