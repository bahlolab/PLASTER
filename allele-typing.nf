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
amplicons_targ = amplicons
    .collect { k, v -> [k, v.target_vcf] }
    .findAll {it[1] }
    .collect { [it[0], path(it[1]), path(it[1] + '.tbi')] }
copy_num = parseCopyNum(manifest, params.copy_num, params.ploidy)

// main workflow
workflow {
    ref = prep_ref(params.ref_fasta, 'fai')

    amp_channel = amplicons
        .collect { k, v -> [k, "$v.chrom:$v.start-$v.end"] }
        .with { Channel.fromList(it) }

    bams = Channel.fromList(manifest) |
        map { (it.values() as ArrayList)[1..4] } |
        prep_bams

    (bams |
        combine(amp_channel, by: 0) |
        map { it[[0, 5, 1]] + ['NA'] + it[[2, 3, 4]] } )
        .with{ gatk_1(it, ref, params.ploidy, params.qd_1) }

    phase(bams, gatk_1.out, copy_num, ref)

    (amp_channel |
        combine(phase.out, by: 0))
        .with { gatk_2(it, ref, 1, params.qd_2) }

    targ_channel = Channel.fromList(amplicons_targ)
    (targ_channel |
        combine(amp_channel, by:0) |
        combine(phase.out, by: 0))
        .with { targeted(it, ref) }

    gatk_2.out |
        map { ['gatk'] + it} |
        mix(targeted.out | map { ['targeted'] + it }) |
        vep
}