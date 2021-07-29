#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// default params
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
include { prep_ref } from './prep_ref'
include { get_pharmvar_vcf } from './typing/get_pharmvar_vcf'
include { prep_bams } from './typing/prep_bams'
include { gatk as gatk_1; gatk as gatk_2 } from './typing/gatk'
include { phase } from './typing/phase'
include { targeted } from './typing/targeted'
include { vep } from './typing/vep'
include { pharmvar_star_allele } from './typing/pharmvar_star_allele'

// main workflow
workflow typing {

    println "\n------- PLASTER: allele-typing -------\n${params.test ? 'Running test dataset' : ''}"

    // check and load inputs
    manifest = parseManifestAT(params.manifest)
    amplicons_json = path(params.amplicons_json)
    (amplicons, fusion, pharmvar) = checkManiAmpsAT(manifest.collect{it.amplicon}.unique(), amplicons_json)
    copy_num = parseCopyNum(manifest, params.copy_num, params.ploidy)

    // run tasks
    ref = prep_ref(params.ref_fasta, 'fai')

    pv_vcf = pharmvar ?
        get_pharmvar_vcf(Channel.fromList(pharmvar), ref) :
        false

    amps = Channel.fromList(amplicons)

    (Channel.fromList(manifest) |
        map { (it.values() as ArrayList)[0..3] })
        .with { prep_bams (it, ref, fusion)}

    gatk_1(prep_bams.out, amps.map{ it.take(2) }, ref,
           qd: params.qd_1, targ: false)

    phase(prep_bams.out, gatk_1.out, channel.fromList(copy_num), ref)

    gatk_2(phase.out, amps, ref,
           ploidy: 1, qd: params.qd_2, targ: pv_vcf) |
        vep

    if (pharmvar) {
        vep.out |
            map { it[0..1] } |
            combine(get_pharmvar_vcf.out, by:0) |
            map { it[0..2] } |
            combine(Channel.fromList(pharmvar), by: 0) |
            pharmvar_star_allele
    }
}