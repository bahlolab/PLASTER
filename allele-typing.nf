#!/usr/bin/env nextflow

nextflow.enable.dsl=2

wf_name = "PLASTER: allele-typing"
println "\n------- $wf_name -------\n"

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
include { prep_ref } from './tasks/prep_ref'
include { get_pharmvar_vcf } from './tasks/allele-typing/get_pharmvar_vcf'
include { prep_bams } from './tasks/allele-typing/prep_bams'
include { gatk as gatk_1; gatk as gatk_2 } from './tasks/allele-typing/gatk'
include { phase } from './tasks/allele-typing/phase'
include { targeted } from './tasks/allele-typing/targeted'
include { vep } from './tasks/allele-typing/vep'

// check and load inputs
manifest = parseManifestAT(params.manifest)
amplicons_json = path(params.amplicons_json)
(amplicons, fusion, pharmvar) = checkManiAmpsAT(manifest.collect{it.amplicon}.unique(), amplicons_json)
copy_num = parseCopyNum(manifest, params.copy_num, params.ploidy)

// main workflow
workflow {
    ref = prep_ref(params.ref_fasta, 'fai')

    pv_vcf = pharmvar ?
        get_pharmvar_vcf(Channel.fromList(pharmvar), ref) :
        false

    amps = Channel.fromList(amplicons)

    (Channel.fromList(manifest) |
        map { (it.values() as ArrayList)[1..4] })
        .with { prep_bams (it, ref, fusion)}

    gatk_1(prep_bams.out, amps.map{ it.take(2) }, ref,
           qd: params.qd_1, targ: false)

    phase(prep_bams.out, gatk_1.out, channel.fromList(copy_num), ref)

    gatk_2(phase.out.bams, amps, ref,
           ploidy: 1, qd: params.qd_2, targ: pv_vcf) |
        vep

    phase.out.smry | view
}