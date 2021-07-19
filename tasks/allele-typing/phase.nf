
include { readTSV } from '../../functions'

workflow phase {
    take:
        bams // am, sm, nr, bam, bai
        vcfs // am, vcf, tbi
        copy_num //am, sm, cn
        ref // fasta, fai, dict
    main:
        out =
            (get_snp_pos(vcfs) |
            combine(bams, by: 0))
            .with { assign_snps(it, ref )} |
                combine(copy_num, by: 0..1 ) |
                amp_phaser
        // TODO - report phased.fail
        phased = amp_phaser.out.phases |
            map { it + [
                readTSV(it[3], ['phase', 'count', 'freq', 'ratio', 'copy_num'])
                    .collect { it.count as int }.min() ]} |
            branch { pass: it[4] >= params.min_reads_phased; fail: true}

        out = phased.pass |
            map { it.take(3) } |
            combine(bams, by: 0..1) |
            split_phases |
            transpose |
            map { it.take(2) +
                [(it[3] =~ '([0-9]+)\\.bam$')[0][1], it[2].toFile().text as int] +
                it.takeRight(2) }

    emit:
        out // am, sm, ps, nr, bam, bai
}


process get_snp_pos {
    label 'S'
    publishDir "progress/get_snp_pos", mode: params.intermediate_pub_mode
    tag { am }

    input:
    tuple val(am), path(vcf), path(tbi)

    output:
    tuple val(am), path(out)

    script:
        out = "${am}.snp_pos.tsv"
        """
        bcftools view $vcf -i 'GT="alt"' -v snps -M2 -HG | cut -f1,2 > $out
        """
}

process assign_snps {
    label 'M_NR'
    publishDir "progress/assign_snps", mode: params.intermediate_pub_mode
    tag { "$sm:$am" }

    input:
        tuple val(am), path(pos), val(sm), val(nr), path(bam), path(bai)
        tuple path(ref), path(fai), path(dict)

    output:
        tuple val(am), val(sm), val(nr), file(vcf)

    script:
    vcf = "SM-${sm}.read_vars.vcf.gz"
    """
    bam_qname_to_rg_sm.py $bam --out tmp.bam
    samtools index tmp.bam
    bcftools mpileup tmp.bam -d 5000 --threads $task.cpus -AIB -a FMT/AD,FMT/DP -Q 0 -f $ref -R $pos -Ou |
        bcftools call -A -m --prior 1e-1 --ploidy 1 -Oz > $vcf
    """
}

process amp_phaser {
    label 'S2_NR'
    publishDir "progress/AmpPhaseR", mode: 'symlink'
    publishDir "output/AmpPhaseR", mode: 'copy', pattern: '*.png'
    publishDir "output/AmpPhaseR", mode: 'copy', pattern: '*_read_phase_smry.tsv.gz'
    tag { "$sm:$am" }

    input:
        tuple val(am), val(sm), val(nr), path(vcf), val(cn)

    output:
        tuple val(am), val(sm), path("${pref}_phased.tsv.gz"), path("${pref}_phase_summary.tsv"), emit: phases
        tuple val(am), val(sm), path("${pref}_phase_plot.png"), path("${pref}_read_phase_smry.tsv.gz"), emit: report

    script:
    pref = "${sm}_${am}"
    """
    AmpPhaseR_wrapper.R $vcf $pref \\
        --copy-num $cn \\
        --max-breakpoints 2 \\
        --max-allele-freq 0.95 \\
        --min-phase-freq 0.05 \\
        --max-freq-abs-delta 0.75 \\
        --max-freq-rel-delta 0.75
    """
}

process split_phases {
    label 'M_NR'
    publishDir "progress/split_phases", mode: params.intermediate_pub_mode
    tag { "$sm:$am" }

    input:
        tuple val(am), val(sm), file(read_phase), val(nr), path(bam), path(bai)

    output:
        tuple val(am), val(sm), path("${pref}_*.bam.count"), path("${pref}_*.bam"), path("${pref}_*.bam.bai")

    script:
        pref = "SM-${sm}.AM-${am}.phase"
        """
        samtools view -u $bam | bam_annotate_phases.py - $read_phase --out ${pref}.bam --update-rg
        samtools split ${pref}.bam -f "%*_%!.%."
        for BAM in ${pref}_*.bam; do samtools index \$BAM; samtools view \$BAM | wc -l > \$BAM.count; done
        """
}

