
process extract_ccs_failed {
    label 'M'
    publishDir "progress/extract_ccs_failed", mode: "$params.intermediate_pub_mode"

    input:
        tuple path(subreads_bam), path(ccs_bam)

    output:
        tuple path("${pref}.nr"), path("${pref}.bam")

    script:
        pref = params.run_id + ".failed_ccs_subreads"
        """
        samtools view -u $subreads_bam |
            extract_pbbam_subset.py - $ccs_bam -u --out - --count ${pref}.nr |
            samtools view -b -o ${pref}.bam
        """
}