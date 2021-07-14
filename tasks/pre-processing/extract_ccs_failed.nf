
process extract_ccs_failed {
    label 'M'
    publishDir "progress/extract_ccs_failed", mode: "$params.intermediate_pub_mode"

    input:
        path subreads_bam
        path ccs_bam

    output:
        path out, emit: ccs_failed_subreads

    script:
        out = params.run_id + ".failed_ccs_subreads.bam"
        """
        samtools view -u $subreads_bam |
            extract_pbbam_subset.py - $ccs_bam -u --out - |
            samtools view -b -o $out
        """
}