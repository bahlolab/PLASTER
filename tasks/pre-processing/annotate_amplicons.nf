

process annotate_amplicons {
    label 'M_NR'
    publishDir "intermediates/annotate_amplicons", mode: "$params.intermediate_pub_mode"
    tag { "$rt:$is_bc" }

    input:
        tuple val(rt), val(is_bc), val(nr), path(bam), path(amplicons), path(manifest)

    output:
        tuple val(rt), val(is_bc), val(nr), path(out), emit: bams

    script:
        out = "$params.run_id.$rt.${is_bc}.sm_am_annot.bam"
        """
        samtools view -u $bam |
            bam_annotate_amplicons.py - --window 500 --max-dist 2 \\
                --manifest $manifest --amplicons $amplicons --out $out
        """
}