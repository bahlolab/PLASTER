
workflow annotate_samples {
    take:
        data
    main:
        data = data.branch { is_bc: it[1] ;  not_bc: true }
        data.is_bc | AS
    emit:
        bams = AS.out.bams |
            mix( data.not_bc.map {it.take(4)} )
}

process AS {
    label 'M_NR'
    publishDir "progress/annotate_samples", mode: "$params.intermediate_pub_mode"
    tag { "$rt:$is_bc" }

    input:
        tuple val(rt), val(is_bc), val(nr), path(bam), path(barcodes_fa)

    output:
        tuple val(rt), val(is_bc), val(nr), path(out), emit: bams

    script:
        out = params.run_id + '.' + rt + '.sm_annot.bam'
        """
        samtools view -u $bam |
            bam_annotate_samples.py - \\
                --barcodes $barcodes_fa \\
                --out $out
        """
}