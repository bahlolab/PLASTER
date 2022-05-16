
include { path } from '../functions'

workflow annotate_samples {
    take:
        data
    main:
        data = data.branch { is_bc: it[1] ;  not_bc: true }

        if (params.barcodes_fasta) {
            barcodes_fasta = path(params.barcodes_fasta)
            bams = data.is_bc |
                combine([barcodes_fasta]) |
                with_fasta |
                mix( data.not_bc.map {it.take(4)} )
        } else {
            sample = params.sample ?: "sample"
            bams = data.is_bc |
                combine([sample]) |
                single_sample
        }

    emit:
        bams
}

process with_fasta {
    label 'M_NR'
    publishDir "progress/annotate_samples", mode: "$params.intermediate_pub_mode"
    tag { "$rt:$is_bc" }

    input:
        tuple val(rt), val(is_bc), val(nr), path(bam), path(barcodes_fa)

    output:
        tuple val(rt), val(is_bc), val(nr), path(out)

    script:
        out = params.run_id + '.' + rt + '.sm_annot.bam'
        """
        samtools view -u $bam |
            bam_annotate_samples.py - \\
                --barcodes $barcodes_fa \\
                --out $out
        """
}

process single_sample {
    label 'M_NR'
    publishDir "progress/annotate_samples", mode: "$params.intermediate_pub_mode"
    tag { "$rt:$is_bc" }

    input:
    tuple val(rt), val(is_bc), val(nr), path(bam), val(sample)

    output:
    tuple val(rt), val(is_bc), val(nr), path(out), emit: bams

    script:
    out = params.run_id + '.' + rt + '.sm_annot.bam'
    """
        samtools view -u $bam |
            bam_annotate_samples.py - \\
                --sample $sample \\
                --out $out
        """
}