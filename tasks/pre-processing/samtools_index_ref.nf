

process samtools_index_ref {
    label 'S'
    publishDir "output/ref", mode: params.output_pub_mode

    input:
        path ref

    output:
        tuple path("${ref}.fai"), path(dict)

    script:
        dict = ref.toString().replaceAll('.fa(sta)?$', '.dict')
        """
        samtools faidx $ref
        samtools dict $ref > $dict
        """
}