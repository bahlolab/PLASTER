
process extract_barcode_set {
    label 'XS'
    publishDir "progress/extract_barcode_set", mode: "$params.intermediate_pub_mode"

    input:
        path manifest
        path barcodes_fa

    output:
        path "bc_order.txt", emit: order
        path "bc_set.fa", emit: fasta

    script:
        """
        samtools faidx $barcodes_fa

        tail -n+2 $manifest | cut -f2 | sort | uniq > bc_order.txt

        cat bc_order.txt | xargs --replace samtools faidx $barcodes_fa {} > bc_set.fa

        BS=\$(grep '>' bc_set.fa | wc -l)
        SB=\$(cat bc_order.txt | wc -l)

        if [[ \${BS} -ne \${SB} ]]; then
            exit 1
        fi
        """
}