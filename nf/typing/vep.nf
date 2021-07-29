

process vep {
    label 'M2'
    publishDir "output", mode: params.output_pub_mode
    tag { "$am" }

    input:
        tuple val(am), path(vcf), file(tbi)

    output:
        tuple val(am), path(out), file("${out}.tbi")

    script:
        out = "${am}.vep.vcf.gz"
        """
        vep --input_file $vcf \\
            --database \\
            --format vcf \\
            --vcf \\
            --everything \\
            --allele_number \\
            --variant_class \\
            --dont_skip \\
            --assembly $params.vep_assembly \\
            --cache_version $params.vep_cache_ver \\
            --allow_non_variant \\
            --pick_allele_gene \\
            --output_file STDOUT |
            bcftools view --no-version -Oz -o $out
        bcftools index -t $out
        """
}