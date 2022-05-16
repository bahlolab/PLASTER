import groovy.json.JsonOutput

process pharmvar_star_allele {
    label 'S2'
    publishDir "output", mode: params.output_pub_mode
    tag { am }

    input:
        tuple val(am), path(sm_vcf), path(pv_vcf), val(pharmvar_meta)


    output:
        tuple val(am), path("${am}.allele_definition.csv"), path("${am}.sample_phase_alleles.csv")

    script:
        json = JsonOutput.toJson(pharmvar_meta)
        """
        pharmvar_star_allele.R $sm_vcf $pv_vcf \\
            --pharmvar '$json' \\
            --out-pref $am
        """
}