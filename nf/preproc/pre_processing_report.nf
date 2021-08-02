

process pre_processing_report {
    label 'S4'
    publishDir "output", mode: params.output_pub_mode

    input:
        tuple path(rmd), path(amplicons), path(barcodes_fa), path(stats), path(lima_smry)

    output:
        tuple path(html), path(amp_counts_smry), emit: report

    script:
        html = "pre_processing_report.html"
        amp_counts_smry = "amp_counts_smry.tsv"
        """
        cp --remove-destination `readlink $rmd` $rmd
        R -e "run_id='$params.run_id';\\
              barcodes='$barcodes_fa';\\
              amplicons_json='$amplicons';\\
              lima_summary='$lima_smry';\\
              rmarkdown::render('$rmd', output_file='$html')"\\
            --slave --vanilla
        """
}