
include { path } from '../functions'

workflow pre_processing_report {
    take:
        data // stats, lima_smry, amplicons
    main:
        rmd = path("${workflow.projectDir}/bin/preproc-report.Rmd")
        if (params.barcodes_fasta) {
            barcodes_fasta = Channel.fromList([path(params.barcodes_fasta)])
        } else {
            sample = params.sample ?: "sample"
            barcodes_fasta = Channel.from(['N']) |
                collectFile(seed: ">$sample", newLine:true, name: 'barcodes.fa')
        }
        out = data |
            map {[rmd] + it } |
            combine(barcodes_fasta) |
            report

    emit:
        out
}


process report {
    label 'S4'
    publishDir "output", mode: params.output_pub_mode

    input:
        tuple path(rmd), path(stats), path(lima_smry), path(amplicons), path(barcodes_fa)

    output:
        tuple path(html), path(amp_counts_smry)

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