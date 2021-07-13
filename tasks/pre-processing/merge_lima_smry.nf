
workflow merge_lima_smry {
    take:
        data
    main:
        data |
            toSortedList() |
            map { [
                (it.find { it[0] == 'CCS' })[2],
                (it.find { it[0] == 'SR' })[2] ] } |
            MLS
    emit:
        smry = MLS.out.smry
}

process MLS {
    label 'XS'
    publishDir "intermediates/merge_lima_smry", mode: "$params.intermediate_pub_mode"

    input:
        tuple path(lima_smry_ccs), path(lima_smry_sr)

    output:
        path out, emit: smry

    script:
        out = params.run_id + ".merged_lima_smry.tsv"
        """
        merge_lima_smry.R $lima_smry_ccs $lima_smry_sr $out
        """
}