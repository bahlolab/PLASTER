
workflow pb_lima {
    take:
        data
    main:
        lima(data)
        lima.out.smry |
            toSortedList() |
            map { [
                (it.find { it[0] == 'CCS' })[2],
                (it.find { it[0] == 'SR' })[2] ] } |
            merge_smry

    emit:
        //  rt, is_bc, nr, bam
        bams = lima.out.bc |
            mix (lima.out.not_bc) |
            map { it.take(2) + [it[2].toFile().text.trim() as int, it[3]] }
        smry = merge_smry.out

}

process lima {
    label 'L'
    publishDir "progress/pb_lima", mode: "$params.intermediate_pub_mode"
    tag { rt }

    input:
        tuple path(bam), val(rt), path(bc_fasta)

    output:
        tuple val(rt), val(true), file('nr_bc'), file(is_bc),  emit: bc
        tuple val(rt), val(false), file('nr_nbc'), file(no_bc), file('nr_nbc'), emit: not_bc
        tuple val(rt), file("${pref}.lima.counts"), file(smry), emit: smry

    script:
        pref = "${params.run_id}.${rt}.lima"
        out = "${pref}.bam"
        smry = "${pref}.lima.summary"
        is_bc = "${pref}.is_bc.bam"
        no_bc = "${pref}.no_bc.bam"
        """
        lima $bam $bc_fasta $out ${rt ==~ /^CCS$/ ? '--ccs' : ''} \\
            --same \\
            --num-threads $task.cpus \\
            --dump-removed
        sed '2q;d' $smry | grep -oP '(?<=: )[0-9]+(?=.+\$)' > nr_bc
        sed '3q;d' $smry | grep -oP '(?<=: )[0-9]+(?=.+\$)' > nr_nbc
        mv $out $is_bc
        mv ${pref}.removed.bam $no_bc
        """
}

process merge_smry {
    label 'XS'
    publishDir "progress/merge_lima_smry", mode: "$params.intermediate_pub_mode"

    input:
        tuple path(lima_smry_ccs), path(lima_smry_sr)

    output:
        path out

    script:
        out = params.run_id + ".merged_lima_smry.tsv"
        """
        merge_lima_smry.R $lima_smry_ccs $lima_smry_sr $out
        """
}