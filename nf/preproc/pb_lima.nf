
include { path } from '../functions'

workflow pb_lima {
    take:
        data //rt, nr, bam
    main:
        if (params.barcodes_fasta) {
            barcodes_fasta = path(params.barcodes_fasta)
            data | combine([barcodes_fasta]) | lima
            lima.out.smry |
                toSortedList() |
                map {
                    [(it.find { it[0] == 'CCS' })[2],
                     (it.find { it[0] == 'SR' })[2]]
                } |
                merge_smry
            bams = lima.out.bc |
                mix(lima.out.not_bc) |
                map { [it[0], it[1], it[2].toFile().text.trim() as int, it[3]] }
                // rt, is_bc, nr, bam
            smry = merge_smry.out
        } else {
            bams = data.map { [it[0], true, it[1], it[2]] } // rt, is_bc, nr, bam
            smry = Channel.from(['not done.']) |
                collectFile(seed: 'note', newLine:true, name: 'lima_smry.tsv')
        }

    emit:
        bams = bams //  rt, is_bc, nr, bam
        smry = smry
}

process lima {
    label 'L_NR'
    publishDir "progress/pb_lima", mode: "$params.intermediate_pub_mode"
    tag { rt }

    input:
        tuple val(rt), val(nr), path(bam), path(bc_fasta)

    output:
        tuple val(rt), val(true), file("${is_bc}.nr"), file(is_bc), emit: bc
        tuple val(rt), val(false), file("${no_bc}.nr"), file(no_bc), emit: not_bc
        tuple val(rt), file("${pref}.lima.counts"), file(smry), emit: smry

    script:
        pref = "${params.run_id}.${rt}.lima"
        is_bc = "${pref}.bam"
        no_bc = "${pref}.removed.bam"
        smry = "${pref}.lima.summary"
        """
        lima $bam $bc_fasta $is_bc ${rt ==~ /^CCS$/ ? '--ccs' : ''} \\
            --same \\
            --num-threads $task.cpus \\
            --dump-removed
        sed '2q;d' $smry | grep -oP '(?<=: )[0-9]+(?=.+\$)' > ${is_bc}.nr
        sed '3q;d' $smry | grep -oP '(?<=: )[0-9]+(?=.+\$)' > ${no_bc}.nr
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