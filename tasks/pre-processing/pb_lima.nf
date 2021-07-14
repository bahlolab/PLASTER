
workflow pb_lima {
    take:
        data
        bc_fasta
    main:
        PBL(data, bc_fasta)
    emit:
        smry = PBL.out.smry
//        // rt, is_bc, nr, bam
        bams = PBL.out.bc.map { it + [true] } |
            mix (PBL.out.not_bc.map { it + [false] }) |
            map { [it[0], it[3], it[2].toFile().text.trim() as int, it[1] ] }
}

process PBL {
    label 'L'
    publishDir "progress/pb_lima", mode: "$params.intermediate_pub_mode"
    tag { rt }

    input:
        tuple val(rt), path(bam)
        path bc_fasta

    output:
        tuple val(rt), file(is_bc), file('nr_bc'), emit: bc
        tuple val(rt), file(no_bc), file('nr_nbc'), emit: not_bc
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