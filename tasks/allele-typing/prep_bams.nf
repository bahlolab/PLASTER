import groovy.json.JsonOutput

workflow prep_bams {
    take:
        data // sm, am, nr, bam
        ref // fasta, fai, dict
        fusion // amplicons map
    main:
        if (fusion) println "checking for fusion: ${fusion.keySet().join('/')}"
        fusion_set = fusion ? fusion.keySet() as ArrayList : []

        bams1 = data |
            groupTuple(by:[0,1]) |
            map { it[2] = it[2].sum(); it } |
            branch {
                multi: it[3].size() > 1
                single: true }

        bams2 = bams1.multi |
            merge |
            mix(bams1.single.map { it[3] = it[3][0]; it }) |
            branch {
                fusion: fusion_set.contains(it[1])
                downsample: it[2] > params.max_reads
                pass_through: true }

//        bams3 =
        (bams2.fusion |
            filter {it[2] > 30} |
            groupTuple(by:0) |
            map { it[0..1] + [it[2].sum(), it[3]] })
            .with { fusion_call(it, ref, fusion) }

        fusion_call.out.bams |
            flatMap { sm, bams -> bams.collect { [sm, it] } } |
//            flatten |
            filter { it[1] ==~ /.+\.clean-[^\/]+\.bam$/ } |
            map { sm, bam ->
                [sm, (bam =~ /.+\.clean-([^\/]+).bam$/)[0][1],  bam]} |
            view
//                too_few: it[2] < params.min_reads

//        bams4 = bams2.downsample |
//            downsample |
//            mix(bams2.pass_through) |
//            index |
//            map { it[[1,0,2,3,4]] }
//
//        // write samples with too few reads to file
//        bams1.too_few |
//            map { it.take(3).collect {it.toString()}.join('\t') } |
//            collectFile(name: 'low_read_count.tsv', storeDir: './output/', newLine: true,
//                seed: ['sample', 'amplicon', 'n_reads'].join('\t')) |
//            map {
//                lines = it.toFile().readLines()
//                if (lines.size() > 1) {
//                    println "Warning: ${lines.size() - 1} sample-amplicons excluded due to low read cout, written to $it"
//                }
//            }
    emit:
        null
        // am, sm, nr, bam, bai
//        bams4
}

process merge {
    label 'S_NR'
    publishDir "progress/merge_bams", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am" }

    input:
        tuple val(sm), val(am), val(nr), path(bams)

    output:
        tuple val(sm), val(am), val(nr), path(merged)

    script:
        merged = "SM-${sm}.AM-${am}.merged.bam"
        """
        for BAM in *.bam; do samtools index \$BAM; done
        samtools merge $merged *.bam
        """
}

process downsample {
    label 'S_NR'
    publishDir "progress/downsample_bams", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am" }

    input:
        tuple val(sm), val(am), val(nr), path(bam)

    output:
        tuple val(sm), val(am), val(params.max_reads), path(out_bam)

    script:
        out_bam = bam.name.replaceAll('.bam', '.sub.bam')
        """
        samtools index $bam
        downsample_bam.py $bam -u --count $params.max_reads | samtools view -b -o $out_bam
        """
}

process index {
    label 'S_NR'
    publishDir "progress/index_bams", mode: "$params.intermediate_pub_mode"
    tag { "$sm:$am" }

    input:
        tuple val(sm), val(am), val(nr), path(bam)

    output:
        tuple val(sm), val(am), val(nr), path(bam), path("${bam}.bai")

    script:
        """
        samtools index $bam
        """
}

process fusion_call {
    label 'M_NR'
    publishDir "output/fusion_call", mode: "$params.output_pub_mode"
    tag { "$sm" }

    input:
        tuple val(sm), val(am), val(nr), path(bam)
        tuple path(ref), path(fai), path(dict)
        val(fusion)

    output:
        tuple val(sm), val(am), file("${prefix}.fus_smry.csv"), file("${prefix}.breakpoints.csv.gz"), emit: results
        tuple val(sm), file("${prefix}.*.bam"), emit: bams

    script:
    json = JsonOutput.toJson(fusion)
    if (am.size() == 2) {
        order = am[0] == fusion.keySet()[0] ? [0,1] : [1,0]
        args = "--primary ${bam[order][0]} --secondary ${bam[order][1]}"
    } else {
        args = am[0] == fusion.keySet()[0] ?
            "--primary ${bam[0]}" :
            "--secondary ${bam[0]}"
    }
    prefix = "SM-${sm}"
    """
    fusion_call.R \\
        $args \\
        --amplicons '$json' \\
        --seq $ref \\
        --max-reads ${params.max_reads * 2} \\
        --min-fusion-reads 5 \\
        --min-total-prop 0.01 \\
        --min-chim-prop 0.25 \\
        --min-score-delta 10 \\
        --fusion-window 10 \\
        --out $prefix
    """
}