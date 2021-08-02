import groovy.json.JsonOutput

include { path } from '../functions'

workflow prep_bams {
    take:
    data // sm, am, nr, bam
    ref // fasta, fai, dict
    fusion // amplicons map
    main:
    if (fusion) println "checking for ${fusion.keySet().join('/')} fusion/chimeric reads"
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
            too_many: it[2] > params.max_reads
            too_few: it[2] < params.min_reads
            pass_through: true }

    if (fusion) {
        (bams2.fusion |
            groupTuple(by:0) |
            map { it[0..1] + [it[2].sum(), it[3]] })
            .with { fusion_call(it, ref, fusion) }

        bams3 = fusion_call.out.bams |
            flatMap { sm, bams, count -> [bams, count].transpose().collect { [sm] + it } } |
            filter { it[1] ==~ /.+\.clean-[^\/]+\.bam$/ } |
            map { sm, bam, count ->
                [sm, (bam =~ /.+\.clean-([^\/]+).bam$/)[0][1], count.toFile().text as int, bam]} |
            branch {
                too_many: it[2] > params.max_reads
                too_few: it[2] < params.min_reads
                pass_through: true
            }

        fusion_call.out.results |
            map { it[2..3] } |
            toSortedList |
            map { it.transpose() } |
            combine([[fusion_set]]) |
            combine([path("${workflow.projectDir}/bin/fusion_report.Rmd")]) |
            fusion_report

    } else {
        bams3 = [too_many: Channel.fromList([]),
                 too_few: Channel.fromList([]),
                 pass_through: Channel.fromList([])]
    }

    // write samples with too few reads to file
    bams2.too_few |
        mix(bams3.too_few) |
        map { it.take(3).collect {it.toString()}.join('\t') } |
        collectFile(name: 'low_read_count.tsv', storeDir: './output/', newLine: true,
            seed: ['sample', 'amplicon', 'n_reads'].join('\t')) |
        map {
            lines = it.toFile().readLines()
            if (lines.size() > 1) {
                println "Note: ${lines.size() - 1} sample-amplicons excluded due to low read count, written to $it"
            }
        }

    bams2.too_many |
        mix(bams3.too_many) |
        downsample |
        mix(bams2.pass_through) |
        mix(bams3.pass_through) |
        index

    emit:
        index.out // am, sm, nr, bam, bai
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
    tuple val(am), val(sm), val(nr), path(bam), path("${bam}.bai")

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
    tuple val(sm), val(am), path("${prefix}.fus_smry.csv"), path("${prefix}.breakpoints.csv.gz"), emit: results
    tuple val(sm), path("${prefix}.*.bam"), path("${prefix}.*.bam.count"), emit: bams

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

process fusion_report {
    label 'M'
    publishDir "output", mode: "$params.output_pub_mode"

    input:
        tuple path(smry), path(bps), val(am), path(rmd)

    output:
        tuple path(html), path(calls)

    script:
    html = "fusion_report.html"
    calls = 'fusion_calls.csv'
    """
    cp --remove-destination `readlink $rmd` $rmd
    R -e "amplicons=\'${am.join(',')}\';\\
          rmarkdown::render('$rmd', output_file='$html')" --slave --vanilla
    """
}