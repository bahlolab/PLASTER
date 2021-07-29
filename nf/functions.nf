import groovy.json.JsonSlurper
import java.nio.file.Path

Path path(String filename) {
    file(filename.replaceAll(/^@PROJECT_DIR@/, workflow.projectDir.toString()),
        checkIfExists: true)
}

ArrayList<Map> readTSV(Object file, List<String> colnames) {
    lines = file instanceof Path ?
        file.toFile().readLines() :
        path(file).toFile().readLines()
    assert colnames == lines[0].split('\t')
    lines.each {assert it.split('\t').size() == colnames.size() }
    lines.drop(1)
        .collect {it.split('\t') }
        .collect {
            [colnames, it].transpose().collectEntries { k, v -> [(k): v] }
        }
}

void checkManiAmps(Path manifest_tsv, Path amplicons_json) {
    // check manifest
    def manifest_lines = manifest_tsv.toFile().readLines().drop(1) as ArrayList
    assert manifest_lines.size() > 0
    manifest_lines.each { assert it.split('\t').size() == 3 }
    def manifest_amp_set = manifest_lines
        .collect { it.split('\t')[2].split(';') }
        .flatten().unique().sort()
    // check amplicons
    def amplicons = (new JsonSlurper().parse(amplicons_json.toFile())) as Map
    assert amplicons.keySet().containsAll(manifest_amp_set)
    amplicons.each { k, v ->
        assert (v as Map).keySet()
            .containsAll(['chrom', 'start', 'end', 'strand', 'fwd_primer', 'rvs_primer'])
    }
}

ArrayList checkManiAmpsAT(ArrayList amplicon_set, Path amplicons_json) {
    // check amplicons
    def amplicons = (new JsonSlurper().parse(amplicons_json.toFile())) as Map
    assert amplicons.keySet().containsAll(amplicon_set)
    amplicons.each { k, v ->
        assert (v as Map).keySet()
            .containsAll(['chrom', 'start', 'end', 'strand'])
    }
    amps = amplicons
        .collect { k, v -> [k, "$v.chrom:$v.start-$v.end"] }
    fusion = amplicons
        .collect { k, v -> [k, v.fusion ] }
        .find { a, f -> a != f & amplicons.containsKey(f) }
        .with { it ? amplicons.subMap(it).collectEntries {
            k, v -> [(k): v.subMap(['chrom', 'start', 'end', 'strand'])]
        } : null }
    pharmvar = amplicons
        .findAll {k, v -> v.pharmvar_gene != null &  v.pharmvar_ver != null }
        .collect { k, v -> [k, v] }

    return [amps, fusion, pharmvar]
}

ArrayList parseManifestPP(String filename) {
    header = ['sample', 'barcode', 'amplicons']
    manifest_path = path(filename)
    manifest = manifest_path.toFile().readLines()
        .with { lines ->
            assert header == lines[0].split('\t')
            lines.each {assert it.split('\t').size() == header.size() }
            lines.drop(1).collect {
                [header, it.split('\t')].transpose().collectEntries { k, v -> [(k): v] } } }
        .collect {
            it.amplicons = it.amplicons.split(';')
            it }
}

ArrayList<Map> parseManifestAT(String filename) {
    readTSV(filename, ['sample', 'amplicon', 'n_reads', 'bam_file'])
        .collect {
            it.n_reads = it.n_reads as Integer
            it.bam_file = path(it.bam_file)
            it }
}

ArrayList parseCopyNum(ArrayList<Map> manifest, String filename, Integer default_cn) {
    sm_am_set = manifest.collect{[it.sample, it.amplicon]}.unique()
    cn = []
    missing = sm_am_set
    if (filename) {
        cn = readTSV(filename, ['sample', 'amplicon', 'copy_num'])
            .collect {
                it.copy_num = it.copy_num as Integer
                it.is_default = false
                it }
            .findAll { sm_am_set.contains([it.sample, it.amplicon]) }
        missing = sm_am_set - cn.collect{[it.sample, it.amplicon]}.unique()
    }
    (cn + missing.collect { s, a ->
        [sample: s, amplicon: a, copy_num: default_cn, is_default: true] })
        .collect { [it.amplicon, it.sample, it.copy_num, it.is_default] }
}