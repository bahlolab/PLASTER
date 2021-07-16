import groovy.json.JsonSlurper
import java.nio.file.Path

Path path(String filename) {
    file(filename.replaceAll(/^@PROJECT_DIR@/, workflow.projectDir.toString()),
        checkIfExists: true)
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

void checkManiAmpsAT(ArrayList amplicon_set, Path amplicons_json) {
    // check amplicons
    def amplicons = (new JsonSlurper().parse(amplicons_json.toFile())) as Map
    assert amplicons.keySet().containsAll(amplicon_set)
    amplicons.each { k, v ->
        assert (v as Map).keySet()
            .containsAll(['chrom', 'start', 'end', 'strand', 'fwd_primer', 'rvs_primer'])
    }
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

ArrayList parseManifestAT(String filename) {
    header = ['run_id', 'sample', 'amplicon', 'n_reads', 'bam_file']
    manifest_path = path(filename)
    manifest = manifest_path.toFile().readLines()
        .with { lines ->
            assert header == lines[0].split('\t')
            lines.each {assert it.split('\t').size() == header.size() }
            lines.drop(1).collect {
                [header, it.split('\t')].transpose().collectEntries { k, v -> [(k): v] } } }
        .collect {
            it.n_reads = it.n_reads as Integer
            it.bam_file = path(it.bam_file)
            it }
}


List refFastaFileMap(String filename) {
    def suffs = ['dict', 'fai', 'ccs.mmi', 'subread.mmi']
    def inclusive=true
    def fn = new File(filename).toPath().toAbsolutePath()
    def basename = fn.fileName.toString()
    def pattern =  /.*\.fa(sta)?(\.gz)?$/
    if (! basename ==~  pattern ){
        errorExit("Error: reference fasta file name must match pattern \"$pattern\".")
    }
    def base = basename.replaceAll(/\.fa(sta)?(\.gz)?$/, '') + '.*\\.'
    def patterns = suffs.collectEntries { [(it): base + it + '$' ] }
    def matches = patterns.collectEntries {
        [ (it.key) : (new FileNameByRegexFinder().getFileNames(fn.parent.toString(), it.value))[0] ]
    }
    def failed = matches.findAll{ it.value == null }.keySet()
    if (failed) {
        errorExit("Error: File(s) with suffix \"${failed.join('", "')}\" not found for file \"$filename\"\n")
    }
    matches = matches.collectEntries { [ (it.key) : (new File(it.value).toPath() ) ] }
    if (inclusive){
        matches.put('fa', fn)
        def gzi = (new FileNameByRegexFinder().getFileNames(fn.parent.toString(), base + 'gzi' + '$'))[0]
        if (gzi) { matches.put('gzi', (new File(gzi).toPath()) )}
    }
    return [ matches.collectEntries { [(it.key) : it.value.fileName.toString()] },
             matches.collect { it.value } ]
}
