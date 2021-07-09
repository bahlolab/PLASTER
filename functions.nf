import groovy.json.JsonSlurper

import java.nio.file.Path

Path path(String filename) {

    file(filename.replaceAll(/^@TESTDATA@/, "${workflow.projectDir}/test/data"),
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
