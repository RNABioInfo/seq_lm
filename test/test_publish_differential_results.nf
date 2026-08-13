#!/usr/bin/env nextflow

nextflow.enable.types = true

include { publish_differential_results } from '../modules/generic_helpers.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-de-publish-')
    def results: Path = root.resolve('batch_7')
    def condition: Path = results.resolve('group_treated_vs_control')
    java.nio.file.Files.createDirectories(condition)
    java.nio.file.Files.writeString(condition.resolve('edgeR_results.tsv'), 'result\n')

    published = publish_differential_results(channel.of(tuple(0, 7, results)))
    published
        .map { result ->
            assert result.batch_index == 0
            assert result.analysis_index == 7
            assert java.nio.file.Files.isRegularFile(
                result.results.resolve('group_treated_vs_control/edgeR_results.tsv')
            )
            'differential publication process passed'
        }
        .view()
}
