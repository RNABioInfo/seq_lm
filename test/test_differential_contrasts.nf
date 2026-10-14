#!/usr/bin/env nextflow
nextflow.enable.types = true

include { Sample ; QuantifiedSampleBatch ; GeneContrast } from '../lib/sample.nf'
include { resolve_gene_contrasts ; gene_contrast_groups ; select_analysis_batch } from '../lib/gene_contrasts.nf'
include { differential_expression } from '../subworkflows/differential_expression.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-explicit-de-')
    def annotation: Path = root.resolve('annotation.gtf')
    java.nio.file.Files.writeString(
        annotation,
        (1..30).collect { index: Integer ->
            "chr1\ttest\ttranscript\t${index * 100}\t${index * 100 + 89}\t.\t+\t.\tgene_id \"g${index}\"; transcript_id \"tx${index}\";\n"
        }.join(''),
    )
    def gene_sets: Path = root.resolve('sets.gmt')
    java.nio.file.Files.writeString(gene_sets, 'UP\tinduced\ttx1\ttx2\ttx3\ttx4\ttx5\nDOWN\trepressed\ttx6\ttx7\ttx8\ttx9\ttx10\n')
    def quantified = ['baseline', 'challenge', 'QC only'].collectMany { group: String ->
        (1..3).collect { replicate: Integer ->
            def sample: Sample = record(
                name: "${group}_${replicate}",
                group: group,
                order: null,
                bam_dir: root,
                is_live: false,
                protocol_run_id: null,
            )
            def counts: Path = root.resolve("${group}_${replicate}.quant")
            java.nio.file.Files.writeString(
                counts,
                'tname\tlen\tnum_reads\n' + (1..30).collect { index: Integer ->
                    def abundance: Integer = group == 'QC only' ? 0 : group == 'challenge' && index <= 5 ? 400 + replicate * 7 : group == 'challenge' && index <= 10 ? 20 + replicate : 100 + ((index + replicate) % 7)
                    "tx${index}\t90\t${abundance}\n"
                }.join(''),
            )
            record(batch_index: 0, sample: sample, counts: counts)
        }
    }
    def csv: Path = root.resolve('contrasts.csv')
    java.nio.file.Files.writeString(csv, 'challenge,baseline\nbaseline,challenge\n')
    def contrasts: List<GeneContrast> = resolve_gene_contrasts(quantified*.sample, csv)
    def batch: QuantifiedSampleBatch = record(batch_index: 0, report_sequence: 0, samples: quantified)
    def selected: QuantifiedSampleBatch = select_analysis_batch(batch, gene_contrast_groups(contrasts))
    assert selected.samples.size() == 6
    outputs = differential_expression(channel.of(selected), contrasts, 0, annotation, gene_sets, true, 1.0, 100, 2)
    outputs.results.view { result ->
        assert result.batch_index == 0 && result.report_sequence == 0
        def definition: String = result.results.resolve('contrasts.tsv').text
        assert definition.contains('contrast_1\tchallenge\tbaseline')
        assert definition.contains('contrast_2\tbaseline\tchallenge')
        assert !result.results.resolve('sample_metadata.tsv').text.contains('QC only')
        [1, 2].each { index: Integer ->
            assert result.results.resolve("contrast_${index}/edgeR_results.tsv").exists()
            assert result.results.resolve("contrast_${index}/fry_results.tsv").exists()
            assert result.results.resolve("contrast_${index}/gsva_limma_results.tsv").exists()
        }
        def forward: Map = result.results
            .resolve('contrast_1/edgeR_results.tsv')
            .splitCsv(header: true, sep: '\t')
            .collectEntries { row -> [(row.feature_id): row.logFC.toDouble()] }
        def reverse: Map = result.results
            .resolve('contrast_2/edgeR_results.tsv')
            .splitCsv(header: true, sep: '\t')
            .collectEntries { row -> [(row.feature_id): row.logFC.toDouble()] }
        assert forward.tx1 > 0 && forward.tx6 < 0
        forward.each { feature: String, value: Double ->
            assert Math.abs(value + (reverse[feature] as Double)) < 1e-8
        }
        def scores: Map = result.results
            .resolve('contrast_1/gsva_limma_results.tsv')
            .splitCsv(header: true, sep: '\t')
            .collectEntries { row -> [(row.gene_set): row.effect_size.toDouble()] }
        assert scores.UP > 0 && scores.DOWN < 0
        return "Explicit DEA/fry/GSVA contrasts passed without control; low-depth QC-only group excluded: ${result.results}"
    }
    outputs.report_batches.view { report ->
        assert report.has_differential_results
        'Selected cohort passed DEA readiness'
    }
}
