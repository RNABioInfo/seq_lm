#!/usr/bin/env nextflow
nextflow.enable.types = true

include { Sample ; GeneContrast ; QuantifiedSampleBatch } from '../lib/sample.nf'
include { resolve_gene_contrasts ; gene_contrast_groups ; resolve_ica_groups ; select_analysis_batch ; gene_contrast_tsv ; parse_delimited_rows } from '../lib/gene_contrasts.nf'
include { stability_config ; apply_run_stability_actions ; read_stability_tsv ; write_stability_tsv } from '../lib/stability.nf'
include { parse_sample_sheet ; sample_sheet_has_order_column } from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-contrasts-')
    def samples: List<Sample> = ['baseline', 'challenge', 'QC only', 'quoted, group'].collectMany { group: String ->
        [1, 2].collect { replicate: Integer ->
            record(
                name: "${group}_${replicate}",
                group: group,
                order: null,
                bam_dir: root,
                is_live: true,
                protocol_run_id: 'shared',
            )
        }
    }
    def csv: Path = root.resolve('contrasts.csv')
    java.nio.file.Files.writeString(csv, 'challenge,baseline\n"quoted, group", challenge\nbaseline,challenge\n')
    def contrasts: List<GeneContrast> = resolve_gene_contrasts(samples, csv)
    assert contrasts*.target_group == ['challenge', 'quoted, group', 'baseline']
    assert contrasts*.reference_group == ['baseline', 'challenge', 'challenge']
    assert contrasts*.contrast_id == ['contrast_1', 'contrast_2', 'contrast_3']
    assert !gene_contrast_groups(contrasts).contains('QC only')
    assert gene_contrast_tsv(contrasts).contains('contrast_2\tquoted, group\tchallenge')
    def batch: QuantifiedSampleBatch = record(
        batch_index: 7,
        report_sequence: 9,
        samples: samples.collect { sample -> record(batch_index: 3, sample: sample, counts: csv) },
    )
    def selected: QuantifiedSampleBatch = select_analysis_batch(batch, gene_contrast_groups(contrasts))
    assert selected.samples.size() == 6 && selected.batch_index == 7 && selected.report_sequence == 9
    assert resolve_ica_groups(samples, 'baseline', null).size() == 4
    assert resolve_ica_groups(samples, 'baseline', '  ').size() == 4
    assert resolve_ica_groups(samples, 'baseline', 'challenge') == ['baseline', 'challenge'].toSet()
    ['unknown', 'challenge,', 'challenge,challenge'].each { value: String ->
        def rejected: Boolean = false
        try {
            resolve_ica_groups(samples, 'baseline', value)
        }
        catch (_exception: Exception) {
            rejected = true
        }
        assert rejected : value
    }
    ['', 'challenge,baseline\nchallenge,baseline\n', 'baseline,baseline\n', 'challenge,unknown\n', 'challenge\n', 'challenge,baseline,extra\n', 'challenge,\n', 'challenge,baseline\n\n', '"challenge,baseline\n'].each { value: String ->
        java.nio.file.Files.writeString(csv, value)
        def rejected: Boolean = false
        try {
            resolve_gene_contrasts(samples, csv)
        }
        catch (_exception: Exception) {
            rejected = true
        }
        assert rejected : "Unexpected accepted CSV: ${value}"
    }
    def fallback_samples: List<Sample> = samples.collect { sample ->
        record(
            name: sample.name,
            group: sample.group == 'baseline' ? 'control' : sample.group,
            order: sample.order,
            bam_dir: sample.bam_dir,
            is_live: sample.is_live,
            protocol_run_id: sample.protocol_run_id,
        )
    }
    def defaults: List<GeneContrast> = resolve_gene_contrasts(fallback_samples, null)
    assert defaults.size() == 3 && defaults*.reference_group.unique() == ['control']
    assert stability_config([contrasts: contrasts], samples) != stability_config([contrasts: defaults], samples)
    def run_rows: List<Map> = [
        [sample: 'selected', group: 'challenge', bam_dir: root.toString(), protocol_run_id: 'shared', effectively_live: true, eligible: true, newly_eligible: true],
        [sample: 'unselected', group: 'QC only', bam_dir: root.toString(), protocol_run_id: 'shared', effectively_live: true, eligible: false, newly_eligible: false],
    ]
    assert apply_run_stability_actions(run_rows, 'log', [:]) == null
    assert !run_rows[0].run_eligible && run_rows[0].action_result == 'waiting_for_run_stability'
    assert run_rows[0].run_blocking_samples == 'QC only/unselected'
    def escaped_sample: Sample = record(
        name: 'quoted',
        group: 'escaped " group',
        order: null,
        bam_dir: root,
        is_live: false,
        protocol_run_id: null,
    )
    java.nio.file.Files.writeString(csv, '  "escaped "" group"  ,  "baseline"  \n')
    def escaped: List<GeneContrast> = resolve_gene_contrasts([samples[0], escaped_sample], csv)
    def canonical: Path = root.resolve('contrasts.tsv')
    java.nio.file.Files.writeString(canonical, gene_contrast_tsv(escaped))
    assert parse_delimited_rows(canonical.text, '\t', true, false)[0].target_group == escaped_sample.group
    write_stability_tsv(canonical, ['group', 'sample'], [[group: escaped_sample.group, sample: 'quoted']])
    assert read_stability_tsv(canonical)[0].group == escaped_sample.group
    assert parse_delimited_rows('group\tsample\nA\t\n', '\t', true, false)[0].sample == ''
    def sample_sheet: Path = root.resolve('samples.csv')
    java.nio.file.Files.writeString(sample_sheet, "alias,group,bam_dir,order\nquoted,\"escaped \"\" group\",${root},\nbaseline,baseline,${root},\n")
    def parsed: List<Sample> = parse_sample_sheet([sample_sheet_path: sample_sheet])
    assert parsed*.group == ['escaped " group', 'baseline']
    assert sample_sheet_has_order_column(sample_sheet)
    assert resolve_gene_contrasts(parsed, csv)*.target_group == ['escaped " group']
    log.info('Validated explicit contrasts, fallback, independent ICA selection, and batch identity.')
}
