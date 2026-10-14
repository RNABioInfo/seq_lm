#!/usr/bin/env nextflow

nextflow.enable.types = true

include { validate_parameters; ica_analysis_enabled } from '../lib/validation.nf'

def assert_validation_error(values: Map, expected_message: String) -> Void {
    try {
        validate_parameters(values)
    }
    catch (exception: Exception) {
        def observed_message: String = exception.message ?: exception.cause?.message ?: "${exception}"
        assert observed_message.contains(expected_message):
            "Expected validation error containing '${expected_message}', observed '${observed_message}'."
        return null
    }
    assert false: "Expected parameter validation to fail with: ${expected_message}"
}

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-parameter-validation-')
    def certificate: Path = root.resolve('client.pem')
    def private_key: Path = root.resolve('key.pem')
    def ca_certificate: Path = root.resolve('ca.crt')
    [certificate, private_key, ca_certificate].each { path: Path ->
        java.nio.file.Files.writeString(path, '')
    }
    def base: Map = [
        strand_specificity: 'forward',
        stopping_method: 'dea',
        ica_analysis: false,
        de_lfc_cutoff: 1.0,
        de_padj_cutoff: 0.05,
        min_read_count: 10000,
        min_replicate_sample_count: 2,
        monitoring_behavior: 'terminate',
        ica_max_fold_change: 1.05,
        ica_min_read_growth_fraction: 0.1,
        differential_expression: true,
        live_analysis: true,
        minknow_host: 'minknow.local',
        minknow_port: 9501,
        minknow_client_certificate: certificate,
        minknow_client_private_key: private_key,
        minknow_ca_certificate: ca_certificate,
        num_stable_batches: 3,
        stability_max_feature_diff_fraction: 0.05,
        stability_min_jaccard_similarity: 0.95,
        stability_max_median_abs_lfc_delta: 0.05,
        stability_min_de_calls_for_fraction_metrics: 20,
        stability_max_small_set_call_changes: 2,
        gene_set_enrichment: false,
        timeline_analysis: false,
        reference_genome: root.resolve('genome.fa'),
        reference_annotation: root.resolve('annotation.gtf'),
        gene_sets: null,
    ]

    validate_parameters(base)
    ['ica_max_fold_change', 'ica_min_read_growth_fraction'].each { name: String ->
        def missing_parameter: Map = new LinkedHashMap(base)
        missing_parameter.remove(name)
        assert_validation_error(missing_parameter, "--${name} must be a finite number.")
    }
    validate_parameters(base + [stopping_method: 'dea', ica_max_fold_change: 1, ica_min_read_growth_fraction: 0.1])
    ['none', 'ICA', '', null, true, 1].each { method ->
        assert_validation_error(base + [stopping_method: method], '--stopping_method must be')
    }
    [-1, 0, 0.95, Double.NaN, Double.POSITIVE_INFINITY, null, '1.05', true].each { value ->
        assert_validation_error(base + [ica_max_fold_change: value], '--ica_max_fold_change must be')
    }
    [-1, 0, Double.NaN, Double.POSITIVE_INFINITY, null, '0.1', true].each { value ->
        assert_validation_error(base + [ica_min_read_growth_fraction: value], '--ica_min_read_growth_fraction must be')
    }
    [0, -1, 1.5, Double.NaN, Double.POSITIVE_INFINITY, null, '3', true].each { value ->
        assert_validation_error(base + [num_stable_batches: value], '--num_stable_batches must be a positive integer')
    }
    assert_validation_error(base + [stopping_method: 'ica'], 'requires enabled --ica_analysis and --ica_matrix')
    assert_validation_error(base + [differential_expression: false], '--stopping_method dea requires --differential_expression')
    ['forward', 'reverse', 'unstranded'].each { strand: String ->
        validate_parameters(base + [strand_specificity: strand])
    }
    ['fw', 'both', 'auto', 'Forward', '', '  ', null, true, 1, ['forward']].each { strand ->
        assert_validation_error(base + [strand_specificity: strand], '--strand_specificity must be')
    }
    validate_parameters(base + [monitoring_behavior: 'disabled', minknow_client_certificate: null, minknow_client_private_key: null, minknow_ca_certificate: null])
    validate_parameters(base + [monitoring_behavior: 'log', minknow_client_certificate: null, minknow_client_private_key: null, minknow_ca_certificate: null])
    validate_parameters(base + [differential_expression: false, monitoring_behavior: 'disabled'])
    validate_parameters(base + [differential_expression: false, monitoring_behavior: 'disabled', reference_genome: null, reference_annotation: null])

    assert_validation_error(base + [minknow_host: '  '], 'non-empty --minknow_host')
    assert_validation_error(base + [minknow_port: 0], '--minknow_port must be between')
    assert_validation_error(base + [minknow_port: 65536], '--minknow_port must be between')
    assert_validation_error(base + [minknow_client_certificate: null], 'requires --minknow_client_certificate')
    assert_validation_error(base + [minknow_client_private_key: null], 'requires --minknow_client_private_key')
    assert_validation_error(base + [minknow_ca_certificate: null], 'requires --minknow_ca_certificate')
    assert_validation_error(
        base + [minknow_client_certificate: root.resolve('missing.pem')],
        '--minknow_client_certificate must point to an existing regular file',
    )
    assert_validation_error(
        base + [differential_expression: false, monitoring_behavior: 'disabled', reference_annotation: null],
        '--reference_genome and --reference_annotation must be supplied together',
    )
    assert_validation_error(
        base + [differential_expression: false, monitoring_behavior: 'disabled', reference_genome: null],
        '--reference_genome and --reference_annotation must be supplied together',
    )
    assert_validation_error(
        base + [timeline_analysis: true],
        '--timeline_analysis requires --gene_set_enrichment',
    )

    def ica: Map = base + [ica_control_group: 'control', ica_analysis: true, monitoring_behavior: 'disabled', ica_matrix: root.resolve('weights.tsv'),
        ica_log_base: 2.0, ica_pseudocount: 1.0, ica_min_gene_coverage: 1.0,
        ica_min_read_count: 10000, ica_padj_cutoff: 0.05]
    assert ica_analysis_enabled(ica)
    assert ica_analysis_enabled(ica + [ica_analysis: true])
    assert !ica_analysis_enabled(ica + [ica_analysis: false])
    assert !ica_analysis_enabled(ica + [ica_matrix: null])
    assert !ica_analysis_enabled(ica + [ica_matrix: '  '])
    validate_parameters(ica + [ica_analysis: false, differential_expression: false,
        reference_genome: null, reference_annotation: null])
    validate_parameters(ica + [ica_analysis: false, ica_matrix: null,
        ica_gene_map: root.resolve('missing map.tsv')])
    validate_parameters(ica)
    validate_parameters(ica + [differential_expression: false])
    validate_parameters(ica + [stopping_method: 'ica', monitoring_behavior: 'log', differential_expression: false])
    validate_parameters(ica + [stopping_method: 'ica', monitoring_behavior: 'terminate', differential_expression: false])
    assert_validation_error(ica + [stopping_method: 'ica', monitoring_behavior: 'log', ica_analysis: false], 'requires enabled --ica_analysis and --ica_matrix')
    assert_validation_error(ica + [ica_matrix: null, ica_gene_map: root.resolve('map.tsv')], '--ica_gene_map requires --ica_matrix')
    assert_validation_error(ica + [reference_genome: null], '--ica_matrix requires both')
    assert_validation_error(ica + [ica_log_base: 1], '--ica_log_base must be greater')
    assert_validation_error(ica + [ica_log_base: Double.NaN], '--ica_log_base must be a finite')
    assert_validation_error(ica + [ica_pseudocount: 0], '--ica_pseudocount must be greater')
    assert_validation_error(ica + [ica_min_gene_coverage: 0], '--ica_min_gene_coverage must be greater')
    assert_validation_error(ica + [ica_min_gene_coverage: 1.1], '--ica_min_gene_coverage must be greater')
    assert_validation_error(ica + [ica_min_read_count: -1], '--ica_min_read_count must be a nonnegative integer')
    assert_validation_error(ica + [ica_min_read_count: 0.5], '--ica_min_read_count must be a nonnegative integer')
    assert_validation_error(ica + [ica_padj_cutoff: 1.1], '--ica_padj_cutoff must be greater')
    assert_validation_error(ica + [ica_padj_cutoff: Double.POSITIVE_INFINITY], '--ica_padj_cutoff must be a finite')

    validate_parameters(base + [detect_minknow_completion: true, minknow_completion_settle_seconds: 30])
    validate_parameters(base + [detect_minknow_completion: false, minknow_completion_settle_seconds: 1])
    [0, -1, 1.5, null, '30', true].each { value ->
        assert_validation_error(base + [minknow_completion_settle_seconds: value], '--minknow_completion_settle_seconds must be a positive integer')
    }
    assert_validation_error(base + [detect_minknow_completion: 'true'], '--detect_minknow_completion must be a Boolean')
    channel.of('MinKNOW termination parameters are conditionally validated').view()
}
