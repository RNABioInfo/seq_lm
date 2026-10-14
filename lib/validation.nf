nextflow.enable.types = true

def valid_strand_specificity(value) -> Boolean {
    return value instanceof String && value in ['forward', 'reverse', 'unstranded']
}

def ica_analysis_enabled(params_map: Map) -> Boolean {
    return params_map.ica_analysis == true && params_map.ica_matrix != null &&
        "${params_map.ica_matrix}".trim() != ''
}

def validate_parameters(params_map: Map) -> Void {
    if (!valid_strand_specificity(params_map.strand_specificity)) {
        error('--strand_specificity must be forward, reverse, or unstranded.')
    }

    if (params_map.ica_analysis == true && params_map.ica_gene_map && !params_map.ica_matrix) {
        error('--ica_gene_map requires --ica_matrix.')
    }

    if (params_map.ica_analysis == true && params_map.ica_imodulon_table && !params_map.ica_matrix) {
        error('--ica_imodulon_table requires --ica_matrix.')
    }

    if (ica_analysis_enabled(params_map)) {
        if (!(params_map.ica_control_group instanceof String) || !params_map.ica_control_group.trim()) {
            error('--ica_control_group must be a nonempty group name.')
        }
        if (!params_map.reference_genome || !params_map.reference_annotation) {
            error('--ica_matrix requires both --reference_genome and --reference_annotation.')
        }
        ['ica_log_base', 'ica_pseudocount', 'ica_min_gene_coverage',
         'ica_min_read_count', 'ica_padj_cutoff'].each { name: String ->
            def value = params_map[name]
            if (!(value instanceof Number) || !Double.isFinite((value as Number).doubleValue())) {
                error("--${name} must be a finite number.")
            }
            if (name == 'ica_log_base' && value <= 1) {
                error('--ica_log_base must be greater than 1.')
            }
            if (name == 'ica_pseudocount' && value <= 0) {
                error('--ica_pseudocount must be greater than 0.')
            }
            if (name in ['ica_min_gene_coverage', 'ica_padj_cutoff'] && (value <= 0 || value > 1)) {
                error("--${name} must be greater than 0 and at most 1.")
            }
            if (name == 'ica_min_read_count' && (value < 0 || value != Math.floor((value as Number).doubleValue()))) {
                error('--ica_min_read_count must be a nonnegative integer.')
            }
        }
    }

    if (params_map.containsKey('detect_minknow_completion') && !(params_map.detect_minknow_completion instanceof Boolean)) {
        error('--detect_minknow_completion must be a Boolean.')
    }
    if (params_map.containsKey('minknow_completion_settle_seconds')) {
        def settle = params_map.minknow_completion_settle_seconds
        if (!(settle instanceof Integer) || settle < 1) {
            error('--minknow_completion_settle_seconds must be a positive integer.')
        }
    }

    if (params_map.de_lfc_cutoff < 0) {
        error('--de_lfc_cutoff must be nonnegative.')
    }

    if (params_map.de_padj_cutoff <= 0 || params_map.de_padj_cutoff > 1) {
        error('--de_padj_cutoff must be greater than 0 and at most 1.')
    }

    if (params_map.min_read_count < 0) {
        error('--min_read_count must be nonnegative.')
    }

    if (params_map.min_replicate_sample_count < 1) {
        error('--min_replicate_sample_count must be at least 1.')
    }

    def stability_behaviors: Set<String> = ['disabled', 'log', 'terminate'].toSet()
    if (!stability_behaviors.contains(params_map.monitoring_behavior as String)) {
        error('--monitoring_behavior must be disabled, log, or terminate.')
    }

    def stopping_method = params_map.stopping_method
    if (!(stopping_method instanceof String) || !(stopping_method in ['dea', 'ica'])) {
        error('--stopping_method must be dea or ica.')
    }

    ['ica_max_fold_change', 'ica_min_read_growth_fraction'].each { name: String ->
        def value = params_map[name]
        if (!(value instanceof Number) || !Double.isFinite((value as Number).doubleValue())) {
            error("--${name} must be a finite number.")
        }
        if (name == 'ica_max_fold_change' && value < 1) {
            error('--ica_max_fold_change must be at least 1.')
        }
        if (name == 'ica_min_read_growth_fraction' && value <= 0) {
            error('--ica_min_read_growth_fraction must be greater than 0.')
        }
    }

    if (params_map.monitoring_behavior != 'disabled') {
        if (stopping_method == 'dea' && !params_map.differential_expression) {
            error('--monitoring_behavior with --stopping_method dea requires --differential_expression.')
        }
        if (stopping_method == 'ica' && !ica_analysis_enabled(params_map)) {
            error('--monitoring_behavior with --stopping_method ica requires enabled --ica_analysis and --ica_matrix.')
        }
    }

    if (params_map.monitoring_behavior == 'terminate' && !params_map.live_analysis) {
        error('--monitoring_behavior terminate requires --live_analysis.')
    }
    if (params_map.monitoring_behavior == 'terminate') {
        if (!params_map.minknow_host || !(params_map.minknow_host as String).trim()) {
            error('--monitoring_behavior terminate requires a non-empty --minknow_host.')
        }
        def minknow_port: Integer = params_map.minknow_port as Integer
        if (minknow_port < 1 || minknow_port > 65535) {
            error('--minknow_port must be between 1 and 65535.')
        }
        def credential_parameters: List<String> = [
            'minknow_client_certificate',
            'minknow_client_private_key',
            'minknow_ca_certificate',
        ]
        credential_parameters.each { parameter_name: String ->
            def value: Object = params_map[parameter_name]
            if (!value || !"${value}".trim()) {
                error("--monitoring_behavior terminate requires --${parameter_name}.")
            }
            def credential_path: Path = file(value)
            if (!credential_path.exists() || !credential_path.isFile()) {
                error("--${parameter_name} must point to an existing regular file: ${credential_path}")
            }
        }
    }
    if (!(params_map.num_stable_batches instanceof Number) ||
        !Double.isFinite((params_map.num_stable_batches as Number).doubleValue()) ||
        params_map.num_stable_batches < 1 || params_map.num_stable_batches != Math.floor((params_map.num_stable_batches as Number).doubleValue())) {
        error('--num_stable_batches must be a positive integer.')
    }

    def fraction_stability_params: Map<String,Float> = [
        stability_max_feature_diff_fraction: params_map.stability_max_feature_diff_fraction,
        stability_min_jaccard_similarity: params_map.stability_min_jaccard_similarity,
    ]
    fraction_stability_params.each { name: String, value: Float ->
        if (value < 0 || value > 1) {
            error("--${name} must be between 0 and 1.")
        }
    }
    if (params_map.stability_max_median_abs_lfc_delta < 0) {
        error('--stability_max_median_abs_lfc_delta must be nonnegative.')
    }
    if (params_map.stability_min_de_calls_for_fraction_metrics < 1) {
        error('--stability_min_de_calls_for_fraction_metrics must be at least 1.')
    }
    if (params_map.stability_max_small_set_call_changes < 0) {
        error('--stability_max_small_set_call_changes must be nonnegative.')
    }

    if (params_map.gene_set_enrichment && !params_map.differential_expression) {
        error('--gene_set_enrichment requires --differential_expression.')
    }
    if (params_map.timeline_analysis && !params_map.gene_set_enrichment) {
        error('--timeline_analysis requires --gene_set_enrichment.')
    }
    def reference_genome_provided: Boolean = params_map.reference_genome != null && "${params_map.reference_genome}".trim() != ''
    def reference_annotation_provided: Boolean = params_map.reference_annotation != null && "${params_map.reference_annotation}".trim() != ''
    if (reference_genome_provided != reference_annotation_provided) {
        error('--reference_genome and --reference_annotation must be supplied together.')
    }
    if (params_map.differential_expression && !params_map.reference_genome) {
        error('Differential expression requires --reference_genome.')
    }
    if (params_map.differential_expression && !params_map.reference_annotation) {
        error('Differential expression requires --reference_annotation.')
    }
    if (params_map.gene_set_enrichment && !params_map.gene_sets) {
        error('Gene-set enrichment requires --gene_sets.')
    }
}
