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
        ['ica_log_base', 'ica_pseudocount', 'ica_min_gene_coverage', 'ica_padj_cutoff'].each { name: String ->
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

    if (!(params_map.min_read_count instanceof Integer) || params_map.min_read_count < 0) {
        error('--min_read_count must be a nonnegative integer.')
    }

    if (!(params_map.min_replicate_sample_count instanceof Integer) || params_map.min_replicate_sample_count < 1) {
        error('--min_replicate_sample_count must be a positive integer.')
    }

    def stability_behaviors: Set<String> = ['disabled', 'log', 'terminate'].toSet()
    if (!stability_behaviors.contains(params_map.monitoring_behavior as String)) {
        error('--monitoring_behavior must be disabled, log, or terminate.')
    }

    def stopping_method = params_map.stopping_method
    if (!(stopping_method instanceof String) || !(stopping_method in ['dea', 'ica'])) {
        error('--stopping_method must be dea or ica.')
    }

    if (params_map.monitoring_behavior != 'disabled') {
        if (stopping_method == 'dea' && !params_map.differential_expression) {
            error('--monitoring_behavior with --stopping_method dea requires --differential_expression.')
        }
        if (stopping_method == 'ica' && !ica_analysis_enabled(params_map)) {
            error('--monitoring_behavior with --stopping_method ica requires enabled --ica_analysis and --ica_matrix.')
        }
        validate_precision_parameters(params_map)
        if (stopping_method == 'dea') {
            validate_stopping_file(params_map.dea_stopping_genes, 'dea_stopping_genes', 'GMT')
        }
        else {
            validate_stopping_file(params_map.ica_stopping_components, 'ica_stopping_components', 'component list')
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
    if (params_map.gene_set_enrichment && !params_map.differential_expression) {
        error('--gene_set_enrichment requires --differential_expression.')
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

def validate_stopping_file(value: Object?, name: String, format: String) -> Void {
    if (value == null) {
        return
    }
    if (!(value instanceof String || value instanceof Path) || !"${value}".trim()) {
        error("--${name} must point to an existing regular ${format} file.")
    }
    def path: Path = file(value)
    if (!path.exists() || !path.isFile()) {
        error("--${name} must point to an existing regular ${format} file: ${path}")
    }
}

def validate_precision_parameters(params_map: Map) -> Void {
    ['stopping_gain_cutoff', 'stopping_information_radius_scale', 'stopping_forecast_minutes', 'stopping_max_processing_lag_seconds', 'stopping_min_points',
     'stopping_min_seconds', 'stopping_confirmations', 'stopping_interruption_fraction', 'stopping_recovery_fraction', 'stopping_recovery_intervals'].each { name: String ->
        def value = params_map[name]
        if (!(value instanceof Number) || !Double.isFinite((value as Number).doubleValue()) || value <= 0) error("--${name} must be a positive finite number.")
        if (name in ['stopping_forecast_minutes', 'stopping_max_processing_lag_seconds', 'stopping_min_points', 'stopping_min_seconds', 'stopping_confirmations', 'stopping_recovery_intervals'] && value != Math.floor((value as Number).doubleValue())) error("--${name} must be an integer.")
        if (name in ['stopping_gain_cutoff', 'stopping_interruption_fraction', 'stopping_recovery_fraction'] && value >= 1) error("--${name} must be less than one.")
    }

    if (params_map.stopping_min_points < 6 || params_map.stopping_min_seconds < 3600) {
        error('R monitoring requires at least six points and at least one hour.')
    }
}
