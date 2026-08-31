nextflow.enable.types = true

include {
    DifferentialExpressionResult ;
    Sample ;
    StabilityAudit ;
    StabilityComparisonInput
} from './sample.nf'
include { optional_file ; safe_name } from '../modules/generic_helpers.nf'

def stability_sample_key(sample) -> String {
    return "${sample.group}\t${sample.name}"
}

def stability_config(params_map: Map, samples: List<Sample>) -> Map {
    return [
        schema_version: 1,
        behavior: params_map.behavior,
        num_stable_batches: params_map.num_stable_batches,
        max_feature_diff_fraction: params_map.max_feature_diff_fraction,
        max_median_abs_lfc_delta: params_map.max_median_abs_lfc_delta,
        min_jaccard_similarity: params_map.min_jaccard_similarity,
        max_call_churn_fraction: params_map.max_call_churn_fraction,
        max_lost_call_fraction: params_map.max_lost_call_fraction,
        max_fdr: params_map.max_fdr,
        min_abs_lfc: params_map.min_abs_lfc,
        min_de_calls_for_fraction_metrics: params_map.min_de_calls_for_fraction_metrics,
        max_small_set_call_changes: params_map.max_small_set_call_changes,
        samples: samples.collect { sample ->
            [name: sample.name, group: sample.group, is_live: sample.is_live]
        }.toSorted { left, right ->
            "${left.group}/${left.name}" <=> "${right.group}/${right.name}"
        },
    ]
}

def discover_stability_state(output_root: Path, first_analysis_index: Integer, config: Map) -> Map {
    if (first_analysis_index < 1) {
        return [previous_results: optional_file(), streaks: [:], eligible: [:]]
    }
    def prior_index: Integer = first_analysis_index - 1
    def previous_results: Path = output_root.resolve('differential_expression').resolve("batch_${prior_index}")
    if (!java.nio.file.Files.isDirectory(previous_results)) {
        log.warn("Cannot restore DE stability baseline because ${previous_results} is missing; stability state starts fresh.")
        return [previous_results: optional_file(), streaks: [:], eligible: [:]]
    }
    def audit_dir: Path = output_root.resolve('stability').resolve("batch_${prior_index}")
    def config_path: Path = audit_dir.resolve('config.json')
    def contrasts_path: Path = audit_dir.resolve('contrast_stability.tsv')
    def samples_path: Path = audit_dir.resolve('sample_stability.tsv')
    if (![config_path, contrasts_path, samples_path].every { path: Path -> java.nio.file.Files.isRegularFile(path) }) {
        log.info("Using DE snapshot ${previous_results} as a stability baseline; no compatible prior audit was found.")
        return [previous_results: file(previous_results), streaks: [:], eligible: [:]]
    }
    def prior_config: Map = new groovy.json.JsonSlurper().parse(config_path.toFile()) as Map
    if (groovy.json.JsonOutput.toJson(prior_config) != groovy.json.JsonOutput.toJson(config)) {
        log.warn('DE stability parameters or sample structure changed; the previous stability streak is reset.')
        return [previous_results: optional_file(), streaks: [:], eligible: [:]]
    }
    def contrast_rows: List<Map> = read_stability_tsv(contrasts_path)
    def sample_rows: List<Map> = read_stability_tsv(samples_path)
    def streaks: Map<String,Integer> = contrast_rows.collectEntries { row: Map ->
        [(row.contrast_id as String): (row.consecutive_stable_batches as String).toInteger()]
    }
    def eligible: Map<String,Boolean> = sample_rows.collectEntries { row: Map ->
        [("${row.group}\t${row.sample}" as String): parse_stability_boolean(row.eligible)]
    }
    log.info("Restored DE stability state from ${audit_dir}.")
    return [previous_results: file(previous_results), streaks: streaks, eligible: eligible]
}

def read_stability_tsv(path: Path) -> List<Map> {
    def lines: List<String> = java.nio.file.Files.readAllLines(path)
    if (lines.empty) {
        return []
    }
    def fields: List<String> = lines[0].split('\t', -1).toList()
    return lines.drop(1).findAll { line: String -> !line.empty }.collect { line: String ->
        def values: List<String> = line.split('\t', -1).toList()
        fields.withIndex().collectEntries { field: String, index: Integer ->
            [(field): index < values.size() ? values[index] : '']
        }
    }
}

def parse_stability_boolean(value: Object) -> Boolean {
    return "${value}".equalsIgnoreCase('true')
}

def expected_contrasts(samples: List<Sample>) -> List<Map> {
    def groups: List<String> = samples*.group.unique()
    if (!groups.contains('control')) {
        error("DE stability requires the reference group 'control'.")
    }
    return groups
        .findAll { group: String -> group != 'control' }
        .toSorted()
        .collect { group: String ->
            [
                contrast_id: "group_${safe_name(group)}_vs_control",
                target_group: group,
                reference_group: 'control',
            ]
        }
}

def tsv_value(value: Object) -> String {
    return value == null ? '' : "${value}".replace('\t', ' ').replace('\n', ' ')
}

def write_stability_tsv(path: Path, columns: List<String>, rows: List<Map>) -> Path {
    def lines: List<String> = [columns.join('\t')]
    lines.addAll(rows.collect { row: Map ->
        columns.collect { column: String -> tsv_value(row[column]) }.join('\t')
    })
    java.nio.file.Files.writeString(path, lines.join('\n') + '\n')
    return path
}

/**
 * The action boundary intentionally owns both logging and STOP creation. A
 * future MinKNOW API call belongs here, before the STOP marker is written.
 */
def apply_stability_action(sample_row: Map, behavior: String) -> String {
    if (!parse_stability_boolean(sample_row.newly_eligible)) {
        return 'none'
    }
    def label: String = "${sample_row.group}/${sample_row.sample}"
    if (behavior == 'log') {
        log.info("DE stability: sample '${label}' reached the configured stability threshold and would be stopped.")
        return 'logged'
    }
    if (behavior != 'terminate') {
        return 'none'
    }
    def stop_path: Path = Path.of(sample_row.bam_dir as String).resolve('STOP')
    try {
        java.nio.file.Files.createFile(stop_path)
        log.info("DE stability: created STOP for sample '${label}' at ${stop_path}.")
        return 'stop_created'
    }
    catch (java.nio.file.FileAlreadyExistsException _ignored) {
        log.info("DE stability: STOP already exists for sample '${label}' at ${stop_path}.")
        return 'stop_exists'
    }
}

workflow differential_stability {
    take:
    differential_results: Channel<DifferentialExpressionResult>
    samples: List<Sample>
    active_samples: List<Sample>
    behavior: String
    settings: Map
    initial_state: Map

    main:
    def previous_result: Path = initial_state.previous_results
    stability_comparison_inputs_ch = differential_results.map { result ->
        def baseline: Boolean = previous_result.name == optional_file().name
        def comparison: StabilityComparisonInput = record(
            batch_index: result.batch_index,
            analysis_index: result.analysis_index,
            baseline: baseline,
            previous_results: previous_result,
            current_results: result.results,
        )
        previous_result = result.results
        return comparison
    }

    assessed_ch = assess_differential_stability(
        stability_comparison_inputs_ch,
        settings.max_feature_diff_fraction,
        settings.max_median_abs_lfc_delta,
        settings.min_jaccard_similarity,
        settings.max_call_churn_fraction,
        settings.max_lost_call_fraction,
        settings.max_fdr,
        settings.min_abs_lfc,
        settings.min_de_calls_for_fraction_metrics,
        settings.max_small_set_call_changes,
    )

    def streaks: Map<String,Integer> = new LinkedHashMap<String,Integer>(initial_state.streaks as Map)
    def previous_eligibility: Map<String,Boolean> = new LinkedHashMap<String,Boolean>(initial_state.eligible as Map)
    def active_keys: Set<String> = active_samples.collect { sample -> stability_sample_key(sample) }.toSet()
    def contrast_definitions: List<Map> = expected_contrasts(samples)
    audits_ch = assessed_ch.map { result ->
        def assessment_rows: List<Map> = read_stability_tsv(result.assessments)
        def assessments_by_id: Map = assessment_rows.collectEntries { row: Map -> [(row.contrast_id): row] }
        def contrast_rows: List<Map> = contrast_definitions.collect { contrast: Map ->
            def assessment: Map = assessments_by_id[contrast.contrast_id] ?: [:]
            def stable: Boolean = !result.baseline && parse_stability_boolean(assessment.stable)
            def streak: Integer = stable ? (streaks[contrast.contrast_id] ?: 0) + 1 : 0
            streaks[contrast.contrast_id] = streak
            return [
                analysis_index: result.analysis_index,
                batch_index: result.batch_index,
                baseline: result.baseline,
                contrast_id: contrast.contrast_id,
                target_group: contrast.target_group,
                reference_group: contrast.reference_group,
            ] + assessment + [
                stable: stable,
                consecutive_stable_batches: streak,
                reason: result.baseline ? 'baseline' : assessment.reason,
            ]
        }
        def sample_rows: List<Map> = samples.toSorted { left, right ->
            "${left.group}/${left.name}" <=> "${right.group}/${right.name}"
        }.collect { sample ->
            def key: String = stability_sample_key(sample)
            def required: List<String> = contrast_definitions
                .findAll { contrast: Map -> sample.group in [contrast.target_group, contrast.reference_group] }
                *.contrast_id
            def consecutive_stable_batches: Integer = required.empty
                ? 0
                : required.collect { contrast_id: String -> streaks[contrast_id] ?: 0 }.min() as Integer
            def effectively_live: Boolean = sample.is_live && active_keys.contains(key)
            def eligible: Boolean = effectively_live && !required.empty && required.every { contrast_id: String ->
                (streaks[contrast_id] ?: 0) >= settings.num_stable_batches
            }
            def newly_eligible: Boolean = eligible && !(previous_eligibility[key] ?: false)
            previous_eligibility[key] = eligible
            return [
                analysis_index: result.analysis_index,
                batch_index: result.batch_index,
                group: sample.group,
                sample: sample.name,
                bam_dir: sample.bam_dir.toAbsolutePath().normalize(),
                effectively_live: effectively_live,
                required_contrasts: required.join(','),
                consecutive_stable_batches: consecutive_stable_batches,
                eligible: eligible,
                newly_eligible: newly_eligible,
                behavior: behavior,
                action_result: 'none',
            ]
        }
        return record(
            batch_index: result.batch_index,
            analysis_index: result.analysis_index,
            behavior: behavior,
            contrast_rows: contrast_rows,
            sample_rows: sample_rows,
            config: settings.config,
        ) as StabilityAudit
    }
    finalize_stability_audit(audits_ch)

    emit:
    finalize_stability_audit.out
}

process assess_differential_stability {
    label 'seq_lm_dea'
    container 'rnabioinfo/seq_lm_report:v1.0.0'
    cpus 1
    maxForks 1
    fair true

    input:
    comparison: StabilityComparisonInput
    max_feature_diff_fraction: Number
    max_median_abs_lfc_delta: Number
    min_jaccard_similarity: Number
    max_call_churn_fraction: Number
    max_lost_call_fraction: Number
    max_fdr: Number
    min_abs_lfc: Number
    min_de_calls_for_fraction_metrics: Integer
    max_small_set_call_changes: Integer

    stage:
    stageAs comparison.current_results, 'current_results'
    stageAs comparison.previous_results, 'previous_results'

    output:
    record(
        batch_index: comparison.batch_index,
        analysis_index: comparison.analysis_index,
        baseline: comparison.baseline,
        current_results: comparison.current_results,
        assessments: file('contrast_assessments.tsv'),
    )

    script:
    if (comparison.baseline) {
        return """
            printf 'contrast_id\ttarget_group\treference_group\tstable\treason\n' > contrast_assessments.tsv
            """
    }
    return """
        stability_de_analysis \\
            --previous-results previous_results \\
            --current-results current_results \\
            --output contrast_assessments.tsv \\
            --max-feature-diff-fraction ${max_feature_diff_fraction} \\
            --max-median-abs-lfc-delta ${max_median_abs_lfc_delta} \\
            --min-jaccard-similarity ${min_jaccard_similarity} \\
            --max-call-churn-fraction ${max_call_churn_fraction} \\
            --max-lost-call-fraction ${max_lost_call_fraction} \\
            --max-fdr ${max_fdr} \\
            --min-abs-lfc ${min_abs_lfc} \\
            --min-de-calls-for-fraction-metrics ${min_de_calls_for_fraction_metrics} \\
            --max-small-set-call-changes ${max_small_set_call_changes}
        """
}

process finalize_stability_audit {
    label 'seq_lm'
    cpus 1
    maxForks 1
    fair true
    cache false

    publishDir "${params.out_dir}/stability", mode: 'copy', overwrite: false, saveAs: { fname -> "batch_${audit.analysis_index}/${fname}" }

    input:
    audit: StabilityAudit

    output:
    record(
        batch_index: audit.batch_index,
        analysis_index: audit.analysis_index,
        contrast_stability: file('contrast_stability.tsv'),
        sample_stability: file('sample_stability.tsv'),
        config: file('config.json'),
    )

    exec:
    audit.sample_rows.each { sample_row: Map ->
        sample_row.action_result = apply_stability_action(sample_row, audit.behavior)
    }
    def contrast_columns: List<String> = [
        'analysis_index', 'batch_index', 'baseline', 'contrast_id', 'target_group', 'reference_group',
        'current_de_call_count', 'previous_de_call_count', 'de_call_union_count', 'changed_de_call_count',
        'added_feature_fraction', 'dropped_feature_fraction', 'median_abs_lfc_delta', 'jaccard_similarity',
        'call_churn_fraction', 'lost_call_fraction', 'feature_identity_stable', 'effect_size_stable',
        'de_calls_stable', 'stable', 'consecutive_stable_batches', 'reason',
    ]
    def sample_columns: List<String> = [
        'analysis_index', 'batch_index', 'group', 'sample', 'bam_dir', 'effectively_live',
        'required_contrasts', 'consecutive_stable_batches', 'eligible', 'newly_eligible',
        'behavior', 'action_result',
    ]
    write_stability_tsv(task.workDir.resolve('contrast_stability.tsv'), contrast_columns, audit.contrast_rows)
    write_stability_tsv(task.workDir.resolve('sample_stability.tsv'), sample_columns, audit.sample_rows)
    java.nio.file.Files.writeString(
        task.workDir.resolve('config.json'),
        new groovy.json.JsonBuilder(audit.config).toPrettyString() + '\n',
    )
}
