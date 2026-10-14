nextflow.enable.types = true

include {
    GeneContrast ;
    DifferentialExpressionResult ;
    Sample ;
    StabilityAudit ;
    StabilityComparisonInput
} from './sample.nf'
include { optional_file } from '../modules/generic_helpers.nf'
include { tsv_field ; parse_delimited_rows } from './gene_contrasts.nf'
include { get_stop_files_in_dir } from './bam_ingress.nf'

def stability_sample_key(sample) -> String {
    return "${sample.group}\t${sample.name}"
}

def stability_config(params_map: Map, samples: List<Sample>) -> Map {
    return [
        schema_version: 3,
        contrasts: params_map.contrasts,
        behavior: params_map.behavior,
        num_stable_batches: params_map.num_stable_batches,
        max_feature_diff_fraction: params_map.max_feature_diff_fraction,
        max_median_abs_lfc_delta: params_map.max_median_abs_lfc_delta,
        min_jaccard_similarity: params_map.min_jaccard_similarity,
        max_fdr: params_map.max_fdr,
        min_abs_lfc: params_map.min_abs_lfc,
        min_de_calls_for_fraction_metrics: params_map.min_de_calls_for_fraction_metrics,
        max_small_set_call_changes: params_map.max_small_set_call_changes,
        samples: samples.collect { sample ->
            [name: sample.name, group: sample.group, is_live: sample.is_live, protocol_run_id: sample.protocol_run_id]
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
        log.info("No complete prior stability audit was found; stability state starts fresh.")
        return [previous_results: optional_file(), streaks: [:], eligible: [:]]
    }

    def prior_config: Map = new groovy.json.JsonSlurper().parse(config_path.toFile()) as Map
    if (groovy.json.JsonOutput.toJson(prior_config) != groovy.json.JsonOutput.toJson(config)) {
        log.warn('DE stability parameters, sample structure, or protocol-run mappings changed; the previous stability streak is reset.')
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
    return parse_delimited_rows(path.text, '\t', true, false)
}

def parse_stability_boolean(value: Object) -> Boolean {
    return "${value}".equalsIgnoreCase('true')
}


def tsv_value(value: Object) -> String {
    return value == null ? '' : "${value}".replace('\t', ' ').replace('\n', ' ')
}

def write_stability_tsv(path: Path, columns: List<String>, rows: List<Map>) -> Path {
    def lines: List<String> = [columns.join('\t')]

    lines.addAll(rows.collect { row: Map ->
        columns.collect { column: String -> tsv_field(tsv_value(row[column])) }.join('\t')
    })

    java.nio.file.Files.writeString(path, lines.join('\n') + '\n')
    return path
}

def run_seq_run_manager_stop(protocol_run_id: String, minknow_connection: Map) -> Boolean {
    def command: List<String> = [
        (minknow_connection.command ?: 'seq-run-manager') as String,
        'stop',
        '--host',
        minknow_connection.host as String,
        '--port',
        "${minknow_connection.port}",
        '--client-certificate-path',
        "${minknow_connection.client_certificate}",
        '--client-private-key-path',
        "${minknow_connection.client_private_key}",
        '--ca-certificate-path',
        "${minknow_connection.ca_certificate}",
        '--run-id',
        protocol_run_id,
    ]

    try {
        def process: Process = new ProcessBuilder(command).redirectErrorStream(true).start()
        def command_output: String = process.inputStream.getText('UTF-8').trim()
        def exit_code: Integer = process.waitFor()
        if (exit_code == 0) {
            if (command_output) {
                log.info("seq-run-manager stop for run '${protocol_run_id}': ${command_output}")
            }
            return true
        }
        log.warn(
            "Automatic MinKNOW run termination failed for run '${protocol_run_id}' (exit code ${exit_code}). " +
                "The workflow will retry when all active live samples are stable in a later batch.${command_output ? " Output: ${command_output}" : ''}"
        )
        return false
    }
    catch (exception: Exception) {
        log.warn(
            "Automatic MinKNOW run termination failed for run '${protocol_run_id}': ${exception.message}. " +
                'The workflow will retry when all active live samples are stable in a later batch.'
        )
        return false
    }
}

def apply_run_stability_actions(sample_rows: List<Map>, behavior: String, minknow_connection: Map) -> String? {
    def candidates: List<Map> = []

    sample_rows.each { row: Map ->
        row.run_eligible = false
        row.run_blocking_samples = ''
        row.action_result = 'none'
        if (!parse_stability_boolean(row.effectively_live)) {
            return
        }
        if (!get_stop_files_in_dir(Path.of(row.bam_dir as String)).empty) {
            row.action_result = 'stop_exists'
            return
        }
        if (!(row.protocol_run_id as String)) {
            if (parse_stability_boolean(row.eligible)) {
                if (behavior == 'terminate') {
                    row.action_result = 'termination_disabled_no_run_id'
                }
                else if (behavior == 'log' && parse_stability_boolean(row.newly_eligible)) {
                    log.info("Sequencing stability: sample '${row.group}/${row.sample}' is eligible; run readiness is unavailable without a protocol run ID.")
                    row.action_result = 'logged'
                }
            }
            return
        }
        candidates.add(row)
    }

    def runs: Map = candidates.groupBy { row: Map -> row.protocol_run_id as String }
    def marker_error: String? = null

    runs.each { protocol_run_id: String, members: List<Map> ->
        if (marker_error != null) {
            return
        }

        def blockers: List<String> = members
            .findAll { row: Map -> !parse_stability_boolean(row.eligible) }
            .collect { row: Map -> "${row.group}/${row.sample}" as String }
            .toSorted()

        def ready: Boolean = blockers.empty
        members.each { row: Map ->
            row.run_eligible = ready
            row.run_blocking_samples = blockers.join(',')
            if (!ready && parse_stability_boolean(row.eligible)) {
                row.action_result = 'waiting_for_run_stability'
            }
        }

        if (!ready) {
            return
        }

        if (behavior == 'log') {
            log.info("Sequencing stability: all ${members.size()} active live sample(s) in MinKNOW run '${protocol_run_id}' are eligible; the run would be stopped.")
            members.each { row: Map -> row.action_result = 'logged' }
            return
        }

        if (behavior != 'terminate') {
            return
        }

        if (!run_seq_run_manager_stop(protocol_run_id, minknow_connection)) {
            members.each { row: Map -> row.action_result = 'stop_failed' }
            return
        }

        def marker_failures: List<String> = []
        members.each { row: Map ->
            def stop_path: Path = Path.of(row.bam_dir as String).resolve('STOP')
            try {
                java.nio.file.Files.createFile(stop_path)
                row.action_result = 'stop_created'
                log.info("Sequencing stability: stopped MinKNOW run '${protocol_run_id}' and created STOP for sample '${row.group}/${row.sample}' at ${stop_path}.")
            }
            catch (java.nio.file.FileAlreadyExistsException _ignored) {
                if (java.nio.file.Files.isRegularFile(stop_path)) {
                    row.action_result = 'stop_exists'
                }
                else {
                    marker_failures.add("${stop_path}: path exists but is not a regular file" as String)
                }
            }
            catch (exception: Exception) {
                marker_failures.add("${stop_path}: ${exception.message}" as String)
            }
        }

        if (!marker_failures.empty) {
            marker_error = "MinKNOW run '${protocol_run_id}' was stopped, but STOP markers could not be created: ${marker_failures.join('; ')}. Repair these directories and create the missing STOP files before resuming."
        }
    }

    return marker_error
}

workflow differential_stability {
    take:
    differential_results: Channel<DifferentialExpressionResult>
    samples: List<Sample>
    active_samples: List<Sample>
    behavior: String
    settings: Map
    initial_state: Map
    minknow_connection: Map

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
        settings.max_fdr,
        settings.min_abs_lfc,
        settings.min_de_calls_for_fraction_metrics,
        settings.max_small_set_call_changes,
    )

    def streaks: Map<String,Integer> = new LinkedHashMap<String,Integer>(initial_state.streaks as Map)
    def previous_eligibility: Map<String,Boolean> = new LinkedHashMap<String,Boolean>(initial_state.eligible as Map)
    def active_keys: Set<String> = active_samples.collect { sample -> stability_sample_key(sample) }.toSet()
    def contrast_definitions: List<GeneContrast> = settings.contrasts

    audits_ch = assessed_ch.map { result ->
        def assessment_rows: List<Map> = read_stability_tsv(result.assessments)
        def assessments_by_id: Map = assessment_rows.collectEntries { row: Map -> [(row.contrast_id): row] }
        def contrast_rows: List<Map> = contrast_definitions.collect { contrast ->
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
                .findAll { contrast -> sample.group in [contrast.target_group, contrast.reference_group] }
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
                protocol_run_id: sample.protocol_run_id,
                effectively_live: effectively_live,
                required_contrasts: required.join(','),
                consecutive_stable_batches: consecutive_stable_batches,
                eligible: eligible,
                newly_eligible: newly_eligible,
                behavior: behavior,
                action_result: 'none',
            ]
        }
        // fair orders process outputs, not task execution. Apply termination in
        // this sequential handler so later snapshots cannot act before earlier ones.
        def marker_error: String? = apply_run_stability_actions(sample_rows, behavior, minknow_connection)
        if (marker_error != null) {
            error(marker_error)
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
    finalized_audits_ch = finalize_stability_audit(audits_ch, minknow_connection)

    emit:
    finalized_audits_ch
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
    minknow_connection: Map

    output:
    record(
        batch_index: audit.batch_index,
        analysis_index: audit.analysis_index,
        contrast_stability: file('contrast_stability.tsv'),
        sample_stability: file('sample_stability.tsv'),
        config: file('config.json'),
    )

    exec:
    def contrast_columns: List<String> = [
        'analysis_index', 'batch_index', 'baseline', 'contrast_id', 'target_group', 'reference_group',
        'current_de_call_count', 'previous_de_call_count', 'de_call_union_count', 'changed_de_call_count',
        'added_feature_fraction', 'dropped_feature_fraction', 'median_abs_lfc_delta', 'jaccard_similarity',
        'feature_identity_stable', 'effect_size_stable', 'de_calls_stable', 'stable',
        'consecutive_stable_batches', 'reason',
    ]
    def sample_columns: List<String> = [
        'analysis_index', 'batch_index', 'group', 'sample', 'bam_dir', 'protocol_run_id', 'effectively_live',
        'required_contrasts', 'consecutive_stable_batches', 'eligible', 'newly_eligible',
        'behavior', 'action_result', 'run_eligible', 'run_blocking_samples',
    ]
    write_stability_tsv(task.workDir.resolve('contrast_stability.tsv'), contrast_columns, audit.contrast_rows)
    write_stability_tsv(task.workDir.resolve('sample_stability.tsv'), sample_columns, audit.sample_rows)
    java.nio.file.Files.writeString(
        task.workDir.resolve('config.json'),
        new groovy.json.JsonBuilder(audit.config).toPrettyString() + '\n',
    )
}
