nextflow.enable.types = true

include { DifferentialExpressionResult ; Sample } from './sample.nf'
include { read_monitoring_tsv ; write_monitoring_tsv ; apply_acquisition_actions } from './acquisition_actions.nf'

/** Ordered state transitions; only this handler can request acquisition stops. */
workflow dea_monitoring {
    take:
    differential_results: Channel<DifferentialExpressionResult>
    samples: List<Sample>
    active_samples: List<Sample>
    behavior: String
    settings: Map
    observations: Map
    minknow_connection: Map

    main:
    def root: Path = Path.of(settings.output_root as String).resolve('termination').resolve('dea')

    def fallback_state: Path = root.resolve('state.json')
    def fallback_started: Number = java.nio.file.Files.exists(fallback_state)
        ? (new groovy.json.JsonSlurper().parse(fallback_state.toFile()) as Map).started_at
        : System.currentTimeMillis() / 1000.0

    def active_keys: Set<String> = active_samples.collect { sample -> "${sample.group}\t${sample.name}" as String }.toSet()

    // Serialize state mutations in completion order; never wait for an unrelated result.
    audits = differential_results.map { result ->
        def source: Map = observations[result.batch_index] ?: [started_at: fallback_started, observed_at: System.currentTimeMillis() / 1000.0, terminal:true]
        def arrivals: Map = result.observation?.samples ?: source.samples ?: [:]
        def observed: Map = [started_at:source.started_at,
            observed_at: arrivals.isEmpty() ? source.observed_at : arrivals.values().collect { value -> value.available_at }.max(),
            terminal:source.terminal,closed_samples:observations._closed ?: source.closed_samples ?: []]

        def destination: Path = root.resolve("batch_${result.analysis_index}")
        java.nio.file.Files.createDirectories(root)
        def staging: Path = java.nio.file.Files.createTempDirectory(root, '.pending-')

        def rows: List<Map> = samples.collect { sample ->
            [
                key: "${sample.group}\t${sample.name}", 
                name: sample.name, group: sample.group,
                bam_dir: sample.bam_dir.toAbsolutePath().normalize().toString(), 
                protocol_run_id: sample.protocol_run_id,
                declared_live: settings.live_analysis && sample.is_live && active_keys.contains("${sample.group}\t${sample.name}" as String),
                live: settings.live_analysis && sample.is_live && active_keys.contains("${sample.group}\t${sample.name}" as String) && !(observed.closed_samples ?: []).contains("${sample.group}\t${sample.name}" as String)
            ]
        }

        def payload: Map = [
            samples: rows, 
            contrasts: settings.contrasts.collect { contrast ->
                [contrast_id: contrast.contrast_id, target_group: contrast.target_group, reference_group: contrast.reference_group]
            }, 
            identities: settings.identities, 
            config: settings.dea_config,
            started_at: observed.started_at, 
            observed_at: observed.observed_at, 
            arrivals: arrivals,
            completed_at: System.currentTimeMillis() / 1000.0, 
            behavior: behavior,
            results: result.results.toAbsolutePath().normalize().toString()
            ]

        def request: Path = staging.resolve('input.json')
        java.nio.file.Files.writeString(request, new groovy.json.JsonBuilder(payload).toPrettyString())

        def command: List<String> = ['python3', settings.command as String, '--input', request.toString(),
            '--state', root.resolve('state.json').toString(), '--output', staging.toString()]
        def process: Process = new ProcessBuilder(command).redirectErrorStream(true).start()
        def output: String = process.inputStream.getText('UTF-8')

        if (process.waitFor() != 0) error("R monitoring failed: ${output}")

        def sample_rows: List<Map> = read_monitoring_tsv(staging.resolve('sample_monitoring.tsv'))
        sample_rows.each { row: Map -> row.batch_index = result.batch_index; row.analysis_index = result.analysis_index }

        // Recheck against the input front immediately before shared acquisition actions.
        def action_now: Number = System.currentTimeMillis() / 1000.0
        def horizon: Number = settings.dea_config.horizon_seconds ?: 1800
        def maximum_lag: Number = settings.dea_config.max_processing_lag_seconds ?: 600

        sample_rows.each { row: Map ->
            def key: String = "${row.group}\t${row.sample}"
            def represented: Map = arrivals[key]
            def latest: Map = observations._latest?.get(key)
            if (parse_guard_blocked(row, represented, latest, action_now, horizon, maximum_lag)) {
                row.eligible = false
                row.newly_eligible = false
                row.reason = 'newer_input_or_stale_result'
            }
        }

        def marker_error: String? = observed.terminal ? null : apply_acquisition_actions(sample_rows, behavior, minknow_connection)
        write_monitoring_tsv(staging.resolve('sample_monitoring.tsv'), sample_rows[0].keySet().toList(), sample_rows)

        if (marker_error != null) error(marker_error)

        def committed: Map = new groovy.json.JsonSlurper().parse(staging.resolve('state.json').toFile()) as Map
        sample_rows.each { row -> committed.sample_eligible["${row.group}\t${row.sample}" as String] = "${row.eligible}" == 'true' }

        java.nio.file.Files.writeString(staging.resolve('state.json'),new groovy.json.JsonBuilder(committed).toString())
        java.nio.file.Files.move(staging, destination, java.nio.file.StandardCopyOption.ATOMIC_MOVE)

        def temporary_state: Path = java.nio.file.Files.createTempFile(root,'.state-', '.json')
        java.nio.file.Files.copy(destination.resolve('state.json'),temporary_state,java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        java.nio.file.Files.move(temporary_state,root.resolve('state.json'),java.nio.file.StandardCopyOption.ATOMIC_MOVE,java.nio.file.StandardCopyOption.REPLACE_EXISTING)

        return record(
            batch_index: result.batch_index, 
            analysis_index: result.analysis_index,
            monitoring_results: destination)
    }
    
    emit:
    audits
}

/** No action on unrepresented files, aged data, or a result delayed in processing. */
def parse_guard_blocked(row: Map, represented: Map?, latest: Map?, now: Number, horizon: Number, maximum_lag: Number) -> Boolean {
    if ("${row.effectively_live}" != 'true') return false
    if (represented == null) return true
    return (latest != null && latest.version > represented.version) ||
        now - (represented.available_at as Number) > horizon ||
        "${row.processing_lag_seconds}".toDouble() + now - "${row.result_checked_at}".toDouble() > maximum_lag
}
