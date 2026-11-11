nextflow.enable.types = true

include { Sample } from './sample.nf'
include { read_monitoring_tsv ; write_monitoring_tsv ; apply_acquisition_actions ; monitoring_sample_key } from './acquisition_actions.nf'
include { sequencing_container_command } from './sequencing_container.nf'

/** Ordered state transitions; only this handler can request acquisition stops. */
workflow precision_monitoring {
    take:
    snapshots: Channel
    samples: List<Sample>
    active_samples: List<Sample>
    behavior: String
    settings: Map
    observations: Map
    minknow_connection: Map

    main:
    def root: Path = Path.of(settings.output_root as String).resolve('termination').resolve(settings.method as String)

    def state_path: Path = root.resolve('state.json')
    def fallback_started: Number = state_path.exists()
        ? read_precision_state(state_path).started_at
        : System.currentTimeMillis() / 1000.0

    def active_keys: Set<String> = active_samples.collect { sample -> monitoring_sample_key(sample) }.toSet()

    // Serialize state mutations in completion order; never wait for an unrelated result.
    audits = snapshots.map { result ->
        def source: Map = observations[result.batch_index] ?: [
            started_at: fallback_started, observed_at: System.currentTimeMillis() / 1000.0, terminal: true,
        ]
        def represented: Map = result.observation != null ? result.observation : source
        def arrivals: Map = represented.samples ?: [:]
        def observed_at: Number = arrivals.isEmpty()
            ? source.observed_at
            : arrivals.values().collect { value -> value.available_at }.max()
        def terminal: Boolean = source.terminal == true

        def destination: Path = root.resolve("batch_${result.analysis_index}")
        java.nio.file.Files.createDirectories(root)
        def staging: Path = java.nio.file.Files.createTempDirectory(root, '.pending-')

        def closed_keys: List = observations._closed ?: source.closed_samples ?: []
        def rows: List<Map> = samples.collect { sample ->
            def key: String = monitoring_sample_key(sample)
            def declared_live: Boolean = settings.live_analysis && sample.is_live && active_keys.contains(key)
            return [
                key: key,
                name: sample.name,
                group: sample.group,
                bam_dir: sample.bam_dir.toAbsolutePath().normalize().toString(),
                protocol_run_id: sample.protocol_run_id,
                declared_live: declared_live,
                live: declared_live && !closed_keys.contains(key),
            ]
        }

        def payload: Map = [
            method: settings.method,
            samples: rows,
            contrasts: settings.contrasts.collect { contrast ->
                [contrast_id: contrast.contrast_id, target_group: contrast.target_group, reference_group: contrast.reference_group]
            },
            identities: settings.identities,
            config: settings.precision_config,
            started_at: source.started_at,
            observed_at: observed_at,
            arrivals: arrivals,
            completed_at: System.currentTimeMillis() / 1000.0,
            behavior: behavior,
            results: result.results.toAbsolutePath().normalize().toString(),
        ]

        evaluate_precision_snapshot(staging, state_path, settings, payload)
        def sample_rows: List<Map> = read_monitoring_tsv(staging.resolve('sample_monitoring.tsv'))
        sample_rows.each { row: Map ->
            row.batch_index = result.batch_index
            row.analysis_index = result.analysis_index
        }

        // Recheck the represented inputs and each contrast's own lag before acting.
        def contrast_rows: List<Map> = read_monitoring_tsv(staging.resolve('contrast_monitoring.tsv'))
        def timing: Map = [
            now: System.currentTimeMillis() / 1000.0,
            checked_at: payload.completed_at,
            horizon: settings.precision_config.horizon_seconds ?: 1800,
            maximum_lag: settings.precision_config.max_processing_lag_seconds ?: 600,
        ]
        guard_precision_actions(sample_rows, contrast_rows, arrivals, observations._latest ?: [:], timing, terminal)

        def marker_error: String? = terminal ? null : apply_acquisition_actions(sample_rows, behavior, minknow_connection)
        write_monitoring_tsv(staging.resolve('sample_monitoring.tsv'), sample_rows[0].keySet().toList(), sample_rows)

        if (marker_error != null) error(marker_error)

        commit_precision_audit(staging, destination, state_path, sample_rows)

        return record(
            batch_index: result.batch_index,
            analysis_index: result.analysis_index,
            monitoring_results: destination)
    }

    emit:
    audits
}

def read_precision_state(path: Path) -> Map {
    return new groovy.json.JsonSlurper().parse(path.toFile()) as Map
}

def evaluate_precision_snapshot(staging: Path, state: Path, settings: Map, payload: Map) -> Void {
    def request: Path = staging.resolve('input.json')
    java.nio.file.Files.writeString(request, new groovy.json.JsonBuilder(payload).toPrettyString())

    def engine: String? = settings.container_engine as String?
    def executable: List<String> = engine ? ['precision-termination'] : ['python3', settings.command as String]
    def arguments: List<String> = executable + ['--input', request.toString(), '--state', state.toString(),
                                                '--output', staging.toString()]
    if (engine) {
        def inputs: List<Path> = [Path.of(payload.results as String)]
        if (state.exists()) inputs.add(state)
        arguments = sequencing_container_command(engine, arguments, inputs, [staging], false)
    }
    def process: Process = new ProcessBuilder(arguments).redirectErrorStream(true).start()
    def output: String = process.inputStream.getText('UTF-8')
    if (process.waitFor() != 0) {
        error("Precision monitoring failed: ${output}")
    }
    return null
}

def commit_precision_audit(staging: Path, destination: Path, state: Path, rows: List<Map>) -> Void {
    def staged_state: Path = staging.resolve('state.json')
    def committed: Map = read_precision_state(staged_state)
    rows.each { row: Map ->
        def key: String = "${row.group}\t${row.sample}"
        committed.sample_eligible[key] = "${row.eligible}" == 'true'
    }

    java.nio.file.Files.writeString(staged_state, new groovy.json.JsonBuilder(committed).toString())
    java.nio.file.Files.move(staging, destination, java.nio.file.StandardCopyOption.ATOMIC_MOVE)

    def temporary_state: Path = java.nio.file.Files.createTempFile(state.parent, '.state-', '.json')
    java.nio.file.Files.copy(destination.resolve('state.json'), temporary_state, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    java.nio.file.Files.move(temporary_state, state, java.nio.file.StandardCopyOption.REPLACE_EXISTING,
        java.nio.file.StandardCopyOption.ATOMIC_MOVE)
    return null
}

/** A newer member invalidates the contrast for every acquisition that depends on it. */
def guard_precision_actions(rows: List<Map>, contrasts: List<Map>, arrivals: Map, latest: Map,
                            timing: Map, terminal: Boolean) -> Void {
    def time_since_check: Number = timing.now - timing.checked_at
    def blocked: Set<String> = contrasts.findAll { contrast: Map ->
        def lag: String = contrast.processing_lag_seconds as String
        lag && lag.toDouble() + time_since_check > timing.maximum_lag
    }.collect { contrast: Map -> contrast.contrast_id as String }.toSet()

    rows.each { row: Map ->
        def key: String = "${row.group}\t${row.sample}"
        def represented: Map? = arrivals[key] as Map?
        def newest: Map? = latest[key] as Map?
        if (represented_input_stale(row, represented, newest, timing.now as Number, timing.horizon as Number)) {
            blocked.addAll((row.required_contrasts as String).tokenize(','))
        }
    }

    rows.each { row: Map ->
        def affected: Boolean = (row.required_contrasts as String).tokenize(',').any { id: String ->
            blocked.contains(id)
        }
        if ("${row.effectively_live}" == 'true' && (terminal || affected)) {
            row.eligible = false
            row.newly_eligible = false
            row.reason = terminal ? 'terminal_snapshot' : 'newer_input_or_stale_result'
        }
    }
    return null
}

/** A newer input invalidates a contrast even if the member has since closed. */
def represented_input_stale(row: Map, represented: Map?, latest: Map?, now: Number, horizon: Number) -> Boolean {
    if (latest != null && (represented == null || latest.version > represented.version)) {
        return true
    }
    if ("${row.effectively_live}" != 'true') {
        return false
    }
    return represented == null || now - (represented.available_at as Number) > horizon
}
