nextflow.enable.types = true

include { Sample } from './sample.nf'
include { read_stability_tsv ; write_stability_tsv ; parse_stability_boolean ; stability_sample_key ; apply_run_stability_actions } from './stability.nf'

def new_ica_monitoring_state() -> Map {
    // Start fresh on every invocation, including resume. Cached analyses must
    // not restore a stale stop decision or count as new sequencing information.
    return [fingerprint: null, histories: [:], eligible: [:]]
}

def finite_ica_number(value: Object) -> Double {
    def number: Double = "${value}".toDouble()
    if (!Double.isFinite(number)) {
        throw new IllegalArgumentException('Nonfinite ICA monitoring input')
    }
    return number
}

/** Assess ALL model components. The window spans N informative updates (N+1 snapshots). */
def assess_ica_monitoring(snapshot: Path, samples: List<Sample>, active_samples: List<Sample>, settings: Map, state: Map) -> Map {
    def components: List<String> = []
    def scales: Map = [:]
    def metadata: Map = [:]
    def activities: Map = [:]
    def log_base: Double = 2.0
    def invalid_reason: String? = null
    try {
        def provenance: Map = new groovy.json.JsonSlurper().parse(snapshot.resolve('provenance.json').toFile()) as Map
        def status: Map = new groovy.json.JsonSlurper().parse(snapshot.resolve('status.json').toFile()) as Map
        components = provenance.model.components as List<String>
        scales = provenance.model.component_max_abs_weights as Map
        log_base = finite_ica_number(provenance.settings.log_base)
        if (components.empty || components.toSet().size() != components.size() || log_base <= 1) {
            throw new IllegalArgumentException('Invalid ICA model or log base')
        }
        components.each { component: String ->
            if (!scales.containsKey(component) || finite_ica_number(scales[component]) <= 0) {
                throw new IllegalArgumentException("Missing or invalid model scale for ${component}")
            }
        }
        def rows: List<Map> = read_stability_tsv(snapshot.resolve('sample_metadata.tsv'))
        rows.each { row: Map ->
            def key: String = "${row.group}\t${row.alias}"
            if (metadata.containsKey(key)) {
                throw new IllegalArgumentException('Duplicate ICA sample identity')
            }
            metadata[key] = row
        }
        def selected_samples: List<Sample> = settings.groups == null ? samples : samples.findAll { sample -> sample.group in settings.groups }
        def expected: Set<String> = selected_samples.collect { sample -> stability_sample_key(sample) }.toSet()
        if (metadata.keySet() != expected || rows*.sample_id.toSet().size() != rows.size()) {
            throw new IllegalArgumentException('ICA sample identities do not match the experiment')
        }
        def fingerprint: String = new groovy.json.JsonBuilder([
            model: provenance.model, settings: provenance.settings,
            software: provenance.software,
            controls: (provenance.control_sample_ids as List).toSorted(),
            identities: rows.collect { row: Map -> [row.group, row.alias, row.sample_id] }.toSorted { left, right -> left.toString() <=> right.toString() },
            monitoring: settings,
        ]).toString()
        if (state.fingerprint != fingerprint) {
            state.histories.clear()
            state.eligible.clear()
            state.fingerprint = fingerprint
        }
        if (status.status != 'ready' || !rows.every { row: Map ->
            parse_stability_boolean(row.ready) && finite_ica_number(row.assigned_abundance) > 0 &&
                finite_ica_number(row.assigned_abundance) >= finite_ica_number(provenance.settings.min_read_count)
        }) {
            invalid_reason = 'snapshot_not_ready'
        }
        else {
            def sample_ids: Set<String> = rows*.sample_id.toSet()
            read_stability_tsv(snapshot.resolve('activities_long.tsv')).each { row: Map ->
                if (!components.contains(row.component_id) || !sample_ids.contains(row.sample_id)) {
                    throw new IllegalArgumentException('Unexpected ICA activity identity')
                }
                def key: String = "${row.sample_id}\t${row.component_id}"
                if (activities.containsKey(key)) {
                    throw new IllegalArgumentException('Duplicate ICA activity identity')
                }
                activities[key] = finite_ica_number(row.activity)
            }
            if (activities.size() != rows.size() * components.size()) {
                throw new IllegalArgumentException('Missing ICA activity')
            }
        }
    }
    catch (exception: Exception) {
        invalid_reason = "invalid_snapshot: ${exception.message}"
    }
    if (invalid_reason != null) {
        state.histories.clear()
    }

    def active_keys: Set<String> = active_samples.collect { sample -> stability_sample_key(sample) }.toSet()
    def component_rows: List<Map> = []
    def sample_rows: List<Map> = samples.collect { sample ->
        def key: String = stability_sample_key(sample)
        def row: Map = metadata[key] ?: [:]
        def history: List<Map> = state.histories[key] ?: []
        def informative: Boolean = false
        def selected: Boolean = metadata.containsKey(key)
        if (invalid_reason == null && selected) {
            def depth: Double = finite_ica_number(row.assigned_abundance)
            def current: Map = [depth: depth, activities: components.collectEntries { component: String ->
                [(component): activities["${row.sample_id}\t${component}" as String]]
            }]
            if (!history.empty && depth < (history.last().depth as Number)) {
                history = []
            }
            informative = history.empty || depth >= (history.last().depth as Number) * (1 + settings.min_read_growth_fraction)
            if (informative) {
                history = (history + [current]).takeRight((settings.num_stable_batches as Integer) + 1)
            }
            state.histories[key] = history
            // Even without enough new reads, a moving control baseline can
            // revoke eligibility. Such snapshots never advance the window.
            def evaluated: List<Map> = informative ? history : history + [current]
            components.each { component: String ->
                def values: List<Number> = evaluated.collect { entry: Map -> entry.activities[component] as Number }
                def span: Double = (values.max() as Double) - (values.min() as Double)
                def expression_span: Double = span * finite_ica_number(scales[component])
                def change: Double = Math.exp(expression_span * Math.log(log_base))
                def stable: Boolean = history.size() > 1 && Double.isFinite(expression_span) &&
                    expression_span <= Math.log((settings.max_fold_change as Number).doubleValue()) / Math.log(log_base)
                component_rows.add([
                    group: sample.group, sample: sample.name, component_id: component,
                    activity_range: span, max_abs_gene_weight: scales[component],
                    max_fold_change: change, stable: stable,
                ])
            }
        }
        else {
            components.each { component: String ->
                component_rows.add([group: sample.group, sample: sample.name, component_id: component, stable: false])
            }
        }
        def assessments: List<Map> = component_rows.findAll { item: Map -> item.group == sample.group && item.sample == sample.name }
        def passing: Integer = assessments.count { item: Map -> item.stable }
        def enough_history: Boolean = history.size() >= (settings.num_stable_batches as Integer) + 1
        def stable: Boolean = selected && invalid_reason == null && !components.empty && passing == components.size()
        def effectively_live: Boolean = sample.is_live && active_keys.contains(key)
        def eligible: Boolean = effectively_live && stable && enough_history
        def newly_eligible: Boolean = eligible && !(state.eligible[key] ?: false)
        state.eligible[key] = eligible
        return [
            group: sample.group, sample: sample.name, bam_dir: sample.bam_dir.toAbsolutePath().normalize(),
            protocol_run_id: sample.protocol_run_id, effectively_live: effectively_live,
            method: 'ica', metric: 'ica_monitoring', required_components: components.size(),
            passing_components: passing, informative: informative,
            blocking_components: assessments.findAll { item: Map -> !item.stable }.collect { item: Map -> item.component_id }.join(','),
            consecutive_stable_batches: stable ? Math.max(0, history.size() - 1) : 0,
            eligible: eligible, newly_eligible: newly_eligible,
            reason: !selected ? 'not_selected_for_ica' : invalid_reason ?: (!enough_history ? 'insufficient_informative_history' : !stable ? 'component_change_exceeds_tolerance' : 'all_components_converged'),
            behavior: settings.behavior,
        ]
    }
    return [sample_rows: sample_rows, component_rows: component_rows]
}

workflow ica_monitoring {
    take:
    snapshots: Channel
    samples: List<Sample>
    active_samples: List<Sample>
    settings: Map
    minknow_connection: Map
    output_root: Path

    main:
    def state: Map = new_ica_monitoring_state()
    audits = snapshots.map { snapshot ->
        def provenance: Map = new groovy.json.JsonSlurper().parse(snapshot.results.resolve('provenance.json').toFile()) as Map
        if (provenance.analysis_index != snapshot.analysis_index || provenance.batch_index != snapshot.batch_index ||
            provenance.report_sequence != snapshot.report_sequence) {
            error('ICA monitoring snapshot indices do not match the source record.')
        }
        def assessment: Map = assess_ica_monitoring(snapshot.results, samples, active_samples, settings, state)
        def marker_error: String? = apply_run_stability_actions(assessment.sample_rows, settings.behavior as String, minknow_connection)
        if (marker_error != null) {
            error(marker_error)
        }
        def root: Path = output_root.resolve('stability').resolve('ica')
        java.nio.file.Files.createDirectories(root)
        def destination: Path = root.resolve("batch_${snapshot.analysis_index}")
        def staging: Path = java.nio.file.Files.createTempDirectory(root, '.pending-')
        def index_fields: Map = [batch_index: snapshot.batch_index, analysis_index: snapshot.analysis_index]
        def sample_rows: List<Map> = assessment.sample_rows.collect { row: Map -> index_fields + row }
        def component_rows: List<Map> = assessment.component_rows.collect { row: Map -> index_fields + row }
        write_stability_tsv(staging.resolve('sample_stability.tsv'), [
            'analysis_index', 'batch_index', 'group', 'sample', 'bam_dir', 'protocol_run_id', 'effectively_live',
            'method', 'metric', 'required_components', 'passing_components', 'blocking_components', 'informative',
            'consecutive_stable_batches', 'eligible', 'newly_eligible', 'reason', 'behavior',
            'action_result', 'run_eligible', 'run_blocking_samples',
        ], sample_rows)
        write_stability_tsv(staging.resolve('component_stability.tsv'), [
            'analysis_index', 'batch_index', 'group', 'sample', 'component_id',
            'activity_range', 'max_abs_gene_weight', 'max_fold_change', 'stable',
        ], component_rows)
        java.nio.file.Files.writeString(staging.resolve('config.json'), new groovy.json.JsonBuilder(settings).toPrettyString() + '\n')
        java.nio.file.Files.move(staging, destination, java.nio.file.StandardCopyOption.ATOMIC_MOVE)
        record(batch_index: snapshot.batch_index, analysis_index: snapshot.analysis_index,
            sample_stability: destination.resolve('sample_stability.tsv'),
            component_stability: destination.resolve('component_stability.tsv'))
    }

    emit:
    audits
}
