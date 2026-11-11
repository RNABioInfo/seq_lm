nextflow.enable.types = true

include { tsv_field ; parse_delimited_rows } from './gene_contrasts.nf'
include { get_stop_files_in_dir } from './bam_ingress.nf'
include { sequencing_container_command } from './sequencing_container.nf'

def monitoring_sample_key(sample) -> String {
    return "${sample.group}\t${sample.name}"
}

def read_monitoring_tsv(path: Path) -> List<Map> {
    return parse_delimited_rows(path.text, '\t', true, false)
}

def parse_monitoring_boolean(value: Object) -> Boolean {
    return "${value}".equalsIgnoreCase('true')
}


def tsv_value(value: Object) -> String {
    return value == null ? '' : "${value}".replace('\t', ' ').replace('\n', ' ')
}

def write_monitoring_tsv(path: Path, columns: List<String>, rows: List<Map>) -> Path {
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
        def engine: String? = minknow_connection.container_engine as String?
        if (engine) {
            def credentials: List<Path> = [
                Path.of(minknow_connection.client_certificate as String),
                Path.of(minknow_connection.client_private_key as String),
                Path.of(minknow_connection.ca_certificate as String),
            ]
            command = sequencing_container_command(engine, command, credentials, [], true)
        }
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

def apply_acquisition_actions(sample_rows: List<Map>, behavior: String, minknow_connection: Map) -> String? {
    def candidates: List<Map> = []

    sample_rows.each { row: Map ->
        row.run_eligible = false
        row.run_blocking_samples = ''
        row.action_result = 'none'
        if (!parse_monitoring_boolean(row.effectively_live)) {
            return
        }
        if (!get_stop_files_in_dir(Path.of(row.bam_dir as String)).empty) {
            row.action_result = 'stop_exists'
            return
        }
        if (!(row.protocol_run_id as String)) {
            if (parse_monitoring_boolean(row.eligible)) {
                if (behavior == 'terminate') {
                    row.action_result = 'termination_disabled_no_run_id'
                }
                else if (behavior == 'log' && parse_monitoring_boolean(row.newly_eligible)) {
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
            .findAll { row: Map -> !parse_monitoring_boolean(row.eligible) }
            .collect { row: Map -> "${row.group}/${row.sample}" as String }
            .toSorted()

        def ready: Boolean = blockers.empty
        members.each { row: Map ->
            row.run_eligible = ready
            row.run_blocking_samples = blockers.join(',')
            if (!ready && parse_monitoring_boolean(row.eligible)) {
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

