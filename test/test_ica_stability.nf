#!/usr/bin/env nextflow
nextflow.enable.types = true

include { Sample } from '../lib/sample.nf'
include { assess_ica_monitoring ; new_ica_monitoring_state ; ica_monitoring } from '../lib/ica_stability.nf'
include { read_stability_tsv } from '../lib/stability.nf'

def snapshot(root: Path, index: Integer, depth: Number, values: Map, scales: Map, extra: Map = [:]) -> Path {
    def path: Path = root.resolve("snapshot_${index}_${System.nanoTime()}")
    java.nio.file.Files.createDirectories(path)
    def provenance: Map = [batch_index: index, analysis_index: index, report_sequence: index,
        model: [components: ['a', 'b'], component_max_abs_weights: scales, hashes: [matrix: 'fixed-model']],
        settings: [log_base: 2, min_read_count: 10, pseudocount: 1],
        software: [version: 'test'], control_sample_ids: ['control-one', 'control-two']]
    provenance.putAll(extra)
    java.nio.file.Files.writeString(path.resolve('provenance.json'), new groovy.json.JsonBuilder(provenance).toString())
    java.nio.file.Files.writeString(path.resolve('status.json'), '{"status":"ready"}')
    java.nio.file.Files.writeString(path.resolve('sample_metadata.tsv'),
        "sample_id\talias\tgroup\tassigned_abundance\tready\ncontrol-one\tone\tcontrol\t${depth}\ttrue\ncontrol-two\ttwo\tcontrol\t${depth}\ttrue\ntreated-one\tone\ttreated\t${depth}\ttrue\n")
    def lines: List<String> = ['sample_id\tcomponent_id\tactivity']
    ['control-one', 'control-two', 'treated-one'].each { id: String ->
        ['a', 'b'].each { component: String ->
            lines.add("${id}\t${component}\t${values["${id}/${component}" as String] ?: 0}")
        }
    }
    java.nio.file.Files.writeString(path.resolve('activities_long.tsv'), lines.join('\n') + '\n')
    return path
}

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-ica-monitoring-')
    def samples: List<Sample> = [
        record(name: 'one', group: 'control', order: null, bam_dir: root.resolve('c1'), is_live: true, protocol_run_id: 'shared'),
        record(name: 'two', group: 'control', order: null, bam_dir: root.resolve('c2'), is_live: true, protocol_run_id: 'shared'),
        record(name: 'one', group: 'treated', order: null, bam_dir: root.resolve('t1'), is_live: true, protocol_run_id: 'shared'),
    ]
    samples.each { sample -> java.nio.file.Files.createDirectories(sample.bam_dir) }
    def settings: Map = [behavior: 'log', num_stable_batches: 2, max_fold_change: 1.05, min_read_growth_fraction: 0.1]
    def state: Map = new_ica_monitoring_state()
    def scales: Map = [a: 0.5, b: 1.0]
    def first: Path = snapshot(root, 0, 100, [:], scales)
    def result: Map = assess_ica_monitoring(first, samples, samples, settings, state)
    assert result.sample_rows.every { row: Map -> !row.eligible }
    result = assess_ica_monitoring(snapshot(root, 1, 120, ['treated-one/a': 0.05], scales), samples, samples, settings, state)
    assert result.sample_rows.every { row: Map -> !row.eligible }
    result = assess_ica_monitoring(snapshot(root, 2, 144, ['treated-one/a': 0.1], scales), samples, samples, settings, state)
    assert result.sample_rows.every { row: Map -> row.eligible }
    def score: Number = result.component_rows.find { row: Map -> row.group == 'treated' && row.component_id == 'a' }.max_fold_change
    assert Math.abs(score.doubleValue() - Math.pow(2, 0.05)) < 1e-12
    // ICA scale/sign changes must not change the convergence score.
    [-10.0, 10.0].each { multiplier: Number ->
        def scaled_state: Map = new_ica_monitoring_state()
        [0, 1, 2].each { index: Integer ->
            result = assess_ica_monitoring(snapshot(root, index, 100 * Math.pow(1.2, index),
                ['treated-one/a': index * 0.05 / multiplier], [a: 0.5 * Math.abs(multiplier.doubleValue()), b: 1]),
                samples, samples, settings, scaled_state)
        }
        assert result.sample_rows.every { row: Map -> row.eligible }
        assert Math.abs(result.component_rows.find { row: Map -> row.group == 'treated' && row.component_id == 'a' }.max_fold_change.doubleValue() - score.doubleValue()) < 1e-12
    }
    // One component in one sample blocks that sample, even if all others pass.
    result = assess_ica_monitoring(snapshot(root, 3, 173, ['treated-one/b': -0.2], scales), samples, samples, settings, state)
    assert !result.sample_rows.find { row: Map -> row.group == 'treated' }.eligible
    assert result.sample_rows.find { row: Map -> row.group == 'treated' }.blocking_components == 'b'
    // Repeating depth cannot establish a window; tiny additions also do not count.
    def repeated_state: Map = new_ica_monitoring_state()
    [100, 100, 101, 109].each { depth: Integer ->
        result = assess_ica_monitoring(snapshot(root, depth, depth, [:], scales), samples, samples, settings, repeated_state)
        assert result.sample_rows.every { row: Map -> !row.eligible }
    }
    // A changing control reference can revoke eligibility without advancing depth.
    result = assess_ica_monitoring(snapshot(root, 4, 173, ['treated-one/b': 1], scales), samples, samples, settings, state)
    assert !result.sample_rows.find { row: Map -> row.group == 'treated' }.eligible
    assert result.sample_rows.every { row: Map -> !row.informative }
    // Pairwise changes can pass while accumulated drift must fail.
    def drift_state: Map = new_ica_monitoring_state()
    [0, 1, 2].each { index: Integer ->
        result = assess_ica_monitoring(snapshot(root, index, 100 * Math.pow(1.2, index), ['treated-one/b': index * 0.06], scales), samples, samples, settings, drift_state)
    }
    assert !result.sample_rows.find { row: Map -> row.group == 'treated' }.eligible
    // Deferred, incomplete and nonfinite data discard the previous window.
    def deferred: Path = snapshot(root, 5, 200, [:], scales)
    java.nio.file.Files.writeString(deferred.resolve('status.json'), '{"status":"deferred"}')
    result = assess_ica_monitoring(deferred, samples, samples, settings, state)
    assert result.sample_rows.every { row: Map -> !row.eligible }
    assert state.histories.isEmpty()
    def missing: Path = snapshot(root, 6, 240, [:], scales)
    java.nio.file.Files.writeString(missing.resolve('activities_long.tsv'), 'sample_id\tcomponent_id\tactivity\ncontrol-one\ta\tNaN\n')
    result = assess_ica_monitoring(missing, samples, samples, settings, state)
    assert result.sample_rows.every { row: Map -> !row.eligible && row.reason.startsWith('invalid_snapshot') }
    // A changed model, control set, transform, or threshold restarts monitoring.
    [0, 1, 2].each { index: Integer ->
        result = assess_ica_monitoring(snapshot(root, index, 100 * Math.pow(1.2, index), [:], scales), samples, samples, settings, state)
    }
    assert result.sample_rows.every { row: Map -> row.eligible }
    result = assess_ica_monitoring(snapshot(root, 7, 200, [:], scales, [software: [version: 'changed']]), samples, samples, settings, state)
    assert result.sample_rows.every { row: Map -> !row.eligible }
    // Analysis selection cannot remove active unselected members from run safety.
    def extra_sample: Sample = record(name: 'omitted', group: 'QC only', order: null,
        bam_dir: root.resolve('omitted'), is_live: true, protocol_run_id: samples[0].protocol_run_id)
    java.nio.file.Files.createDirectories(extra_sample.bam_dir)
    def selected_state: Map = new_ica_monitoring_state()
    [0, 1, 2].each { index: Integer ->
        result = assess_ica_monitoring(snapshot(root, index, 100 * Math.pow(1.2, index), [:], scales),
            samples + [extra_sample], samples + [extra_sample], settings + [groups: ['control', 'treated']], selected_state)
    }
    assert result.sample_rows.find { row: Map -> row.group == 'QC only' }.reason == 'not_selected_for_ica'
    assert !result.sample_rows.find { row: Map -> row.group == 'QC only' }.eligible
    assert result.sample_rows.findAll { row: Map -> row.group != 'QC only' }.every { row: Map -> row.eligible }
    // Integration: publish immutable audits and stop a shared run only when all pass.
    def manager: Path = root.resolve('fake-manager')
    def calls: Path = root.resolve('calls')
    java.nio.file.Files.writeString(manager, "#!/bin/sh\nprintf 'stop\\n' >> '${calls}'\n")
    manager.toFile().setExecutable(true)
    def inputs: List = (0..3).collect { index: Integer ->
        record(batch_index: index, analysis_index: index, report_sequence: index,
            results: snapshot(root, index, 100 * Math.pow(1.2, index), ['treated-one/b': index == 0 ? 1 : 0], scales))
    }
    ica_audits = ica_monitoring(channel.fromList(inputs), samples, samples, settings + [behavior: 'terminate'],
        [command: manager.toString(), host: 'test', port: 9501, client_certificate: root, client_private_key: root, ca_certificate: root], root)
    ica_audits.collect().view { audits ->
        audits = audits.toSorted { left, right -> left.analysis_index <=> right.analysis_index }
        assert audits.size() == 4
        assert read_stability_tsv(audits[2].sample_stability).find { row: Map -> row.group == 'treated' }.action_result == 'none'
        assert read_stability_tsv(audits[3].sample_stability).every { row: Map -> row.action_result == 'stop_created' }
        assert java.nio.file.Files.readAllLines(calls).size() == 1
        assert samples.every { sample -> java.nio.file.Files.isRegularFile(sample.bam_dir.resolve('STOP')) }
        assert read_stability_tsv(audits[3].component_stability).size() == 6
        'Verified ICA metric, all-component gating, history safeguards, and shared-run termination'
    }
}
