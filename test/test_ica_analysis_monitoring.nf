nextflow.enable.types = true

include { Sample ; QuantifiedSampleBatch ; GeneContrast } from '../lib/sample.nf'
include { imodulon_analysis } from '../subworkflows/imodulon_analysis.nf'
include { precision_monitoring } from '../lib/precision_monitoring.nf'
include { read_monitoring_tsv } from '../lib/acquisition_actions.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-ica-precision-')
    def matrix: Path = root.resolve('matrix.csv')
    def selection: Path = root.resolve('selected.txt')
    def annotation: Path = root.resolve('annotation.gtf')
    java.nio.file.Files.writeString(matrix, ',a,b\ng1,1,0\ng2,0,-2\n')
    java.nio.file.Files.writeString(selection, ' a \n\n')
    java.nio.file.Files.writeString(annotation,
        'chr\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' +
        'chr\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n')

    def samples: List<Sample> = ['control', 'early', 'late'].collectMany { group: String ->
        def replicate_count: Integer = group == 'early' ? 3 : 2
        (0..<replicate_count).collect { replicate: Integer ->
            def directory: Path = root.resolve("${group}_${replicate}")
            java.nio.file.Files.createDirectories(directory)
            record(name: "rep${replicate}" as String, group: group, order: null, bam_dir: directory,
                is_live: group != 'control', protocol_run_id: 'shared')
        }
    }
    def contrasts: List<GeneContrast> = ['early', 'late'].collect { group: String ->
        record(contrast_id: group, target_group: group, reference_group: 'control')
    }
    def observations: Map = [:]
    def now: Number = System.currentTimeMillis() / 1000.0
    def batches: List<QuantifiedSampleBatch> = (0..2).collect { index: Integer ->
        def time: Number = now - (2 - index) * 600
        def arrivals: Map = samples.findAll { sample -> sample.is_live }.collectEntries { sample ->
            [("${sample.group}\t${sample.name}" as String): [available_at: time, version: index]]
        }
        observations[index] = [started_at: now - 1200, observed_at: time, samples: arrivals, terminal: false]
        record(batch_index: index, report_sequence: index, observation: [samples: arrivals],
            samples: samples.withIndex().collect { sample, position: Integer ->
                def counts: Path = root.resolve("counts_${position}_${index}.quant")
                def scale: Integer = sample.group == 'control' ? 1 : index + 1
                if (sample.group == 'late' && sample.name == 'rep1' && index < 2) {
                    scale = 0
                }
                def depth_factor: Integer = sample.group == 'early' && sample.name == 'rep2' ? 20 : 1500
                java.nio.file.Files.writeString(counts,
                    "tname\tlen\tnum_reads\nt1\t100\t${(position + 1) * depth_factor * scale}\nt2\t100\t${(7 - position) * depth_factor * scale}\n")
                record(sample: sample, batch_index: index, counts: counts)
            })
    }

    def snapshots = imodulon_analysis(channel.fromList(batches), matrix, annotation, matrix, matrix,
        [control_group: 'control', has_gene_map: false, has_imodulon_table: false, min_gene_coverage: 1.0,
         log_base: 2.0, pseudocount: 1.0, min_read_count: 10000, min_replicate_sample_count: 2, padj_cutoff: 0.05,
         monitoring: true, stopping_components: selection], 0, root)
    def settings: Map = [method: 'ica', live_analysis: true, output_root: root.toString(),
        command: projectDir.parent.resolve('bin/precision-termination').toString(),
        container_engine: workflow.containerEngine,
        contrasts: contrasts, identities: [:], precision_config: [max_processing_lag_seconds: 3600]]
    def audits = precision_monitoring(snapshots, samples, samples, 'log', settings, observations, [:])
    audits.collect().view { collected ->
        assert collected.size() == 3
        def first: Path = root.resolve('termination/ica/batch_0')
        def final_batch: Path = root.resolve('termination/ica/batch_2')
        def first_metrics: List<Map> = read_monitoring_tsv(first.resolve('contrast_monitoring.tsv'))
        assert first_metrics.find { row: Map -> row.contrast_id == 'early' }.uncertainty_feature_count == '1'
        assert first_metrics.find { row: Map -> row.contrast_id == 'late' }.reason == 'results_not_ready'
        def partial: List<Map> = read_monitoring_tsv(root.resolve('ica/batch_0/monitoring_differential_activity.tsv'))
        assert partial.size() == 2
        assert partial.every { row: Map -> row.target_group == 'early' && row.target_n == '3' }
        assert !root.resolve('ica/batch_0/differential_activity.tsv').exists()
        assert read_monitoring_tsv(root.resolve('ica/batch_2/differential_activity.tsv')).size() == 4

        def metrics: List<Map> = read_monitoring_tsv(final_batch.resolve('contrast_monitoring.tsv'))
        assert metrics.every { row: Map -> row.uncertainty_scope == 'selected_features' && row.uncertainty_feature_count == '1' }
        def rows: List<Map> = read_monitoring_tsv(final_batch.resolve('sample_monitoring.tsv'))
        assert rows.every { row: Map -> row.method == 'ica' && row.metric == 'R_information' && row.eligible == 'false' }
        assert rows.findAll { row: Map -> row.group == 'control' }.every { row: Map -> row.reason == 'not_monitored' }
        assert samples.every { sample -> !sample.bam_dir.resolve('STOP').exists() }
        'Verified ICA projection, independent readiness, selected precision and audit publication without DEA'
    }
}
