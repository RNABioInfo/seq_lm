nextflow.enable.types = true

include { Sample ; QuantifiedSampleBatch } from '../lib/sample.nf'
include { imodulon_analysis } from '../subworkflows/imodulon_analysis.nf'
include { ica_monitoring } from '../lib/ica_stability.nf'
include { read_stability_tsv } from '../lib/stability.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-ica-analysis-monitoring-')
    def matrix: Path = root.resolve('matrix.csv')
    def annotation: Path = root.resolve('annotation.gtf')
    java.nio.file.Files.writeString(matrix, ',a,b\ng1,1,0\ng2,0,-2\n')
    java.nio.file.Files.writeString(annotation,
        'chr\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' +
        'chr\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n')
    def samples: List<Sample> = (0..3).collect { index: Integer ->
        record(name: "rep${index % 2}", group: index < 2 ? 'control' : 'treated', order: null,
            bam_dir: root.resolve("sample_${index}"), is_live: true, protocol_run_id: 'shared')
    }
    samples.each { sample -> java.nio.file.Files.createDirectories(sample.bam_dir) }
    def batches: List<QuantifiedSampleBatch> = (0..2).collect { index: Integer ->
        record(batch_index: index, report_sequence: index,
            samples: samples.withIndex().collect { sample, position: Integer ->
                def counts: Path = root.resolve("counts_${position}_${index}.quant")
                java.nio.file.Files.writeString(counts, "tname\tlen\tnum_reads\nt1\t100\t${(position + 1) * 20 * (index + 1)}\nt2\t100\t${(5 - position) * 20 * (index + 1)}\n")
                record(sample: sample, batch_index: index, counts: counts)
            })
    }
    ica_snapshots_ch = imodulon_analysis(channel.fromList(batches), matrix, annotation, matrix, matrix,
        [has_gene_map: false, has_imodulon_table: false, min_gene_coverage: 1.0,
         log_base: 2.0, pseudocount: 1.0, min_read_count: 10, padj_cutoff: 0.05], 0, root)
    ica_audits = ica_monitoring(ica_snapshots_ch, samples, samples,
        [behavior: 'log', num_stable_batches: 2, max_fold_change: 1.05, min_read_growth_fraction: 0.1], [:], root)
    ica_audits.collect().view { audits ->
        audits = audits.toSorted { left, right -> left.analysis_index <=> right.analysis_index }
        assert audits.size() == 3
        def rows: List<Map> = read_stability_tsv(audits.last().sample_stability)
        assert rows.every { row: Map -> row.eligible == 'true' && row.action_result == 'logged' }
        assert rows.every { row: Map -> row.required_components == '2' && row.passing_components == '2' }
        assert samples.every { sample -> !java.nio.file.Files.exists(sample.bam_dir.resolve('STOP')) }
        'Verified actual ICA preparation/projection through all-component monitoring without DEA'
    }
}
