nextflow.enable.types = true

include { imodulon_analysis } from '../subworkflows/imodulon_analysis.nf'
include { next_ica_snapshot_index ; publish_ica_snapshot } from '../lib/ica_publication.nf'
include { Sample ; QuantifiedSampleBatch } from '../lib/sample.nf'

workflow {
    def root: Path = file(params.ica_test_fixture)
    def out_root: Path = file(params.ica_test_output)
    def first: Integer = next_ica_snapshot_index(out_root)
    def samples: List<Sample> = (0..3).collect { index: Integer ->
        record(name: "rep${index % 2}", group: index < 2 ? 'control' : 'treated',
            order: index < 2 ? 0 : 10, bam_dir: root, is_live: index == 0, protocol_run_id: null)
    }
    def batches: List<QuantifiedSampleBatch> = (0..1).collect { index: Integer ->
        record(batch_index: 4 + index * 4, report_sequence: index,
            samples: samples.withIndex().collect { sample, sample_index: Integer ->
                record(sample: sample, batch_index: sample_index == 0 ? index : 0,
                    counts: root.resolve(index == 0 && sample_index == 0 ? 'zero.quant' : "q${sample_index}.quant"))
            })
    }
    imodulon_analysis(channel.fromList(batches), root.resolve('matrix.csv'),
        root.resolve('annotation.gtf'), root.resolve('matrix.csv'),
        [has_gene_map: false, min_gene_coverage: 1.0, log_base: 2.0,
         pseudocount: 1.0, min_read_count: 0, padj_cutoff: 0.05], first, out_root)
    imodulon_analysis.out.snapshots.collect().view { results ->
        assert results*.analysis_index == [first, first + 1]
        assert results*.batch_index == [4, 8]
        def latest = new groovy.json.JsonSlurper().parse(out_root.resolve('ica/latest.json').toFile())
        assert latest.analysis_index == first + 1
        assert latest.status == 'ready'
        assert latest.path == "batch_${first + 1}"
        def initial = new groovy.json.JsonSlurper().parse(results[0].results.resolve('status.json').toFile())
        assert initial.status == 'deferred'
        assert !results[0].results.resolve('activities.tsv').exists()
        assert results[1].results.resolve('activities.tsv').exists()
        assert next_ica_snapshot_index(out_root) == first + 2
        def rejected: Boolean = false
        try {
            publish_ica_snapshot(results[1], out_root)
        }
        catch (exception: Exception) {
            def message: String = exception.message ?: exception.cause?.message ?: "${exception}"
            assert message.contains('refusing to overwrite'): message
            rejected = true
        }
        assert rejected
        'ICA deferred/ready, stopped samples, ordered publication and restart indexing passed'
    }
}
