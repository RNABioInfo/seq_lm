#!/usr/bin/env nextflow
nextflow.enable.types = true
include { precision_monitoring } from '../lib/precision_monitoring.nf'
include { DifferentialExpressionResult ; Sample ; GeneContrast } from '../lib/sample.nf'
params {
    test_stopping_genes: Boolean = false
}
workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seqlm-r-integration-')
    def samples: List<Sample> = ['control','target'].collectMany { group: String ->
        (1..3).collect { index: Integer ->
            def directory: Path = root.resolve("${group}_${index}")
            java.nio.file.Files.createDirectories(directory)
            return record(name:"${index}" as String,group:group,order:null,bam_dir:directory,is_live:group=='target',protocol_run_id:null,minknow_run_dir:null)
        }
    }
    def contrasts: List<GeneContrast> = [record(contrast_id:'pair',target_group:'target',reference_group:'control')]
    def observations: Map = [:]
    def now: Number = System.currentTimeMillis()/1000.0
    def results: List<DifferentialExpressionResult> = (0..7).collect { index: Integer ->
        def dir: Path = root.resolve("results_${index}/pair")
        java.nio.file.Files.createDirectories(dir)
        def header: String = 'feature_id\tR\tci_lower\tci_upper\n'
        java.nio.file.Files.writeString(dir.resolve('precision.tsv'),header+'g1\t1\t-1\t1\n'+(params.test_stopping_genes ? 'g2\t10\t-10\t10\n' : ''))
        if (params.test_stopping_genes) {
            java.nio.file.Files.writeString(dir.parent.resolve('dea_stopping_features.tsv'), 'feature_id\ng1\n')
        }
        java.nio.file.Files.writeString(dir.resolve('precision_samples.tsv'),'key\tgroup\tassigned_reads\n'+samples.collect { sample ->
            "\"${sample.group}\t${sample.name}\"\t${sample.group}\t${sample.group=='control' ? 100000 : (index+1)*10000}\n"
        }.join(''))
        observations[index]=[started_at:now-4200,observed_at:now-4200+index*600,terminal:false]
        return record(batch_index:index,analysis_index:index,report_sequence:index,results:dir.parent)
    }
    def settings: Map = [method:'dea',live_analysis:true,output_root:root.toString(),command:projectDir.parent.resolve('bin/precision-termination').toString(),contrasts:contrasts,
        container_engine: workflow.containerEngine,
        identities:params.test_stopping_genes ? [dea_stopping_genes:'fixture-gmt-hash'] : [:],precision_config:[:]]
    def audits = precision_monitoring(channel.fromList(results),samples,samples.findAll { sample -> sample.is_live },'log',settings,observations,[:])
    audits.collect().map { rows ->
        assert rows.size()==8
        assert root.resolve('termination/dea/state.json').exists()
        assert rows[-1].monitoring_results.resolve('sample_monitoring.tsv').exists()
        assert rows[-1].monitoring_results.resolve('forecast_history.tsv').exists()
        def audit: Map = new groovy.json.JsonSlurper().parse(rows[-1].monitoring_results.resolve('audit.json').toFile()) as Map
        assert audit.contrasts[0].R_information == 1
        assert audit.contrasts[0].uncertainty_feature_count == 1
        assert audit.contrasts[0].uncertainty_scope == (params.test_stopping_genes ? 'selected_features' : 'all_features')
        return 'R snapshot ingestion, serialized state, explicit times and audit publication passed'
    }.view()
}
