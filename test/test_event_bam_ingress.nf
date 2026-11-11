#!/usr/bin/env nextflow
nextflow.enable.types = true
include { bam_ingress } from '../lib/bam_ingress.nf'
workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seqlm-event-ingress-')
    def sheet: Path = root.resolve('samples.csv')
    ['fast','slow'].each { name: String -> java.nio.file.Files.createDirectories(root.resolve(name)) }
    java.nio.file.Files.writeString(root.resolve('slow/start.bam'),'initial')
    java.nio.file.Files.writeString(sheet,"alias,group,bam_dir,is_live\nfast,control,${root}/fast,true\nslow,control,${root}/slow,true\n")
    Thread.startDaemon('event-test-writer') {
        Thread.sleep(400)
        java.nio.file.Files.writeString(root.resolve('fast/a.bam'),'partial')
        Thread.sleep(50)
        java.nio.file.Files.writeString(root.resolve('fast/a.bam'),'complete')
        Thread.sleep(450)
        java.nio.file.Files.writeString(root.resolve('fast/b.bam'),'complete')
        Thread.sleep(450)
        java.nio.file.Files.writeString(root.resolve('fast/STOP'),'')
        java.nio.file.Files.writeString(root.resolve('slow/STOP'),'')
    }
    def args: Map = [live_analysis:true,sample_sheet_path:sheet,
        bam_poll_interval_ms:50,bam_stability_polls:3,termination_requested:false,detect_minknow_completion:false,
        precision_event_mode:true,precision_state_directory:root.resolve('state').toString()]
    bam_ingress(args).collect().map { batches ->
        assert batches*.batch_index.toSorted() == (0..<batches.size()).toList()
        def fast: List = batches.findAll { batch -> batch.chunks.any { chunk -> chunk.sample.name=='fast' } }.toSorted { left, right -> left.observation.samples['control\tfast'].version <=> right.observation.samples['control\tfast'].version }
        assert fast.size()==2
        assert fast[0].observation.samples['control\tfast'].version==0
        assert fast[1].observation.samples['control\tfast'].version==1
        assert fast[1].observation.observed_at-fast[0].observation.observed_at < 10
        assert batches.findAll { batch -> batch.chunks.empty }.size()==2
        assert fast[0].chunks[0].bam_paths[0].text=='complete'
        def clock: Map = new groovy.json.JsonSlurper().parse(root.resolve('state/clock.json').toFile()) as Map
        assert clock.schema_version==2
        assert clock.next_versions['control\tfast']==2
        'Existing BAM trigger emits independent timestamped updates, stable files, and terminal drain notifications'
    }.view()
}
