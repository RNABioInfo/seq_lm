#!/usr/bin/env nextflow
nextflow.enable.types = true
include { apply_acquisition_actions } from '../lib/acquisition_actions.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seqlm-acquisition-actions-')
    def fake: Path = root.resolve('manager')
    def calls: Path = root.resolve('calls')
    java.nio.file.Files.writeString(fake, "#!/bin/sh\nprintf '%s\\n' \"\$*\" >> '${calls}'\nexit 0\n")
    fake.toFile().setExecutable(true)
    def connection: Map = [command: fake.toString(),host:'unused',port:9501,client_certificate:'unused',client_private_key:'unused',ca_certificate:'unused']
    def rows: List<Map> = ['a','b','fixed'].collect { name: String ->
        def directory: Path = root.resolve(name)
        java.nio.file.Files.createDirectories(directory)
        return [sample:name,group:'target',bam_dir:directory.toString(),protocol_run_id:'shared',effectively_live:name!='fixed',eligible:name=='a',newly_eligible:true]
    }
    assert apply_acquisition_actions(rows,'terminate',connection)==null
    assert !calls.exists() && rows[0].run_blocking_samples=='target/b'
    rows[1].eligible=true
    assert apply_acquisition_actions(rows,'log',connection)==null
    assert !calls.exists() && rows[0].action_result=='logged'
    assert apply_acquisition_actions(rows,'terminate',connection)==null
    assert java.nio.file.Files.readAllLines(calls).size()==1
    assert root.resolve('a/STOP').exists() && root.resolve('b/STOP').exists() && !root.resolve('fixed/STOP').exists()
    assert apply_acquisition_actions(rows,'terminate',connection)==null
    assert java.nio.file.Files.readAllLines(calls).size()==1
    def retry: Path = root.resolve('retry')
    java.nio.file.Files.createDirectories(retry)
    def retry_rows: List<Map> = [[sample:'retry',group:'target',bam_dir:retry.toString(),protocol_run_id:'retry',effectively_live:true,eligible:true,newly_eligible:true]]
    java.nio.file.Files.writeString(fake,'#!/bin/sh\nexit 1\n')
    assert apply_acquisition_actions(retry_rows,'terminate',connection)==null
    assert retry_rows[0].action_result=='stop_failed' && !retry.resolve('STOP').exists()
    java.nio.file.Files.writeString(fake,'#!/bin/sh\nexit 0\n')
    assert apply_acquisition_actions(retry_rows,'terminate',connection)==null && retry.resolve('STOP').exists()
    def broken: Path = root.resolve('missing')
    def broken_rows: List<Map> = [[sample:'broken',group:'target',bam_dir:broken.toString(),protocol_run_id:'broken',effectively_live:true,eligible:true,newly_eligible:true]]
    assert apply_acquisition_actions(broken_rows,'terminate',connection).contains('STOP markers could not be created')
    channel.of('Shared acquisition, fixed control, log mode, retry, idempotence and marker-failure checks passed').view()
}
