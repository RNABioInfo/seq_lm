nextflow.enable.types = true

include { ChunkQCResult ; QuantifiedSample ; Sample } from './sample.nf'
include { safe_name } from '../modules/generic_helpers.nf'

def sample_checkpoint_schema() -> Integer {
    return 1
}

def sample_final_marker() -> String {
    return 'FINAL'
}

def sample_checkpoint_tools() -> Map<String,String> {
    return [samtools: 'rnabioinfo/seq_lm_samtools:v1.0.0', quality_control: 'rnabioinfo/seq_lm_quality_control:v1.0.0', oarfish: 'rnabioinfo/seq_lm_oarfish:v1.0.0']
}

def sample_checkpoint_key(sample: Sample) -> String {
    return "${sample.group}\t${sample.name}"
}

def sample_output_dir(output_root: Path, sample: Sample) -> Path {
    return output_root.resolve(safe_name(sample.group)).resolve(safe_name(sample.name))
}

def file_identity(path: Path) -> Map {
    def canonical: Path = path.toAbsolutePath().normalize()
    return [path: canonical.toString(), size: java.nio.file.Files.size(canonical), mtime_ms: java.nio.file.Files.getLastModifiedTime(canonical).toMillis()]
}

def sample_bam_inventory(sample: Sample) -> List<Map> {
    def bams: List<Path> = []
    java.nio.file.Files
        .walk(sample.bam_dir)
        .withCloseable { paths ->
            paths
                .filter { path: Path ->
                    java.nio.file.Files.isRegularFile(path) && path.fileName.toString().endsWith('.bam')
                }
                .forEach { path: Path -> bams.add(path) }
        }
    return bams
        .collect { bam: Path -> file_identity(bam) }
        .toSorted { left, right -> left.path <=> right.path }
}

def sha256_file(path: Path) -> String {
    def digest: java.security.MessageDigest = java.security.MessageDigest.getInstance('SHA-256')
    path
        .toFile()
        .withInputStream { input ->
            def stream: java.security.DigestInputStream = new java.security.DigestInputStream(input, digest)
            stream.transferTo(java.io.OutputStream.nullOutputStream())
        }
    return digest.digest().encodeHex().toString()
}

def next_analysis_snapshot_index(output_root: Path) -> Integer {
    def differential_root: Path = output_root.resolve('differential_expression')
    if (!java.nio.file.Files.isDirectory(differential_root)) {
        return 0
    }
    def indices: List<Integer> = []
    java.nio.file.Files
        .newDirectoryStream(differential_root)
        .withCloseable { entries ->
            entries.each { entry: Path ->
                def matcher = entry.fileName.toString() =~ /^batch_(\d+)$/
                if (matcher.matches() && java.nio.file.Files.isDirectory(entry)) {
                    indices.add(matcher.group(1).toInteger())
                }
            }
        }
    return indices.empty ? 0 : indices.max() + 1
}

def discover_sample_checkpoints(samples: List<Sample>, output_root: Path, genome: Path, annotation: Path) -> Map {
    def restored: List = []
    def active: List = []
    def safe_paths: Set<String> = [] as Set<String>

    samples.each { sample ->
        def safe_path: String = "${safe_name(sample.group)}/${safe_name(sample.name)}"
        if (!safe_paths.add(safe_path)) {
            error("Samples resolve to the same output path '${safe_path}'. Rename the colliding group or alias.")
        }

        def sample_dir: Path = sample_output_dir(output_root, sample)
        def marker: Path = sample_dir.resolve(sample_final_marker())
        if (!java.nio.file.Files.exists(marker)) {
            active.add(sample)
            return null
        }
        if (!java.nio.file.Files.isRegularFile(marker)) {
            checkpoint_error(sample, sample_dir, "${sample_final_marker()} is not a regular file")
        }

        def manifest: Map
        try {
            manifest = new groovy.json.JsonSlurper().parse(marker.toFile()) as Map
        }
        catch (exception: Exception) {
            checkpoint_error(sample, sample_dir, "cannot parse ${sample_final_marker()}: ${exception.message}")
        }

        validate_checkpoint_manifest(manifest, sample, sample_dir, genome, annotation)
        def qc_results: List<ChunkQCResult> = (manifest.qc as List).collect { qc: Map ->
            def nanoplot: Path = sample_dir.resolve(qc.nanoplot as String).normalize()
            def flagstat: Path = sample_dir.resolve(qc.flagstat as String).normalize()
            record(
                batch_index: (qc.batch_index as Number).intValue(),
                sample: sample,
                bam: sample.bam_dir,
                nanoplot_data: file(nanoplot),
                flagstat: file(flagstat),
            )
        }
        def counts: Path = sample_dir.resolve(manifest.quantification.path as String).normalize()
        restored.add(
            [sample: sample, quantification: record(batch_index: 0, sample: sample, counts: file(counts)), qc_results: qc_results]
        )
        log.info("Restoring finalized sample '${sample.group}/${sample.name}' from ${sample_dir}.")
    }
    return [restored: restored, active: active]
}

def validate_checkpoint_manifest(manifest: Map, sample: Sample, sample_dir: Path, genome: Path, annotation: Path) {
    def problems: List<String> = []
    if (manifest.schema_version != sample_checkpoint_schema()) {
        problems.add("unsupported schema_version '${manifest.schema_version}'")
    }
    if (manifest.sample?.name != sample.name || manifest.sample?.group != sample.group) {
        problems.add('sample identity does not match the sample sheet')
    }
    if (manifest.bams != sample_bam_inventory(sample)) {
        problems.add('BAM inventory changed after finalization')
    }
    if (manifest.reference?.genome != file_identity(genome)) {
        problems.add('reference genome changed')
    }
    if (manifest.reference?.annotation != file_identity(annotation)) {
        problems.add('reference annotation changed')
    }
    if (manifest.tools != sample_checkpoint_tools()) {
        problems.add('upstream tool versions changed')
    }

    def quant: Map = manifest.quantification as Map
    validate_checkpoint_artifact(sample_dir, quant, 'quantification', problems)
    def qc_entries: List = manifest.qc instanceof List ? manifest.qc as List : []
    if (qc_entries.empty) {
        problems.add('QC artifact inventory is empty')
    }
    else {
        qc_entries.eachWithIndex { qc: Map, index: Integer ->
            validate_checkpoint_artifact(
                sample_dir,
                [path: qc.nanoplot, sha256: qc.nanoplot_sha256],
                "QC NanoPlot artifact ${index}",
                problems,
            )
            validate_checkpoint_artifact(
                sample_dir,
                [path: qc.flagstat, sha256: qc.flagstat_sha256],
                "QC flagstat artifact ${index}",
                problems,
            )
        }
    }
    if (!problems.empty) {
        checkpoint_error(sample, sample_dir, problems.join('; '))
    }
}

def validate_checkpoint_artifact(sample_dir: Path, artifact: Map, label: String, problems: List<String>) {
    if (!artifact?.path || !artifact?.sha256) {
        problems.add("${label} metadata is incomplete")
        return null
    }
    def path: Path = sample_dir.resolve(artifact.path as String).normalize()
    if (!path.startsWith(sample_dir.normalize())) {
        problems.add("${label} path escapes the sample output directory")
    }
    else if (!java.nio.file.Files.isRegularFile(path)) {
        problems.add("${label} is missing")
    }
    else if (sha256_file(path) != artifact.sha256) {
        problems.add("${label} checksum does not match")
    }
}

def checkpoint_error(sample: Sample, sample_dir: Path, reason: String) {
    error(
        "Finalized sample checkpoint '${sample.group}/${sample.name}' is invalid: ${reason}. " + "Delete ${sample_dir} (or the complete output directory) and rerun to recompute it."
    )
}

def checkpoint_static_manifest(sample: Sample, genome: Path, annotation: Path) -> Map {
    return checkpoint_static_manifest(
        sample,
        [genome: file_identity(genome), annotation: file_identity(annotation)],
    )
}

def checkpoint_static_manifest(sample: Sample, reference_identity: Map) -> Map {
    return [schema_version: sample_checkpoint_schema(), sample: [name: sample.name, group: sample.group], bams: sample_bam_inventory(sample), reference: reference_identity, tools: sample_checkpoint_tools()]
}

/**
 * Materialize one durable, self-validating sample result bundle. The FINAL
 * manifest is created only after every artifact has been copied and hashed.
 */
process write_sample_checkpoints {
    label 'seq_lm_qc'
    container 'rnabioinfo/seq_lm_quality_control:v1.0.0'
    cpus 1

    publishDir params.out_dir, mode: 'copy', overwrite: true, saveAs: { path ->
        def value: String = "${path}"
        def prefix: String = 'checkpoint_output/'
        def offset: Integer = value.lastIndexOf(prefix)
        return offset >= 0 ? value.substring(offset + prefix.size()) : value
    }

    input:
    quantifications: List
    qc_results: List
    reference_identity: Map

    stage:
    stageAs quantifications*.counts, 'checkpoint_inputs/quant/input?.quant'
    stageAs qc_results*.nanoplot_data, 'checkpoint_inputs/nanoplot/input?.tsv.gz'
    stageAs qc_results*.flagstat, 'checkpoint_inputs/flagstat/input?.tsv'

    output:
    files('checkpoint_output/*/*', arity: '1..*')

    script:
    def indexed_quantifications: List = quantifications
        .toList()
        .withIndex()
        .collect { quantified_sample, index ->
            [value: quantified_sample, input_index: index + 1]
        }
    def indexed_qc_results: List = qc_results
        .toList()
        .withIndex()
        .collect { qc_result, index ->
            [value: qc_result, input_index: index + 1]
        }
    def quantifications_by_sample: Map<String,List> = indexed_quantifications.groupBy { indexed_quantification ->
        sample_checkpoint_key(indexed_quantification.value.sample)
    }
    def qc_by_sample: Map<String,List> = indexed_qc_results.groupBy { indexed_qc_result ->
        sample_checkpoint_key(indexed_qc_result.value.sample)
    }
    def commands: List<String> = []
    quantifications_by_sample.each { key: String, sample_quantifications: List ->
        def final_quantification = sample_quantifications.max { indexed_result ->
            indexed_result.value.batch_index
        }
        def sample_qc: List = qc_by_sample.containsKey(key) ? qc_by_sample[key] : []
        if (sample_qc.empty) {
            error("Cannot finalize sample '${key}': QC outputs are missing.")
        }
        def sample: Sample = final_quantification.value.sample
        def sample_dir: String = "checkpoint_output/${safe_name(sample.group)}/${safe_name(sample.name)}"
        def quant_index: Integer = final_quantification.input_index
        commands.add("mkdir -p ${sample_dir}/quantification ${sample_dir}/qc/nanoplot ${sample_dir}/qc/flagstat")
        def sorted_qc: List = sample_qc.toSorted { left, right ->
            left.value.batch_index <=> right.value.batch_index
        }
        def qc_manifest: List = sorted_qc.collect { indexed_qc ->
            def qc = indexed_qc.value
            def qc_index: Integer = indexed_qc.input_index
            def nanoplot_path: String = "qc/nanoplot/chunk_${qc.batch_index}.tsv.gz"
            def flagstat_path: String = "qc/flagstat/chunk_${qc.batch_index}.tsv"
            commands.add(
                "cp checkpoint_inputs/nanoplot/input${qc_index}.tsv.gz ${sample_dir}/${nanoplot_path}"
            )
            commands.add(
                "cp checkpoint_inputs/flagstat/input${qc_index}.tsv ${sample_dir}/${flagstat_path}"
            )
            [batch_index: qc.batch_index, nanoplot: nanoplot_path, flagstat: flagstat_path]
        }
        def manifest: Map = checkpoint_static_manifest(sample, reference_identity)
        manifest.quantification = [path: 'quantification/final.quant']
        manifest.qc = qc_manifest
        def manifest_json: String = groovy.json.JsonOutput.toJson(manifest)
        def quoted_manifest: String = "'" + manifest_json.replace("'", "'\"'\"'") + "'"
        commands.add(
            "cp checkpoint_inputs/quant/input${quant_index}.quant ${sample_dir}/quantification/final.quant"
        )
        commands.add("printf '%s\\n' ${quoted_manifest} > ${sample_dir}/manifest.in.json")
    }
    """
        ${commands.join('\n')}

        python3 - <<'PY'
import hashlib
import json
import pathlib

def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()

for manifest_path in pathlib.Path("checkpoint_output").glob("*/*/manifest.in.json"):
    sample_dir = manifest_path.parent
    manifest = json.loads(manifest_path.read_text())
    quantification = manifest["quantification"]
    quantification["sha256"] = digest(sample_dir / quantification["path"])
    for qc in manifest["qc"]:
        qc["nanoplot_sha256"] = digest(sample_dir / qc["nanoplot"])
        qc["flagstat_sha256"] = digest(sample_dir / qc["flagstat"])
    (sample_dir / "FINAL").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\\n"
    )
    manifest_path.unlink()
PY
        """
}
