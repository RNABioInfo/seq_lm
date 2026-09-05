nextflow.enable.types = true

/** Allocate ICA indices independently of differential-expression snapshots. */
def next_ica_snapshot_index(output_root: Path) -> Integer {
    def root: Path = output_root.resolve('ica')
    if (!java.nio.file.Files.isDirectory(root)) {
        return 0
    }
    def indices: List<Integer> = []
    java.nio.file.Files.newDirectoryStream(root).withCloseable { entries ->
        entries.each { entry: Path ->
            def match = entry.fileName.toString() =~ /^batch_(\d+)$/
            if (match.matches() && java.nio.file.Files.isDirectory(entry)) {
                indices.add(match.group(1).toInteger())
            }
        }
    }
    return indices.empty ? 0 : indices.max() + 1
}

/** Publish a complete immutable snapshot before atomically advancing its pointer. */
def publish_ica_snapshot(result, output_root: Path) -> Path {
    def root: Path = output_root.resolve('ica')
    def source_root: Path = result.results.toRealPath()
    java.nio.file.Files.createDirectories(root)
    def destination: Path = root.resolve("batch_${result.analysis_index}")
    if (java.nio.file.Files.exists(destination)) {
        error("ICA snapshot already exists: ${destination}; refusing to overwrite it.")
    }
    def staging: Path = java.nio.file.Files.createTempDirectory(root, '.pending-')
    java.nio.file.Files.walk(source_root).withCloseable { paths ->
        paths.forEach { source: Path ->
            def target: Path = staging.resolve(source_root.relativize(source).toString())
            if (java.nio.file.Files.isDirectory(source)) {
                java.nio.file.Files.createDirectories(target)
            }
            else {
                java.nio.file.Files.copy(source, target)
            }
        }
    }
    // Both paths are on the output filesystem. Directory rename exposes complete data.
    java.nio.file.Files.move(staging, destination, java.nio.file.StandardCopyOption.ATOMIC_MOVE)
    def status = new groovy.json.JsonSlurper().parse(destination.resolve('status.json').toFile())
    def pointer: Map = [schema_version: 1, analysis_index: result.analysis_index,
        batch_index: result.batch_index, report_sequence: result.report_sequence,
        path: destination.fileName.toString(), status: status.status]
    def temporary: Path = java.nio.file.Files.createTempFile(root, '.latest-', '.json')
    java.nio.file.Files.writeString(temporary, new groovy.json.JsonBuilder(pointer).toPrettyString() + '\n')
    java.nio.file.Files.move(temporary, root.resolve('latest.json'),
        java.nio.file.StandardCopyOption.ATOMIC_MOVE, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    return destination
}
