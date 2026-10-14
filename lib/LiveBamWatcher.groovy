import java.nio.file.Files
import java.nio.file.Path

import nextflow.Global
import nextflow.dataflow.ChannelImpl
import nextflow.extension.CH

/**
 * Closeable per-sample polling source used by live BAM ingress.
 *
 * Windows applications do not reliably produce Linux file notification events
 * for files below WSL DrvFs mounts such as /mnt/c. Polling also lets ingress
 * wait for a BAM to stop changing before downstream tools open it.
 */
class LiveBamWatcher {

    private static final log = org.slf4j.LoggerFactory.getLogger(LiveBamWatcher)

    static ChannelImpl poll(
        Path root,
        Collection<Path> initialBams,
        long pollIntervalMillis,
        int requiredStablePolls,
        Map completion,
        Closure onActive,
        Closure onError
    ) {
        def target = CH.create()
        Path pollRoot = root.toAbsolutePath().normalize()
        Set<String> initiallyAccepted = initialBams.collect { Path path -> pathKey(path) } as Set<String>

        Global.session.addIgniter {
            Thread pollerThread = Thread.startDaemon("bam-ingress-poller-${pollRoot.fileName}") {
                runPoller(
                    pollRoot,
                    initiallyAccepted,
                    pollIntervalMillis,
                    requiredStablePolls,
                    completion,
                    target,
                    onActive,
                    onError
                )
            }
            Global.onCleanup({ pollerThread.interrupt() })
        }

        return new ChannelImpl(target)
    }

    private static void runPoller(
        Path root,
        Set<String> initiallyAccepted,
        long pollIntervalMillis,
        int requiredStablePolls,
        Map completion,
        def target,
        Closure onActive,
        Closure onError
    ) {
        try {
            Set<String> acceptedBams = new LinkedHashSet<>(initiallyAccepted)
            Map<String, Map<String, Long>> observations = [:]
            Path runDirectory = completion.enabled ? MinknowCompletion.resolveRunDirectory(root, completion.run_dir as Path) : null
            Map completionState = [:]
            String lastReason = null
            onActive.call()

            while (!Thread.currentThread().isInterrupted()) {
                List<Path> stopFiles = findStopFiles(root)
                List<Path> currentBams = findBamFiles(root)
                Set<String> currentKeys = currentBams.collect { Path path -> pathKey(path) } as Set<String>
                observations.keySet().removeIf { String key -> !currentKeys.contains(key) }

                Map<String, List<Long>> inventory = [:]
                currentBams.each { Path bamPath ->
                    String key = pathKey(bamPath)
                    // A producer may remove a file between discovery and stat.
                    def attributes
                    try {
                        attributes = Files.readAttributes(bamPath, java.nio.file.attribute.BasicFileAttributes)
                    }
                    catch (java.nio.file.NoSuchFileException ignored) {
                        observations.remove(key)
                        return
                    }
                    long size = attributes.size()
                    long modified = attributes.lastModifiedTime().toMillis()
                    inventory[key] = [size, modified]
                    if (acceptedBams.contains(key)) {
                        return
                    }
                    Map<String, Long> previous = observations[key]
                    long stablePolls = previous != null && previous.size == size && previous.modified == modified
                        ? previous.stablePolls + 1L
                        : 1L
                    observations[key] = [size: size, modified: modified, stablePolls: stablePolls]

                    if (stablePolls >= requiredStablePolls) {
                        acceptedBams.add(key)
                        observations.remove(key)
                        target.bind([source: "poll", path: bamPath])
                    }
                }

                boolean hasPendingBams = currentKeys.any { String key -> !acceptedBams.contains(key) }
                if (!stopFiles.isEmpty() && !hasPendingBams) {
                    target.bind([source: "poll", path: stopFiles[0]])
                    return
                }

                if (completion.enabled) {
                    Map evidence = MinknowCompletion.inspect(runDirectory, completion.protocol_run_id as String)
                    if (evidence.reason != null && evidence.reason != lastReason) {
                        if (evidence.reason == 'waiting for a final summary') log.info("MinKNOW completion for ${root}: ${evidence.reason}")
                        else log.warn("MinKNOW completion for ${root}: ${evidence.reason}; continuing to watch BAMs")
                    }
                    lastReason = evidence.reason
                    if (MinknowCompletion.settled(completionState, evidence, inventory, hasPendingBams,
                            requiredStablePolls, (completion.settle_seconds as long) * 1000L, System.nanoTime())) {
                        log.info("MinKNOW output completed for ${root}: ${evidence.path}; BAM output settled")
                        target.bind([source: 'minknow-final-summary', path: evidence.path])
                        return
                    }
                }

                Thread.sleep(pollIntervalMillis)
            }
        }
        catch (InterruptedException ignored) {
            Thread.currentThread().interrupt()
        }
        catch (Throwable error) {
            onError.call(error)
        }
        finally {
            target.bind(CH.stop())
        }
    }

    private static List<Path> findBamFiles(Path root) {
        List<Path> bamFiles = []
        Files.walk(root).withCloseable { stream ->
            def iterator = stream.iterator()
            while (iterator.hasNext()) {
                Path path = iterator.next()
                if (Files.isRegularFile(path) && path.fileName.toString().endsWith('.bam')) {
                    bamFiles.add(path.toAbsolutePath().normalize())
                }
            }
        }
        return bamFiles.sort { left, right -> left.toString() <=> right.toString() }
    }

    private static List<Path> findStopFiles(Path root) {
        List<Path> stopFiles = []
        Files.walk(root).withCloseable { stream ->
            def iterator = stream.iterator()
            while (iterator.hasNext()) {
                Path path = iterator.next()
                if (Files.isRegularFile(path) && isStopPath(path)) {
                    stopFiles.add(path.toAbsolutePath().normalize())
                }
            }
        }
        return stopFiles.sort { left, right -> left.toString() <=> right.toString() }
    }

    private static String pathKey(Path path) {
        return path.toAbsolutePath().normalize().toString()
    }

    private static boolean isStopPath(Path path) {
        return path.fileName.toString().startsWith("STOP")
    }
}
