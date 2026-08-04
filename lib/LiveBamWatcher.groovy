import java.nio.file.ClosedWatchServiceException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardWatchEventKinds
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService

import nextflow.Global
import nextflow.dataflow.ChannelImpl
import nextflow.extension.CH

/**
 * Closeable per-sample file watcher used by live BAM ingress.
 *
 * Nextflow's watchPath source is intentionally open-ended. This helper owns the
 * underlying WatchService so the sample stream can bind Channel.STOP as soon as
 * that sample's STOP file is observed.
 */
class LiveBamWatcher {

    static ChannelImpl watch(Path root, Closure onActive, Closure onError) {
        def target = CH.create()
        Path watchRoot = root.toAbsolutePath().normalize()

        Global.session.addIgniter {
            Thread watcherThread = Thread.startDaemon("bam-ingress-${watchRoot.fileName}") {
                runWatcher(watchRoot, target, onActive, onError)
            }
            Global.onCleanup({ watcherThread.interrupt() })
        }

        return new ChannelImpl(target)
    }

    private static void runWatcher(Path root, def target, Closure onActive, Closure onError) {
        WatchService watchService = null

        try {
            watchService = root.fileSystem.newWatchService()
            Map<WatchKey, Path> watchKeys = [:]
            registerWatchDirs(root, watchService, watchKeys)
            onActive.call(watchKeys.size())

            List<Path> existingStops = findStopFiles(root)
            if (!existingStops.isEmpty()) {
                target.bind([source: "existing", path: existingStops[0]])
                return
            }

            while (!Thread.currentThread().isInterrupted()) {
                WatchKey watchKey = watchService.take()
                Path eventDir = watchKeys[watchKey]
                if (eventDir == null) {
                    watchKey.reset()
                    continue
                }

                boolean shouldStop = false
                for (WatchEvent<?> watchEvent : watchKey.pollEvents()) {
                    if (watchEvent.kind() == StandardWatchEventKinds.OVERFLOW) {
                        continue
                    }

                    String source = eventSource(watchEvent.kind())
                    Path eventPath = eventDir.resolve((Path) watchEvent.context()).toAbsolutePath().normalize()

                    if (source == "create" && Files.isDirectory(eventPath)) {
                        registerWatchDirs(eventPath, watchService, watchKeys)
                    }

                    target.bind([source: source, path: eventPath])

                    if (isStopPath(eventPath)) {
                        shouldStop = true
                        break
                    }
                }

                if (shouldStop) {
                    return
                }

                if (!watchKey.reset()) {
                    watchKeys.remove(watchKey)
                    if (watchKeys.isEmpty()) {
                        return
                    }
                }
            }
        }
        catch (InterruptedException ignored) {
            Thread.currentThread().interrupt()
        }
        catch (ClosedWatchServiceException ignored) {
            // The watch service is closed during normal cleanup.
        }
        catch (Throwable error) {
            onError.call(error)
        }
        finally {
            try {
                if (watchService != null) {
                    watchService.close()
                }
            }
            catch (Throwable ignored) {
                // Ignore cleanup failures after the sample stream has closed.
            }
            target.bind(CH.stop())
        }
    }

    private static void registerWatchDirs(Path root, WatchService watchService, Map<WatchKey, Path> watchKeys) {
        Files.walk(root).withCloseable { stream ->
            def iterator = stream.iterator()
            while (iterator.hasNext()) {
                Path dir = iterator.next()
                if (Files.isDirectory(dir)) {
                    registerWatchDir(dir.toAbsolutePath().normalize(), watchService, watchKeys)
                }
            }
        }
    }

    private static void registerWatchDir(Path dir, WatchService watchService, Map<WatchKey, Path> watchKeys) {
        if (watchKeys.values().contains(dir)) {
            return
        }

        WatchKey watchKey = dir.register(
            watchService,
            StandardWatchEventKinds.ENTRY_CREATE,
            StandardWatchEventKinds.ENTRY_MODIFY
        )
        watchKeys[watchKey] = dir
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

    private static String eventSource(WatchEvent.Kind<?> kind) {
        if (kind == StandardWatchEventKinds.ENTRY_CREATE) {
            return "create"
        }
        if (kind == StandardWatchEventKinds.ENTRY_MODIFY) {
            return "modify"
        }
        return kind.name().replace("ENTRY_", "").toLowerCase()
    }

    private static boolean isStopPath(Path path) {
        return path.fileName.toString().startsWith("STOP")
    }
}
