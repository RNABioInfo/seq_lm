import java.nio.file.Files
import java.nio.file.Path
import java.time.OffsetDateTime

/** Filesystem completion evidence, restricted to one dedicated MinKNOW run. */
class MinknowCompletion {
    static String directoryRunId(Path directory) {
        def match = directory.fileName.toString() =~ /^\d{8}_\d{4}_.+_([a-fA-F0-9]{8})$/
        return match.matches() ? match.group(1).toLowerCase() : null
    }

    static Path resolveRunDirectory(Path bamDirectory, Path override = null) {
        if (override != null) return override.toAbsolutePath().normalize()
        Path candidate = bamDirectory.toAbsolutePath().normalize()
        if (candidate.fileName.toString() ==~ /barcode\d+|unclassified/) {
            candidate = candidate.parent
            if (candidate == null || !(candidate.fileName.toString() in ['bam_pass', 'bam_fail'])) return null
        }
        if (candidate != null && candidate.fileName.toString() in ['bam_pass', 'bam_fail']) candidate = candidate.parent
        return candidate != null && directoryRunId(candidate) != null ? candidate : null
    }

    static Map inspect(Path directory, String expectedProtocolId = null) {
        if (directory == null) return [reason: 'unrecognized MinKNOW layout; use minknow_run_dir or manual STOP']
        try {
            List<Path> summaries = []
            Files.newDirectoryStream(directory).withCloseable { entries ->
                entries.each { Path entry ->
                    if (Files.isRegularFile(entry) && entry.fileName.toString() ==~ /final_summary(?:_.*)?\.txt/) summaries.add(entry)
                }
            }
            if (summaries.empty) return [reason: 'waiting for a final summary']
            if (summaries.size() != 1) return [reason: 'multiple final summaries; pooled/reused directories are unsupported']
            Path summary = summaries[0]
            if (Files.size(summary) == 0) return [reason: 'empty final summary']
            // Read once: content and file metadata both participate in settling.
            String content = Files.readString(summary)
            Map<String, String> fields = [:]
            for (String line : content.readLines()) {
                if (!line.trim()) continue
                int split = line.indexOf('=')
                if (split < 1) return [reason: 'malformed final summary']
                String key = line.substring(0, split).trim()
                if (!key) return [reason: 'malformed final summary key']
                if (fields.containsKey(key)) return [reason: 'duplicate final summary field']
                fields[key] = line.substring(split + 1).trim()
            }
            for (String key : ['protocol_run_id', 'acquisition_run_id', 'acquisition_stopped', 'processing_stopped']) {
                if (!fields[key]) return [reason: "incomplete final summary: missing ${key}"]
            }
            OffsetDateTime acquisitionEnd = OffsetDateTime.parse(fields.acquisition_stopped)
            OffsetDateTime processingEnd = OffsetDateTime.parse(fields.processing_stopped)
            if (processingEnd.isBefore(acquisitionEnd)) return [reason: 'processing stopped before acquisition']
            if (fields.started && acquisitionEnd.isBefore(OffsetDateTime.parse(fields.started))) return [reason: 'acquisition stopped before it started']
            if (expectedProtocolId && fields.protocol_run_id != expectedProtocolId) return [reason: 'final summary protocol_run_id does not match sample metadata']
            String shortId = directoryRunId(directory)
            if (shortId && !fields.protocol_run_id.replace('-', '').toLowerCase().startsWith(shortId)) return [reason: 'final summary protocol_run_id does not match run directory']
            return [path: summary, signature: [summary.toString(), Files.size(summary), Files.getLastModifiedTime(summary).toMillis(), content], reason: null]
        }
        catch (Exception exception) {
            return [reason: "cannot validate final summary: ${exception.message}"]
        }
    }

    /** Independent of wall-clock changes; reset on any summary or BAM change. */
    static boolean settled(Map state, Map evidence, Map inventory, boolean pendingBams,
                           int stablePolls, long settleMillis, long nowNanos) {
        if (evidence.reason != null) {
            state.clear()
            return false
        }
        def signature = [evidence.signature, inventory]
        if (state.signature != signature) {
            state.signature = signature
            state.since = nowNanos
            state.polls = 1
            return false
        }
        state.polls = (state.polls as int) + 1
        return !pendingBams && state.polls >= stablePolls && nowNanos - (state.since as long) >= settleMillis * 1_000_000L
    }
}
