nextflow.enable.types = true

include { GeneContrast ; Sample ; QuantifiedSampleBatch } from './sample.nf'

/** Resolve a headerless CSV once, or construct the default control contrasts. */
def resolve_gene_contrasts(samples: List<Sample>, contrast_file: Path?) -> List<GeneContrast> {
    def groups: Set<String> = samples*.group.toSet()
    def pairs: List<List<String>> = []
    if (contrast_file != null) {
        if (!contrast_file.exists() || !contrast_file.isFile()) {
            error("--gene_contrast_file must be an existing regular file: ${contrast_file}")
        }
        // One contrast per physical line; the CSV parser otherwise skips blank lines.
        contrast_file
            .readLines()
            .each { line: String ->
                if (!line.matches(/\s*(?:"(?:[^"]|"")*"|[^",\r\n]*)\s*,\s*(?:"(?:[^"]|"")*"|[^",\r\n]*)\s*/)) {
                    error('Malformed gene contrast CSV: each line must contain exactly two CSV fields.')
                }
            }
        pairs = parse_delimited_rows(contrast_file.text, ',', false, true)
            .collect { row ->
                if (row.size() != 2 || row.any { field -> !field || field.trim() == '' }) {
                    error("Each gene contrast CSV row must contain exactly two nonempty group names: ${row}")
                }
                return row.collect { field -> field.trim() }
            }
    }
    else {
        if (samples.count { sample -> sample.group == 'control' } < 2) {
            error('Sample sheet requires at least two control samples in the group column.')
        }
        pairs = groups
            .findAll { group: String -> group != 'control' }
            .toSorted()
            .collect { group: String -> [group, 'control'] }
    }
    if (pairs.empty) {
        error('At least one gene contrast is required.')
    }
    def seen: Set<List<String>> = [] as Set<List<String>>
    return pairs
        .withIndex()
        .collect { pair: List<String>, index: Integer ->
            if (pair.any { group: String -> !groups.contains(group) }) {
                error("Gene contrast references an unknown sample-sheet group: ${pair}")
            }
            if (pair[0] == pair[1] || !seen.add(pair)) {
                error("Self-comparisons and duplicate ordered gene contrasts are invalid: ${pair}")
            }
            if (pair.any { group: String -> group.contains('\t') || group.contains('\n') || group.contains('\r') }) {
                error('Gene contrast group names cannot contain tabs or newlines.')
            }
            return record(contrast_id: "contrast_${index + 1}", target_group: pair[0], reference_group: pair[1])
        }
}

/** Preserve standard doubled-quote escapes while using Nextflow's CSV parser. */
def parse_delimited_rows(text: String, separator: String, header: Boolean, strip: Boolean) -> List {
    // Nextflow treats the first closing quote as the field terminator. Transcode
    // quote delimiters to an unused character so literal double quotes survive.
    def quote_code: Integer? = (57344..63743).find { code: Integer -> !text.contains(Character.toString(code)) }
    if (quote_code == null) {
        error('CSV contains every private-use quote delimiter; cannot parse the file.')
    }
    def quote: String = Character.toString(quote_code)
    def normalized: String = text.replaceAll(/"((?:[^"]|"")*)"/) { _match: String, content: String ->
        quote + content.replace('""', '"') + quote
    }
    if (strip) {
        normalized = normalized.replaceAll(" *${quote}([^${quote}]*)${quote} *") { _match: String, content: String ->
            quote + content + quote
        }
    }
    def parser: nextflow.util.CsvParser = new nextflow.util.CsvParser().setQuote(quote).setSeparator(separator).setStrip(strip)
    def rows: List<List<String>> = normalized.split(/\r?\n/).findAll { line: String -> line.trim() != '' }.collect { line: String -> parser.parse(line) }
    if (!header || rows.empty) {
        return rows
    }
    def columns: List<String> = rows[0]
    return rows.drop(1).collect { row: List<String> ->
        columns.withIndex().collectEntries { column: String, index: Integer -> [(column): index < row.size() ? row[index] ?: '' : ''] }
    }
}

def gene_contrast_groups(contrasts: List<GeneContrast>) -> Set<String> {
    return contrasts.collectMany { contrast -> [contrast.target_group, contrast.reference_group] }.toSet()
}

/** ICA conditions are independent of gene contrasts; always include its control. */
def resolve_ica_groups(samples: List<Sample>, control_group: String, condition_groups: String?) -> Set<String> {
    def groups: Set<String> = samples*.group.toSet()
    def control: String = control_group.trim()
    if (!control || samples.count { sample -> sample.group == control } < 2) {
        error("ICA requires at least two samples in --ica_control_group '${control}'.")
    }
    if (condition_groups == null || condition_groups.trim() == '') {
        return groups
    }
    def conditions: List<String> = condition_groups.split(',', -1).collect { group: String -> group.trim() }
    if (conditions.any { group: String -> !group || !groups.contains(group) } || conditions.toSet().size() != conditions.size()) {
        error('--ica_condition_groups must contain distinct, nonempty sample-sheet group names.')
    }
    return (conditions + [control]).toSet()
}

/** Retain batch identity and select only the analysis cohort. */
def select_analysis_batch(batch, groups: Set<String>) -> QuantifiedSampleBatch {
    return record(
        batch_index: batch.batch_index,
        report_sequence: batch.report_sequence,
        samples: batch.samples.findAll { sample -> groups.contains(sample.sample.group) },
    )
}

def gene_contrast_tsv(contrasts: List<GeneContrast>) -> String {
    return 'contrast_id\ttarget_group\treference_group\n' + contrasts.collect { contrast ->
        [contrast.contrast_id, contrast.target_group, contrast.reference_group].collect { field: String -> tsv_field(field) }.join('\t')
    }.join('\n') + '\n'
}

/** Quote embedded double quotes for the canonical TSV readers. */
def tsv_field(value: String) -> String {
    def text: String = "${value}".replaceAll(/[\t\r\n]+/, ' ')
    return text.contains('"') ? '"' + text.replace('"', '""') + '"' : text
}
