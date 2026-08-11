#!/usr/bin/env Rscript

library(optparse)
library(edgeR)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
    stop("Unable to determine the edgeR analysis script directory.")
}
script_path = sub("^--file=", "", script_arg[[1]])
script_dir = dirname(normalizePath(script_path))

source(file.path(script_dir, "read_gmt.R"))

options = list(
    make_option(
        c("--quant_manifest"),
        type = "character"
    ),
    make_option(
        c("--annotation"),
        type = "character"
    ),
    make_option(
        c("--output_dir"),
        type = "character"
    ),
    make_option(
        c("--gene_sets"),
        type = "character",
        help = "Gene sets in Gene Matrix Transposed (GMT) format."
    ),
    make_option(
        c("--control_group"),
        type = "character",
        default = "control",
        help = "Reference group used for pairwise contrasts [default: %default]."
    ),
    make_option(
        c("--lfc"),
        type = "double",
        default = 1,
        help = "Minimum absolute log2 fold change for glmTreat [default: %default]."
    ),
    make_option(
        c("--quant_id_column"),
        type = "character",
        default = NULL,
        help = paste0(
            "Feature-ID column in quantification files. Common names are ",
            "auto-detected when omitted."
        )
    ),
    make_option(
        c("--quant_count_column"),
        type = "character",
        default = NULL,
        help = paste0(
            "Count column in quantification files. Common names are ",
            "auto-detected when omitted."
        )
    ),
    make_option(
        c("--annotation_id_attributes"),
        type = "character",
        default = NULL,
        help = paste0(
            "Comma-separated GTF/GFF attribute names usable as identifiers. ",
            "Identifier-like attributes are auto-detected when omitted."
        )
    ),
    make_option(
        c("--strip_id_versions"),
        action = "store_true",
        default = FALSE,
        help = paste0(
            "Also match identifiers after removing a terminal numeric version, ",
            "for example ENSG000001.2 -> ENSG000001."
        )
    ),
    make_option(
        c("--plot_top_n"),
        type = "integer",
        default = 30L,
        help = "Maximum number of gene sets in each fry plot [default: %default]."
    )
)

parser = OptionParser(
    description = "edgeR differential expression and fry gene-set analysis",
    option_list = options
)

args = parse_args(parser)

if (is.null(args$quant_manifest) || is.null(args$output_dir)) {
    stop("--quant_manifest and --output_dir are required.")
}
if (args$plot_top_n < 1L) {
    stop("--plot_top_n must be at least 1.")
}

output_dir = args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

quant_manifest_path = args$quant_manifest
quant_manifest = read.delim(quant_manifest_path, sep = "\t")

required_manifest_columns = c("name", "group", "count_file")
missing_manifest_columns = setdiff(required_manifest_columns, colnames(quant_manifest))
if (length(missing_manifest_columns) > 0) {
    stop(
        paste0(
            "Quantification manifest is missing required column(s): ",
            paste(missing_manifest_columns, collapse = ", ")
        )
    )
}
if (anyDuplicated(quant_manifest$name)) {
    stop("Sample names in --quant_manifest must be unique.")
}
if (!(args$control_group %in% quant_manifest$group)) {
    stop(sprintf(
        "Control group '%s' is absent from --quant_manifest.",
        args$control_group
    ))
}

choose_column = function(
    available,
    requested,
    candidates,
    role,
    file
) {
    if (!is.null(requested)) {
        if (!(requested %in% available)) {
            stop(sprintf(
                "Requested %s column '%s' is absent from %s. Available columns: %s",
                role,
                requested,
                file,
                paste(available, collapse = ", ")
            ))
        }
        return(requested)
    }

    available_lower = tolower(available)
    for (candidate in candidates) {
        position = match(tolower(candidate), available_lower)
        if (!is.na(position)) {
            return(available[[position]])
        }
    }
    stop(sprintf(
        paste0(
            "Could not auto-detect the %s column in %s. Available columns: %s. ",
            "Use the corresponding command-line override."
        ),
        role,
        file,
        paste(available, collapse = ", ")
    ))
}

read_quantification = function(
    file,
    sample_name,
    requested_id_column,
    requested_count_column
) {
    if (!file.exists(file)) {
        stop(sprintf("Quantification file does not exist: %s", file))
    }
    first_line = readLines(file, n = 1L, warn = FALSE)
    delimiter = if (grepl("\t", first_line, fixed = TRUE)) {
        "\t"
    } else if (grepl(",", first_line, fixed = TRUE)) {
        ","
    } else {
        stop(sprintf(
            "Could not detect a tab or comma delimiter in %s.",
            file
        ))
    }
    quantification = read_delim(
        file,
        delim = delimiter,
        show_col_types = FALSE,
        name_repair = "minimal",
        progress = FALSE
    )
    id_column = choose_column(
        colnames(quantification),
        requested_id_column,
        c(
            "tname",
            "target_id",
            "transcript_id",
            "gene_id",
            "feature_id",
            "name",
            "id"
        ),
        "feature-ID",
        file
    )
    count_column = choose_column(
        colnames(quantification),
        requested_count_column,
        c(
            "num_reads",
            "expected_count",
            "est_counts",
            "count",
            "counts",
            "read_count"
        ),
        "count",
        file
    )

    feature_ids = as.character(quantification[[id_column]])
    raw_counts = quantification[[count_column]]
    numeric_counts = suppressWarnings(as.numeric(raw_counts))
    invalid_counts = is.na(numeric_counts) & !is.na(raw_counts)
    if (any(invalid_counts)) {
        stop(sprintf(
            "Count column '%s' in %s contains non-numeric values.",
            count_column,
            file
        ))
    }
    if (any(is.na(feature_ids)) || any(!nzchar(feature_ids))) {
        stop(sprintf(
            "Feature-ID column '%s' in %s contains missing or empty values.",
            id_column,
            file
        ))
    }
    if (any(numeric_counts < 0, na.rm = TRUE)) {
        stop(sprintf(
            "Count column '%s' in %s contains negative values.",
            count_column,
            file
        ))
    }

    message(sprintf(
        "Reading %s: feature IDs from '%s', counts from '%s'.",
        basename(file),
        id_column,
        count_column
    ))
    tibble(
        feature_id = feature_ids,
        sample = sample_name,
        count = replace_na(numeric_counts, 0)
    ) %>%
        group_by(feature_id, sample) %>%
        summarise(count = sum(count), .groups = "drop")
}

resolve_manifest_paths = function(paths, manifest_path) {
    manifest_dir = dirname(normalizePath(manifest_path))
    vapply(
        paths,
        function(path) {
            if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
                path
            } else {
                file.path(manifest_dir, path)
            }
        },
        character(1L),
        USE.NAMES = FALSE
    )
}

collect_counts = function(
    quant_manifest,
    quant_manifest_path,
    requested_id_column,
    requested_count_column
) {
    count_files = resolve_manifest_paths(
        quant_manifest$count_file,
        quant_manifest_path
    )
    map2_dfr(
        count_files,
        quant_manifest$name,
        \(file, sample_name) {
            read_quantification(
                file,
                sample_name,
                requested_id_column,
                requested_count_column
            )
        }
    ) %>%
    pivot_wider(
        names_from = sample,
        values_from = count,
        values_fill = 0
    ) %>%
    select(feature_id, all_of(quant_manifest$name)) %>%
    column_to_rownames("feature_id") %>%
    arrange(row.names(.))
}

extract_annotation_attribute = function(attributes, attribute) {
    escaped_attribute = paste0("\\Q", attribute, "\\E")
    extract_matches = function(pattern) {
        matches = regexec(pattern, attributes, perl = TRUE)
        values = regmatches(attributes, matches)
        vapply(
            values,
            function(value) {
                if (length(value) >= 2L) value[[2L]] else NA_character_
            },
            character(1L)
        )
    }

    gtf_values = extract_matches(paste0(
        "(?:^|;[[:space:]]*)",
        escaped_attribute,
        '[[:space:]]+"([^"]*)"'
    ))
    gff_values = extract_matches(paste0(
        "(?:^|;[[:space:]]*)",
        escaped_attribute,
        "=([^;]*)"
    ))
    values = ifelse(is.na(gtf_values), gff_values, gtf_values)
    present = which(!is.na(values) & nzchar(values))
    if (length(present) == 0L) {
        return(tibble())
    }

    split_values = lapply(
        utils::URLdecode(values[present]),
        function(value) trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
    )
    tibble(
        record_id = rep(present, lengths(split_values)),
        identifier_type = attribute,
        identifier = unlist(split_values, use.names = FALSE)
    ) %>%
        filter(nzchar(identifier))
}

identifier_aliases = function(identifier, strip_versions = FALSE) {
    decoded = utils::URLdecode(identifier)
    aliases = c(identifier, decoded)
    aliases = c(
        aliases,
        sub(
            "^(?:gene|transcript|rna|mrna|protein|cds|feature):",
            "",
            aliases,
            ignore.case = TRUE,
            perl = TRUE
        )
    )
    if (strip_versions) {
        aliases = c(aliases, sub("\\.[0-9]+$", "", aliases, perl = TRUE))
    }
    unique(aliases[nzchar(aliases)])
}

add_identifier_aliases = function(table, strip_versions = FALSE) {
    table %>%
        mutate(
            match_key = map(
                identifier,
                identifier_aliases,
                strip_versions = strip_versions
            )
        ) %>%
        unnest_longer(match_key) %>%
        distinct()
}

read_annotation = function(
    annotation_path,
    requested_attributes = NULL,
    strip_versions = FALSE
) {
    if (is.null(annotation_path)) {
        return(tibble(
            record_id = integer(),
            identifier_type = character(),
            identifier = character(),
            match_key = character()
        ))
    }
    if (!file.exists(annotation_path)) {
        stop(sprintf("Annotation file does not exist: %s", annotation_path))
    }

    annotation = read_tsv(
        annotation_path,
        comment = "#",
        col_names = FALSE,
        quote = "",
        show_col_types = FALSE,
        progress = FALSE
    )
    if (ncol(annotation) < 9L) {
        stop("Annotation must be a nine-column GTF or GFF3 file.")
    }
    colnames(annotation)[seq_len(9L)] = c(
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes"
    )

    if (!is.null(requested_attributes)) {
        selected_attributes = trimws(strsplit(
            requested_attributes,
            ",",
            fixed = TRUE
        )[[1L]])
    } else {
        selected_attributes = c(
            "gene_id",
            "transcript_id",
            "gene",
            "gene_name",
            "locus_tag",
            "old_locus_tag",
            "protein_id",
            "feature_id",
            "ID",
            "Name",
            "Alias",
            "Parent",
            "geneID",
            "transcriptID",
            "proteinID"
        )
    }
    annotation_long = map_dfr(
        selected_attributes,
        function(attribute) {
            extract_annotation_attribute(annotation$attributes, attribute)
        }
    )
    if (!is.null(requested_attributes)) {
        missing_attributes = setdiff(
            selected_attributes,
            unique(annotation_long$identifier_type)
        )
        if (length(missing_attributes) > 0L) {
            stop(sprintf(
                "Requested annotation identifier attribute(s) not found: %s",
                paste(missing_attributes, collapse = ", ")
            ))
        }
    }
    if (nrow(annotation_long) == 0L) {
        stop(
            paste0(
                "No identifier-like attributes were found in the annotation. ",
                "Use --annotation_id_attributes to select custom attributes."
            )
        )
    }

    add_identifier_aliases(annotation_long, strip_versions)
}

build_identifier_map = function(
    annotation,
    count_ids,
    strip_versions = FALSE
) {
    count_aliases = tibble(
        identifier = count_ids,
        feature_id = count_ids
    ) %>%
        add_identifier_aliases(strip_versions) %>%
        select(feature_id, match_key)

    direct_map = count_aliases %>%
        transmute(
            identifier = feature_id,
            identifier_type = "count_id",
            match_key,
            feature_id,
            mapping_source = "direct"
        )
    if (nrow(annotation) == 0L) {
        return(direct_map)
    }

    anchor_rows = annotation %>%
        inner_join(
            count_aliases,
            by = "match_key",
            relationship = "many-to-many"
        ) %>%
        distinct(record_id, feature_id)
    if (nrow(anchor_rows) == 0L) {
        return(direct_map)
    }

    strong_identifier = grepl(
        "(^|_)(id|parent|locus_tag)$",
        tolower(annotation$identifier_type),
        perl = TRUE
    ) |
        tolower(annotation$identifier_type) %in% c(
            "id",
            "parent",
            "locus_tag",
            "old_locus_tag"
        )
    strong_links = annotation[strong_identifier, ] %>%
        select(record_id, match_key) %>%
        distinct()
    known = strong_links %>%
        inner_join(
            anchor_rows,
            by = "record_id",
            relationship = "many-to-many"
        ) %>%
        select(match_key, feature_id) %>%
        distinct()

    for (iteration in seq_len(10L)) {
        reached_rows = strong_links %>%
            inner_join(
                known,
                by = "match_key",
                relationship = "many-to-many"
            ) %>%
            select(record_id, feature_id) %>%
            distinct()
        expanded = strong_links %>%
            inner_join(
                reached_rows,
                by = "record_id",
                relationship = "many-to-many"
            ) %>%
            select(match_key, feature_id) %>%
            distinct()
        updated = bind_rows(known, expanded) %>% distinct()
        if (nrow(updated) == nrow(known)) {
            break
        }
        known = updated
    }

    reached_rows = bind_rows(
        anchor_rows,
        strong_links %>%
            inner_join(
                known,
                by = "match_key",
                relationship = "many-to-many"
            ) %>%
            select(record_id, feature_id)
    ) %>%
        distinct()
    annotation_map = annotation %>%
        inner_join(
            reached_rows,
            by = "record_id",
            relationship = "many-to-many"
        ) %>%
        transmute(
            identifier,
            identifier_type,
            match_key,
            feature_id,
            mapping_source = "annotation"
        ) %>%
        distinct()

    bind_rows(direct_map, annotation_map) %>% distinct()
}

resolve_gene_sets = function(
    gmt,
    identifier_map,
    strip_versions = FALSE
) {
    gmt_members = gmt_to_long(gmt) %>%
        transmute(
            gene_set = set_name,
            description,
            gmt_member = gene
        ) %>%
        distinct()
    member_aliases = gmt_members %>%
        filter(!is.na(gmt_member)) %>%
        distinct(identifier = gmt_member) %>%
        add_identifier_aliases(strip_versions) %>%
        select(gmt_member = identifier, match_key)

    matched = gmt_members %>%
        select(gene_set, gmt_member) %>%
        inner_join(
            member_aliases,
            by = "gmt_member",
            relationship = "many-to-many"
        ) %>%
        inner_join(
            identifier_map,
            by = "match_key",
            relationship = "many-to-many"
        ) %>%
        group_by(gene_set, gmt_member, feature_id) %>%
        summarise(
            match_method = if_else(
                any(mapping_source == "direct"),
                "direct",
                paste(sort(unique(identifier_type)), collapse = ",")
            ),
            .groups = "drop"
        )
    resolution = gmt_members %>%
        left_join(matched, by = c("gene_set", "gmt_member"))
    resolved_sets = lapply(names(gmt), function(set_name) {
        resolution %>%
            filter(gene_set == set_name, !is.na(feature_id)) %>%
            pull(feature_id) %>%
            unique()
    })
    names(resolved_sets) = names(gmt)

    list(sets = resolved_sets, table = resolution)
}

make_gene_set_coverage = function(resolution, retained_ids) {
    resolution %>%
        group_by(gene_set, description) %>%
        summarise(
            gmt_members = n_distinct(gmt_member, na.rm = TRUE),
            matched_gmt_members = n_distinct(
                gmt_member[!is.na(feature_id)]
            ),
            count_matrix_members = n_distinct(
                feature_id[!is.na(feature_id)]
            ),
            tested_members = n_distinct(
                feature_id[
                    !is.na(feature_id) & feature_id %in% retained_ids
                ]
            ),
            tested_gmt_members = n_distinct(
                gmt_member[
                    !is.na(feature_id) & feature_id %in% retained_ids
                ]
            ),
            .groups = "drop"
        ) %>%
        mutate(
            count_matrix_coverage = if_else(
                gmt_members > 0,
                matched_gmt_members / gmt_members,
                NA_real_
            ),
            tested_coverage = if_else(
                gmt_members > 0,
                tested_gmt_members / gmt_members,
                NA_real_
            )
        )
}

make_count_annotation = function(identifier_map) {
    identifier_map %>%
        filter(mapping_source == "annotation") %>%
        select(feature_id, identifier_type, identifier) %>%
        distinct() %>%
        group_by(feature_id, identifier_type) %>%
        summarise(
            identifier = paste(sort(unique(identifier)), collapse = "|"),
            .groups = "drop"
        ) %>%
        pivot_wider(
            names_from = identifier_type,
            values_from = identifier
        )
}

safe_name = function(value) {
    gsub("[^A-Za-z0-9._-]", "_", value)
}

write_fry_plot = function(fry_table, output_path, title, top_n) {
    if (nrow(fry_table) == 0L) {
        return(invisible(NULL))
    }

    plot_table = fry_table %>%
        arrange(FDR) %>%
        slice_head(n = top_n) %>%
        mutate(
            display_name = vapply(
                gene_set,
                function(label) {
                    paste(strwrap(label, width = 48L), collapse = "\n")
                },
                character(1L)
            ),
            signed_significance = if_else(
                Direction == "Down",
                -1,
                1
            ) * -log10(pmax(FDR, .Machine$double.xmin))
        )

    plot_height = max(1200L, 90L * nrow(plot_table) + 350L)
    png(output_path, width = 2200, height = plot_height, res = 180)
    old_par = par(no.readonly = TRUE)
    on.exit({
        par(old_par)
        dev.off()
    }, add = TRUE)
    par(mar = c(5, 22, 4, 2))
    colors = ifelse(
        plot_table$Direction == "Up",
        "#B2182B",
        ifelse(plot_table$Direction == "Down", "#2166AC", "#666666")
    )
    limits = range(
        c(plot_table$signed_significance, -log10(0.05), log10(0.05)),
        finite = TRUE
    )
    bars = barplot(
        plot_table$signed_significance,
        names.arg = plot_table$display_name,
        horiz = TRUE,
        las = 1,
        col = colors,
        border = NA,
        xlab = "Signed -log10 directional FDR",
        main = title,
        xlim = limits
    )
    abline(v = c(-1, 1) * -log10(0.05), lty = 2, col = "#444444")
    abline(v = 0, col = "#222222")
    invisible(bars)
}

counts = collect_counts(
    quant_manifest,
    quant_manifest_path,
    args$quant_id_column,
    args$quant_count_column
)
annotation = read_annotation(
    args$annotation,
    args$annotation_id_attributes,
    args$strip_id_versions
)
identifier_map = build_identifier_map(
    annotation,
    rownames(counts),
    args$strip_id_versions
)
count_annotation = make_count_annotation(identifier_map)

sample_groups = factor(quant_manifest$group)
sample_groups = relevel(sample_groups, ref = args$control_group)

sample_metadata = quant_manifest %>%
    transmute(
        sample = as.character(name),
        group = as.character(group)
    )
write_tsv(
    sample_metadata,
    file.path(output_dir, "sample_metadata.tsv")
)

dge_list = DGEList(counts = counts, group = sample_groups)
keep = filterByExpr(dge_list)
dge_list = dge_list[keep, , keep.lib.sizes = FALSE]
dge_list = normLibSizes(dge_list)

mds = plotMDS(dge_list, top = 500L, plot = FALSE)
mds_data = tibble(
    sample = colnames(dge_list),
    group = as.character(sample_groups),
    dimension_1 = mds$x,
    dimension_2 = mds$y,
    dimension_1_variance = mds$var.explained[[mds$dim.plot[[1L]]]],
    dimension_2_variance = mds$var.explained[[mds$dim.plot[[2L]]]],
    axis_label = mds$axislabel,
    top_features = mds$top,
    gene_selection = mds$gene.selection
)
write_tsv(
    mds_data,
    file.path(output_dir, "edgeR_mds_data.tsv")
)

design = model.matrix(~0 + sample_groups)
colnames(design) = levels(sample_groups)
dge_list = estimateDisp(dge_list, design, robust = TRUE)

bcv_data = tibble(
    feature_id = rownames(dge_list),
    average_log_cpm = dge_list$AveLogCPM,
    tagwise_dispersion = dge_list$tagwise.dispersion,
    tagwise_bcv = sqrt(dge_list$tagwise.dispersion),
    trended_dispersion = dge_list$trended.dispersion,
    trended_bcv = sqrt(dge_list$trended.dispersion),
    common_dispersion = dge_list$common.dispersion,
    common_bcv = sqrt(dge_list$common.dispersion)
)
write_tsv(
    bcv_data,
    file.path(output_dir, "edgeR_bcv_data.tsv")
)

fit = glmQLFit(dge_list, design, robust = TRUE)

norm_feature_counts = as.data.frame(cpm(dge_list)) %>%
    rownames_to_column("feature_id")
    
write_tsv(
    norm_feature_counts,
    file.path(output_dir, "feature_counts.tsv")
)

gene_sets_path = args$gene_sets
gmt = NULL
resolved_gene_sets = NULL
coverage = NULL
if (!is.null(gene_sets_path)) {
    gmt = read_gmt(gene_sets_path)
    resolution = resolve_gene_sets(
        gmt,
        identifier_map,
        args$strip_id_versions
    )
    resolved_gene_sets = resolution$sets
    coverage = make_gene_set_coverage(
        resolution$table,
        rownames(dge_list)
    )
    write_tsv(
        resolution$table,
        file.path(output_dir, "gene_set_resolution.tsv")
    )
    write_tsv(
        coverage,
        file.path(output_dir, "gene_set_coverage.tsv")
    )
    resolved_members = sum(coverage$matched_gmt_members)
    if (resolved_members == 0L) {
        stop(
            paste0(
                "None of the GMT members resolved to count-matrix feature IDs. ",
                "Use matching IDs, supply --annotation, select custom attributes ",
                "with --annotation_id_attributes, or enable ",
                "--strip_id_versions. See gene_set_resolution.tsv."
            )
        )
    }
    message(sprintf(
        "Resolved %d of %d GMT members to the count matrix.",
        resolved_members,
        sum(coverage$gmt_members)
    ))
}

target_groups = setdiff(levels(sample_groups), args$control_group)
if (length(target_groups) == 0L) {
    stop("At least one non-control group is required.")
}

for (target_group in target_groups) {
    contrast = numeric(ncol(design))
    names(contrast) = colnames(design)
    contrast[[target_group]] = 1
    contrast[[args$control_group]] = -1

    contrast_name = paste0(
        "group_",
        safe_name(target_group),
        "_vs_",
        safe_name(args$control_group)
    )
    contrast_dir = file.path(output_dir, contrast_name)
    dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)

    treated = glmTreat(fit, contrast = contrast, lfc = args$lfc)
    de_table = topTags(
        treated,
        n = Inf,
        sort.by = "none"
    )$table %>%
        rownames_to_column("feature_id")
    if (nrow(count_annotation) > 0L) {
        de_table = de_table %>%
            left_join(count_annotation, by = "feature_id")
    }
    write_tsv(
        de_table,
        file.path(contrast_dir, "edgeR_results.tsv")
    )

    if (!is.null(gmt)) {
        tested_sets = lapply(
            resolved_gene_sets,
            function(genes) {
                match(
                    intersect(genes, rownames(dge_list)),
                    rownames(dge_list)
                )
            }
        )
        tested_sets = tested_sets[
            vapply(tested_sets, length, integer(1L)) >= 2L
        ]

        if (length(tested_sets) == 0L) {
            warning(sprintf(
                "No gene sets contain at least two retained genes for %s.",
                contrast_name
            ))
        } else {
            fry_result = fry(
                dge_list,
                index = tested_sets,
                design = design,
                contrast = contrast,
                sort = "none"
            ) %>%
                as.data.frame() %>%
                rownames_to_column("gene_set") %>%
                left_join(
                    coverage %>%
                        select(
                            gene_set,
                            description,
                            gmt_members,
                            matched_gmt_members,
                            count_matrix_members,
                            tested_members,
                            tested_gmt_members,
                            count_matrix_coverage,
                            tested_coverage
                        ),
                    by = "gene_set"
                ) %>%
                arrange(FDR, PValue)

            write_tsv(
                fry_result,
                file.path(contrast_dir, "fry_results.tsv")
            )
            write_fry_plot(
                fry_result,
                file.path(contrast_dir, "fry_signed_significance.png"),
                paste0("fry: ", target_group, " vs ", args$control_group),
                args$plot_top_n
            )
        }
    }
}
