#!/usr/bin/env Rscript

library(optparse)
library(GSVA)
library(limma)
library(readr)
library(tibble)

options = list(
    make_option(c("--feature_counts"), type = "character"),
    make_option(c("--sample_metadata"), type = "character"),
    make_option(c("--gene_set_resolution"), type = "character"),
    make_option(c("--output_dir"), type = "character"),
    make_option(
        c("--control_group"),
        type = "character",
        default = "control",
        help = "Reference group used for pairwise contrasts [default: %default]."
    ),
    make_option(
        c("--min_size"),
        type = "integer",
        default = 2L,
        help = "Minimum number of variable features in a scored set [default: %default]."
    )
)

parser = OptionParser(
    description = paste0(
        "GSVA scoring and limma testing from edgeR-normalized CPM and ",
        "resolved gene-set memberships"
    ),
    option_list = options
)
args = parse_args(parser)

required_options = c(
    "feature_counts",
    "sample_metadata",
    "gene_set_resolution",
    "output_dir"
)
missing_options = required_options[vapply(
    required_options,
    function(name) is.null(args[[name]]) || !nzchar(args[[name]]),
    logical(1L)
)]
if (length(missing_options) > 0L) {
    stop(
        paste0(
            "Missing required option(s): --",
            paste(missing_options, collapse = ", --")
        )
    )
}
if (args$min_size < 2L) {
    stop("--min_size must be at least 2.")
}

require_file = function(path, label) {
    if (!file.exists(path)) {
        stop(sprintf("%s does not exist: %s", label, path))
    }
}

require_columns = function(table, required, label) {
    missing = setdiff(required, colnames(table))
    if (length(missing) > 0L) {
        stop(sprintf(
            "%s is missing required column(s): %s",
            label,
            paste(missing, collapse = ", ")
        ))
    }
}

nonempty_text = function(values, label) {
    values = as.character(values)
    if (any(is.na(values)) || any(!nzchar(values))) {
        stop(sprintf("%s contains missing or empty values.", label))
    }
    values
}

safe_name = function(value) {
    gsub("[^A-Za-z0-9._-]", "_", value)
}

read_expression_matrix = function(path) {
    table = read_tsv(
        path,
        col_types = cols(.default = col_character()),
        show_col_types = FALSE,
        progress = FALSE
    )
    require_columns(table, "feature_id", "Feature-count table")
    if (ncol(table) < 3L) {
        stop("Feature-count table must contain at least two sample columns.")
    }
    feature_ids = nonempty_text(table$feature_id, "Feature-count feature_id column")
    if (anyDuplicated(feature_ids)) {
        duplicates = unique(feature_ids[duplicated(feature_ids)])
        stop(sprintf(
            "Feature-count table contains duplicate feature_id value(s): %s",
            paste(duplicates, collapse = ", ")
        ))
    }

    sample_names = setdiff(colnames(table), "feature_id")
    if (any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
        stop("Feature-count sample column names must be nonempty and unique.")
    }
    raw_values = as.matrix(table[, sample_names, drop = FALSE])
    numeric_values = suppressWarnings(matrix(
        as.numeric(raw_values),
        nrow = nrow(raw_values),
        ncol = ncol(raw_values),
        dimnames = list(feature_ids, sample_names)
    ))
    if (any(is.na(numeric_values))) {
        stop("Feature-count table contains missing or nonnumeric CPM values.")
    }
    if (!all(is.finite(numeric_values))) {
        stop("Feature-count table contains infinite CPM values.")
    }
    if (any(numeric_values < 0)) {
        stop("Feature-count table contains negative CPM values.")
    }
    numeric_values
}

read_sample_metadata = function(path, sample_names, control_group) {
    metadata = read_tsv(
        path,
        col_types = cols(.default = col_character()),
        show_col_types = FALSE,
        progress = FALSE
    )
    require_columns(metadata, c("sample", "group"), "Sample metadata")
    metadata = metadata[, c("sample", "group")]
    metadata$sample = nonempty_text(metadata$sample, "Sample metadata sample column")
    metadata$group = nonempty_text(metadata$group, "Sample metadata group column")
    if (anyDuplicated(metadata$sample)) {
        stop("Sample metadata contains duplicate sample values.")
    }
    missing_metadata = setdiff(sample_names, metadata$sample)
    extra_metadata = setdiff(metadata$sample, sample_names)
    if (length(missing_metadata) > 0L || length(extra_metadata) > 0L) {
        stop(sprintf(
            paste0(
                "Feature-count samples and sample metadata do not match ",
                "(missing metadata: %s; extra metadata: %s)."
            ),
            ifelse(length(missing_metadata), paste(missing_metadata, collapse = ", "), "none"),
            ifelse(length(extra_metadata), paste(extra_metadata, collapse = ", "), "none")
        ))
    }
    metadata = metadata[match(sample_names, metadata$sample), , drop = FALSE]
    if (!(control_group %in% metadata$group)) {
        stop(sprintf("Control group '%s' is absent from sample metadata.", control_group))
    }
    if (length(setdiff(unique(metadata$group), control_group)) == 0L) {
        stop("At least one non-control group is required.")
    }
    metadata
}

read_resolution = function(path) {
    resolution = read_tsv(
        path,
        col_types = cols(.default = col_character()),
        show_col_types = FALSE,
        progress = FALSE
    )
    require_columns(
        resolution,
        c("gene_set", "description", "feature_id"),
        "Gene-set resolution table"
    )
    resolution$gene_set = nonempty_text(
        resolution$gene_set,
        "Gene-set resolution gene_set column"
    )
    resolution$description[is.na(resolution$description)] = ""
    descriptions_per_set = tapply(
        resolution$description,
        resolution$gene_set,
        function(values) length(unique(values))
    )
    if (any(descriptions_per_set > 1L)) {
        stop("Gene-set resolution contains conflicting descriptions for a gene set.")
    }
    resolution
}

prepare_gene_sets = function(resolution, retained_ids, variable_ids, min_size) {
    gene_set_order = unique(resolution$gene_set)
    descriptions = vapply(
        gene_set_order,
        function(gene_set) {
            values = resolution$description[resolution$gene_set == gene_set]
            if (length(values) > 0L) values[[1L]] else ""
        },
        character(1L)
    )
    names(descriptions) = gene_set_order

    resolved_sets = lapply(gene_set_order, function(gene_set) {
        feature_ids = resolution$feature_id[resolution$gene_set == gene_set]
        unique(feature_ids[!is.na(feature_ids) & nzchar(feature_ids)])
    })
    names(resolved_sets) = gene_set_order
    retained_sets = lapply(resolved_sets, intersect, y = retained_ids)
    variable_sets = lapply(resolved_sets, intersect, y = variable_ids)
    scoreable = lengths(variable_sets) >= min_size
    status = ifelse(scoreable, "scored", "below_min_size")

    coverage = tibble(
        gene_set = gene_set_order,
        description = unname(descriptions),
        resolved_members = lengths(resolved_sets),
        retained_members = lengths(retained_sets),
        variable_members = lengths(variable_sets),
        scored_members = ifelse(scoreable, lengths(variable_sets), 0L),
        status = status
    )
    list(
        sets = variable_sets[scoreable],
        coverage = coverage,
        descriptions = descriptions
    )
}

write_heatmap = function(scores, metadata, output_path) {
    scaled_scores = t(scale(t(scores)))
    scaled_scores[!is.finite(scaled_scores)] = 0
    group_levels = unique(metadata$group)
    group_colors = setNames(
        grDevices::hcl.colors(length(group_levels), palette = "Dark 3"),
        group_levels
    )
    column_colors = unname(group_colors[metadata$group])
    longest_gene_set = max(nchar(rownames(scores), type = "width"))
    right_margin = max(24L, ceiling(longest_gene_set * 0.72))
    width = max(2200L, 160L * ncol(scores) + 50L * longest_gene_set + 700L)
    height = max(1400L, 90L * nrow(scores) + 500L)

    png(output_path, width = width, height = height, res = 180)
    old_par = par(no.readonly = TRUE)
    on.exit({
        par(old_par)
        dev.off()
    }, add = TRUE)
    par(xpd = NA)
    stats::heatmap(
        scaled_scores,
        Rowv = if (nrow(scores) > 1L) TRUE else NA,
        Colv = if (ncol(scores) > 1L) TRUE else NA,
        scale = "none",
        col = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(256L),
        ColSideColors = column_colors,
        margins = c(10, right_margin),
        xlab = "Samples",
        ylab = "Gene sets",
        main = "GSVA scores (row-standardized for display)"
    )
    legend(
        "topright",
        legend = group_levels,
        fill = unname(group_colors),
        title = "Group",
        bty = "n",
        inset = c(-0.02, -0.02)
    )
}

write_boxplots = function(scores, metadata, output_path) {
    group_levels = unique(metadata$group)
    group_colors = setNames(
        grDevices::hcl.colors(length(group_levels), palette = "Dark 3"),
        group_levels
    )
    group_factor = factor(metadata$group, levels = group_levels)

    pdf(output_path, width = 11, height = 8.5, onefile = TRUE)
    old_par = par(no.readonly = TRUE)
    on.exit({
        par(old_par)
        dev.off()
    }, add = TRUE)
    set.seed(1L)
    for (gene_set in rownames(scores)) {
        values = as.numeric(scores[gene_set, ])
        par(mar = c(8, 5, 6, 2))
        boxplot(
            values ~ group_factor,
            col = unname(group_colors[group_levels]),
            border = "#444444",
            outline = FALSE,
            las = 2,
            ylab = "GSVA score",
            xlab = "",
            main = gene_set
        )
        stripchart(
            values ~ group_factor,
            vertical = TRUE,
            method = "jitter",
            jitter = 0.08,
            pch = 21,
            bg = "white",
            col = "#222222",
            add = TRUE
        )
    }
}

require_file(args$feature_counts, "Feature-count table")
require_file(args$sample_metadata, "Sample metadata")
require_file(args$gene_set_resolution, "Gene-set resolution table")
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

expression_cpm = read_expression_matrix(args$feature_counts)
metadata = read_sample_metadata(
    args$sample_metadata,
    colnames(expression_cpm),
    args$control_group
)
resolution = read_resolution(args$gene_set_resolution)

log_expression = log2(expression_cpm + 1)
variable_features = apply(
    log_expression,
    1L,
    function(values) max(values) > min(values)
)
variable_expression = log_expression[variable_features, , drop = FALSE]
if (nrow(variable_expression) == 0L) {
    stop("No variable expression features remain for GSVA.")
}

prepared = prepare_gene_sets(
    resolution,
    rownames(log_expression),
    rownames(variable_expression),
    args$min_size
)
write_tsv(
    prepared$coverage,
    file.path(args$output_dir, "gsva_gene_set_coverage.tsv")
)
if (length(prepared$sets) == 0L) {
    stop(sprintf(
        "No gene sets contain at least %d variable retained features.",
        args$min_size
    ))
}

parameter = GSVA::gsvaParam(
    variable_expression,
    prepared$sets,
    kcdf = "Gaussian",
    minSize = args$min_size,
    maxSize = Inf,
    tau = 1,
    maxDiff = TRUE,
    absRanking = FALSE,
    checkNA = "yes",
    use = "all.obs"
)
scores = GSVA::gsva(
    parameter,
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
)
scores = as.matrix(scores)
expected_gene_sets = names(prepared$sets)
if (!setequal(rownames(scores), expected_gene_sets)) {
    stop("GSVA returned a different set of gene-set identifiers than expected.")
}
scores = scores[expected_gene_sets, colnames(expression_cpm), drop = FALSE]
if (any(!is.finite(scores))) {
    stop("GSVA returned missing or infinite scores.")
}

coverage_by_set = prepared$coverage[match(rownames(scores), prepared$coverage$gene_set), ]
gene_counts = setNames(coverage_by_set$scored_members, coverage_by_set$gene_set)
descriptions = prepared$descriptions[rownames(scores)]

wide_scores = data.frame(
    gene_set = rownames(scores),
    description = unname(descriptions),
    n_genes = unname(gene_counts[rownames(scores)]),
    scores,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
write_tsv(wide_scores, file.path(args$output_dir, "gsva_scores.tsv"))

long_scores = do.call(
    rbind,
    lapply(seq_len(nrow(scores)), function(index) {
        data.frame(
            gene_set = rownames(scores)[[index]],
            description = unname(descriptions[[index]]),
            n_genes = unname(gene_counts[[rownames(scores)[[index]]]]),
            sample = colnames(scores),
            group = metadata$group,
            score = as.numeric(scores[index, ]),
            stringsAsFactors = FALSE
        )
    })
)
write_tsv(long_scores, file.path(args$output_dir, "gsva_scores_long.tsv"))

parameters = tibble(
    parameter = c(
        "input",
        "transformation",
        "method",
        "kcdf",
        "min_size",
        "max_size",
        "tau",
        "max_diff",
        "abs_ranking",
        "parallel_backend"
    ),
    value = c(
        "edgeR-filtered TMM-normalized CPM",
        "log2(CPM + 1)",
        "GSVA",
        "Gaussian",
        as.character(args$min_size),
        "Inf",
        "1",
        "TRUE",
        "FALSE",
        "SerialParam"
    )
)
write_tsv(parameters, file.path(args$output_dir, "gsva_parameters.tsv"))

write_heatmap(
    scores,
    metadata,
    file.path(args$output_dir, "gsva_score_heatmap.png")
)
write_boxplots(
    scores,
    metadata,
    file.path(args$output_dir, "gsva_group_boxplots.pdf")
)

group_levels = c(
    args$control_group,
    setdiff(unique(metadata$group), args$control_group)
)
group_factor = factor(metadata$group, levels = group_levels)
design = model.matrix(~0 + group_factor)
colnames(design) = group_levels
fit = limma::lmFit(scores, design)

for (target_group in setdiff(group_levels, args$control_group)) {
    contrast = numeric(ncol(design))
    names(contrast) = colnames(design)
    contrast[[target_group]] = 1
    contrast[[args$control_group]] = -1
    contrast_fit = limma::contrasts.fit(fit, contrasts = contrast)
    contrast_fit = limma::eBayes(contrast_fit)
    statistics = limma::topTable(
        contrast_fit,
        number = Inf,
        adjust.method = "BH",
        sort.by = "none"
    )
    statistics = rownames_to_column(statistics, "gene_set")

    results = data.frame(
        gene_set = statistics$gene_set,
        description = unname(prepared$descriptions[statistics$gene_set]),
        n_genes = unname(gene_counts[statistics$gene_set]),
        target_group = target_group,
        control_group = args$control_group,
        effect_size = statistics$logFC,
        average_score = statistics$AveExpr,
        t_statistic = statistics$t,
        p_value = statistics$P.Value,
        adjusted_p_value = statistics$adj.P.Val,
        log_odds = statistics$B,
        stringsAsFactors = FALSE
    )
    contrast_dir = file.path(
        args$output_dir,
        paste0(
            "group_",
            safe_name(target_group),
            "_vs_",
            safe_name(args$control_group)
        )
    )
    dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
    write_tsv(results, file.path(contrast_dir, "gsva_limma_results.tsv"))
}
