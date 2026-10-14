# Contrast definitions shared by edgeR and GSVA.
read_contrasts = function(path, groups) {
    if (is.null(path) || !file.exists(path) || dir.exists(path)) {
        stop("--contrasts does not match a valid file.")
    }

    definitions = readr::read_tsv(
        path, col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE, progress = FALSE
    )

    columns = c("contrast_id", "target_group", "reference_group")
    if (!identical(names(definitions), columns) || nrow(definitions) < 1L) {
        stop("The contrast TSV requires contrast_id, target_group, reference_group and at least one contrast row.")
    }

    if (anyNA(definitions) || any(as.matrix(definitions) == "") ||
        anyDuplicated(definitions$contrast_id) ||
        anyDuplicated(definitions[, c("target_group", "reference_group")]) ||
        any(definitions$target_group == definitions$reference_group) ||
        any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", definitions$contrast_id))) {
        stop("Check if contrasts TSV contains invalid IDs, empty fields, self-comparisons, or duplicates.")
    }

    selected = unique(c(definitions$target_group, definitions$reference_group))
    if (!setequal(groups, selected)) {
        stop("Analysis groups must exactly match groups referenced by the contrast TSV.")
    }

    definitions
}
