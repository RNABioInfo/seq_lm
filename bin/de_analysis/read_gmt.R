#!/usr/bin/env Rscript

# Read Gene Matrix Transposed (GMT) files using base R only.
#
# A GMT file is tab-delimited, with one gene set per line:
#   set_name<TAB>description<TAB>gene_1<TAB>gene_2<TAB>...
#
# Example:
#   HALLMARK_APOPTOSIS<TAB>Apoptosis genes<TAB>CASP3<TAB>BAX<TAB>TP53
#
# Usage from another R script:
#   source("read_gmt.R")
#   pathways <- read_gmt("pathways.gmt")
#   pathways[["HALLMARK_APOPTOSIS"]]$genes
#
# Command-line usage:
#   Rscript read_gmt.R pathways.gmt

read_gmt <- function(file,
                     keep_description = TRUE,
                     skip_blank = TRUE,
                     comment_char = "",
                     duplicate_names = c("error", "make_unique", "merge"),
                     encoding = "UTF-8") {
  duplicate_names <- match.arg(duplicate_names)

  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop("'file' must be a single, non-missing character string.", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop(sprintf("GMT file does not exist: %s", file), call. = FALSE)
  }

  lines <- readLines(
    con = file,
    warn = FALSE,
    encoding = encoding,
    skipNul = TRUE
  )

  # Remove a possible UTF-8 byte-order mark from the first field.
  if (length(lines) > 0L) {
    lines[1L] <- sub("^\ufeff", "", lines[1L])
  }

  line_numbers <- seq_along(lines)

  if (skip_blank) {
    keep <- nzchar(trimws(lines))
    lines <- lines[keep]
    line_numbers <- line_numbers[keep]
  }

  if (nzchar(comment_char)) {
    keep <- !startsWith(trimws(lines), comment_char)
    lines <- lines[keep]
    line_numbers <- line_numbers[keep]
  }

  if (length(lines) == 0L) {
    return(structure(
      list(),
      class = c("gmt", "list"),
      source = normalizePath(file, mustWork = FALSE)
    ))
  }

  records <- lapply(seq_along(lines), function(i) {
    fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1L]]

    if (length(fields) < 2L) {
      stop(
        sprintf(
          "Invalid GMT record on line %d: expected at least a name and description.",
          line_numbers[i]
        ),
        call. = FALSE
      )
    }

    set_name <- fields[1L]
    description <- fields[2L]
    genes <- if (length(fields) > 2L) fields[-c(1L, 2L)] else character()

    if (!nzchar(set_name)) {
      stop(
        sprintf("Invalid GMT record on line %d: gene-set name is empty.", line_numbers[i]),
        call. = FALSE
      )
    }

    # Empty trailing or intermediate gene fields are not gene identifiers.
    genes <- genes[nzchar(genes)]

    record <- list(genes = genes)
    if (keep_description) {
      record <- c(list(description = description), record)
    }
    record
  })

  set_names <- vapply(
    lines,
    function(line) strsplit(line, "\t", fixed = TRUE)[[1L]][1L],
    character(1L),
    USE.NAMES = FALSE
  )

  duplicated_names <- unique(set_names[duplicated(set_names)])
  if (length(duplicated_names) > 0L) {
    if (duplicate_names == "error") {
      stop(
        sprintf(
          "Duplicate gene-set name(s): %s",
          paste(duplicated_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (duplicate_names == "make_unique") {
      set_names <- make.unique(set_names, sep = "_")
    }

    if (duplicate_names == "merge") {
      groups <- split(seq_along(records), factor(set_names, levels = unique(set_names)))
      records <- lapply(groups, function(indices) {
        merged_genes <- unique(unlist(
          lapply(records[indices], `[[`, "genes"),
          use.names = FALSE
        ))

        record <- list(genes = merged_genes)
        if (keep_description) {
          descriptions <- unique(vapply(
            records[indices],
            `[[`,
            character(1L),
            "description"
          ))
          record <- c(
            list(description = paste(descriptions[nzchar(descriptions)], collapse = "; ")),
            record
          )
        }
        record
      })
      set_names <- names(groups)
    }
  }

  names(records) <- set_names
  structure(
    records,
    class = c("gmt", "list"),
    source = normalizePath(file, mustWork = FALSE)
  )
}

# Convert the nested GMT representation into a long data frame.
gmt_to_long <- function(gmt) {
  if (!is.list(gmt)) {
    stop("'gmt' must be the list returned by read_gmt().", call. = FALSE)
  }

  if (length(gmt) == 0L) {
    return(data.frame(
      set_name = character(),
      description = character(),
      gene = character(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(seq_along(gmt), function(i) {
    record <- gmt[[i]]
    genes <- record$genes
    description <- if (!is.null(record$description)) record$description else NA_character_

    if (length(genes) == 0L) {
      return(data.frame(
        set_name = names(gmt)[i],
        description = description,
        gene = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      set_name = rep(names(gmt)[i], length(genes)),
      description = rep(description, length(genes)),
      gene = genes,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

print.gmt <- function(x, ...) {
  gene_counts <- vapply(x, function(record) length(record$genes), integer(1L))
  cat(sprintf(
    "<GMT collection: %d gene set%s, %d gene-set memberships>\n",
    length(x),
    if (length(x) == 1L) "" else "s",
    sum(gene_counts)
  ))

  if (length(x) > 0L) {
    preview_count <- min(length(x), 6L)
    preview <- data.frame(
      gene_set = names(x)[seq_len(preview_count)],
      genes = gene_counts[seq_len(preview_count)],
      row.names = NULL,
      check.names = FALSE
    )
    print(preview, row.names = FALSE)
    if (length(x) > preview_count) {
      cat(sprintf("... and %d more\n", length(x) - preview_count))
    }
  }

  invisible(x)
}

# Run a small command-line interface only when this file is executed directly.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 1L || args %in% c("-h", "--help")) {
    cat(
      "Usage: Rscript read_gmt.R FILE.gmt\n",
      "\n",
      "Reads and summarizes a Gene Matrix Transposed file using base R only.\n",
      sep = ""
    )
    quit(status = if (length(args) == 1L) 0L else 1L)
  }

  result <- read_gmt(args[1L])
  print(result)
}