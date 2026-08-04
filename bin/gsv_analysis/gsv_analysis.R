#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(GSVA)

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
    stop("Unable to determine the edgeR analysis script directory.")
}
script_path = sub("^--file=", "", script_arg[[1]])
script_dir = dirname(normalizePath(script_path))

options = list(
    make_option(
        c("--norm_counts"),
        type = "character"
    ),
    make_option(
        c("")
    )
)