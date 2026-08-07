#!/usr/bin/env bash
set -euo pipefail

test_dir="$(mktemp -d "${TMPDIR:-/tmp}/seq_lm_gsva_generic.XXXXXX")"
trap 'rm -rf -- "${test_dir}"' EXIT

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
analysis_script="${script_dir}/gsva_analysis.R"

Rscript - "${test_dir}" <<'RSCRIPT'
args = commandArgs(trailingOnly = TRUE)
test_dir = args[[1L]]

sample_names = c("control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3")
feature_ids = sprintf("tx%03d", seq_len(30L))
counts = matrix(100, nrow = length(feature_ids), ncol = length(sample_names))
rownames(counts) = feature_ids
colnames(counts) = sample_names
for (sample_index in seq_along(sample_names)) {
    counts[, sample_index] = counts[, sample_index] + ((seq_len(30L) + sample_index) %% 11L)
}
counts[seq_len(5L), 4:6] = counts[seq_len(5L), 4:6] + 350
counts[6:10, 4:6] = counts[6:10, 4:6] * 0.12

write.table(
    data.frame(feature_id = feature_ids, counts, check.names = FALSE),
    file.path(test_dir, "feature_counts.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
write.table(
    data.frame(
        sample = sample_names,
        group = rep(c("baseline", "challenge"), each = 3L)
    ),
    file.path(test_dir, "sample_metadata.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
resolution = rbind(
    data.frame(
        gene_set = "UP_SET",
        description = "synthetic induced program",
        gmt_member = sprintf("GENE%03d", 1:5),
        feature_id = feature_ids[1:5],
        match_method = "synthetic"
    ),
    data.frame(
        gene_set = "DOWN_SET",
        description = "synthetic repressed program",
        gmt_member = sprintf("GENE%03d", 6:10),
        feature_id = feature_ids[6:10],
        match_method = "synthetic"
    )
)
write.table(
    resolution,
    file.path(test_dir, "gene_set_resolution.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
RSCRIPT

Rscript "${analysis_script}" \
    --feature_counts "${test_dir}/feature_counts.tsv" \
    --sample_metadata "${test_dir}/sample_metadata.tsv" \
    --gene_set_resolution "${test_dir}/gene_set_resolution.tsv" \
    --output_dir "${test_dir}/results" \
    --control_group baseline

for output in \
    gsva_scores.tsv \
    gsva_scores_long.tsv \
    gsva_gene_set_coverage.tsv \
    gsva_parameters.tsv \
    gsva_score_heatmap.png \
    gsva_group_boxplots.pdf \
    group_challenge_vs_baseline/gsva_limma_results.tsv
do
    test -s "${test_dir}/results/${output}"
done

Rscript - "${test_dir}/results" <<'RSCRIPT'
args = commandArgs(trailingOnly = TRUE)
results_dir = args[[1L]]
wide = read.delim(
    file.path(results_dir, "gsva_scores.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)
long = read.delim(
    file.path(results_dir, "gsva_scores_long.tsv"),
    stringsAsFactors = FALSE
)
coverage = read.delim(
    file.path(results_dir, "gsva_gene_set_coverage.tsv"),
    stringsAsFactors = FALSE
)
limma = read.delim(
    file.path(results_dir, "group_challenge_vs_baseline", "gsva_limma_results.tsv"),
    stringsAsFactors = FALSE
)

stopifnot(nrow(wide) == 2L)
stopifnot(nrow(long) == 12L)
stopifnot(identical(unique(long$sample), names(wide)[4:9]))
stopifnot(all(coverage$status == "scored"))
stopifnot(all(coverage$resolved_members == 5L))
stopifnot(all(coverage$variable_members == 5L))
stopifnot(all(limma$adjusted_p_value >= 0 & limma$adjusted_p_value <= 1))
stopifnot(limma$effect_size[limma$gene_set == "UP_SET"] > 0)
stopifnot(limma$effect_size[limma$gene_set == "DOWN_SET"] < 0)

for (row_index in seq_len(nrow(wide))) {
    gene_set = wide$gene_set[[row_index]]
    for (sample in names(wide)[4:9]) {
        long_value = long$score[long$gene_set == gene_set & long$sample == sample]
        stopifnot(length(long_value) == 1L)
        stopifnot(isTRUE(all.equal(wide[[sample]][[row_index]], long_value)))
    }
}
RSCRIPT

expect_failure() {
    local label=$1
    shift
    if "$@" >"${test_dir}/${label}.stdout" 2>"${test_dir}/${label}.stderr"; then
        printf 'Expected failure for %s\n' "${label}" >&2
        exit 1
    fi
}

{
    head -n 2 "${test_dir}/feature_counts.tsv"
    tail -n +2 "${test_dir}/feature_counts.tsv"
} > "${test_dir}/duplicate_feature_counts.tsv"
expect_failure duplicate_ids Rscript "${analysis_script}" \
    --feature_counts "${test_dir}/duplicate_feature_counts.tsv" \
    --sample_metadata "${test_dir}/sample_metadata.tsv" \
    --gene_set_resolution "${test_dir}/gene_set_resolution.tsv" \
    --output_dir "${test_dir}/duplicate_results" \
    --control_group baseline

sed '$d' "${test_dir}/sample_metadata.tsv" > "${test_dir}/mismatched_metadata.tsv"
expect_failure mismatched_samples Rscript "${analysis_script}" \
    --feature_counts "${test_dir}/feature_counts.tsv" \
    --sample_metadata "${test_dir}/mismatched_metadata.tsv" \
    --gene_set_resolution "${test_dir}/gene_set_resolution.tsv" \
    --output_dir "${test_dir}/mismatch_results" \
    --control_group baseline

awk -F '\t' 'BEGIN { OFS = "\t" } NR == 2 { $2 = "Inf" } { print }' \
    "${test_dir}/feature_counts.tsv" > "${test_dir}/infinite_feature_counts.tsv"
expect_failure nonfinite_expression Rscript "${analysis_script}" \
    --feature_counts "${test_dir}/infinite_feature_counts.tsv" \
    --sample_metadata "${test_dir}/sample_metadata.tsv" \
    --gene_set_resolution "${test_dir}/gene_set_resolution.tsv" \
    --output_dir "${test_dir}/nonfinite_results" \
    --control_group baseline

expect_failure absent_control Rscript "${analysis_script}" \
    --feature_counts "${test_dir}/feature_counts.tsv" \
    --sample_metadata "${test_dir}/sample_metadata.tsv" \
    --gene_set_resolution "${test_dir}/gene_set_resolution.tsv" \
    --output_dir "${test_dir}/absent_control_results" \
    --control_group absent

head -n 2 "${test_dir}/gene_set_resolution.tsv" > "${test_dir}/undersized_resolution.tsv"
expect_failure no_scoreable_sets Rscript "${analysis_script}" \
    --feature_counts "${test_dir}/feature_counts.tsv" \
    --sample_metadata "${test_dir}/sample_metadata.tsv" \
    --gene_set_resolution "${test_dir}/undersized_resolution.tsv" \
    --output_dir "${test_dir}/undersized_results" \
    --control_group baseline

printf 'Generic GSVA smoke test passed.\n'
