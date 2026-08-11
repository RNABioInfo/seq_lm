#!/usr/bin/env bash
set -euo pipefail

test_dir="$(mktemp -d "${TMPDIR:-/tmp}/seq_lm_edgeR_generic.XXXXXX")"
trap 'rm -rf -- "${test_dir}"' EXIT

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "${test_dir}/counts" "${test_dir}/custom_counts"

{
    printf '##gff-version 3\n'
    for index in $(seq 1 20); do
        printf 'chr1\ttest\tgene\t%d\t%d\t.\t+\t.\tID=gene:GENE%03d;Name=Symbol%03d\n' \
            "$((index * 100))" "$((index * 100 + 89))" "${index}" "${index}"
        printf 'chr1\ttest\tmRNA\t%d\t%d\t.\t+\t.\tID=transcript:tx%03d;Parent=gene:GENE%03d;gene_id=GENE%03d\n' \
            "$((index * 100))" "$((index * 100 + 89))" \
            "${index}" "${index}" "${index}"
    done
} > "${test_dir}/annotation.gff3"

{
    printf 'UP_SET\tunrelated synthetic pathway'
    for index in $(seq 1 5); do
        printf '\tGENE%03d' "${index}"
    done
    printf '\nDOWN_SET\tsecond synthetic pathway'
    for index in $(seq 6 10); do
        printf '\tGENE%03d' "${index}"
    done
    printf '\n'
} > "${test_dir}/pathways.gmt"

{
    printf 'DIRECT_UP\tmembers already use count-table IDs'
    for index in $(seq 1 5); do
        printf '\ttx%03d.1' "${index}"
    done
    printf '\nDIRECT_DOWN\tsecond direct-ID pathway'
    for index in $(seq 6 10); do
        printf '\ttx%03d.1' "${index}"
    done
    printf '\n'
} > "${test_dir}/direct_pathways.gmt"

printf 'name\tgroup\tcount_file\n' > "${test_dir}/manifest.tsv"
printf 'name\tgroup\tcount_file\n' > "${test_dir}/custom_manifest.tsv"
for sample_index in $(seq 1 6); do
    if ((sample_index <= 3)); then
        group='baseline'
        sample_name="control_${sample_index}"
    else
        group='challenge'
        sample_name="treated_$((sample_index - 3))"
    fi
    count_file="counts/${sample_name}.csv"
    printf '%s\t%s\t%s\n' "${sample_name}" "${group}" "${count_file}" \
        >> "${test_dir}/manifest.tsv"
    printf '%s\t%s\tcustom_counts/%s.csv\n' \
        "${sample_name}" "${group}" "${sample_name}" \
        >> "${test_dir}/custom_manifest.tsv"

    awk -v sample="${sample_index}" 'BEGIN {
        print "target_id,length,est_counts"
        for (gene_index = 1; gene_index <= 20; gene_index++) {
            count = 100 + ((gene_index + sample) % 7)
            if (sample > 3 && gene_index <= 5) {
                count = 380 + 7 * sample + gene_index
            } else if (sample > 3 && gene_index >= 6 && gene_index <= 10) {
                count = 18 + sample + (gene_index % 3)
            }
            printf "tx%03d.1,900,%d\n", gene_index, count
        }
    }' > "${test_dir}/${count_file}"

    awk -F ',' 'BEGIN { OFS = "," }
        NR == 1 { print "custom_feature", "custom_abundance"; next }
        { print $1, $3 }
    ' "${test_dir}/${count_file}" \
        > "${test_dir}/custom_counts/${sample_name}.csv"
done

Rscript "${script_dir}/edgeR_analysis.R" \
    --quant_manifest "${test_dir}/manifest.tsv" \
    --annotation "${test_dir}/annotation.gff3" \
    --gene_sets "${test_dir}/pathways.gmt" \
    --output_dir "${test_dir}/results" \
    --control_group baseline \
    --strip_id_versions

awk -F '\t' '
    NR > 1 && ($4 != $3 || $8 != 1) {
        print "Incomplete GMT resolution for " $1 > "/dev/stderr"
        exit 1
    }
    END {
        if (NR != 3) {
            print "Expected two gene sets in coverage output" > "/dev/stderr"
            exit 1
        }
    }
' "${test_dir}/results/gene_set_coverage.tsv"

test -s "${test_dir}/results/gene_set_resolution.tsv"
test -s "${test_dir}/results/sample_metadata.tsv"
test -s "${test_dir}/results/edgeR_bcv_data.tsv"
awk -F '\t' '
    NR == 1 && $0 != "feature_id\taverage_log_cpm\ttagwise_dispersion\ttagwise_bcv\ttrended_dispersion\ttrended_bcv\tcommon_dispersion\tcommon_bcv" {
        print "Unexpected BCV export columns" > "/dev/stderr"
        exit 1
    }
    NR > 1 && (NF != 8 || $1 == "" || $3 < 0 || $4 < 0 || $5 < 0 || $6 < 0 || $7 < 0 || $8 < 0) {
        print "Invalid BCV export row" > "/dev/stderr"
        exit 1
    }
    END {
        if (NR < 2) {
            print "BCV export contains no retained features" > "/dev/stderr"
            exit 1
        }
    }
' "${test_dir}/results/edgeR_bcv_data.tsv"
test -s "${test_dir}/results/group_challenge_vs_baseline/fry_results.tsv"
test -s "${test_dir}/results/group_challenge_vs_baseline/edgeR_results.tsv"
test -s "${test_dir}/results/group_challenge_vs_baseline/fry_signed_significance.png"

Rscript "${script_dir}/edgeR_analysis.R" \
    --quant_manifest "${test_dir}/custom_manifest.tsv" \
    --gene_sets "${test_dir}/direct_pathways.gmt" \
    --output_dir "${test_dir}/direct_results" \
    --control_group baseline \
    --quant_id_column custom_feature \
    --quant_count_column custom_abundance

awk -F '\t' '
    NR > 1 && ($4 != $3 || $8 != 1) {
        print "Incomplete direct GMT resolution for " $1 > "/dev/stderr"
        exit 1
    }
' "${test_dir}/direct_results/gene_set_coverage.tsv"
test -s "${test_dir}/direct_results/sample_metadata.tsv"
test -s "${test_dir}/direct_results/edgeR_bcv_data.tsv"
test -s "${test_dir}/direct_results/group_challenge_vs_baseline/fry_results.tsv"

printf 'Generic edgeR/GMT smoke test passed.\n'
