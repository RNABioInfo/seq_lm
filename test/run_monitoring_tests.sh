#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d "${TMPDIR:-/tmp}/seq-lm-monitoring-tests.XXXXXX")"
export NXF_OFFLINE=true
export PATH="${repo_dir}/bin:$PATH"
cd "$test_root"

run_case() {
    local label="$1" script="$2"
    shift 2
    local log="$test_root/${label}.log"
    if ! nextflow -log "$log" run -lib "$repo_dir/lib" "$repo_dir/test/$script" \
        -c "$repo_dir/test/monitoring_test.config" --out_dir "$test_root/output-$label" \
        -work-dir "$test_root/work-$label" -ansi-log false "$@" >"$log.stdout" 2>&1; then
        cat "$log.stdout"
        exit 1
    fi
    printf 'Verified %s\n' "$label"
}

run_case dea-monitoring test_dea_monitoring.nf
run_case dea-selected-monitoring test_dea_monitoring.nf --test_stopping_genes true
run_case event-ingress test_event_bam_ingress.nf
run_case independent-collation test_independent_collation.nf
run_case stable-bam-trigger test_bam_ingress_stable_polling.nf
run_case report-join test_report_batch_join.nf
run_case ica-report-join test_ica_report_batch_join.nf
run_case dea-action-guards test_dea_action_guards.nf
run_case acquisition-actions test_acquisition_actions.nf
run_case ica test_ica_analysis_monitoring.nf
run_case gene-contrasts test_gene_contrasts.nf
run_case validation test_parameter_validation.nf
run_case checkpoint-resume test_sample_checkpoints.nf
run_case ingress-stop-drain test_bam_ingress_unequal_stops.nf
printf 'Monitoring integration artifacts: %s\n' "$test_root"
