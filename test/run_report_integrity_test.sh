#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d "${TMPDIR:-/tmp}/seq-lm-report-integrity.XXXXXX")"
export NXF_OFFLINE=true
cd "$test_root"
run_case() {
    local kind="$1" scenario="$2" expected="$3"
    local log="$test_root/${kind}_${scenario}.log"
    if nextflow -log "$log" run -lib "$repo_dir/lib" "$repo_dir/test/test_report_join_integrity.nf" \
        --kind "$kind" --scenario "$scenario" -work-dir "$test_root/work" -ansi-log false >"$log.stdout" 2>&1; then
        if [[ "$expected" != "pass" ]]; then
            cat "$log.stdout"
            printf 'Unexpected successful %s/%s join\n' "$kind" "$scenario" >&2
            exit 1
        fi
    else
        if [[ "$expected" == "pass" ]] || ! rg -qi "$expected" "$log.stdout"; then
            cat "$log.stdout"
            exit 1
        fi
    fi
    printf 'Verified %s/%s\n' "$kind" "$scenario"
}
for kind in ica qc; do
    run_case "$kind" valid pass
    run_case "$kind" missing 'join mismatch'
    run_case "$kind" missing_report 'join mismatch'
    run_case "$kind" duplicate_report 'duplicate'
    run_case "$kind" duplicate_snapshot 'duplicate'
done
run_case ica sequence_mismatch 'Inconsistent ICA report sequence'
run_case qc sequence_gap 'Missing report sequence'
run_case qc duplicate_sequence 'duplicate report sequence'
nextflow -log "$test_root/publication.log" run -lib "$repo_dir/lib" "$repo_dir/test/test_publish_qc_report_snapshot.nf" -work-dir "$test_root/work" -ansi-log false
if nextflow -log "$test_root/stale.log" run -lib "$repo_dir/lib" "$repo_dir/test/test_publish_qc_report_snapshot.nf" --stale true -work-dir "$test_root/work" -ansi-log false >"$test_root/stale.stdout" 2>&1; then
    printf 'Unexpected successful stale report publication\n' >&2
    exit 1
fi
rg -q 'stale or invalid QC report revision' "$test_root/stale.stdout"
printf 'Report integrity checks passed: %s\n' "$test_root"
