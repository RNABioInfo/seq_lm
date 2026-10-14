#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d "${TMPDIR:-/tmp}/seq-lm-stability-tests.XXXXXX")"
export NXF_OFFLINE=true
export PATH="${repo_dir}/bin:$PATH"
cd "$test_root"

run_case() {
    local label="$1" script="$2"
    shift 2
    local log="$test_root/${label}.log"
    if ! nextflow -log "$log" run -lib "$repo_dir/lib" "$repo_dir/test/$script" \
        -c "$repo_dir/test/stability_test.config" --out_dir "$test_root/output-$label" \
        -work-dir "$test_root/work-$label" -ansi-log false "$@" >"$log.stdout" 2>&1; then
        cat "$log.stdout"
        exit 1
    fi
    printf 'Verified %s\n' "$label"
}

run_case log test_stability_workflow.nf --test_behavior log
run_case reset-log test_stability_workflow.nf --test_behavior log --test_reset true
run_case terminate test_stability_workflow.nf --test_behavior terminate
run_case missing-id test_stability_workflow.nf --test_behavior terminate --test_stop_mode missing
run_case retry test_stability_workflow.nf --test_behavior terminate --test_stop_mode retry
run_case reset-terminate test_stability_workflow.nf --test_behavior terminate --test_reset true
run_case ica test_ica_stability.nf
for scenario in shared separate excluded retry log; do
    run_case "multiplexed-$scenario" test_multiplexed_stability.nf --test_scenario "$scenario"
done

failure_log="$test_root/marker-failure.log"
if nextflow -log "$failure_log" run -lib "$repo_dir/lib" "$repo_dir/test/test_multiplexed_stability.nf" \
    -c "$repo_dir/test/stability_test.config" --out_dir "$test_root/output-marker-failure" \
    -work-dir "$test_root/work-marker-failure" -ansi-log false --test_scenario marker_failure >"$failure_log.stdout" 2>&1; then
    printf 'Unexpected successful STOP marker failure test\n' >&2
    exit 1
fi
if ! rg -q "MinKNOW run 'shared-run' was stopped, but STOP markers could not be created:.*a2/STOP" "$failure_log.stdout"; then
    cat "$failure_log.stdout"
    exit 1
fi
printf 'Verified marker-failure diagnostic\n'
run_case ingress-stop-drain test_bam_ingress_unequal_stops.nf
printf 'Stability integration artifacts: %s\n' "$test_root"
