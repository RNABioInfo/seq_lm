#!/usr/bin/env bash
# Native orchestration test: requires samtools and the report Python dependencies.
# Uses a synthetic quantifier, not a biological quantification benchmark.
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d "${TMPDIR:-/tmp}/seq-lm-event-quantification.XXXXXX")"
mkdir -p "$test_root/tools"
cp "$repo_dir/test/fixtures/synthetic_oarfish.py" "$test_root/tools/oarfish"
chmod +x "$test_root/tools/oarfish"
export PATH="$test_root/tools:$repo_dir/bin:$PATH"
export NXF_OFFLINE=true
cd "$test_root"
nextflow -log "$test_root/nextflow.log" run -lib "$repo_dir/lib" \
    "$repo_dir/test/test_event_quantification.nf" -c "$repo_dir/test/monitoring_test.config" \
    --stopping_method dea --monitoring_behavior log --live_analysis true \
    -work-dir "$test_root/work" -ansi-log false
printf 'Event quantification test artifacts: %s\n' "$test_root"
