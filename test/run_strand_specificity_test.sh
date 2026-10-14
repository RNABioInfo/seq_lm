#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d "${TMPDIR:-/tmp}/seq-lm-strand-tests.XXXXXX")"
export NXF_OFFLINE=true
mkdir -p "$test_root/bin"
cat > "$test_root/bin/oarfish" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
args=("$@")
output=''
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) output="$2" ;;
        --genome-alignments|--annotation|--genome-fasta) test -f "$2" ;;
        --filter-group) test "$2" = no-filters ;;
        -d) [[ "$2" == fw || "$2" == rc || "$2" == both ]] ;;
        -j) test "$2" = 4 ;;
        *) printf 'Unexpected Oarfish argument: %s\n' "$1" >&2; exit 1 ;;
    esac
    shift 2
done
test -n "$output"
printf '%s\n' "${args[@]}" > "$output.args"
printf 'tname\tlen\tnum_reads\ntx1\t4\t10\n' > "$output.quant"
MOCK
chmod +x "$test_root/bin/oarfish"
export PATH="$test_root/bin:$PATH"
cd "$test_root"
run_case() {
    local name="$1" script="$2"
    shift 2
    nextflow -log "$test_root/$name.log" run -lib "$repo_dir/lib" \
        "$repo_dir/test/$script" -c "$repo_dir/test/strand_specificity_test.config" \
        -work-dir "$test_root/work" -ansi-log false "$@"
}
run_case parameters test_parameter_validation.nf
run_case checkpoints test_sample_checkpoints.nf
if run_case stale test_sample_checkpoints.nf --stale true > "$test_root/stale.stdout" 2>&1; then
    printf 'Unexpected successful corrupt checkpoint restoration\n' >&2
    exit 1
fi
rg -q 'quantification checksum does not match' "$test_root/stale.stdout"
run_case oarfish_default test_oarfish_strand_specificity.nf
run_case writer_default test_write_sample_checkpoints.nf
for strand in forward reverse unstranded; do
    run_case "oarfish_$strand" test_oarfish_strand_specificity.nf --strand_specificity "$strand"
    run_case "writer_$strand" test_write_sample_checkpoints.nf --strand_specificity "$strand"
done
printf "params.validate_params = false\nparams.disable_ping = true\n" > "$test_root/no-schema.config"
if nextflow -log "$test_root/no-schema.log" run -lib "$repo_dir/lib" "$repo_dir/main.nf" \
    -c "$test_root/no-schema.config" --strand_specificity invalid -ansi-log false \
    > "$test_root/no-schema.stdout" 2>&1; then
    printf 'Unexpected successful invalid strand configuration without schema validation\n' >&2
    exit 1
fi
rg -q -- '--strand_specificity must be forward, reverse, or unstranded' "$test_root/no-schema.stdout"
printf 'alias,group,bam_dir\nsample,control,%s\n' "$test_root" > "$test_root/samples.csv"
if nextflow -log "$test_root/schema.log" run -lib "$repo_dir/lib" "$repo_dir/main.nf" \
    --sample_sheet "$test_root/samples.csv" --strand_specificity invalid -ansi-log false \
    > "$test_root/schema.stdout" 2>&1; then
    printf 'Unexpected successful invalid strand configuration with schema validation\n' >&2
    exit 1
fi
rg -q -- '--strand_specificity' "$test_root/schema.stdout"
rg -q 'Validation of pipeline parameters failed' "$test_root/schema.stdout"
nextflow -log "$test_root/help.log" run -lib "$repo_dir/lib" "$repo_dir/main.nf" \
    --help -ansi-log false > "$test_root/help.stdout" 2>&1
rg -q -- '--strand_specificity' "$test_root/help.stdout"
printf 'Strand specificity checks passed: %s\n' "$test_root"
