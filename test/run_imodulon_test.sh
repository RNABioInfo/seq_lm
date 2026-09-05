#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d "${TMPDIR:-/tmp}/seq-lm-ica-test.XXXXXX")"
export PATH="${repo_dir}/bin:$PATH"
python3 - "$repo_dir" "$test_root" <<'PY'
from pathlib import Path
import sys
sys.path.insert(0,str(Path(sys.argv[1])/'bin/imodulon_analysis/tests'))
from smoke import fixture
root=Path(sys.argv[2])/'fixture with spaces'
root.mkdir()
fixture(root)
(root/'zero.quant').write_text('tname\tnum_reads\nt1\t0\nt1b\t0\nt2\t0\n')
PY
cd "$repo_dir"
for invocation in 1 2; do
    NXF_OFFLINE=true nextflow run test/test_imodulon_analysis.nf \
        -c test/imodulon_test.config \
        --ica_test_fixture "$test_root/fixture with spaces" \
        --ica_test_output "$test_root/output" \
        -work-dir "$test_root/work" -ansi-log false
done
printf 'ICA integration artifacts: %s\n' "$test_root"
