#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
docker build -t rnabioinfo/seq_lm_ica:v1.0.0 "${repo_dir}/../docker_containers/seq_lm_ica"
docker run --rm \
    -v "${repo_dir}/bin:/opt/seq_lm/bin:ro" \
    -e PATH=/opt/seq_lm/bin:/usr/local/bin:/usr/bin:/bin \
    rnabioinfo/seq_lm_ica:v1.0.0 \
    python /opt/seq_lm/bin/imodulon_analysis/tests/smoke.py
