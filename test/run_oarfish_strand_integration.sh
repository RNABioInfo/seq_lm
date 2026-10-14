#!/usr/bin/env bash
# Requires a running Docker daemon, the workflow's Oarfish image, samtools, and Python 3.
set -euo pipefail
image='rnabioinfo/seq_lm_oarfish:v1.0.0'
docker image inspect "$image" >/dev/null
test_root="$(mktemp -d "${TMPDIR:-/tmp}/seq-lm-oarfish-strands.XXXXXX")"
docker run --rm --entrypoint oarfish "$image" --version > "$test_root/version.txt"
rg -q '0\.10\.0' "$test_root/version.txt"
docker run --rm --entrypoint oarfish "$image" --help > "$test_root/help.txt"
rg -q -- '--strand-filter' "$test_root/help.txt"
python3 - "$test_root" <<'PY'
import pathlib
import random
import sys

root = pathlib.Path(sys.argv[1])
rng = random.Random(42)
genome = ''.join(rng.choice('ACGT') for _ in range(700))
(root / 'genome.fa').write_text(f'>chr1\n{genome}\n')
gtf = []
for transcript, start, end, strand in [('plus', 50, 249, '+'), ('minus', 350, 549, '-')]:
    attrs = f'gene_id "{transcript}"; transcript_id "{transcript}";'
    for feature in ['transcript', 'exon']:
        gtf.append(f'chr1\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}')
(root / 'annotation.gtf').write_text('\n'.join(gtf) + '\n')
sam = ['@HD\tVN:1.6\tSO:unsorted', '@SQ\tSN:chr1\tLN:700']
# Forward means sense relative to the transcript, including minus-strand genes.
for transcript, pos, flag, count in [('plus', 70, 0, 2), ('plus', 70, 16, 3),
                                     ('minus', 380, 16, 7), ('minus', 380, 0, 11)]:
    sequence = genome[pos - 1:pos + 99]
    for index in range(count):
        sam.append(f'{transcript}_{flag}_{index}\t{flag}\tchr1\t{pos}\t60\t100=\t*\t0\t0\t'
                   f'{sequence}\t{"I" * 100}\tAS:i:100\tNM:i:0\tMD:Z:100')
(root / 'reads.sam').write_text('\n'.join(sam) + '\n')
PY
samtools view -b "$test_root/reads.sam" | samtools collate -O -u - > "$test_root/reads.bam"
for strand in fw rc both; do
    docker run --rm -v "$test_root:/data" --entrypoint oarfish "$image" \
        -j 1 --genome-alignments /data/reads.bam --annotation /data/annotation.gtf \
        --genome-fasta /data/genome.fa --filter-group no-filters -d "$strand" \
        --output "/data/$strand"
done
python3 - "$test_root" <<'PY'
import csv
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
for strand, expected in [('fw', {'plus': 2, 'minus': 7}),
                         ('rc', {'plus': 3, 'minus': 11}),
                         ('both', {'plus': 5, 'minus': 18})]:
    with (root / f'{strand}.quant').open() as handle:
        observed = {row['tname']: float(row['num_reads']) for row in csv.DictReader(handle, delimiter='\t')}
    assert observed.keys() == expected.keys(), (strand, observed)
    for transcript, count in expected.items():
        assert math.isclose(observed[transcript], count, abs_tol=1e-6), (strand, observed, expected)
    print(f'Oarfish {strand} transcript-relative counts passed: {observed}')
PY
printf 'Oarfish container integration passed: %s\n' "$test_root"
