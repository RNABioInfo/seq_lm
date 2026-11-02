#!/usr/bin/env python3
# Synthetic quantifier for orchestration testing, never a biological analysis.
import subprocess
import sys
import time
from pathlib import Path

args = sys.argv
out = args[args.index("--output") + 1]
bam = args[args.index("--genome-alignments") + 1]
if out == "fast_0":
    time.sleep(4)
n = int(subprocess.check_output(["samtools", "view", "-c", bam], text=True))
Path(out + ".quant").write_text("tname\tlen\tnum_reads\ntx1\t100\t" + str(n) + "\n")
