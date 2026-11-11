# Live analysis completion

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

For BAM files written directly by MinKNOW, seq_lm can detect run completion
from `final_summary.txt` or `final_summary_*.txt`. In the `disabled`
(default) and `log` monitoring modes, seq_lm checks for these files at
startup and during analysis. Completion detection requires access to the
files only.

By default, seq_lm waits until the final summary and BAM files have remained
unchanged for at least 30 seconds. Set `minknow_completion_settle_seconds`
to change this interval. Other runs can continue after a sample completes.
The final report is available after all live inputs close and analysis ends.

## Run directory

For automatic detection, set `bam_dir` to the dated MinKNOW run directory,
its `bam_pass` or `bam_fail` subdirectory, or a barcode or `unclassified`
directory directly within either subdirectory. For other layouts, set
`minknow_run_dir` in the sample sheet to the directory for that run.

```csv
alias,group,bam_dir,is_live,minknow_run_dir
sample1,control,/data/run/bams/barcode01,true,/data/run
```

## Manual completion

For externally aligned, reprocessed, or copied BAMs, set
`detect_minknow_completion = false`. After all BAM files are available,
create a `STOP` file in the sample's BAM directory. Use the same procedure
if the final summary is missing or completion is not detected.

```bash
touch /data/run/bams/barcode01/STOP
```

The marker signals seq_lm to finish processing the sample's BAM files.
To configure seq_lm to stop the acquisition in MinKNOW, see the
[automatic sequencing termination instructions](dea_termination.md#stop-sequencing-automatically).
