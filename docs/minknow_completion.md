# MinKNOW completion without gRPC

Live samples in `monitoring_behavior = "disabled"` or `"log"` automatically close
their input stream after detecting MinKNOW's final summary. This feature reads
files only: it does not contact MinKNOW, stop sequencing, or write a persistent
`STOP` file. Automated termination (`"terminate"`) retains its existing behavior.
Manual `STOP*` files remain available in every live mode.

## Layout and identity

Use direct live MinKNOW BAM output in a dedicated directory for one run.
Automatically recognized `bam_dir` values are the standard dated run directory
(`YYYYMMDD_HHMM_device_flowcell_8hexprotocolid`), its `bam_pass` or `bam_fail`
folder, or a `barcodeNN`/`unclassified` directory immediately below those folders.
The watcher never searches arbitrary ancestors or sibling runs. For a different
layout, specify the optional `minknow_run_dir` column in the seqLM sample sheet:

```csv
alias,group,bam_dir,is_live,minknow_run_dir
sample1,control,/data/run/bams/barcode01,true,/data/run
```

Blank or missing values use automatic recognition. The override must exist and
asserts that its directory belongs to one run. Pooled or reused directories with
multiple runs are outside the supported scope.

Only `final_summary.txt` and `final_summary_*.txt` in the resolved run directory
are considered. Sequencing summaries and HTML reports never signal completion.
The summary must be nonempty, stable key-value text containing `protocol_run_id`,
`acquisition_run_id`, and valid `acquisition_stopped` and `processing_stopped`
timestamps. Processing must stop at or after acquisition. If a `started` timestamp
is present, it is also checked. The protocol identifier must match available
MinKNOW sample-sheet metadata and the short identifier in a standard run-directory
name. Conflicting identity, malformed summaries, and multiple candidate summaries
are logged and leave the watcher open.

## Draining and finalization

Set `detect_minknow_completion = false` to opt out. The default is `true`.
`minknow_completion_settle_seconds` is a positive integer, default `30`.
Once a valid summary is observed, its contents and file metadata and the complete
BAM inventory must remain unchanged for that interval and at least
`bam_stability_polls` observations. New BAMs still follow the normal stability
rules. Any appearing, disappearing, or changing BAM, or any summary change,
restarts settling; invalid or disappearing summary evidence cancels settling.
Completion is never inferred from inactivity alone.

Each sample emits one internal stop event (`minknow-final-summary`), including
the summary path in the log. Samples sharing a run summary independently drain
their own BAM directories. Samples from other runs continue receiving data and
completed samples leave subsequent batch barriers. A valid summary already
present at startup follows the same settling rules.

Input streams close independently. Persistent `FINAL` checkpoints keep their
existing timing: they are written when the collected workflow streams close,
not immediately when one sample's input closes. Final report/checkpoint writing
therefore still waits for the remaining live inputs and downstream work.

## Reliability

[ONT's output specification](https://software-docs.nanoporetech.com/output-specifications/latest/minknow/support/)
describes final summaries as written after acquisition and file writing complete.
This provides completion evidence for direct MinKNOW BAM output, supplemented by
settling to drain files visible to the watcher. Missing summaries, crashes before
summary creation, and incompatible summary formats leave inputs open; use manual
`STOP` as the fallback.

Disable detection for delayed external alignment, post-run reprocessing, or
asynchronously copied BAMs. MinKNOW's summary cannot guarantee those producers
have finished; arrivals after input closure cannot be accepted. Increasing the
settling interval is not a completeness guarantee for those setups.

Synthetic tests cover directory resolution, identity/format checks, settling,
late BAMs, shared barcode runs, independent completion, manual STOP, and opt-out.
They cannot establish compatibility across MinKNOW versions. Before relying on
this feature with a particular installation, compare a real completed run's final
summary with its final BAM inventory and verify all expected BAMs appear in seqLM's
accepted input logs. No completed-run MinKNOW dataset was available for this
implementation's validation.
