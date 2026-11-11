# Automated stopping for live sequencing

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

seq_lm estimates the reduction in uncertainty expected from additional
sequencing. The stopping rule can use differential expression (`dea`) or
iModulon activity (`ica`). Enable the corresponding analysis in the run
configuration.
[Differential expression](differential_expression.md), [iModulon analysis](imodulon_analysis.md)

## Monitoring settings

Add the following settings to the `params` block of the run configuration.
In `log` mode, the report indicates when the stopping criteria are met
while sequencing continues.

```groovy
live_analysis = true
stopping_method = "dea"
monitoring_behavior = "log"
```

Monitoring is disabled by default (`monitoring_behavior = "disabled"`).
With the default settings, both the predicted reduction in overall
uncertainty over the next 30 minutes and the recent observed reduction
must be below 1%. Individual estimates can remain uncertain even when
the stopping criteria are met.

Use `stopping_gain_cutoff` to change the threshold and
`stopping_forecast_minutes` to change the forecast period. Set the forecast
period to at least the typical interval between BAM files. For hourly BAM
files, use at least 60 minutes.

seq_lm requires at least six informative updates and one hour of data. With
the default settings, three consecutive updates must meet the stopping
criteria over a span of at least 20 minutes. Processing delays, outdated
inputs, or interruptions in BAM delivery can prevent a stop.

## Feature selection

By default, the stopping rule uses all features in the selected analysis.
Use the following settings to limit the stopping rule to a subset.

- For `dea`, set `dea_stopping_genes` to a GMT file. The decision uses the
  union of matched features across all gene sets in the file. This setting
  is independent of gene-set enrichment.
- For `ica`, set `ica_stopping_components` to a text file with one exact
  model component ID per line.

The stopping rule uses the selected subset. The report presents the full
analysis results.

## Stop sequencing automatically

Run the certificate setup on the computer running MinKNOW with the seq_lm
Conda environment active. Specify the MinKNOW client certificate directory.

```bash
bin/seq-run-manager cert --minknow-client-certs-directory /path/to/minknow/conf/rpc-client-certs
```

In the `params` block, enter the paths that the command reports for the
client certificate, private key, and MinKNOW CA certificate. Set
`monitoring_behavior` to `terminate`.

```groovy
monitoring_behavior = "terminate"
minknow_client_certificate = "/path/to/minknow_cert.pem"
minknow_client_private_key = "/path/to/minknow_key.pem"
minknow_ca_certificate = "/path/to/minknow_cert.crt"
```

The workflow connects to `host.docker.internal` on port `9501` by default.
Set `minknow_host` and `minknow_port` if MinKNOW is at a different address.
Use an address reachable from your container.

Each live sample's BAM directory must have exactly one MinKNOW
`*sample_sheet*.csv` in its parent directory. The sheet must have one row
whose `sample_id` matches the seq_lm `alias` and whose `protocol_run_id`
identifies the acquisition. All live samples in the same acquisition must
meet the stopping criteria before seq_lm can stop it. After MinKNOW confirms
the stop, seq_lm writes `STOP` markers and processes the remaining BAM files.

## Monitoring results

The **DEA Termination** or **ICA Termination** tab presents the forecast,
the status of the stopping criteria, and any action taken. Detailed results
are available in `<out_dir>/termination/dea/` or
`<out_dir>/termination/ica/`.

If seq_lm reports an unsupported monitoring state on restart, use a new
output directory.
