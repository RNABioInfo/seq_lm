"""Current-format R monitoring tables and plots; no legacy DEA audit reader."""

import json
from pathlib import Path

import pandas as pd
from bokeh.models import Span
from dominate.tags import p as paragraph
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import BokehPlot


def load_dea_monitoring(path: str | None, metadata: pd.DataFrame, behavior: str, report_samples=None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(
            {
                "Sample": metadata.index,
                "Status": "R monitoring disabled or not yet available",
            }
        )
    folder = Path(path)
    methods = json.loads((folder / "audit.json").read_text())
    if methods.get("forecast_method") != "depth_curve":
        raise ValueError("DEA monitoring requires current curve-only audit data")
    audit = pd.read_csv(
        folder / "sample_monitoring.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    required = {
        "group",
        "sample",
        "method",
        "metric",
        "eligible",
        "consecutive_confirmations",
        "reason",
        "behavior",
        "action_result",
    }
    if not required.issubset(audit):
        raise ValueError("R monitoring audit is missing required columns")
    if set(audit.method) != {"dea"} or set(audit.metric) != {"R_information"}:
        raise ValueError("Unexpected R monitoring method")
    if set(audit.behavior) != {behavior}:
        raise ValueError("R monitoring behavior mismatch")
    pairs = list(zip(audit["group"], audit["sample"]))
    if len(set(pairs)) != len(pairs):
        raise ValueError("Duplicate R monitoring sample identities")
    if report_samples is not None and set(pairs) != set(
        zip(report_samples["Group"], report_samples["Name"])
    ):
        raise ValueError("R monitoring samples do not match report")
    if not set(audit.eligible).issubset({"true", "false"}):
        raise ValueError("Invalid R monitoring eligibility")
    result = audit.rename(
        columns={
            "sample": "Sample",
            "group": "Group",
            "consecutive_confirmations": "Consecutive confirmations",
            "eligible": "Eligible",
            "reason": "Reason",
            "action_result": "Action",
        }
    )
    result = result[
        ["Sample", "Group", "Consecutive confirmations", "Eligible", "Reason", "Action"]
    ]
    result.attrs["dea_monitoring_directory"] = str(folder)
    return result


def add_dea_monitoring(table):
    paragraph(
        "R summarizes model-based uncertainty in log2 fold changes. Smaller values indicate narrower intervals. The aggregate and forecast are custom heuristics; they do not guarantee accuracy for every gene or sequential confidence coverage."
    )
    DataTable.from_pandas(table, use_index=False)
    directory = table.attrs.get("dea_monitoring_directory")
    if not directory:
        return
    root = Path(directory)
    config = (
        json.loads((root / "audit.json").read_text())["config"]
        if (root / "audit.json").exists()
        else {"gain_cutoff": 0.01}
    )
    horizon_minutes = config.get("horizon_seconds", 1800) / 60
    paragraph(
        f"Stopping cutoff: {100 * config['gain_cutoff']:g}% expected reduction in aggregate R over {horizon_minutes:g} minutes. The prediction uses the fitted depth curve plus its calibrated error margin. Solid forecast lines show their sum; dotted lines show the raw curve forecast; dashed lines show later observed improvement where available. Circles mark eligibility, which does not imply that a shared acquisition was stopped."
    )
    history = pd.read_csv(root / "forecast_history.tsv", sep="\t")
    rates = pd.read_csv(root / "read_rate_history.tsv", sep="\t")
    if history.empty:
        return
    plots = [
        (
            "elapsed_seconds",
            "R_information",
            "Uncertainty over observation time",
            "Observation time (hours)",
            3600.0,
        ),
        (
            "assigned_monitored_reads",
            "R_information",
            "Uncertainty versus read depth",
            "Assigned reads in monitored libraries (millions)",
            1e6,
        ),
        (
            "elapsed_seconds",
            "predicted_conservative_gain",
            f"Next {horizon_minutes:g}-minute improvement",
            "Observation time (hours)",
            3600.0,
        ),
        (
            "elapsed_seconds",
            "fraction_R_le_05",
            "Fraction of genes with R ≤ 0.5",
            "Observation time (hours)",
            3600.0,
        ),
        (
            "elapsed_seconds",
            "fraction_R_le_1",
            "Fraction of genes with R ≤ 1",
            "Observation time (hours)",
            3600.0,
        ),
    ]
    plots.extend(
        [
            (
                "elapsed_seconds",
                "maximum_sample_data_age_seconds",
                "Age of represented data",
                "Observation time (hours)",
                3600.0,
            ),
            (
                "elapsed_seconds",
                "processing_lag_seconds",
                "Processing lag (seconds)",
                "Observation time (hours)",
                3600.0,
            ),
            (
                "elapsed_seconds",
                "confirmation_span_seconds",
                "Confirmation span (seconds)",
                "Observation time (hours)",
                3600.0,
            ),
        ]
    )
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    for x, y, title, xlabel, scale in plots:
        if y not in history:
            continue
        p = BokehPlot(
            title=title,
            x_axis_label=xlabel,
            y_axis_label="Relative improvement (%)" if "gain" in y else y,
        )
        if y == "R_information":
            p._fig.y_range.start = 0
        for i, (cid, rows) in enumerate(history.groupby("contrast_id", sort=False)):
            xx = pd.to_numeric(rows[x], errors="coerce") / scale
            yy = pd.to_numeric(rows[y], errors="coerce")
            factor = 100 if "gain" in y else 1
            p._fig.line(
                xx,
                yy * factor,
                color=colors[i % len(colors)],
                legend_label=str(cid),
                line_width=2,
            )
            if y == "R_information" and "fit_R" in rows:
                p._fig.line(
                    xx,
                    pd.to_numeric(rows.fit_R, errors="coerce"),
                    color=colors[i % len(colors)],
                    line_dash="dashed",
                )
            if "gain" in y and "actual_gain_next_horizon" in rows:
                p._fig.line(
                    xx,
                    pd.to_numeric(rows.actual_gain_next_horizon, errors="coerce") * 100,
                    color=colors[i % len(colors)],
                    line_dash="dashed",
                )
            if y == "predicted_conservative_gain":
                p._fig.line(
                    xx,
                    pd.to_numeric(rows.curve_gain, errors="coerce") * 100,
                    color=colors[i % len(colors)],
                    line_dash="dotted",
                )
            if "eligible" in rows:
                selected = rows.eligible.astype(str).str.lower().eq("true")
                p._fig.scatter(
                    xx[selected],
                    yy[selected] * factor,
                    marker="circle",
                    size=8,
                    color=colors[i % len(colors)],
                )
        if "gain" in y:
            p._fig.add_layout(
                Span(
                    location=100 * config["gain_cutoff"],
                    dimension="width",
                    line_dash="dotted",
                    line_color="#444444",
                )
            )
        EZChart(p, "epi2melabs", height="400px")
        
    if not rates.empty and "predicted_increment" in rates:
        p = BokehPlot(
            title=f"Forecast reads in the next {horizon_minutes:g} minutes",
            x_axis_label="Observation time (hours)",
            y_axis_label="Assigned reads",
        )
        for i, ((cid, sample), rows) in enumerate(
            rates.groupby(["contrast_id", "sample"], sort=False)
        ):
            p._fig.line(
                rows.elapsed_seconds / 3600,
                pd.to_numeric(rows.predicted_increment, errors="coerce"),
                color=colors[i % len(colors)],
                legend_label=f"{cid}: {sample}",
            )

        EZChart(p, "epi2melabs", height="400px")

    for field, title, ylabel in [
        (
            "last_interval_reads_per_minute",
            "Observed throughput",
            "Assigned reads per minute",
        ),
        ("last_interval_seconds", "Read-arrival observation gaps", "Seconds"),
        ("rate_coverage_seconds", "Forecast-window coverage", "Seconds"),
    ]:
        if rates.empty or field not in rates:
            continue
        plot = BokehPlot(
            title=title, x_axis_label="Observation time (hours)", y_axis_label=ylabel
        )
        for i, ((cid, sample), rows) in enumerate(
            rates.groupby(["contrast_id", "sample"], sort=False)
        ):
            plot._fig.line(
                rows.elapsed_seconds / 3600,
                pd.to_numeric(rows[field], errors="coerce"),
                color=colors[i % len(colors)],
                legend_label=f"{cid}: {sample}",
            )

        EZChart(plot, "epi2melabs", height="400px")
