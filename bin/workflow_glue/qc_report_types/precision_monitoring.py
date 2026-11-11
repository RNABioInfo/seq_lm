"""Shared precision stopping tables and plots for DEA and ICA."""

import json
from pathlib import Path

import pandas as pd
from bokeh.models import Span
from dominate.tags import p as paragraph
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import BokehPlot
from precision_termination.engine import FORECAST_METHOD
from precision_termination.forecast import FIT_HALF_LIFE_SECONDS


def load_precision_monitoring(
    path: str | Path, behavior: str,
    report_samples=None, *, method="dea", expected_snapshot=None,
) -> pd.DataFrame:
    folder = Path(path)
    methods = json.loads((folder / "audit.json").read_text())
    if methods.get("forecast_method") not in {"depth_curve", FORECAST_METHOD}:
        raise ValueError("Unrecognized precision forecast in monitoring audit")
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
        raise ValueError("Stopping audit is missing required columns")
    if set(audit.method) != {method} or set(audit.metric) != {"R_information"}:
        raise ValueError("Unexpected stopping audit method or metric")
    if expected_snapshot is not None:
        for column, expected in zip(("batch_index", "analysis_index"), expected_snapshot, strict=True):
            if column not in audit or set(audit[column]) != {str(expected)}:
                raise ValueError("ICA termination audit does not match the analysis snapshot")
    if set(audit.behavior) != {behavior}:
        raise ValueError("Stopping audit behavior mismatch")
    pairs = list(zip(audit["group"], audit["sample"]))
    if len(set(pairs)) != len(pairs):
        raise ValueError("Duplicate stopping audit sample identities")
    if report_samples is not None and set(pairs) != set(
        zip(report_samples["Group"], report_samples["Name"])
    ):
        raise ValueError("Stopping audit samples do not match report")
    if not set(audit.eligible).issubset({"true", "false"}):
        raise ValueError("Invalid stopping audit eligibility")
    result = audit.rename(
        columns={
            "sample": "Sample",
            "group": "Group",
            "consecutive_confirmations": "Checks passed in a row",
            "eligible": "Stopping rule met",
            "reason": "Reason",
            "action_result": "Action",
        }
    )
    result = result[
        ["Sample", "Group", "Checks passed in a row", "Stopping rule met", "Reason", "Action"]
    ].copy()
    result["Stopping rule met"] = result["Stopping rule met"].map({"true": "Yes", "false": "No"})
    result.attrs["precision_monitoring_directory"] = str(folder)
    return result


def add_precision_monitoring(table: pd.DataFrame | None, *, method: str, behavior: str) -> None:
    ica = method == "ica"
    quantity = "iModulon activity differences" if ica else "gene expression differences"
    features = "components" if ica else "features"
    paragraph(
        "This tab shows how much collecting more reads is expected to improve the precision of the results, "
        "whether the stopping rule has been met, and what action was taken. "
        f"The uncertainty score summarizes uncertainty in {quantity} between groups. "
        "A lower score means a more precise estimate."
    )
    if behavior == "log":
        paragraph("This run records stopping checks for review.")
    else:
        paragraph("SeqLM can request a stop once all active samples sharing the acquisition meet the stopping rule.")

    if table is None:
        paragraph("No stopping checks have been recorded yet.")
        return

    DataTable.from_pandas(table, use_index=False)
    root = Path(table.attrs["precision_monitoring_directory"])
    audit = json.loads((root / "audit.json").read_text())
    config = audit["config"]
    contrasts = audit.get("contrasts", [])
    selected = any(row.get("uncertainty_scope") == "selected_features" for row in contrasts)
    labels = {
        row["contrast_id"]: f"{row['target_group']} vs {row['reference_group']}"
        for row in contrasts if row.get("target_group") is not None
    }
    counts = ", ".join(
        f"{labels.get(row['contrast_id'], row['contrast_id'])}: {row['uncertainty_feature_count']} {features}"
        for row in contrasts if row.get("uncertainty_feature_count") is not None
    )
    all_features = "all model components" if ica else "all quantified features"
    selected_features = "selected iModulon components" if ica else "features in the selected gene sets"
    paragraph(
        "The uncertainty score uses "
        + (selected_features if selected else all_features)
        + (f" ({counts})" if counts else "")
        + f". The total also includes {features} with unavailable estimates."
    )
    horizon_minutes = config.get("horizon_seconds", 1800) / 60
    if audit.get("forecast_method") == FORECAST_METHOD:
        paragraph(
            "The curve uses raw uncertainty scores. It gives newer checkpoints more weight, "
            f"with weights halving every {FIT_HALF_LIFE_SECONDS / 3600:g} hours. "
            "Large misses have less influence on the fit."
        )
    paragraph(
        f"The stopping cutoff is a {100 * config['gain_cutoff']:g}% drop in uncertainty over "
        f"{horizon_minutes:g} minutes. The rule looks for improvements smaller than this cutoff "
        "in both the previous window and the forecast. The horizontal line marks the cutoff. "
        "The solid line shows the forecast used for the stopping check, including an allowance for forecast error. "
        "The dotted line shows the forecast before that allowance. "
        "The dashed line shows the improvement measured once later data are available. "
        "Circles mark checks that met the stopping rule. See the Action column for what happened to each sample."
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
            "Uncertainty score",
            3600.0,
        ),
        (
            "assigned_monitored_reads",
            "R_information",
            "Uncertainty versus read depth",
            "Assigned reads in monitored libraries (millions)",
            "Uncertainty score",
            1e6,
        ),
        (
            "elapsed_seconds",
            "predicted_conservative_gain",
            f"Next {horizon_minutes:g}-minute improvement",
            "Observation time (hours)",
            "Relative improvement (%)",
            3600.0,
        ),
        (
            "elapsed_seconds",
            "fraction_R_le_05",
            f"{features.capitalize()} with uncertainty score of 0.5 or less",
            "Observation time (hours)",
            f"Fraction of {features}",
            3600.0,
        ),
        (
            "elapsed_seconds",
            "fraction_R_le_1",
            f"{features.capitalize()} with uncertainty score of 1 or less",
            "Observation time (hours)",
            f"Fraction of {features}",
            3600.0,
        ),
        (
            "elapsed_seconds",
            "maximum_sample_data_age_seconds",
            "Age of the data used in the check",
            "Observation time (hours)",
            "Age (seconds)",
            3600.0,
        ),
        (
            "elapsed_seconds",
            "processing_lag_seconds",
            "Delay from newest data to results",
            "Observation time (hours)",
            "Processing lag (seconds)",
            3600.0,
        ),
        (
            "elapsed_seconds",
            "confirmation_span_seconds",
            "Time covered by consecutive passing checks",
            "Observation time (hours)",
            "Time span (seconds)",
            3600.0,
        ),
    ]
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    for x, y, title, xlabel, ylabel, scale in plots:
        if y not in history:
            continue
        p = BokehPlot(
            title=title,
            x_axis_label=xlabel,
            y_axis_label=ylabel,
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
                legend_label=labels.get(cid, str(cid)),
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
                legend_label=f"{labels.get(cid, cid)}: {sample}",
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
                legend_label=f"{labels.get(cid, cid)}: {sample}",
            )

        EZChart(plot, "epi2melabs", height="400px")
