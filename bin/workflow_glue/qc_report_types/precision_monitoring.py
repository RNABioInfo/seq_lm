"""Shared precision stopping tables and plots for DEA and ICA."""

import json
from pathlib import Path

import pandas as pd
from bokeh.models import Span
from dominate.tags import p as paragraph
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import BokehPlot
from precision_termination.engine import (
    DEFAULTS,
    FORECAST_METHOD,
    minimum_confirmation_span_seconds,
)
from precision_termination.reasons import audit_reason

REASON_LABELS = {
    "precision_met": "Precision check passed",
    "improvement_continues": "More reads may improve precision",
    "stability_not_confirmed": "Waiting for a stable uncertainty trend",
    "waiting_for_reads": "Waiting for more sequencing data",
    "results_not_ready": "Waiting for usable analysis results",
    "not_monitored": "No stop check applies",
}


def reason_label(reason: str) -> str:
    return ". ".join(
        REASON_LABELS.get(code, code.replace("_", " ").capitalize())
        for code in audit_reason(reason).split(",")
    )


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
    result["Reason"] = result["Reason"].map(reason_label)
    result.attrs["precision_monitoring_directory"] = str(folder)
    return result


def add_precision_monitoring(table: pd.DataFrame | None, *, method: str, behavior: str) -> None:
    ica = method == "ica"
    features = "components" if ica else "features"
    analysis = "iModulon activity differences (ICA)" if ica else "gene expression differences (DEA)"
    paragraph(f"Stopping checks use {analysis}.")
    paragraph(
        """
        Automatic termination evaluates whether continued sequencing yields improvement in analysis uncertainty.
        This table shows the current evaluation status of the different comparisons.
        A lower score means a more precise analysis result.
        """
    )
    if behavior == "log":
        paragraph("Behavior is log only.")
    else:
        paragraph("Behavior is automatic stopping once all active samples sharing the acquisition meet the cutoff.")

    if table is None:
        paragraph("No stopping checks have been recorded yet.")
        return

    DataTable.from_pandas(table, use_index=False)
    root = Path(table.attrs["precision_monitoring_directory"])
    audit = json.loads((root / "audit.json").read_text())
    config = DEFAULTS | audit["config"]
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
    )
    horizon_minutes = config["horizon_seconds"] / 60
    gain_cutoff = 100 * config["gain_cutoff"]
    required_span = minimum_confirmation_span_seconds(config["horizon_seconds"]) / 60
    paragraph(
        f"A precision check passes when both the decrease in uncertainty over the past {horizon_minutes:g} minutes "
        f"and the forecast with its error allowance for the next {horizon_minutes:g} minutes are below {gain_cutoff:g}%. "
        f"The stopping rule requires at least {config['confirmations']} passing checks in a row "
        f"spanning at least {required_span:g} minutes. See Action for sample outcomes."
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
            "confirmation_span_seconds",
            "Time covered by consecutive passing checks",
            "Observation time (hours)",
            "Time span (minutes)",
            3600.0,
        ),
    ]
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
    line_styles = {
        "R_information": [("fit_R", "dashed", "Fitted uncertainty")],
        "predicted_conservative_gain": [
            ("curve_gain", "dotted", "Forecast without margin"),
            ("actual_gain_next_horizon", "dashed", "Observed improvement"),
        ],
    }
    cutoffs = {
        "predicted_conservative_gain": (gain_cutoff, f"Stopping cutoff ({gain_cutoff:g}%)"),
        "confirmation_span_seconds": (required_span, f"Required passing span ({required_span:g} min)"),
    }
    for x, y, title, xlabel, ylabel, scale in plots:
        if y not in history:
            continue
        p = BokehPlot(
            title=title,
            x_axis_label=xlabel,
            y_axis_label=ylabel,
        )
        if y in {"R_information", "confirmation_span_seconds"}:
            p._fig.y_range.start = 0
        factor = {"predicted_conservative_gain": 100, "confirmation_span_seconds": 1 / 60}.get(y, 1)
        extra_lines = line_styles.get(y, [])

        for i, (cid, rows) in enumerate(history.groupby("contrast_id", sort=False)):
            xx = pd.to_numeric(rows[x], errors="coerce") / scale
            yy = pd.to_numeric(rows[y], errors="coerce") * factor
            color = colors[i % len(colors)]
            p._fig.line(
                xx,
                yy,
                color=color,
                legend_label=labels.get(cid, str(cid)),
                line_width=2,
            )
            for field, dash, _ in extra_lines:
                if field not in rows:
                    continue
                p._fig.line(
                    xx,
                    pd.to_numeric(rows[field], errors="coerce") * factor,
                    color=color,
                    line_dash=dash,
                )
            if "eligible" in rows:
                eligible = rows.eligible.astype(str).str.lower().eq("true")
                p._fig.scatter(
                    xx[eligible],
                    yy[eligible],
                    marker="circle",
                    size=8,
                    color=color,
                )

        # Bokeh needs a point for each legend swatch. NaNs keep the swatches off the plot.
        legend_point = [float("nan")]
        if y == "predicted_conservative_gain":
            p._fig.line(
                legend_point, legend_point, color="#333333", line_width=2,
                legend_label="Forecast with margin",
            )
        for _, dash, label in extra_lines:
            p._fig.line(
                legend_point, legend_point, color="#333333", line_dash=dash,
                legend_label=label,
            )
        if y in cutoffs:
            cutoff, label = cutoffs[y]
            cutoff_style = dict(line_dash=[8, 4], line_color="#777777")
            p._fig.add_layout(
                Span(location=cutoff, dimension="width", **cutoff_style)
            )
            p._fig.line(legend_point, legend_point, legend_label=label, **cutoff_style)
            # Spans do not affect the automatic range. Include the cutoff even at the first check.
            elapsed = history[x] / scale
            p._fig.scatter([elapsed.min(), elapsed.max()], [cutoff, cutoff], alpha=0)

        EZChart(p, "epi2melabs", height="400px")

    if rates.empty or "last_interval_reads_per_minute" not in rates:
        return

    plot = BokehPlot(
        title="Observed throughput",
        x_axis_label="Observation time (hours)",
        y_axis_label="Assigned reads per minute",
    )
    for i, ((cid, sample), rows) in enumerate(
        rates.groupby(["contrast_id", "sample"], sort=False)
    ):
        plot._fig.line(
            rows.elapsed_seconds / 3600,
            pd.to_numeric(rows.last_interval_reads_per_minute, errors="coerce"),
            color=colors[i % len(colors)],
            legend_label=f"{labels.get(cid, cid)}: {sample}",
        )

    EZChart(plot, "epi2melabs", height="400px")
