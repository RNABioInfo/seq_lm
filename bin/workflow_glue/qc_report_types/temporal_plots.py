"""Validate temporal metadata and add descriptive gene-set time-course plots."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bokeh.layouts import column
from bokeh.models import (
    ColorBar,
    ColumnDataSource,
    Div,
    HoverTool,
    LinearColorMapper,
)
import colorcet as cc
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Tabs
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd

from .differential_plots import (
    DifferentialResult,
    _display_labels,
    cluster_heatmap_rows,
)
from .gsva_plots import GSVAResult


GENE_SET_RESOLUTION_FILE = "gene_set_resolution.tsv"
TEMPORAL_CAVEAT = "Descriptive only. Heatmap colors are gene-wise z-scores."


@dataclass
class TemporalResult:
    """Validated inputs for one temporal trajectory across all sample groups."""

    metadata: pd.DataFrame
    scores_long: pd.DataFrame
    log_cpm: pd.DataFrame
    members: dict[str, tuple[str, ...]]
    feature_labels: dict[str, str]
    gene_set_order: list[str]


def _temporal_metadata(
    samples_df: pd.DataFrame,
    differential: DifferentialResult,
) -> pd.DataFrame:
    """Validate sample times and return metadata in elapsed-minute order."""
    if {"Name", "Group", "Time (min)"}.issubset(samples_df.columns):
        metadata = samples_df[["Name", "Group", "Time (min)"]].rename(
            columns={
                "Name": "sample",
                "Group": "group",
                "Time (min)": "time_minutes",
            }
        )
    elif {"name", "group", "order"}.issubset(samples_df.columns):
        metadata = samples_df[["name", "group", "order"]].rename(
            columns={"name": "sample", "order": "time_minutes"}
        )
    else:
        raise ValueError(
            "Temporal analysis requires sample, group, and order metadata."
        )

    metadata = metadata.copy()
    for column in ("sample", "group"):
        if metadata[column].isna().any():
            raise ValueError(f"Temporal analysis contains missing {column} values.")
        metadata[column] = metadata[column].astype(str).str.strip()
        if metadata[column].eq("").any():
            raise ValueError(f"Temporal analysis contains empty {column} values.")
    if metadata["sample"].duplicated().any():
        raise ValueError("Temporal analysis requires globally unique sample names.")

    time_text = metadata["time_minutes"].fillna("").astype(str).str.strip()
    invalid_time = ~time_text.str.fullmatch(r"[+-]?\d+")
    if invalid_time.any():
        invalid_samples = metadata.loc[invalid_time, "sample"].tolist()
        raise ValueError(
            "Temporal analysis requires signed integer order values in elapsed "
            "minutes; invalid sample(s): " + ", ".join(invalid_samples)
        )
    metadata["time_minutes"] = time_text.astype(int)

    expected_samples = differential.sample_metadata.index.tolist()
    if set(metadata["sample"]) != set(expected_samples):
        raise ValueError(
            "Temporal sample metadata does not match differential metadata."
        )
    metadata = metadata.set_index("sample", drop=False).loc[expected_samples]
    expected_groups = differential.sample_metadata.loc[expected_samples, "group"]
    if metadata["group"].tolist() != expected_groups.tolist():
        raise ValueError(
            "Temporal sample groups do not match differential metadata."
        )

    inconsistent_groups = (
        metadata.groupby("group", sort=False)["time_minutes"]
        .nunique()
        .loc[lambda values: values.ne(1)]
        .index.tolist()
    )
    if inconsistent_groups:
        raise ValueError(
            "Temporal analysis requires every group to use one elapsed minute; "
            "inconsistent group(s): " + ", ".join(inconsistent_groups)
        )
    ambiguous_times = (
        metadata.groupby("time_minutes", sort=False)["group"]
        .nunique()
        .loc[lambda values: values.ne(1)]
        .index.tolist()
    )
    if ambiguous_times:
        raise ValueError(
            "Temporal analysis requires every elapsed minute to identify one "
            "group; ambiguous minute(s): "
            + ", ".join(str(value) for value in sorted(ambiguous_times))
        )
    if metadata["time_minutes"].nunique() < 2:
        raise ValueError(
            "Temporal analysis requires at least two distinct elapsed minutes."
        )

    metadata["sample_order"] = np.arange(len(metadata))
    return metadata.sort_values(
        ["time_minutes", "sample_order"],
        kind="mergesort",
    )


def _read_scored_members(
    results_dir: Path,
    differential: DifferentialResult,
    gsva: GSVAResult,
) -> tuple[pd.DataFrame, dict[str, tuple[str, ...]]]:
    """Reconstruct the retained, variable members used by GSVA."""
    resolution_path = results_dir / GENE_SET_RESOLUTION_FILE
    if not resolution_path.is_file():
        raise ValueError(
            f"Missing temporal gene-set resolution table: {resolution_path}"
        )
    resolution = pd.read_csv(resolution_path, sep="\t")
    missing = [
        column for column in ("gene_set", "feature_id") if column not in resolution
    ]
    if missing:
        raise ValueError(
            f"{resolution_path} is missing columns: " + ", ".join(missing)
        )
    if resolution["gene_set"].isna().any():
        raise ValueError(f"{resolution_path} contains missing gene_set values.")
    resolution["gene_set"] = resolution["gene_set"].astype(str).str.strip()
    if resolution["gene_set"].eq("").any():
        raise ValueError(f"{resolution_path} contains empty gene_set values.")
    resolution["feature_id"] = resolution["feature_id"].fillna("").astype(str).str.strip()
    resolution = resolution.loc[resolution["feature_id"].ne("")].drop_duplicates(
        ["gene_set", "feature_id"]
    )

    log_cpm = np.log2(differential.feature_counts + 1)
    variable_features = set(
        log_cpm.index[log_cpm.max(axis=1).gt(log_cpm.min(axis=1))]
    )
    score_counts = (
        gsva.scores_long.drop_duplicates("gene_set")
        .set_index("gene_set")["n_genes"]
        .astype(int)
    )
    coverage_counts = gsva.coverage.set_index("gene_set")["scored_members"]
    members = {}
    for gene_set in gsva.gene_set_order:
        resolved = resolution.loc[
            resolution["gene_set"].eq(gene_set), "feature_id"
        ].tolist()
        scored = tuple(
            feature_id for feature_id in resolved if feature_id in variable_features
        )
        expected = int(score_counts.at[gene_set])
        coverage_expected = int(coverage_counts.at[gene_set])
        if len(scored) != expected or len(scored) != coverage_expected:
            raise ValueError(
                f"Gene set '{gene_set}' has {len(scored)} reconstructed temporal "
                f"members, but GSVA reports {expected} scored members and coverage "
                f"reports {coverage_expected}."
            )
        members[gene_set] = scored
    return log_cpm, members


def load_temporal_results(
    results_dir: str | Path,
    samples_df: pd.DataFrame,
    differential: DifferentialResult,
    gsva: GSVAResult,
) -> TemporalResult:
    """Load and cross-validate one elapsed-minute gene-set trajectory."""
    metadata = _temporal_metadata(samples_df, differential)
    log_cpm, members = _read_scored_members(
        Path(results_dir),
        differential,
        gsva,
    )
    labels = _display_labels(differential.contrasts[0].results)
    feature_labels = dict(
        zip(
            differential.contrasts[0].results["feature_id"].astype(str),
            labels.astype(str),
            strict=True,
        )
    )
    scores = gsva.scores_long.merge(
        metadata[["sample", "time_minutes", "sample_order"]].reset_index(
            drop=True
        ),
        on="sample",
        how="left",
        validate="many_to_one",
    ).sort_values(
        ["gene_set", "time_minutes", "sample_order"],
        kind="mergesort",
    )
    return TemporalResult(
        metadata=metadata,
        scores_long=scores,
        log_cpm=log_cpm,
        members=members,
        feature_labels=feature_labels,
        gene_set_order=gsva.gene_set_order,
    )


def prepare_temporal_scores(
    data: TemporalResult,
    gene_set: str,
    condition_colors: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return time-point score summaries and deterministically jittered samples."""
    subset = data.scores_long.loc[
        data.scores_long["gene_set"].eq(gene_set)
    ].copy()
    if subset.empty:
        raise ValueError(f"Unknown temporal GSVA gene set: {gene_set}")
    summary = (
        subset.groupby(["time_minutes", "group"], sort=True)["score"]
        .agg(mean_score="mean", sd_score="std", n="size")
        .reset_index()
        .sort_values("time_minutes", kind="mergesort")
    )
    summary["lower"] = summary["mean_score"] - summary["sd_score"]
    summary["upper"] = summary["mean_score"] + summary["sd_score"]
    summary["sd_label"] = summary["sd_score"].map(
        lambda value: "unavailable (n=1)" if pd.isna(value) else f"{value:.3f}"
    )
    summary["plot_color"] = summary["group"].map(condition_colors)
    if summary["plot_color"].isna().any():
        raise ValueError("Temporal GSVA groups are missing condition colors.")

    raw = subset.merge(
        summary[
            ["time_minutes", "group", "mean_score", "sd_score", "sd_label", "n"]
        ],
        on=["time_minutes", "group"],
        how="left",
        validate="many_to_one",
    )
    times = np.sort(summary["time_minutes"].unique())
    minimum_gap = float(np.diff(times).min())
    jitter_width = minimum_gap * 0.06
    raw["plot_time"] = raw["time_minutes"].astype(float)
    for time_minutes, indices in raw.groupby("time_minutes", sort=True).groups.items():
        ordered = raw.loc[indices].sort_values("sample_order", kind="mergesort").index
        offsets = (
            np.linspace(-jitter_width, jitter_width, len(ordered))
            if len(ordered) > 1
            else np.array([0.0])
        )
        raw.loc[ordered, "plot_time"] = float(time_minutes) + offsets
    raw["plot_color"] = raw["group"].map(condition_colors)
    return summary, raw


def _add_error_bars(plot, source: ColumnDataSource, color: str):
    """Draw SD stems and caps for rows with a finite interval."""
    stem = plot.segment(
        x0="time_minutes",
        x1="time_minutes",
        y0="lower",
        y1="upper",
        source=source,
        line_color=color,
        line_alpha=0.55,
    )
    lower_cap = plot.segment(
        x0="cap_left",
        x1="cap_right",
        y0="lower",
        y1="lower",
        source=source,
        line_color=color,
        line_alpha=0.55,
    )
    upper_cap = plot.segment(
        x0="cap_left",
        x1="cap_right",
        y0="upper",
        y1="upper",
        source=source,
        line_color=color,
        line_alpha=0.55,
    )
    return stem, lower_cap, upper_cap


def create_temporal_score_plot(
    data: TemporalResult,
    gene_set: str,
    condition_colors: dict[str, str],
) -> BokehPlot:
    """Create mean GSVA score, SD whiskers, and raw samples over time."""
    summary, raw = prepare_temporal_scores(data, gene_set, condition_colors)
    plot = BokehPlot(
        title=f"GSVA score over time — {gene_set}",
        x_axis_label="Elapsed time (min)",
        y_axis_label="GSVA score",
        height=470,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    summary_source = ColumnDataSource(summary)
    plot._fig.line(
        "time_minutes",
        "mean_score",
        source=summary_source,
        line_color="#333333",
        line_width=2.5,
    )
    means = plot._fig.scatter(
        "time_minutes",
        "mean_score",
        source=summary_source,
        size=10,
        fill_color="plot_color",
        line_color="#333333",
    )
    finite = summary.loc[summary["sd_score"].notna()].copy()
    times = np.sort(summary["time_minutes"].unique())
    cap_width = float(np.diff(times).min()) * 0.025
    finite["cap_left"] = finite["time_minutes"] - cap_width
    finite["cap_right"] = finite["time_minutes"] + cap_width
    _add_error_bars(
        plot._fig,
        ColumnDataSource(finite),
        color="#333333",
    )

    raw_source = ColumnDataSource(raw)
    samples = plot._fig.scatter(
        "plot_time",
        "score",
        source=raw_source,
        marker="circle",
        size=8,
        fill_color="plot_color",
        fill_alpha=0.75,
        line_color="white",
        line_width=1,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[means],
            tooltips=[
                ("Group", "@group"),
                ("Time", "@time_minutes min"),
                ("Mean score", "@mean_score{0.000}"),
                ("SD", "@sd_label"),
                ("Replicates", "@n"),
            ],
        ),
        HoverTool(
            renderers=[samples],
            tooltips=[
                ("Sample", "@sample"),
                ("Group", "@group"),
                ("Time", "@time_minutes min"),
                ("Score", "@score{0.000}"),
                ("Mean score", "@mean_score{0.000}"),
                ("SD", "@sd_label"),
                ("Replicates", "@n"),
            ],
        ),
    )
    plot._fig.grid.grid_line_alpha = 0.15
    plot.report_height = 500
    return plot


def prepare_temporal_gene_expression(
    data: TemporalResult,
    gene_set: str,
) -> pd.DataFrame:
    """Summarize and row-standardize each member's temporal log-CPM profile."""
    if gene_set not in data.members:
        raise ValueError(f"Unknown temporal gene set: {gene_set}")
    rows = []
    for feature_id in data.members[gene_set]:
        for time_minutes, time_rows in data.metadata.groupby(
            "time_minutes", sort=True
        ):
            samples = time_rows["sample"].tolist()
            values = data.log_cpm.loc[feature_id, samples].to_numpy(dtype=float)
            standard_deviation = (
                float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
            )
            rows.append(
                {
                    "feature_id": feature_id,
                    "display_label": data.feature_labels.get(feature_id, feature_id),
                    "time_minutes": int(time_minutes),
                    "group": time_rows["group"].iloc[0],
                    "mean_log_cpm": float(np.mean(values)),
                    "sd_log_cpm": standard_deviation,
                    "n": len(values),
                }
            )
    summary = pd.DataFrame(rows)
    summary["lower"] = summary["mean_log_cpm"] - summary["sd_log_cpm"]
    summary["upper"] = summary["mean_log_cpm"] + summary["sd_log_cpm"]
    summary["sd_label"] = summary["sd_log_cpm"].map(
        lambda value: "unavailable (n=1)" if pd.isna(value) else f"{value:.3f}"
    )
    gene_means = summary.groupby("feature_id", sort=False)["mean_log_cpm"].transform(
        "mean"
    )
    gene_standard_deviations = summary.groupby("feature_id", sort=False)[
        "mean_log_cpm"
    ].transform("std")
    summary["z_score"] = (
        summary["mean_log_cpm"]
        .sub(gene_means)
        .div(gene_standard_deviations.replace(0, np.nan))
        .fillna(0.0)
    )
    summary["time_label"] = summary["time_minutes"].astype(str)
    return summary


def prepare_temporal_gene_heatmap(
    data: TemporalResult,
    gene_set: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Prepare hover data and cluster genes by standardized temporal profile."""
    summary = prepare_temporal_gene_expression(data, gene_set)
    members = list(data.members[gene_set])
    times = sorted(summary["time_minutes"].unique())
    z_scores = (
        summary.pivot(index="feature_id", columns="time_minutes", values="z_score")
        .loc[members, times]
        .fillna(0.0)
    )
    feature_order = cluster_heatmap_rows(z_scores)

    labels = summary.drop_duplicates("feature_id").set_index("feature_id")[
        "display_label"
    ]
    duplicate_labels = labels.duplicated(keep=False)
    row_labels = labels.astype(str).to_dict()
    for feature_id in labels.index[duplicate_labels]:
        row_labels[feature_id] = f"{labels.at[feature_id]} ({feature_id})"
    summary["row_label"] = summary["feature_id"].map(row_labels)
    summary["cluster_order"] = summary["feature_id"].map(
        {feature_id: index for index, feature_id in enumerate(feature_order)}
    )
    summary = summary.sort_values(
        ["cluster_order", "time_minutes"], kind="mergesort"
    ).reset_index(drop=True)
    return summary, feature_order


def create_temporal_gene_plot(
    data: TemporalResult,
    gene_set: str,
) -> BokehPlot:
    """Create a clustered gene-by-time heatmap of temporal expression profiles."""
    summary, feature_order = prepare_temporal_gene_heatmap(data, gene_set)
    row_labels = (
        summary.drop_duplicates("feature_id")
        .set_index("feature_id")["row_label"]
        .to_dict()
    )
    time_order = [str(value) for value in sorted(summary["time_minutes"].unique())]
    z_limit = max(float(summary["z_score"].abs().max()), 1.0)
    mapper = LinearColorMapper(
        palette=cc.b_diverging_bwr_20_95_c54,
        low=-z_limit,
        high=z_limit,
    )
    plot_height = max(420, min(1200, 20 * len(feature_order) + 120))
    plot = BokehPlot(
        title=f"Gene expression over time — {gene_set}",
        x_range=time_order,
        y_range=list(
            reversed([row_labels[feature_id] for feature_id in feature_order])
        ),
        x_axis_location="above",
        x_axis_label="Elapsed time (min)",
        y_axis_label="Gene",
        height=plot_height,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    source = ColumnDataSource(summary)
    cells = plot._fig.rect(
        x="time_label",
        y="row_label",
        width=0.98,
        height=0.98,
        source=source,
        fill_color={"field": "z_score", "transform": mapper},
        line_color=None,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[cells],
            tooltips=[
                ("Gene", "@display_label"),
                ("Feature ID", "@feature_id"),
                ("Group", "@group"),
                ("Time", "@time_minutes min"),
                ("Gene z-score", "@z_score{0.000}"),
                ("Mean log CPM", "@mean_log_cpm{0.000}"),
                ("SD", "@sd_label"),
                ("Replicates", "@n"),
            ],
        )
    )
    plot._fig.add_layout(ColorBar(color_mapper=mapper, title="Gene z-score"), "right")
    plot._fig.grid.grid_line_color = None
    plot.report_height = plot_height
    return plot


def create_temporal_view(
    data: TemporalResult,
    gene_set: str,
    condition_colors: dict[str, str],
) -> BokehPlot:
    """Stack the pathway- and gene-level temporal figures with interpretation."""
    score_plot = create_temporal_score_plot(data, gene_set, condition_colors)
    gene_plot = create_temporal_gene_plot(data, gene_set)
    note = Div(
        text=f"<p>{TEMPORAL_CAVEAT}</p>",
        sizing_mode="stretch_width",
    )
    combined = BokehPlot()
    combined._fig = column(
        note,
        score_plot._fig,
        gene_plot._fig,
        sizing_mode="stretch_width",
    )
    combined.report_height = score_plot.report_height + gene_plot.report_height + 80
    return combined


def add_temporal_analysis(
    data: TemporalResult,
    condition_colors: dict[str, str],
) -> None:
    """Add synchronized temporal figures behind a gene-set dropdown."""
    tabs = Tabs()
    with tabs.add_dropdown_menu("Gene set", change_header=True):  # type: ignore
        for gene_set in data.gene_set_order:
            with tabs.add_dropdown_tab(gene_set):  # type: ignore
                temporal_view = create_temporal_view(
                    data,
                    gene_set,
                    condition_colors,
                )
                EZChart(
                    temporal_view,
                    "epi2melabs",
                    height=f"{temporal_view.report_height}px",
                )
