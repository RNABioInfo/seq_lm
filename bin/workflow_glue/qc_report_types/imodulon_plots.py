"""Validate fixed-matrix ICA snapshots and add interactive report views."""

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from html import escape
import json
import math
from pathlib import Path

from bokeh.layouts import column
from bokeh.models import (
    ColorBar,
    ColumnDataSource,
    Div,
    FixedTicker,
    HoverTool,
    LinearColorMapper,
    Span,
)
from bokeh.palettes import Category20
import colorcet as cc
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd


READY_FILES = {
    "activities": "activities_long.tsv",
    "summary": "activity_summary.tsv",
    "differential": "differential_activity.tsv",
    "coverage": "component_coverage.tsv",
    "qc": "projection_qc.tsv",
    "mapping": "gene_mapping.tsv",
}


@dataclass
class IModulonResult:
    """One validated immutable ICA snapshot."""

    root: Path
    status: dict
    provenance: dict
    samples: pd.DataFrame
    activities: pd.DataFrame | None = None
    summary: pd.DataFrame | None = None
    differential: pd.DataFrame | None = None
    coverage: pd.DataFrame | None = None
    qc: pd.DataFrame | None = None
    mapping: pd.DataFrame | None = None
    components: tuple[str, ...] = ()
    groups: tuple[str, ...] = ()
    cutoff: float = 0.05

    @property
    def ready(self) -> bool:
        return self.status["status"] == "ready"


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise ValueError(f"ICA report input is missing: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _require(frame: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise ValueError(
            f"{path} is missing required ICA column(s): {', '.join(missing)}"
        )


def _numbers(frame: pd.DataFrame, columns: list[str], path: Path, optional=()) -> None:
    optional = set(optional)
    for name in columns:
        raw = frame[name]
        converted = pd.to_numeric(raw.mask(raw.eq("")), errors="coerce")
        if (converted.isna() & raw.ne("")).any() or (
            converted.isna().any() and name not in optional
        ):
            raise ValueError(
                f"{path} contains invalid or missing numeric values in {name}"
            )
        if converted.dropna().map(math.isfinite).eq(False).any():
            raise ValueError(f"{path} contains nonfinite values in {name}")
        frame[name] = converted


def load_imodulon_results(
    results_dir,
    expected_batch=None,
    expected_sequence=None,
    expected_analysis_index=None,
    expected_samples=None,
) -> IModulonResult:
    """Load a ready or deferred snapshot without recomputing analysis."""
    root = Path(results_dir)
    try:
        status = json.loads((root / "status.json").read_text())
        provenance = json.loads((root / "provenance.json").read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid ICA snapshot metadata in {root}: {exc}") from exc
    if status.get("status") not in {"ready", "deferred"}:
        raise ValueError("ICA status must be ready or deferred")
    if provenance.get("schema_version") != 1:
        raise ValueError("Unsupported ICA provenance schema version")
    for field, expected in (
        ("batch_index", expected_batch),
        ("report_sequence", expected_sequence),
        ("analysis_index", expected_analysis_index),
    ):
        if expected is not None and int(provenance.get(field, -1)) != int(expected):
            raise ValueError(f"ICA snapshot {field} does not match its report record")

    sample_path = root / "sample_metadata.tsv"
    samples = _read_tsv(sample_path)
    _require(
        samples,
        {"sample_id", "alias", "group", "order", "assigned_abundance", "ready"},
        sample_path,
    )
    if samples["sample_id"].eq("").any() or samples["sample_id"].duplicated().any():
        raise ValueError("ICA sample IDs must be unique and nonempty")
    if (
        samples.empty
        or samples[["alias", "group"]].eq("").any().any()
        or samples.duplicated(["group", "alias"]).any()
    ):
        raise ValueError("ICA sample group/alias pairs must be unique and nonempty")
    if expected_samples is not None:
        columns = (
            ("Name", "Group")
            if {"Name", "Group"}.issubset(expected_samples.columns)
            else ("name", "group")
        )
        if not set(columns).issubset(expected_samples.columns):
            raise ValueError("QC sample metadata cannot be matched to ICA samples")
        expected_pairs = Counter(
            zip(
                expected_samples[columns[0]].astype(str),
                expected_samples[columns[1]].astype(str),
            )
        )
        observed_pairs = Counter(zip(samples["alias"], samples["group"]))
        if expected_pairs != observed_pairs:
            raise ValueError("ICA samples do not match the QC report sample cohort")
    _numbers(samples, ["assigned_abundance", "order"], sample_path, optional=["order"])
    samples["ready"] = samples["ready"].str.lower().map({"true": True, "false": False})
    if samples["ready"].isna().any():
        raise ValueError("ICA sample readiness must be true or false")
    threshold = float(provenance.get("settings", {}).get("min_read_count", -1))
    if (
        not math.isfinite(threshold)
        or threshold < 0
        or samples["assigned_abundance"].lt(0).any()
    ):
        raise ValueError(
            "ICA sample readiness requires nonnegative abundance and a valid threshold"
        )
    expected_ready = samples["assigned_abundance"].gt(0) & samples[
        "assigned_abundance"
    ].ge(threshold)
    if not samples["ready"].eq(expected_ready).all() or (
        status["status"] == "ready"
    ) != bool(expected_ready.all()):
        raise ValueError(
            "ICA snapshot status/readiness does not match sample abundance"
        )
    control_ids = provenance.get("control_sample_ids", [])
    control_samples = samples.loc[samples["group"].str.lower().eq("control")]
    if (
        len(control_samples) < 2
        or control_samples["group"].nunique() != 1
        or len(control_ids) != len(set(control_ids))
        or set(control_ids) != set(control_samples["sample_id"])
    ):
        raise ValueError(
            "ICA control sample identities must identify the complete control group with at least two samples"
        )
    result = IModulonResult(root, status, provenance, samples)
    if not result.ready:
        if set(status.get("sample_ids", [])) != set(
            samples.loc[~expected_ready, "sample_id"]
        ):
            raise ValueError("ICA deferred sample identities do not match readiness")
        return result

    loaded = {
        name: _read_tsv(root / filename) for name, filename in READY_FILES.items()
    }
    activities, summary, differential = (
        loaded["activities"],
        loaded["summary"],
        loaded["differential"],
    )
    _require(
        activities,
        {"component_id", "sample_id", "alias", "group", "order", "activity"},
        root / READY_FILES["activities"],
    )
    _require(
        summary,
        {"component_id", "group", "mean", "sd", "n"},
        root / READY_FILES["summary"],
    )
    _require(
        differential,
        {
            "component_id",
            "target_group",
            "control_group",
            "activity_difference",
            "target_mean",
            "control_mean",
            "target_sd",
            "control_sd",
            "target_n",
            "control_n",
            "standard_error",
            "degrees_of_freedom",
            "t_statistic",
            "ci_lower",
            "ci_upper",
            "p_value",
            "adjusted_p_value",
            "significant",
            "status",
        },
        root / READY_FILES["differential"],
    )
    _numbers(
        activities,
        ["activity", "order"],
        root / READY_FILES["activities"],
        optional=["order"],
    )
    _numbers(
        summary, ["mean", "sd", "n"], root / READY_FILES["summary"], optional=["sd"]
    )
    statistic_columns = [
        "activity_difference",
        "target_mean",
        "control_mean",
        "target_sd",
        "control_sd",
        "target_n",
        "control_n",
        "standard_error",
        "degrees_of_freedom",
        "t_statistic",
        "ci_lower",
        "ci_upper",
        "p_value",
        "adjusted_p_value",
    ]
    _numbers(
        differential,
        statistic_columns,
        root / READY_FILES["differential"],
        optional=statistic_columns,
    )
    parsed_significance = (
        differential["significant"]
        .str.lower()
        .map({"true": True, "false": False, "": pd.NA})
    )
    if parsed_significance.isna().ne(differential["significant"].eq("")).any():
        raise ValueError(
            "ICA differential significance must be true, false, or unavailable"
        )
    differential["is_significant"] = parsed_significance.eq(True)
    differential["significant"] = parsed_significance.map(
        {True: "true", False: "false"}
    ).fillna("unavailable")
    if activities.duplicated(["component_id", "sample_id"]).any():
        raise ValueError("ICA activities contain duplicate component/sample rows")
    sample_ids = samples["sample_id"].tolist()
    if set(activities["sample_id"]) != set(sample_ids):
        raise ValueError("ICA activity samples do not match sample metadata")
    model = provenance.get("model", {})
    components = tuple(str(value) for value in model.get("components", []))
    if (
        not components
        or "" in components
        or len(set(components)) != len(components)
        or set(activities["component_id"]) != set(components)
    ):
        raise ValueError("ICA activity components do not match provenance")
    if len(activities) != len(components) * len(sample_ids):
        raise ValueError("ICA activities must contain every component/sample pair")
    for sample in samples.itertuples(index=False):
        rows = activities.loc[activities["sample_id"].eq(sample.sample_id)]
        if (
            not rows["alias"].eq(sample.alias).all()
            or not rows["group"].eq(sample.group).all()
        ):
            raise ValueError("ICA activity metadata does not match sample metadata")
        if not (
            rows["order"].isna().all()
            if pd.isna(sample.order)
            else rows["order"].eq(sample.order).all()
        ):
            raise ValueError("ICA activity order does not match sample metadata")
    settings = provenance.get("settings", {})
    cutoff = float(settings.get("padj_cutoff", 0.05))
    if not 0 < cutoff <= 1:
        raise ValueError("ICA provenance contains an invalid adjusted p-value cutoff")
    result.activities, result.summary, result.differential = (
        activities,
        summary,
        differential,
    )
    result.coverage, result.qc, result.mapping = (
        loaded["coverage"],
        loaded["qc"],
        loaded["mapping"],
    )
    result.components = components
    result.groups = tuple(samples["group"].drop_duplicates())
    result.cutoff = cutoff
    coverage_path, qc_path, mapping_path = (
        root / READY_FILES[name] for name in ("coverage", "qc", "mapping")
    )
    _require(
        result.coverage,
        {"component_id", "gene_coverage", "retained_squared_weight_fraction"},
        coverage_path,
    )
    _require(
        result.qc,
        {
            "sample_id",
            "residual_sum_squares",
            "centered_sum_squares",
            "residual_rmse",
            "normalized_residual",
        },
        qc_path,
    )
    _require(
        result.mapping,
        {"model_gene_id", "gene_id", "transcript_id", "method", "status"},
        mapping_path,
    )
    _numbers(
        result.coverage,
        ["gene_coverage", "retained_squared_weight_fraction"],
        coverage_path,
    )
    _numbers(
        result.qc,
        [
            "residual_sum_squares",
            "centered_sum_squares",
            "residual_rmse",
            "normalized_residual",
        ],
        qc_path,
        optional=["normalized_residual"],
    )
    if (
        result.coverage["component_id"].duplicated().any()
        or tuple(result.coverage["component_id"]) != components
    ):
        raise ValueError(
            "ICA component coverage does not match provenance component order"
        )
    if result.qc["sample_id"].duplicated().any() or set(result.qc["sample_id"]) != set(
        sample_ids
    ):
        raise ValueError("ICA projection QC samples do not match sample metadata")
    expected_summary = {
        (component, group) for component in components for group in result.groups
    }
    if (
        summary.duplicated(["component_id", "group"]).any()
        or set(zip(summary["component_id"], summary["group"])) != expected_summary
    ):
        raise ValueError(
            "ICA activity summary does not contain every component/group pair"
        )
    observed_summary = activities.groupby(["component_id", "group"])["activity"].agg(
        ["mean", "std", "count"]
    )
    indexed_summary = summary.set_index(["component_id", "group"]).loc[
        observed_summary.index
    ]
    for published, observed in (("mean", "mean"), ("sd", "std"), ("n", "count")):
        if not np.allclose(
            indexed_summary[published],
            observed_summary[observed],
            rtol=1e-7,
            atol=1e-10,
            equal_nan=True,
        ):
            raise ValueError(
                f"ICA activity summary {published} does not match sample activities"
            )
    valid_test_statuses = {"tested", "insufficient_replicates", "zero_variance"}
    if not set(differential["status"]).issubset(valid_test_statuses):
        raise ValueError("ICA differential activity contains an unknown test status")
    control_ids = provenance.get("control_sample_ids", [])
    if not control_ids or not set(control_ids).issubset(set(sample_ids)):
        raise ValueError("ICA provenance contains invalid control sample identities")
    control_groups = set(samples.loc[samples["sample_id"].isin(control_ids), "group"])
    if len(control_groups) != 1:
        raise ValueError("ICA control sample identities do not identify one group")
    expected_contrasts = {
        (component, group, next(iter(control_groups)))
        for component in components
        for group in result.groups
        if group not in control_groups
    }
    if (
        set(
            zip(
                differential["component_id"],
                differential["target_group"],
                differential["control_group"],
            )
        )
        != expected_contrasts
    ):
        raise ValueError(
            "ICA differential activity must contain every component/contrast pair"
        )
    if not differential.empty:
        if differential.duplicated(
            ["component_id", "target_group", "control_group"]
        ).any():
            raise ValueError(
                "ICA differential activity contains duplicate component/contrast rows"
            )
        if (
            not set(differential["component_id"]).issubset(set(components))
            or not set(differential["target_group"]).issubset(set(result.groups))
            or set(differential["control_group"]) != control_groups
        ):
            raise ValueError(
                "ICA differential activity identities do not match snapshot metadata"
            )
        for _contrast, rows in differential.groupby(
            ["target_group", "control_group"], sort=False
        ):
            if tuple(rows["component_id"]) != components:
                raise ValueError(
                    "ICA differential activity must preserve every component per contrast"
                )
    _validate_statistics(result)
    return result


def _validate_statistics(data: IModulonResult) -> None:
    """Check published statistics and availability without rerunning inference."""
    rows = data.differential
    summary = data.summary.set_index(["component_id", "group"])
    inference = [
        "standard_error",
        "degrees_of_freedom",
        "t_statistic",
        "ci_lower",
        "ci_upper",
        "p_value",
        "adjusted_p_value",
    ]
    for row in rows.itertuples(index=False):
        for prefix, group in (
            ("target", row.target_group),
            ("control", row.control_group),
        ):
            expected = summary.loc[(row.component_id, group)]
            count, mean, sd = (
                getattr(row, f"{prefix}_{name}") for name in ("n", "mean", "sd")
            )
            if count != expected["n"] or not np.isclose(
                mean, expected["mean"], rtol=1e-7, atol=1e-10
            ):
                raise ValueError(
                    "ICA differential replicate counts/means do not match activity summaries"
                )
            if (count < 2 and not pd.isna(sd)) or (
                count >= 2 and (pd.isna(sd) or sd < 0)
            ):
                raise ValueError(
                    "ICA differential SD availability does not match replication"
                )
        if not np.isclose(
            row.activity_difference,
            row.target_mean - row.control_mean,
            rtol=1e-7,
            atol=1e-10,
        ):
            raise ValueError(
                "ICA activity difference does not match target-minus-control means"
            )
        if row.status == "tested":
            if min(row.target_n, row.control_n) < 2 or any(
                pd.isna(getattr(row, name)) for name in inference
            ):
                raise ValueError(
                    "ICA tested rows require replication and complete inferential statistics"
                )
            if (
                row.standard_error <= 0
                or row.degrees_of_freedom <= 0
                or not row.ci_lower <= row.activity_difference <= row.ci_upper
                or not 0 <= row.p_value <= 1
                or not 0 <= row.adjusted_p_value <= 1
                or (
                    row.adjusted_p_value < row.p_value
                    and not np.isclose(
                        row.adjusted_p_value, row.p_value, rtol=1e-12, atol=0
                    )
                )
                or (row.target_sd == 0 and row.control_sd == 0)
            ):
                raise ValueError(
                    "ICA tested row contains invalid inferential statistics"
                )
            if row.significant != (
                "true" if row.adjusted_p_value <= data.cutoff else "false"
            ):
                raise ValueError(
                    "ICA significance does not match the provenance cutoff"
                )
        else:
            if row.significant != "unavailable" or any(
                not pd.isna(getattr(row, name)) for name in inference
            ):
                raise ValueError(
                    "ICA untestable rows must have unavailable inference and significance"
                )
            if (
                row.status == "insufficient_replicates"
                and min(row.target_n, row.control_n) >= 2
            ):
                raise ValueError(
                    "ICA insufficient-replicates status contradicts sample counts"
                )
            if row.status == "zero_variance" and (
                min(row.target_n, row.control_n) < 2
                or row.target_sd != 0
                or row.control_sd != 0
            ):
                raise ValueError(
                    "ICA zero-variance status contradicts sample statistics"
                )
    tested = int(rows["status"].eq("tested").sum())
    availability = (
        "complete"
        if tested and tested == len(rows)
        else "partial"
        if tested
        else "unavailable"
    )
    if (
        data.status.get("tested_count") != tested
        or data.status.get("untested_count") != len(rows) - tested
        or data.status.get("statistical_availability") != availability
    ):
        raise ValueError(
            "ICA statistical availability does not match differential activity rows"
        )


def validate_imodulon_timecourse(data: IModulonResult) -> None:
    """Validate the independent-replicate, one-group-per-time contract."""
    if data.samples["order"].isna().any():
        raise ValueError(
            "ICA time-course reporting requires an order value for every sample"
        )
    if (data.samples["order"] % 1).ne(0).any():
        raise ValueError(
            "ICA time-course order values must be signed integer elapsed minutes"
        )
    group_times = data.samples.groupby("group", sort=False)["order"].nunique()
    if group_times.ne(1).any():
        raise ValueError(
            "ICA time-course reporting requires one elapsed minute per group"
        )
    time_groups = data.samples.groupby("order", sort=False)["group"].nunique()
    if time_groups.ne(1).any():
        raise ValueError(
            "ICA time-course reporting requires one group per elapsed minute"
        )
    if data.samples["order"].nunique() < 2:
        raise ValueError(
            "ICA time-course reporting requires at least two elapsed minutes"
        )


def _colors(groups):
    palette = Category20[20]
    return {group: palette[index % len(palette)] for index, group in enumerate(groups)}


def _wrap(plot: BokehPlot, height=None) -> None:
    EZChart(
        plot, "epi2melabs", height=f"{height or getattr(plot, 'report_height', 500)}px"
    )


def _empty(title, message, height=250) -> BokehPlot:
    plot = BokehPlot()
    plot._fig = column(
        Div(text=f"<h4>{escape(title)}</h4><p>{escape(message)}</p>"),
        sizing_mode="stretch_width",
    )
    plot.report_height = height
    return plot


def create_activity_heatmap(data: IModulonResult, selected: list[str]) -> BokehPlot:
    rows = data.activities.loc[data.activities["component_id"].isin(selected)].copy()
    labels = (
        data.samples.set_index("sample_id")
        .apply(lambda row: f"{row['group']} / {row['alias']}", axis=1)
        .to_dict()
    )
    rows["sample_label"] = rows["sample_id"].map(labels)
    limit = max(float(data.activities["activity"].abs().max()), np.finfo(float).eps)
    mapper = LinearColorMapper(
        palette=cc.b_diverging_bwr_20_95_c54, low=-limit, high=limit
    )
    plot = BokehPlot(
        title="Control-centered iModulon activities",
        x_range=data.samples["sample_id"].tolist(),
        y_range=list(reversed(selected)),
        x_axis_location="above",
        x_axis_label="Biological sample",
        y_axis_label="Component",
        height=max(430, len(selected) * 20 + 150),
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    plot._fig.xaxis.major_label_overrides = labels
    glyph = plot._fig.rect(
        "sample_id",
        "component_id",
        0.98,
        0.98,
        source=ColumnDataSource(rows),
        fill_color={"field": "activity", "transform": mapper},
        line_color=None,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[glyph],
            tooltips=[
                ("Component", "@component_id"),
                ("Sample", "@alias"),
                ("Group", "@group"),
                ("Activity", "@activity{0.0000}"),
            ],
        )
    )
    plot._fig.add_layout(ColorBar(color_mapper=mapper, title="Activity"), "right")
    plot._fig.xaxis.major_label_orientation = np.pi / 4
    plot._fig.grid.grid_line_color = None
    plot.report_height = plot._fig.height
    return plot


def create_component_distribution(data: IModulonResult, component: str) -> BokehPlot:
    rows = data.activities.loc[data.activities["component_id"].eq(component)].copy()
    groups, colors = list(data.groups), _colors(data.groups)
    plot = BokehPlot(
        title=f"Biological sample activities — {component}",
        x_range=(0.5, len(groups) + 0.5),
        x_axis_label="Group",
        y_axis_label="Control-centered activity",
        height=470,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    plot._fig.xaxis.ticker = FixedTicker(ticks=list(range(1, len(groups) + 1)))
    plot._fig.xaxis.major_label_overrides = {
        float(i): group for i, group in enumerate(groups, 1)
    }
    points = []
    for index, group in enumerate(groups, 1):
        subset = rows.loc[rows["group"].eq(group)]
        values = subset["activity"].to_numpy(float)
        mean = float(values.mean())
        sd = float(values.std(ddof=1)) if len(values) > 1 else None
        if sd is not None:
            plot._fig.segment(
                [index],
                [mean - sd],
                [index],
                [mean + sd],
                line_color=colors[group],
                line_width=2,
            )
        plot._fig.scatter(
            [index], [mean], marker="diamond", size=12, color=colors[group]
        )
        offsets = np.linspace(-0.16, 0.16, len(subset)) if len(subset) > 1 else [0]
        for offset, row in zip(offsets, subset.itertuples(index=False), strict=True):
            points.append(
                {
                    "x": index + float(offset),
                    "activity": row.activity,
                    "alias": row.alias,
                    "group": group,
                    "color": colors[group],
                }
            )
    source = ColumnDataSource(pd.DataFrame(points))
    glyph = plot._fig.scatter(
        "x", "activity", source=source, size=9, color="color", line_color="white"
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[glyph],
            tooltips=[
                ("Sample", "@alias"),
                ("Group", "@group"),
                ("Activity", "@activity{0.0000}"),
            ],
        )
    )
    return plot


def create_volcano(rows: pd.DataFrame, label: str, cutoff: float) -> BokehPlot:
    tested = rows.loc[
        rows["status"].eq("tested") & rows["adjusted_p_value"].notna()
    ].copy()
    if tested.empty:
        return _empty(
            f"Differential activity — {label}",
            "No components have an available Welch test for this contrast.",
        )
    # Preserve positive values exactly; zero is a numerical underflow sentinel.
    plotted_p = tested["adjusted_p_value"].mask(
        tested["adjusted_p_value"].eq(0), np.nextafter(0.0, 1.0)
    )
    tested["minus_log10_adjusted_p"] = -np.log10(plotted_p)
    tested["plot_color"] = np.where(
        tested["adjusted_p_value"].le(cutoff),
        np.where(tested["activity_difference"].ge(0), "#B2182B", "#2166AC"),
        "#888888",
    )
    plot = BokehPlot(
        title=f"Differential activity — {label}",
        x_axis_label="Activity difference (target − control)",
        y_axis_label="−log10 adjusted p-value",
        height=500,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    source = ColumnDataSource(tested)
    plot._fig.add_layout(
        Span(
            location=-math.log10(cutoff),
            dimension="width",
            line_dash="dashed",
            line_color="#555555",
        )
    )
    plot._fig.add_layout(
        Span(location=0, dimension="height", line_dash="dotted", line_color="#555555")
    )
    glyph = plot._fig.scatter(
        "activity_difference",
        "minus_log10_adjusted_p",
        source=source,
        color="plot_color",
        size=9,
        alpha=0.85,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[glyph],
            tooltips=[
                ("Component", "@component_id"),
                ("Difference", "@activity_difference{0.0000}"),
                ("Adjusted p-value", "@adjusted_p_value{0.000e}"),
                ("Status", "@status"),
            ],
        )
    )
    return plot


def create_timecourse(data: IModulonResult, component: str) -> BokehPlot:
    rows = data.activities.loc[data.activities["component_id"].eq(component)].copy()
    colors = _colors(data.groups)
    points, summaries = [], []
    for group, subset in rows.groupby("group", sort=False):
        time = int(subset["order"].iloc[0])
        mean, count = float(subset["activity"].mean()), len(subset)
        sd = float(subset["activity"].std(ddof=1)) if count > 1 else np.nan
        span = sorted(data.samples["order"].unique())
        gap = min(np.diff(span)) if len(span) > 1 else 1
        offsets = np.linspace(-0.025 * gap, 0.025 * gap, count) if count > 1 else [0]
        summaries.append(
            {
                "time": time,
                "mean": mean,
                "sd": sd,
                "lower": mean - sd,
                "upper": mean + sd,
                "group": group,
                "n": count,
                "color": colors[group],
            }
        )
        for offset, row in zip(offsets, subset.itertuples(index=False), strict=True):
            points.append(
                {
                    "plot_time": time + float(offset),
                    "time": time,
                    "activity": row.activity,
                    "alias": row.alias,
                    "group": group,
                    "color": colors[group],
                }
            )
    summary = pd.DataFrame(summaries).sort_values("time")
    plot = BokehPlot(
        title=f"Activity over elapsed time — {component}",
        x_axis_label="Elapsed time (min)",
        y_axis_label="Control-centered activity",
        height=480,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    source = ColumnDataSource(summary)
    plot._fig.line("time", "mean", source=source, line_color="#333333", line_width=2.5)
    finite = summary.loc[summary["sd"].notna()]
    if not finite.empty:
        plot._fig.segment(
            "time",
            "lower",
            "time",
            "upper",
            source=ColumnDataSource(finite),
            line_color="color",
            line_width=2,
        )
    means = plot._fig.scatter(
        "time", "mean", source=source, marker="diamond", size=12, color="color"
    )
    point_source = ColumnDataSource(pd.DataFrame(points))
    samples = plot._fig.scatter(
        "plot_time",
        "activity",
        source=point_source,
        size=9,
        color="color",
        line_color="white",
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[samples],
            tooltips=[
                ("Sample", "@alias"),
                ("Group", "@group"),
                ("Time", "@time min"),
                ("Activity", "@activity{0.0000}"),
            ],
        ),
        HoverTool(
            renderers=[means],
            tooltips=[
                ("Group", "@group"),
                ("Time", "@time min"),
                ("Mean", "@mean{0.0000}"),
                ("SD", "@sd{0.0000}"),
                ("n", "@n"),
            ],
        ),
    )
    return plot


def create_time_effects(data: IModulonResult, component: str) -> BokehPlot:
    rows = data.differential.loc[data.differential["component_id"].eq(component)].copy()
    group_time = data.samples.groupby("group", sort=False)["order"].first().to_dict()
    rows["time"] = rows["target_group"].map(group_time)
    rows = rows.sort_values("time")
    if rows.empty:
        return _empty(
            f"Activity differences over time — {component}",
            "No non-control groups are available.",
        )
    rows["color"] = np.where(rows["is_significant"], "#B2182B", "#777777")
    plot = BokehPlot(
        title=f"Target-minus-control activity difference — {component}",
        x_axis_label="Target elapsed time (min)",
        y_axis_label="Activity difference",
        height=420,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    tested = rows.loc[rows["ci_lower"].notna() & rows["ci_upper"].notna()]
    if not tested.empty:
        plot._fig.segment(
            "time",
            "ci_lower",
            "time",
            "ci_upper",
            source=ColumnDataSource(tested),
            line_color="color",
            line_width=2,
        )
    source = ColumnDataSource(rows)
    glyph = plot._fig.scatter(
        "time", "activity_difference", source=source, size=10, color="color"
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[glyph],
            tooltips=[
                ("Group", "@target_group"),
                ("Time", "@time min"),
                ("Difference", "@activity_difference{0.0000}"),
                ("95% CI", "@ci_lower{0.0000} to @ci_upper{0.0000}"),
                ("Adjusted p-value", "@adjusted_p_value{0.000e}"),
                ("Status", "@status"),
            ],
        )
    )
    return plot


def create_component_effects(data: IModulonResult, component: str) -> BokehPlot:
    """Show all target-minus-control differences and nominal confidence intervals."""
    rows = data.differential.loc[data.differential["component_id"].eq(component)].copy()
    if rows.empty:
        return _empty(
            f"Activity differences — {component}",
            "No non-control groups are available.",
        )
    rows["x"] = np.arange(1, len(rows) + 1)
    rows["color"] = np.where(rows["is_significant"], "#B2182B", "#777777")
    plot = BokehPlot(
        title=f"Target-minus-control activity differences — {component}",
        x_range=(0.5, len(rows) + 0.5),
        x_axis_label="Target group",
        y_axis_label="Activity difference",
        height=390,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    plot._fig.xaxis.ticker = FixedTicker(ticks=rows["x"].tolist())
    plot._fig.xaxis.major_label_overrides = {
        float(row.x): row.target_group for row in rows.itertuples()
    }
    tested = rows.loc[rows["ci_lower"].notna() & rows["ci_upper"].notna()]
    if not tested.empty:
        plot._fig.segment(
            "x",
            "ci_lower",
            "x",
            "ci_upper",
            source=ColumnDataSource(tested),
            line_color="color",
            line_width=2,
        )
    source = ColumnDataSource(rows)
    glyph = plot._fig.scatter(
        "x", "activity_difference", source=source, size=10, color="color"
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[glyph],
            tooltips=[
                ("Target", "@target_group"),
                ("Control", "@control_group"),
                ("Difference", "@activity_difference{0.0000}"),
                ("95% CI", "@ci_lower{0.0000} to @ci_upper{0.0000}"),
                ("Adjusted p-value", "@adjusted_p_value{0.000e}"),
                ("Status", "@status"),
            ],
        )
    )
    return plot


def create_multi_contrast_effects(data: IModulonResult) -> BokehPlot:
    """Show signed activity differences and significance across all contrasts."""
    rows = data.differential.copy()
    rows["contrast"] = rows["target_group"] + " vs " + rows["control_group"]
    rows["significance"] = np.where(rows["is_significant"], "●", "")
    limit = max(float(rows["activity_difference"].abs().max()), np.finfo(float).eps)
    mapper = LinearColorMapper(
        palette=cc.b_diverging_bwr_20_95_c54, low=-limit, high=limit
    )
    plot = BokehPlot(
        title="Activity differences across contrasts",
        x_range=rows["contrast"].drop_duplicates().tolist(),
        y_range=list(reversed(data.components)),
        x_axis_label="Contrast",
        y_axis_label="Component",
        height=max(440, 20 * len(data.components) + 140),
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    source = ColumnDataSource(rows)
    glyph = plot._fig.rect(
        "contrast",
        "component_id",
        0.98,
        0.98,
        source=source,
        fill_color={"field": "activity_difference", "transform": mapper},
        line_color=None,
    )
    plot._fig.text(
        "contrast",
        "component_id",
        text="significance",
        source=source,
        text_align="center",
        text_baseline="middle",
        text_color="black",
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[glyph],
            tooltips=[
                ("Contrast", "@contrast"),
                ("Component", "@component_id"),
                ("Difference", "@activity_difference{0.0000}"),
                ("Adjusted p-value", "@adjusted_p_value{0.000e}"),
                ("Status", "@status"),
            ],
        )
    )
    plot._fig.add_layout(
        ColorBar(color_mapper=mapper, title="Activity difference"), "right"
    )
    plot._fig.xaxis.major_label_orientation = np.pi / 4
    plot._fig.grid.grid_line_color = None
    plot.report_height = plot._fig.height
    return plot


def create_time_heatmap(data: IModulonResult) -> BokehPlot:
    rows = data.summary.merge(
        data.samples[["group", "order"]].drop_duplicates(),
        on="group",
        validate="many_to_one",
    ).copy()
    times = sorted(rows["order"].unique())
    gap = float(np.diff(times).min()) if len(times) > 1 else 1.0
    limit = max(float(rows["mean"].abs().max()), np.finfo(float).eps)
    mapper = LinearColorMapper(
        palette=cc.b_diverging_bwr_20_95_c54, low=-limit, high=limit
    )
    plot = BokehPlot(
        title="Mean iModulon activity over elapsed time",
        y_range=list(reversed(data.components)),
        x_axis_label="Elapsed time (min)",
        y_axis_label="Component",
        height=max(440, 20 * len(data.components) + 140),
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    glyph = plot._fig.rect(
        "order",
        "component_id",
        gap * 0.8,
        0.98,
        source=ColumnDataSource(rows),
        fill_color={"field": "mean", "transform": mapper},
        line_color=None,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[glyph],
            tooltips=[
                ("Component", "@component_id"),
                ("Group", "@group"),
                ("Time", "@order min"),
                ("Mean activity", "@mean{0.0000}"),
                ("SD", "@sd{0.0000}"),
                ("n", "@n"),
            ],
        )
    )
    plot._fig.add_layout(ColorBar(color_mapper=mapper, title="Mean activity"), "right")
    plot._fig.grid.grid_line_color = None
    plot.report_height = plot._fig.height
    return plot


def _overview(data: IModulonResult) -> pd.DataFrame:
    diagnostics = data.provenance.get("model", {}).get("diagnostics", {})
    settings = data.provenance.get("settings", {})
    values = {
        "Status": data.status["status"],
        "Batch index": data.provenance.get("batch_index"),
        "ICA analysis index": data.provenance.get("analysis_index"),
        "Samples": len(data.samples),
        "Components": len(data.components) if data.ready else "unavailable",
        "Statistical availability": data.status.get(
            "statistical_availability", "unavailable"
        ),
        "Control samples": ", ".join(data.provenance.get("control_sample_ids", [])),
        "Gene coverage": diagnostics.get("gene_coverage"),
        "Log base": settings.get("log_base"),
        "Pseudocount": settings.get("pseudocount"),
        "Minimum assigned abundance": settings.get("min_read_count"),
        "Adjusted p-value cutoff": settings.get("padj_cutoff"),
    }
    return pd.DataFrame({"Metric": values.keys(), "Value": values.values()})


def imodulon_timecourse_enabled(data: IModulonResult) -> bool:
    """Infer temporal display from complete order metadata."""
    present = data.samples["order"].notna()
    if present.any() and not present.all():
        raise ValueError(
            "ICA time-course reporting requires order for every sample when order is used"
        )
    return bool(present.all())


def add_imodulon_analysis(data: IModulonResult) -> None:
    """Render ready and deferred ICA snapshots as one primary report tab."""
    timecourse = imodulon_timecourse_enabled(data)
    if timecourse:
        validate_imodulon_timecourse(data)
    tabs = Tabs()
    with tabs.add_tab("Overview"):
        DataTable.from_pandas(_overview(data), use_index=False)
        readiness = data.samples[
            ["sample_id", "alias", "group", "order", "assigned_abundance", "ready"]
        ]
        DataTable.from_pandas(readiness, use_index=False)
    if not data.ready:
        return
    with tabs.add_tab("Activities"):
        pages = Tabs()
        for start in range(0, len(data.components), 50):
            selected = list(data.components[start : start + 50])
            label = f"{start + 1}–{start + len(selected)}"
            with pages.add_tab(label):
                _wrap(create_activity_heatmap(data, selected))
    with tabs.add_tab("Component details"):
        selector = Tabs()
        with selector.add_dropdown_menu("Component", change_header=True):  # type: ignore
            for component in data.components:
                with selector.add_dropdown_tab(component):  # type: ignore
                    combined = BokehPlot()
                    combined._fig = column(
                        create_component_distribution(data, component)._fig,
                        create_component_effects(data, component)._fig,
                        sizing_mode="stretch_width",
                    )
                    combined.report_height = 900
                    _wrap(combined, 900)
    with tabs.add_tab("Differential activity"):
        contrasts = list(
            data.differential.groupby(["target_group", "control_group"], sort=False)
        )
        if not contrasts:
            _wrap(
                _empty(
                    "Differential activity",
                    "This snapshot contains controls only; no contrasts are available.",
                )
            )
            DataTable.from_pandas(
                data.differential.drop(columns="is_significant"), use_index=False
            )
        elif len(contrasts) == 1:
            (target, control), rows = contrasts[0]
            _wrap(create_volcano(rows, f"{target} vs {control}", data.cutoff), 530)
            DataTable.from_pandas(rows.drop(columns="is_significant"), use_index=False)
        else:
            _wrap(create_multi_contrast_effects(data))
            contrast_tabs = Tabs()
            for (target, control), rows in contrasts:
                with contrast_tabs.add_tab(f"{target} vs {control}"):
                    _wrap(
                        create_volcano(rows, f"{target} vs {control}", data.cutoff), 530
                    )
                    DataTable.from_pandas(
                        rows.drop(columns="is_significant"), use_index=False
                    )
    if timecourse:
        with tabs.add_tab("Time course"):
            _wrap(create_time_heatmap(data))
            selector = Tabs()
            with selector.add_dropdown_menu("Component", change_header=True):  # type: ignore
                for component in data.components:
                    with selector.add_dropdown_tab(component):  # type: ignore
                        combined = BokehPlot()
                        combined._fig = column(
                            create_timecourse(data, component)._fig,
                            create_time_effects(data, component)._fig,
                            sizing_mode="stretch_width",
                        )
                        combined.report_height = 950
                        _wrap(combined, 950)
    with tabs.add_tab("Diagnostics"):
        diagnostics = data.provenance["model"]["diagnostics"]
        singular = pd.DataFrame(
            {
                "singular_value_index": range(
                    1, len(diagnostics["singular_values"]) + 1
                ),
                "singular_value": diagnostics["singular_values"],
            }
        )
        DataTable.from_pandas(
            pd.DataFrame(
                {
                    "Metric": [
                        "Rank",
                        "Rank tolerance",
                        "Condition number",
                        "Shared genes",
                        "Model genes",
                        "Gene coverage",
                    ],
                    "Value": [
                        diagnostics.get("rank"),
                        diagnostics.get("rank_tolerance"),
                        diagnostics.get("condition_number"),
                        diagnostics.get("shared_gene_count"),
                        diagnostics.get("model_gene_count"),
                        diagnostics.get("gene_coverage"),
                    ],
                }
            ),
            use_index=False,
        )
        DataTable.from_pandas(singular, use_index=False)
        DataTable.from_pandas(data.coverage, use_index=False)
        DataTable.from_pandas(
            data.qc.merge(
                data.samples[["sample_id", "alias", "group"]],
                on="sample_id",
                validate="one_to_one",
            ),
            use_index=False,
        )
        DataTable.from_pandas(data.mapping, use_index=False)
    with tabs.add_tab("Methods"):
        methods = _empty(
            "Interpretation",
            "Activities are projections of log-transformed, million-scaled abundance onto the supplied fixed matrix and are centered on the shared controls. Effects are activity differences, not log fold changes. Each component uses an independent two-sided Welch test; Benjamini–Hochberg correction is applied separately within each contrast. Confidence intervals are nominal 95% intervals. Missing tests remain unavailable. Component signs and scales belong to the supplied model, so magnitudes are not directly comparable between components. Coverage and reconstruction measures are technical diagnostics, not biological confidence estimates. Live-batch tests are snapshot analyses and do not provide sequential error control.",
            300,
        )
        _wrap(methods, 300)
