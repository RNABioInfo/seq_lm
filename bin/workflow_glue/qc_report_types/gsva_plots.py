"""Parse GSVA outputs and add score- and limma-level report plots."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bokeh.layouts import gridplot
from bokeh.models import (
    ColorBar,
    ColumnDataSource,
    FixedTicker,
    HoverTool,
    LinearColorMapper,
    Span,
)
from bokeh.palettes import RdBu11
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd

from .differential_plots import (
    _empty_plot,
    ContrastResult,
    DifferentialResult,
    NONSIGNIFICANT_COLOR,
)


GSVA_SCORES_LONG_FILE = "gsva_scores_long.tsv"
GSVA_COVERAGE_FILE = "gsva_gene_set_coverage.tsv"
GSVA_LIMMA_FILE = "gsva_limma_results.tsv"
REQUIRED_SCORE_COLUMNS = (
    "gene_set",
    "description",
    "n_genes",
    "sample",
    "group",
    "score",
)
REQUIRED_COVERAGE_COLUMNS = (
    "gene_set",
    "description",
    "resolved_members",
    "retained_members",
    "variable_members",
    "scored_members",
    "status",
)
REQUIRED_LIMMA_COLUMNS = (
    "gene_set",
    "description",
    "n_genes",
    "target_group",
    "control_group",
    "effect_size",
    "average_score",
    "t_statistic",
    "p_value",
    "adjusted_p_value",
    "log_odds",
)
MAX_HEATMAP_GENE_SETS = 50
MAX_DIFFERENTIAL_GENE_SETS = 30


@dataclass
class GSVAContrastResult:
    """One validated limma contrast over GSVA scores."""

    contrast: ContrastResult
    results: pd.DataFrame

    @property
    def label(self) -> str:
        """Human-readable contrast label."""
        return self.contrast.label


@dataclass
class GSVAResult:
    """Validated sample-level GSVA scores and pathway-level limma results."""

    scores_long: pd.DataFrame
    coverage: pd.DataFrame
    contrasts: list[GSVAContrastResult]

    @property
    def gene_set_order(self) -> list[str]:
        """Return the stable score-table gene-set order."""
        return self.scores_long["gene_set"].drop_duplicates().tolist()


def _require_columns(table: pd.DataFrame, required, path: Path) -> None:
    missing = [column for column in required if column not in table]
    if missing:
        raise ValueError(f"{path} is missing columns: " + ", ".join(missing))


def _validate_identifiers(table: pd.DataFrame, columns, path: Path) -> None:
    for column in columns:
        if table[column].isna().any():
            raise ValueError(f"{path} contains missing {column} values.")
        table[column] = table[column].astype(str).str.strip()
        if table[column].eq("").any():
            raise ValueError(f"{path} contains empty {column} values.")


def _validate_finite(table: pd.DataFrame, columns, path: Path) -> None:
    for column in columns:
        values = pd.to_numeric(table[column], errors="coerce")
        if values.isna().any() or not np.isfinite(values.to_numpy(dtype=float)).all():
            raise ValueError(f"{path} contains nonfinite {column} values.")
        table[column] = values


def _validate_probabilities(table: pd.DataFrame, columns, path: Path) -> None:
    _validate_finite(table, columns, path)
    for column in columns:
        if table[column].lt(0).any() or table[column].gt(1).any():
            raise ValueError(f"{path} contains invalid {column} values.")


def _read_scores(results_dir: Path, differential: DifferentialResult) -> pd.DataFrame:
    path = results_dir / GSVA_SCORES_LONG_FILE
    if not path.is_file():
        raise ValueError(f"Missing GSVA score table: {path}")
    scores = pd.read_csv(path, sep="\t")
    _require_columns(scores, REQUIRED_SCORE_COLUMNS, path)
    _validate_identifiers(scores, ("gene_set", "sample", "group"), path)
    _validate_finite(scores, ("n_genes", "score"), path)
    if (
        scores["n_genes"].lt(2).any()
        or not np.equal(scores["n_genes"], np.floor(scores["n_genes"])).all()
    ):
        raise ValueError(f"{path} contains invalid n_genes values.")
    scores["n_genes"] = scores["n_genes"].astype(int)
    if scores.duplicated(["gene_set", "sample"]).any():
        raise ValueError(f"{path} contains duplicate gene_set/sample rows.")

    expected_samples = differential.sample_metadata.index.tolist()
    observed_samples = scores["sample"].drop_duplicates().tolist()
    if set(observed_samples) != set(expected_samples):
        raise ValueError(f"{path} samples do not match differential metadata.")
    expected_pairs = {
        (sample, differential.sample_metadata.at[sample, "group"])
        for sample in expected_samples
    }
    observed_pairs = set(zip(scores["sample"], scores["group"], strict=True))
    if observed_pairs != expected_pairs:
        raise ValueError(f"{path} sample groups do not match differential metadata.")

    gene_set_order = scores["gene_set"].drop_duplicates().tolist()
    expected_grid = pd.MultiIndex.from_product(
        [gene_set_order, expected_samples],
        names=["gene_set", "sample"],
    )
    observed_grid = pd.MultiIndex.from_frame(scores[["gene_set", "sample"]])
    if set(observed_grid) != set(expected_grid):
        raise ValueError(f"{path} does not contain every gene-set/sample combination.")

    for column in ("description", "n_genes"):
        if scores.groupby("gene_set", sort=False)[column].nunique(dropna=False).gt(1).any():
            raise ValueError(f"{path} contains inconsistent {column} values.")
    scores["description"] = scores["description"].fillna("").astype(str)
    order = pd.MultiIndex.from_product([gene_set_order, expected_samples])
    scores = scores.set_index(["gene_set", "sample"]).loc[order].reset_index()
    return scores


def _read_coverage(results_dir: Path, scored_gene_sets: list[str]) -> pd.DataFrame:
    path = results_dir / GSVA_COVERAGE_FILE
    if not path.is_file():
        raise ValueError(f"Missing GSVA coverage table: {path}")
    coverage = pd.read_csv(path, sep="\t")
    _require_columns(coverage, REQUIRED_COVERAGE_COLUMNS, path)
    _validate_identifiers(coverage, ("gene_set", "status"), path)
    if coverage["gene_set"].duplicated().any():
        raise ValueError(f"{path} contains duplicate gene_set values.")
    count_columns = (
        "resolved_members",
        "retained_members",
        "variable_members",
        "scored_members",
    )
    _validate_finite(coverage, count_columns, path)
    for column in count_columns:
        if (
            coverage[column].lt(0).any()
            or not np.equal(coverage[column], np.floor(coverage[column])).all()
        ):
            raise ValueError(f"{path} contains invalid {column} values.")
        coverage[column] = coverage[column].astype(int)
    invalid_status = set(coverage["status"]) - {"scored", "below_min_size"}
    if invalid_status:
        raise ValueError(f"{path} contains invalid status values.")
    observed_scored = coverage.loc[coverage["status"].eq("scored"), "gene_set"].tolist()
    if observed_scored != scored_gene_sets:
        raise ValueError(f"{path} scored gene sets do not match {GSVA_SCORES_LONG_FILE}.")
    coverage["description"] = coverage["description"].fillna("").astype(str)
    return coverage


def _read_limma(
    path: Path,
    contrast: ContrastResult,
    scores: pd.DataFrame,
) -> pd.DataFrame:
    if not path.is_file():
        raise ValueError(f"Missing GSVA limma table: {path}")
    results = pd.read_csv(path, sep="\t")
    _require_columns(results, REQUIRED_LIMMA_COLUMNS, path)
    _validate_identifiers(
        results,
        ("gene_set", "target_group", "control_group"),
        path,
    )
    if results["gene_set"].duplicated().any():
        raise ValueError(f"{path} contains duplicate gene_set values.")
    expected_gene_sets = scores["gene_set"].drop_duplicates().tolist()
    if results["gene_set"].tolist() != expected_gene_sets:
        raise ValueError(f"{path} gene-set order does not match the GSVA score table.")
    if set(results["target_group"]) != {contrast.target_group}:
        raise ValueError(f"{path} has an unexpected target_group.")
    if set(results["control_group"]) != {contrast.reference_group}:
        raise ValueError(f"{path} has an unexpected control_group.")
    _validate_finite(
        results,
        ("n_genes", "effect_size", "average_score", "t_statistic", "log_odds"),
        path,
    )
    _validate_probabilities(results, ("p_value", "adjusted_p_value"), path)
    score_metadata = scores.drop_duplicates("gene_set").set_index("gene_set")
    for row in results.itertuples(index=False):
        if int(row.n_genes) != int(score_metadata.at[row.gene_set, "n_genes"]):
            raise ValueError(f"{path} n_genes does not match the GSVA score table.")
    results["n_genes"] = results["n_genes"].astype(int)
    results["description"] = results["description"].fillna("").astype(str)
    results["display_label"] = results["gene_set"]
    return results


def load_gsva_results(
    results_dir: str | Path,
    differential: DifferentialResult,
) -> GSVAResult:
    """Load and cross-validate all GSVA report inputs."""
    results_dir = Path(results_dir)
    scores = _read_scores(results_dir, differential)
    gene_set_order = scores["gene_set"].drop_duplicates().tolist()
    coverage = _read_coverage(results_dir, gene_set_order)
    contrasts = [
        GSVAContrastResult(
            contrast=contrast,
            results=_read_limma(
                results_dir / contrast.contrast_id / GSVA_LIMMA_FILE,
                contrast,
                scores,
            ),
        )
        for contrast in differential.contrasts
    ]
    return GSVAResult(scores_long=scores, coverage=coverage, contrasts=contrasts)


def _score_matrix(data: GSVAResult) -> pd.DataFrame:
    sample_order = data.scores_long["sample"].drop_duplicates().tolist()
    return data.scores_long.pivot(
        index="gene_set",
        columns="sample",
        values="score",
    ).loc[data.gene_set_order, sample_order]


def _display_labels(data: GSVAResult) -> dict[str, str]:
    return {gene_set: gene_set for gene_set in data.gene_set_order}


def create_score_heatmap(
    data: GSVAResult,
    condition_colors: dict[str, str],
    max_gene_sets: int = MAX_HEATMAP_GENE_SETS,
) -> BokehPlot:
    """Create a row-standardized heatmap of the most variable GSVA scores."""
    matrix = _score_matrix(data)
    variability = matrix.var(axis=1, ddof=1).fillna(0)
    selected = variability.sort_values(ascending=False, kind="mergesort").head(
        max_gene_sets
    ).index.tolist()
    matrix = matrix.loc[selected]
    means = matrix.mean(axis=1)
    standard_deviations = matrix.std(axis=1).replace(0, np.nan)
    z_scores = matrix.sub(means, axis=0).div(standard_deviations, axis=0).fillna(0)
    labels = _display_labels(data)
    melted = (
        z_scores.rename_axis("gene_set")
        .reset_index()
        .melt(id_vars="gene_set", var_name="sample", value_name="z_score")
    )
    melted["display_label"] = melted["gene_set"].map(labels)
    metadata = data.scores_long.drop_duplicates("sample").set_index("sample")
    sample_order = matrix.columns.tolist()
    limit = max(float(melted["z_score"].abs().max()), 1.0)
    mapper = LinearColorMapper(palette=RdBu11, low=-limit, high=limit)
    suffix = f" (top {len(selected)} by variance)" if len(data.gene_set_order) > len(selected) else ""
    heatmap = BokehPlot(
        title=f"GSVA scores across samples{suffix}",
        x_range=sample_order,
        y_range=list(reversed([labels[gene_set] for gene_set in selected])),
        x_axis_location="above",
        x_axis_label="Sample",
        y_axis_label="Gene set",
        height=max(440, 22 * len(selected) + 140),
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    rectangles = heatmap._fig.rect(
        x="sample",
        y="display_label",
        width=0.98,
        height=0.98,
        source=ColumnDataSource(melted),
        fill_color={"field": "z_score", "transform": mapper},
        line_color=None,
    )
    heatmap._fig.add_tools(
        HoverTool(
            renderers=[rectangles],
            tooltips=[
                ("Gene set", "@gene_set"),
                ("Sample", "@sample"),
                ("Row z-score", "@z_score{0.000}"),
            ],
        )
    )
    heatmap._fig.add_layout(ColorBar(color_mapper=mapper, title="Row z-score"), "right")
    heatmap._fig.xaxis.major_label_orientation = np.pi / 4
    heatmap._fig.grid.grid_line_color = None

    annotation = BokehPlot(
        x_range=heatmap._fig.x_range,
        y_range=(0, 1),
        height=90,
        sizing_mode="stretch_width",
        tools="",
    )
    annotation._fig.grid.visible = False
    annotation._fig.xaxis.visible = False
    annotation._fig.yaxis.visible = False
    for group in metadata["group"].drop_duplicates():
        samples = [sample for sample in sample_order if metadata.at[sample, "group"] == group]
        annotation._fig.rect(
            x=samples,
            y=[0.5] * len(samples),
            width=0.98,
            height=0.7,
            fill_color=condition_colors[group],
            line_color=None,
            legend_label=group,
        )
    annotation._fig.legend.title = "Condition"
    annotation._fig.legend.orientation = "horizontal"
    annotation._fig.legend.location = "center"

    combined = BokehPlot()
    combined._fig = gridplot(
        [[annotation._fig], [heatmap._fig]],
        toolbar_location="right",
        merge_tools=True,
        sizing_mode="stretch_width",
    )
    combined.report_height = heatmap._fig.height + annotation._fig.height + 40
    return combined


def create_score_distribution(
    data: GSVAResult,
    gene_set: str,
    condition_colors: dict[str, str],
) -> BokehPlot:
    """Create a box-and-points view of raw scores for one gene set."""
    subset = data.scores_long.loc[data.scores_long["gene_set"].eq(gene_set)].copy()
    if subset.empty:
        raise ValueError(f"Unknown GSVA gene set: {gene_set}")
    groups = subset["group"].drop_duplicates().tolist()
    labels = _display_labels(data)
    plot = BokehPlot(
        title=f"Raw GSVA scores — {labels[gene_set]}",
        x_axis_label="Group",
        y_axis_label="GSVA score",
        x_range=(0.5, len(groups) + 0.5),
        height=470,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    plot._fig.xaxis.ticker = FixedTicker(ticks=list(range(1, len(groups) + 1)))
    plot._fig.xaxis.major_label_overrides = {
        float(index): group for index, group in enumerate(groups, start=1)
    }
    point_rows = []
    for group_index, group in enumerate(groups, start=1):
        group_rows = subset.loc[subset["group"].eq(group)].copy()
        values = group_rows["score"].to_numpy(dtype=float)
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        iqr = q3 - q1
        lower = max(float(values.min()), float(q1 - 1.5 * iqr))
        upper = min(float(values.max()), float(q3 + 1.5 * iqr))
        color = condition_colors[group]
        plot._fig.segment([group_index], [lower], [group_index], [upper], line_color=color)
        plot._fig.vbar(
            x=[group_index],
            width=0.5,
            bottom=[q1],
            top=[q3],
            fill_color=color,
            fill_alpha=0.25,
            line_color=color,
        )
        plot._fig.segment(
            [group_index - 0.25],
            [median],
            [group_index + 0.25],
            [median],
            line_color=color,
            line_width=2,
        )
        offsets = np.linspace(-0.16, 0.16, len(group_rows)) if len(group_rows) > 1 else [0]
        for offset, row in zip(offsets, group_rows.itertuples(index=False), strict=True):
            point_rows.append(
                {
                    "x": group_index + float(offset),
                    "sample": row.sample,
                    "group": row.group,
                    "score": row.score,
                    "plot_color": color,
                }
            )
    source = ColumnDataSource(pd.DataFrame(point_rows))
    points = plot._fig.scatter(
        x="x",
        y="score",
        source=source,
        marker="circle",
        size=9,
        fill_color="plot_color",
        fill_alpha=0.9,
        line_color="white",
        line_width=1,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[points],
            tooltips=[("Sample", "@sample"), ("Group", "@group"), ("Score", "@score{0.000}")],
        )
    )
    plot._fig.grid.grid_line_alpha = 0.15
    return plot


def prepare_limma_volcano(
    analysis: GSVAContrastResult,
    condition_colors: dict[str, str],
    padj_cutoff: float,
) -> pd.DataFrame:
    """Prepare pathway-score differences and adjusted significance."""
    results = analysis.results.copy()
    positive = results["adjusted_p_value"].gt(0)
    floor = (
        max(float(results.loc[positive, "adjusted_p_value"].min()), np.finfo(float).tiny)
        if positive.any()
        else np.finfo(float).tiny
    )
    results["neg_log10_adjusted_p_value"] = -np.log10(
        results["adjusted_p_value"].clip(lower=floor)
    )
    results["significant"] = results["adjusted_p_value"].le(padj_cutoff)
    results["plot_color"] = NONSIGNIFICANT_COLOR
    results.loc[
        results["significant"] & results["effect_size"].gt(0), "plot_color"
    ] = condition_colors[analysis.contrast.target_group]
    results.loc[
        results["significant"] & results["effect_size"].lt(0), "plot_color"
    ] = condition_colors[analysis.contrast.reference_group]
    return results


def create_limma_volcano(
    analysis: GSVAContrastResult,
    condition_colors: dict[str, str],
    padj_cutoff: float,
) -> BokehPlot:
    """Create a GSVA-score-difference versus adjusted-significance plot."""
    prepared = prepare_limma_volcano(analysis, condition_colors, padj_cutoff)
    if prepared.empty:
        return _empty_plot(
            f"GSVA limma volcano — {analysis.label}",
            "No gene sets were tested by limma for this contrast.",
            height=460,
        )
    plot = BokehPlot(
        title=f"GSVA limma volcano — {analysis.label}",
        x_axis_label=(
            f"GSVA score difference ({analysis.contrast.target_group} − "
            f"{analysis.contrast.reference_group})"
        ),
        y_axis_label="-log10(BH-adjusted p-value)",
        height=500,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    source = ColumnDataSource(prepared)
    points = plot._fig.scatter(
        x="effect_size",
        y="neg_log10_adjusted_p_value",
        source=source,
        marker="circle",
        size=9,
        fill_color="plot_color",
        fill_alpha=0.85,
        line_color=None,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[points],
            tooltips=[
                ("Gene set", "@gene_set"),
                ("N genes", "@n_genes"),
                ("Score difference", "@effect_size{0.000}"),
                ("Adjusted p-value", "@adjusted_p_value{0.000e}"),
                ("P-value", "@p_value{0.000e}"),
                ("t statistic", "@t_statistic{0.000}"),
            ],
        )
    )
    plot._fig.add_layout(
        Span(location=0, dimension="height", line_color="#777777", line_width=1)
    )
    plot._fig.add_layout(
        Span(
            location=-np.log10(padj_cutoff),
            dimension="width",
            line_color="#555555",
            line_dash="dashed",
            line_width=1,
        )
    )
    return plot


def _select_limma_heatmap_gene_sets(
    analysis: GSVAContrastResult,
    padj_cutoff: float,
    top_n: int = 20,
) -> list[str]:
    significant = analysis.results.loc[
        analysis.results["adjusted_p_value"].le(padj_cutoff)
    ]
    positive = (
        significant.loc[significant["effect_size"].gt(0)]
        .sort_values(["effect_size", "adjusted_p_value"], ascending=[False, True])
        .head(top_n)
    )
    negative = (
        significant.loc[significant["effect_size"].lt(0)]
        .sort_values(["effect_size", "adjusted_p_value"], ascending=[True, True])
        .head(top_n)
    )
    return pd.concat([negative, positive])["gene_set"].tolist()


def create_limma_heatmap(
    data: GSVAResult,
    analysis: GSVAContrastResult,
    condition_colors: dict[str, str],
    padj_cutoff: float,
) -> BokehPlot:
    """Create a raw-score heatmap for limma-significant gene sets."""
    selected = _select_limma_heatmap_gene_sets(analysis, padj_cutoff)
    if not selected:
        return _empty_plot(
            f"Differential GSVA scores — {analysis.label}",
            f"No gene sets pass BH-adjusted p-value ≤ {padj_cutoff:g}.",
            height=440,
        )
    matrix = _score_matrix(data)
    metadata = data.scores_long.drop_duplicates("sample").set_index("sample")
    contrast_groups = {
        analysis.contrast.reference_group,
        analysis.contrast.target_group,
    }
    sample_order = metadata.index[metadata["group"].isin(contrast_groups)].tolist()
    matrix = matrix.loc[selected, sample_order]
    z_scores = matrix.sub(matrix.mean(axis=1), axis=0).div(
        matrix.std(axis=1).replace(0, np.nan), axis=0
    ).fillna(0)
    labels = _display_labels(data)
    melted = (
        z_scores.rename_axis("gene_set")
        .reset_index()
        .melt(id_vars="gene_set", var_name="sample", value_name="z_score")
    )
    melted["display_label"] = melted["gene_set"].map(labels)
    limit = max(float(melted["z_score"].abs().max()), 1.0)
    mapper = LinearColorMapper(palette=RdBu11, low=-limit, high=limit)
    heatmap = BokehPlot(
        title=f"Differential GSVA scores — {analysis.label}",
        x_range=sample_order,
        y_range=list(reversed([labels[gene_set] for gene_set in selected])),
        x_axis_location="above",
        x_axis_label="Sample",
        y_axis_label="Gene set",
        height=max(440, 22 * len(selected) + 140),
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    rectangles = heatmap._fig.rect(
        x="sample",
        y="display_label",
        width=0.98,
        height=0.98,
        source=ColumnDataSource(melted),
        fill_color={"field": "z_score", "transform": mapper},
        line_color=None,
    )
    heatmap._fig.add_tools(
        HoverTool(
            renderers=[rectangles],
            tooltips=[
                ("Gene set", "@gene_set"),
                ("Sample", "@sample"),
                ("Row z-score", "@z_score{0.000}"),
            ],
        )
    )
    heatmap._fig.add_layout(ColorBar(color_mapper=mapper, title="Row z-score"), "right")
    heatmap._fig.xaxis.major_label_orientation = np.pi / 4
    heatmap._fig.grid.grid_line_color = None
    heatmap.report_height = heatmap._fig.height
    return heatmap


def create_multi_contrast_dot_plot(
    data: GSVAResult,
    padj_cutoff: float,
    top_n: int = MAX_DIFFERENTIAL_GENE_SETS,
) -> BokehPlot:
    """Create an effect/significance summary across all limma contrasts."""
    combined = pd.concat(
        [
            analysis.results.assign(contrast=analysis.label)
            for analysis in data.contrasts
        ],
        ignore_index=True,
    )
    if combined.empty:
        return _empty_plot(
            "GSVA limma across contrasts",
            "No gene sets were tested by limma.",
            height=440,
        )
    ranking = (
        combined.groupby("gene_set", sort=False)["adjusted_p_value"]
        .min()
        .sort_values(kind="mergesort")
        .head(top_n)
    )
    prepared = combined.loc[combined["gene_set"].isin(ranking.index)].copy()
    order = ranking.index.tolist()
    label_map = (
        prepared.drop_duplicates("gene_set").set_index("gene_set")["display_label"].to_dict()
    )
    prepared["display_label"] = prepared["gene_set"].map(label_map)
    positive = prepared["adjusted_p_value"].gt(0)
    floor = (
        max(float(prepared.loc[positive, "adjusted_p_value"].min()), np.finfo(float).tiny)
        if positive.any()
        else np.finfo(float).tiny
    )
    significance = -np.log10(prepared["adjusted_p_value"].clip(lower=floor))
    prepared["point_size"] = 7 + 11 * significance.clip(upper=20) / 20
    prepared["significance_label"] = significance
    limit = max(float(prepared["effect_size"].abs().max()), 0.1)
    mapper = LinearColorMapper(palette=RdBu11, low=-limit, high=limit)
    plot_height = max(460, 24 * len(order) + 160)
    plot = BokehPlot(
        title=(
            "GSVA limma across contrasts"
            + (f" (top {len(order)} by adjusted p-value)" if len(ranking) < combined["gene_set"].nunique() else "")
        ),
        x_range=[analysis.label for analysis in data.contrasts],
        y_range=list(reversed([label_map[gene_set] for gene_set in order])),
        x_axis_label="Contrast",
        y_axis_label="Gene set",
        height=plot_height,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    points = plot._fig.scatter(
        x="contrast",
        y="display_label",
        source=ColumnDataSource(prepared),
        marker="circle",
        size="point_size",
        fill_color={"field": "effect_size", "transform": mapper},
        fill_alpha=0.9,
        line_color=None,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[points],
            tooltips=[
                ("Contrast", "@contrast"),
                ("Gene set", "@gene_set"),
                ("Score difference", "@effect_size{0.000}"),
                ("Adjusted p-value", "@adjusted_p_value{0.000e}"),
            ],
        )
    )
    plot._fig.add_layout(ColorBar(color_mapper=mapper, title="Score difference"), "right")
    plot._fig.xaxis.major_label_orientation = np.pi / 4
    plot._fig.grid.grid_line_alpha = 0.15
    plot.report_height = plot_height
    return plot


def add_gsva_scores(
    data: GSVAResult,
    condition_colors: dict[str, str],
) -> None:
    """Add raw-score plot-type tabs and the long gene-set dropdown."""
    tabs = Tabs()
    with tabs.add_tab("Heatmap"):
        heatmap = create_score_heatmap(data, condition_colors)
        EZChart(
            heatmap,
            "epi2melabs",
            height=f"{getattr(heatmap, 'report_height', 720)}px",
        )
    with tabs.add_tab("Distributions"):
        distribution_tabs = Tabs()
        with distribution_tabs.add_dropdown_menu("Gene set", change_header=True):  # type: ignore
            for gene_set in data.gene_set_order:
                with distribution_tabs.add_dropdown_tab(gene_set):  # type: ignore
                    EZChart(
                        create_score_distribution(data, gene_set, condition_colors),
                        "epi2melabs",
                        height="500px",
                    )
    with tabs.add_tab("Coverage"):
        DataTable.from_pandas(
            data.coverage.drop(columns="description"),
            use_index=False,
        )


def add_gsva_differential(
    data: GSVAResult,
    condition_colors: dict[str, str],
    padj_cutoff: float,
) -> None:
    """Add across-contrast, contrast, and plot-type tabs for GSVA limma."""
    tabs = Tabs()
    with tabs.add_tab("Across contrasts"):
        summary = create_multi_contrast_dot_plot(data, padj_cutoff)
        EZChart(
            summary,
            "epi2melabs",
            height=f"{getattr(summary, 'report_height', 600)}px",
        )
    with tabs.add_tab("Contrasts"):
        contrast_tabs = Tabs()
        for analysis in data.contrasts:
            with contrast_tabs.add_tab(analysis.label):
                plot_tabs = Tabs()
                with plot_tabs.add_tab("Volcano"):
                    EZChart(
                        create_limma_volcano(
                            analysis,
                            condition_colors,
                            padj_cutoff,
                        ),
                        "epi2melabs",
                        height="530px",
                    )
                with plot_tabs.add_tab("Score heatmap"):
                    heatmap = create_limma_heatmap(
                        data,
                        analysis,
                        condition_colors,
                        padj_cutoff,
                    )
                    EZChart(
                        heatmap,
                        "epi2melabs",
                        height=f"{getattr(heatmap, 'report_height', 620)}px",
                    )
