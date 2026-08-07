"""Parse fry outputs and add gene-set enrichment plots to the report."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

from bokeh.layouts import column
from bokeh.models import (
    BoxAnnotation,
    ColumnDataSource,
    Div,
    FixedTicker,
    HoverTool,
    Range1d,
    Span,
)
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Tabs
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd

from .differential_plots import (
    _display_labels,
    _empty_plot,
    ContrastResult,
    DifferentialResult,
    NONSIGNIFICANT_COLOR,
)


FRY_RESULTS_FILE = "fry_results.tsv"
GENE_SET_RESOLUTION_FILE = "gene_set_resolution.tsv"
REQUIRED_FRY_COLUMNS = (
    "gene_set",
    "NGenes",
    "Direction",
    "PValue",
    "FDR",
    "PValue.Mixed",
    "FDR.Mixed",
)
REQUIRED_RESOLUTION_COLUMNS = ("gene_set", "feature_id")
COUNT_COLUMNS = (
    "gmt_members",
    "matched_gmt_members",
    "count_matrix_members",
    "tested_members",
    "tested_gmt_members",
)
COVERAGE_COLUMNS = ("count_matrix_coverage", "tested_coverage")
MAX_SIGNED_GENE_SETS = 30
WORM_SPAN = 0.45
WORM_COLOR = "#333333"


@dataclass
class GeneSetContrastResult:
    """Validated fry results and retained members for one contrast."""

    contrast: ContrastResult
    fry_results: pd.DataFrame
    members: dict[str, tuple[str, ...]]

    @property
    def label(self) -> str:
        """Human-readable contrast label."""
        return self.contrast.label


def _nonempty_string(value) -> str:
    """Return a stripped string, treating common missing values as empty."""
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"", "na", "nan"}:
        return ""
    return value


def _gene_set_display_labels(fry_results: pd.DataFrame) -> pd.Series:
    """Prefer descriptions and disambiguate duplicate display labels."""
    if "description" in fry_results:
        labels = fry_results["description"].map(_nonempty_string)
    else:
        labels = pd.Series("", index=fry_results.index, dtype=object)
    labels = labels.where(labels.ne(""), fry_results["gene_set"])
    duplicates = labels.duplicated(keep=False)
    labels.loc[duplicates] = (
        labels.loc[duplicates]
        + " ["
        + fry_results.loc[duplicates, "gene_set"].astype(str)
        + "]"
    )
    return labels


def _gene_set_selector_labels(fry_results: pd.DataFrame) -> dict[str, str]:
    """Use the GMT first-column identifier as each dropdown label."""
    gene_sets = fry_results["gene_set"].astype(str)
    return dict(zip(gene_sets, gene_sets, strict=True))


def _validate_probability_columns(
    table: pd.DataFrame,
    columns: tuple[str, ...],
    path: Path,
) -> None:
    """Require finite probabilities in the closed unit interval."""
    for column_name in columns:
        table[column_name] = pd.to_numeric(table[column_name], errors="coerce")
        values = table[column_name].to_numpy(dtype=float)
        if (
            np.isnan(values).any()
            or not np.isfinite(values).all()
            or (values < 0).any()
            or (values > 1).any()
        ):
            raise ValueError(
                f"{path} contains invalid {column_name} values; expected finite "
                "numbers in [0, 1]."
            )


def _read_resolution(results_dir: Path) -> pd.DataFrame:
    """Read the analysis-level gene-set identifier resolution table."""
    resolution_path = results_dir / GENE_SET_RESOLUTION_FILE
    if not resolution_path.is_file():
        raise ValueError(
            f"Missing edgeR gene-set resolution table: {resolution_path}"
        )
    resolution = pd.read_csv(resolution_path, sep="\t")
    missing = [
        column
        for column in REQUIRED_RESOLUTION_COLUMNS
        if column not in resolution
    ]
    if missing:
        raise ValueError(
            f"{resolution_path} is missing columns: " + ", ".join(missing)
        )
    for column_name in REQUIRED_RESOLUTION_COLUMNS:
        resolution[column_name] = resolution[column_name].map(_nonempty_string)
        if resolution[column_name].eq("").any():
            raise ValueError(
                f"{resolution_path} contains empty {column_name} values."
            )
    return resolution.drop_duplicates(["gene_set", "feature_id"])


def _read_fry_results(path: Path) -> pd.DataFrame:
    """Read and validate a single fry result table."""
    if not path.is_file():
        raise ValueError(f"Missing fry result table: {path}")
    fry_results = pd.read_csv(path, sep="\t")
    missing = [
        column for column in REQUIRED_FRY_COLUMNS if column not in fry_results
    ]
    if missing:
        raise ValueError(f"{path} is missing columns: " + ", ".join(missing))

    fry_results["gene_set"] = fry_results["gene_set"].map(_nonempty_string)
    if fry_results["gene_set"].eq("").any():
        raise ValueError(f"{path} contains empty gene_set values.")
    if fry_results["gene_set"].duplicated().any():
        duplicates = sorted(
            fry_results.loc[
                fry_results["gene_set"].duplicated(keep=False),
                "gene_set",
            ].unique()
        )
        raise ValueError(
            f"{path} contains duplicate gene_set values: "
            + ", ".join(duplicates)
        )

    invalid_directions = sorted(
        set(fry_results["Direction"].dropna().astype(str)) - {"Up", "Down"}
    )
    if fry_results["Direction"].isna().any() or invalid_directions:
        details = ", ".join(invalid_directions) or "missing values"
        raise ValueError(
            f"{path} contains invalid Direction values: {details}; expected "
            "Up or Down."
        )

    _validate_probability_columns(
        fry_results,
        ("PValue", "FDR", "PValue.Mixed", "FDR.Mixed"),
        path,
    )

    n_genes = pd.to_numeric(fry_results["NGenes"], errors="coerce")
    if (
        n_genes.isna().any()
        or not np.isfinite(n_genes.to_numpy(dtype=float)).all()
        or n_genes.lt(2).any()
        or not np.equal(n_genes, np.floor(n_genes)).all()
    ):
        raise ValueError(
            f"{path} contains invalid NGenes values; expected integers of at "
            "least two."
        )
    fry_results["NGenes"] = n_genes.astype(int)

    for column_name in COUNT_COLUMNS:
        if column_name not in fry_results:
            fry_results[column_name] = np.nan
            continue
        values = pd.to_numeric(fry_results[column_name], errors="coerce")
        invalid = values.notna() & (
            ~np.isfinite(values)
            | values.lt(0)
            | ~np.equal(values, np.floor(values))
        )
        if invalid.any():
            raise ValueError(
                f"{path} contains invalid {column_name} values; expected "
                "nonnegative integers."
            )
        fry_results[column_name] = values

    for column_name in COVERAGE_COLUMNS:
        if column_name not in fry_results:
            fry_results[column_name] = np.nan
            continue
        values = pd.to_numeric(fry_results[column_name], errors="coerce")
        invalid = values.notna() & (
            ~np.isfinite(values) | values.lt(0) | values.gt(1)
        )
        if invalid.any():
            raise ValueError(
                f"{path} contains invalid {column_name} values; expected "
                "numbers in [0, 1]."
            )
        fry_results[column_name] = values

    if "description" not in fry_results:
        fry_results["description"] = ""
    fry_results["display_label"] = _gene_set_display_labels(fry_results)
    return fry_results


def load_gene_set_results(
    results_dir: str | Path,
    differential: DifferentialResult,
) -> list[GeneSetContrastResult]:
    """Load fry tables and validate their retained feature memberships."""
    results_dir = Path(results_dir)
    resolution = _read_resolution(results_dir)
    analyses = []

    for contrast in differential.contrasts:
        fry_path = results_dir / contrast.contrast_id / FRY_RESULTS_FILE
        fry_results = _read_fry_results(fry_path)
        retained_features = set(contrast.results["feature_id"])
        members = {}
        for row in fry_results.itertuples(index=False):
            resolved = set(
                resolution.loc[
                    resolution["gene_set"].eq(row.gene_set),
                    "feature_id",
                ]
            )
            tested = tuple(sorted(resolved & retained_features))
            if len(tested) != row.NGenes:
                raise ValueError(
                    f"{fry_path}: gene set '{row.gene_set}' reports "
                    f"NGenes={row.NGenes}, but {len(tested)} unique resolved "
                    "features occur in the corresponding edgeR result."
                )
            members[row.gene_set] = tested

        analyses.append(
            GeneSetContrastResult(
                contrast=contrast,
                fry_results=fry_results,
                members=members,
            )
        )
    return analyses


def prepare_signed_gene_sets(
    analysis: GeneSetContrastResult,
    condition_colors: dict[str, str],
    padj_cutoff: float,
    top_n: int = MAX_SIGNED_GENE_SETS,
) -> pd.DataFrame:
    """Prepare signed directional-fry values and report colors."""
    results = (
        analysis.fry_results.sort_values(
            ["FDR", "gene_set"],
            kind="mergesort",
        )
        .head(top_n)
        .copy()
    )
    if results.empty:
        results["signed_significance"] = pd.Series(dtype=float)
        results["plot_color"] = pd.Series(dtype=object)
        return results

    positive_fdr = analysis.fry_results.loc[
        analysis.fry_results["FDR"].gt(0),
        "FDR",
    ]
    fdr_floor = (
        max(float(positive_fdr.min()), np.finfo(float).tiny)
        if not positive_fdr.empty
        else np.finfo(float).tiny
    )
    magnitude = -np.log10(results["FDR"].clip(lower=fdr_floor))
    results["signed_significance"] = magnitude.where(
        results["Direction"].eq("Up"),
        -magnitude,
    )
    significant = results["FDR"].le(padj_cutoff)
    results["plot_color"] = NONSIGNIFICANT_COLOR
    results.loc[
        significant & results["Direction"].eq("Up"),
        "plot_color",
    ] = condition_colors[analysis.contrast.target_group]
    results.loc[
        significant & results["Direction"].eq("Down"),
        "plot_color",
    ] = condition_colors[analysis.contrast.reference_group]
    results["direction_label"] = np.where(
        results["Direction"].eq("Up"),
        f"Up in {analysis.contrast.target_group}",
        f"Up in {analysis.contrast.reference_group}",
    )
    return results


def create_signed_significance_plot(
    analysis: GeneSetContrastResult,
    condition_colors: dict[str, str],
    padj_cutoff: float,
    top_n: int = MAX_SIGNED_GENE_SETS,
) -> BokehPlot:
    """Create a signed -log10 directional-FDR summary."""
    prepared = prepare_signed_gene_sets(
        analysis,
        condition_colors,
        padj_cutoff,
        top_n,
    )
    if prepared.empty:
        return _empty_plot(
            f"fry directional enrichment — {analysis.label}",
            "No gene sets were tested for this contrast.",
            height=420,
        )

    truncated = len(analysis.fry_results) > top_n
    suffix = (
        f" (top {top_n} of {len(analysis.fry_results)})"
        if truncated
        else ""
    )
    threshold = -np.log10(padj_cutoff)
    extent = max(
        threshold,
        float(prepared["signed_significance"].abs().max()),
        1.0,
    ) * 1.12
    plot_height = max(420, 27 * len(prepared) + 150)
    plot = BokehPlot(
        title=f"fry signed directional significance — {analysis.label}{suffix}",
        x_axis_label="Signed -log10(directional FDR)",
        y_axis_label="Gene set",
        x_range=(-extent, extent),
        y_range=list(reversed(prepared["gene_set"].tolist())),
        height=plot_height,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    source = ColumnDataSource(prepared)
    bars = plot._fig.hbar(
        y="gene_set",
        right="signed_significance",
        left=0,
        height=0.72,
        fill_color="plot_color",
        fill_alpha=0.85,
        line_color=None,
        source=source,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[bars],
            tooltips=[
                ("Gene set", "@gene_set"),
                ("Description", "@description"),
                ("Direction", "@direction_label"),
                ("N genes", "@NGenes"),
                ("P-value", "@PValue{0.000e}"),
                ("Directional FDR", "@FDR{0.000e}"),
                ("Mixed P-value", "@{PValue.Mixed}{0.000e}"),
                ("Mixed FDR", "@{FDR.Mixed}{0.000e}"),
                ("GMT members", "@gmt_members{0,0}"),
                ("Retained members", "@tested_members{0,0}"),
                ("Count-matrix coverage", "@count_matrix_coverage{0.0%}"),
                ("Tested coverage", "@tested_coverage{0.0%}"),
            ],
        )
    )
    for location in (-threshold, threshold):
        plot._fig.add_layout(
            Span(
                location=location,
                dimension="height",
                line_color="#555555",
                line_dash="dashed",
                line_width=1,
            )
        )
    plot._fig.add_layout(
        Span(location=0, dimension="height", line_color="#888888", line_width=1)
    )
    plot.report_height = plot_height
    return plot


def tricube_moving_average(
    values: np.ndarray | list[float],
    span: float = WORM_SPAN,
    power: float = 3,
) -> np.ndarray:
    """Port limma::tricubeMovingAverage, including its edge correction."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise ValueError("tricube_moving_average expects a one-dimensional vector.")
    if not np.isfinite(values).all():
        raise ValueError("tricube_moving_average values must be finite.")
    if values.size == 0:
        return values.copy()
    span = min(float(span), 1.0)
    if span <= 0:
        return values.copy()
    power = max(float(power), 0.0)

    width = span * len(values)
    half_width = int(width // 2)
    width = 2 * half_width + 1
    if width > len(values):
        width -= 2
        half_width -= 1
    if half_width <= 0:
        return values.copy()

    positions = np.linspace(-1, 1, width) * width / (width + 1)
    weights = (1 - np.abs(positions) ** 3) ** power
    weights /= weights.sum()
    padded = np.pad(values, half_width, mode="constant")
    smoothed = np.convolve(padded, weights, mode="valid")

    cumulative = np.cumsum(weights)
    leading = cumulative[width - half_width - 1:width - 1]
    smoothed[:half_width] /= leading
    smoothed[-half_width:] /= leading[::-1]
    return smoothed


def _rank_contrast_features(
    analysis: GeneSetContrastResult,
    condition_colors: dict[str, str],
) -> pd.DataFrame:
    """Rank finite gene statistics deterministically by logFC then feature ID."""
    ranked = analysis.contrast.results.copy()
    ranked["display_label"] = _display_labels(ranked)
    ranked = (
        ranked.loc[np.isfinite(ranked["logFC"])]
        .sort_values(["logFC", "feature_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    ranked["plot_color"] = NONSIGNIFICANT_COLOR
    ranked.loc[
        ranked["logFC"].lt(0),
        "plot_color",
    ] = condition_colors[analysis.contrast.reference_group]
    ranked.loc[
        ranked["logFC"].gt(0),
        "plot_color",
    ] = condition_colors[analysis.contrast.target_group]
    return ranked


def _format_optional_number(value, format_spec: str = ".3g") -> str:
    """Format optional numeric metadata for the self-contained summary div."""
    if pd.isna(value):
        return "not available"
    return format(float(value), format_spec)


def _gene_set_summary_html(row: pd.Series) -> str:
    """Build escaped HTML describing one selected fry result."""
    return (
        "<details class=\"mb-3\">"
        "<summary><strong>Gene-set details</strong></summary>"
        f"<h4 class=\"mt-2\">{escape(str(row['display_label']))}</h4>"
        f"<p><code>{escape(str(row['gene_set']))}</code><br>"
        f"Direction: <strong>{escape(str(row['Direction']))}</strong>; "
        f"NGenes: <strong>{int(row['NGenes'])}</strong>; "
        f"directional FDR: <strong>{row['FDR']:.3g}</strong>; "
        f"mixed FDR: <strong>{row['FDR.Mixed']:.3g}</strong>; "
        "tested coverage: "
        f"<strong>{escape(_format_optional_number(row['tested_coverage'], '.1%'))}"
        "</strong>.</p>"
        "</details>"
    )


def _barcode_payload(
    analysis: GeneSetContrastResult,
    ranked: pd.DataFrame,
) -> dict[str, dict]:
    """Build all dropdown states for one self-contained barcode plot."""
    payload = {}
    feature_ids = ranked["feature_id"]
    for _, row in analysis.fry_results.sort_values(
        ["FDR", "gene_set"],
        kind="mergesort",
    ).iterrows():
        member_mask = feature_ids.isin(analysis.members[row["gene_set"]])
        member_rows = ranked.loc[member_mask]
        if member_rows.empty:
            continue
        indicator = member_mask.to_numpy(dtype=float)
        mean_membership = indicator.mean()
        worm = tricube_moving_average(indicator) / mean_membership
        y_end = max(2.1, float(worm.max()) * 1.08)
        payload[row["gene_set"]] = {
            "bars": {
                "rank": member_rows["rank"].astype(int).tolist(),
                "bottom": [0.0] * len(member_rows),
                "top": [0.42] * len(member_rows),
                "feature_id": member_rows["feature_id"].astype(str).tolist(),
                "display_label": member_rows["display_label"].astype(str).tolist(),
                "logFC": member_rows["logFC"].astype(float).tolist(),
                "logCPM": member_rows["logCPM"].astype(float).tolist(),
                "PValue": member_rows["PValue"].astype(float).tolist(),
                "FDR": member_rows["FDR"].astype(float).tolist(),
                "plot_color": member_rows["plot_color"].astype(str).tolist(),
            },
            "worm": {
                "rank": ranked["rank"].astype(int).tolist(),
                "relative_enrichment": worm.astype(float).tolist(),
            },
            "summary": _gene_set_summary_html(row),
            "y_end": y_end,
        }
    return payload


def _create_barcode_plot(
    analysis: GeneSetContrastResult,
    condition_colors: dict[str, str],
    ranked: pd.DataFrame,
    gene_set: str,
    state: dict,
) -> BokehPlot:
    """Create one limma-style barcode/worm plot for a selected gene set."""
    row = analysis.fry_results.set_index("gene_set").loc[gene_set]
    details = Div(text=state["summary"], sizing_mode="stretch_width")
    bar_source = ColumnDataSource(state["bars"])
    worm_source = ColumnDataSource(state["worm"])

    n_features = len(ranked)
    y_range = Range1d(0, state["y_end"])
    plot = BokehPlot(
        title=(
            "Gene-set barcode and enrichment worm — "
            f"{row['display_label']} — {analysis.label}"
        ),
        x_axis_label=(
            f"{analysis.contrast.reference_group} (negative logFC) \u2190 "
            "genes ranked by edgeR log2 fold change \u2192 "
            f"{analysis.contrast.target_group} (positive logFC)"
        ),
        y_axis_label="Relative local enrichment",
        x_range=Range1d(0.5, n_features + 0.5),
        y_range=y_range,
        height=430,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )

    zero_location = (
        int(np.searchsorted(ranked["logFC"].to_numpy(dtype=float), 0, side="left"))
        + 0.5
    )
    plot._fig.add_layout(
        BoxAnnotation(
            left=0.5,
            right=zero_location,
            fill_color=condition_colors[analysis.contrast.reference_group],
            fill_alpha=0.07,
            line_alpha=0,
        )
    )
    plot._fig.add_layout(
        BoxAnnotation(
            left=zero_location,
            right=n_features + 0.5,
            fill_color=condition_colors[analysis.contrast.target_group],
            fill_alpha=0.07,
            line_alpha=0,
        )
    )
    bars = plot._fig.vbar(
        x="rank",
        width=max(0.8, n_features / 1200),
        bottom="bottom",
        top="top",
        fill_color="plot_color",
        fill_alpha=0.9,
        line_color=None,
        source=bar_source,
    )
    plot._fig.line(
        x="rank",
        y="relative_enrichment",
        source=worm_source,
        line_color=WORM_COLOR,
        line_width=2.2,
    )
    plot._fig.add_tools(
        HoverTool(
            renderers=[bars],
            tooltips=[
                ("Feature", "@display_label"),
                ("Feature ID", "@feature_id"),
                ("Rank", "@rank{0,0}"),
                ("logFC", "@logFC{0.000}"),
                ("logCPM", "@logCPM{0.000}"),
                ("P-value", "@PValue{0.000e}"),
                ("FDR", "@FDR{0.000e}"),
            ],
        )
    )
    plot._fig.add_layout(
        Span(
            location=1,
            dimension="width",
            line_color="#666666",
            line_dash="dashed",
            line_width=1,
        )
    )
    plot._fig.add_layout(
        Span(
            location=zero_location,
            dimension="height",
            line_color="#888888",
            line_dash="dotted",
            line_width=1,
        )
    )
    tick_indices = np.unique(
        np.rint(np.linspace(0, n_features - 1, min(5, n_features))).astype(int)
    )
    tick_ranks = [int(index + 1) for index in tick_indices]
    plot._fig.xaxis.ticker = FixedTicker(ticks=tick_ranks)
    plot._fig.xaxis.major_label_overrides = {
        float(rank): f"{ranked.iloc[index]['logFC']:.2g}"
        for rank, index in zip(tick_ranks, tick_indices, strict=True)
    }
    plot._fig.grid.grid_line_alpha = 0.15

    combined = BokehPlot()
    combined._fig = column(
        details,
        plot._fig,
        sizing_mode="stretch_width",
    )
    combined.report_height = 570
    return combined


def _barcode_plot_inputs(
    analysis: GeneSetContrastResult,
    condition_colors: dict[str, str],
):
    """Prepare shared ranked data and ordered gene-set states."""
    if analysis.fry_results.empty:
        return None, None, []

    ranked = _rank_contrast_features(analysis, condition_colors)
    if ranked.empty:
        return ranked, None, []
    payload = _barcode_payload(analysis, ranked)
    ordered_sets = [
        gene_set
        for gene_set in analysis.fry_results.sort_values(
            ["FDR", "gene_set"],
            kind="mergesort",
        )["gene_set"]
        if gene_set in payload
    ]
    return ranked, payload, ordered_sets


def create_barcode_plot(
    analysis: GeneSetContrastResult,
    condition_colors: dict[str, str],
    gene_set: str | None = None,
) -> BokehPlot:
    """Create a barcode/worm plot for one gene set."""
    ranked, payload, ordered_sets = _barcode_plot_inputs(
        analysis,
        condition_colors,
    )
    if analysis.fry_results.empty:
        return _empty_plot(
            f"Gene-set barcode — {analysis.label}",
            "No gene sets were tested for this contrast.",
            height=480,
        )
    if ranked is not None and ranked.empty:
        return _empty_plot(
            f"Gene-set barcode — {analysis.label}",
            "No finite gene-level logFC values are available for ranking.",
            height=480,
        )
    if not payload or not ordered_sets:
        return _empty_plot(
            f"Gene-set barcode — {analysis.label}",
            "No resolved gene-set members are available for the ranked genes.",
            height=480,
        )

    gene_set = gene_set or ordered_sets[0]
    if gene_set not in payload:
        raise ValueError(f"Unknown or unresolved gene set: {gene_set}")
    return _create_barcode_plot(
        analysis,
        condition_colors,
        ranked,
        gene_set,
        payload[gene_set],
    )


def add_gene_set_barcode_dropdown(
    analysis: GeneSetContrastResult,
    condition_colors: dict[str, str],
) -> None:
    """Add barcode plots behind the same dropdown tabs used for QC samples."""
    ranked, payload, ordered_sets = _barcode_plot_inputs(
        analysis,
        condition_colors,
    )
    if not payload or not ordered_sets:
        barcode = create_barcode_plot(analysis, condition_colors)
        EZChart(
            barcode,
            "epi2melabs",
            height=f"{getattr(barcode, 'report_height', 590)}px",
        )
        return

    selector_labels = _gene_set_selector_labels(analysis.fry_results)
    tabs = Tabs()
    with tabs.add_dropdown_menu("Gene set", change_header=True):  # type: ignore
        for gene_set in ordered_sets:
            with tabs.add_dropdown_tab(selector_labels[gene_set]):  # type: ignore
                barcode = _create_barcode_plot(
                    analysis,
                    condition_colors,
                    ranked,
                    gene_set,
                    payload[gene_set],
                )
                EZChart(
                    barcode,
                    "epi2melabs",
                    height=f"{getattr(barcode, 'report_height', 590)}px",
                )


def add_gene_set_enrichment(
    analyses: list[GeneSetContrastResult],
    condition_colors: dict[str, str],
    padj_cutoff: float,
) -> None:
    """Add contrast and plot-type tabs for directional fry enrichment."""
    contrast_tabs = Tabs()
    for analysis in analyses:
        with contrast_tabs.add_tab(analysis.label):
            plot_tabs = Tabs()
            with plot_tabs.add_tab("Directional summary"):
                signed = create_signed_significance_plot(
                    analysis,
                    condition_colors,
                    padj_cutoff,
                )
                EZChart(
                    signed,
                    "epi2melabs",
                    height=f"{getattr(signed, 'report_height', 520)}px",
                )
            with plot_tabs.add_tab("Barcode"):
                add_gene_set_barcode_dropdown(analysis, condition_colors)
