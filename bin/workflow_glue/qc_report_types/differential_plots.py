"""Parse edgeR outputs and add differential-analysis plots to the report."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from bokeh.layouts import gridplot
from bokeh.models import (
    ColorBar,
    ColumnDataSource,
    HoverTool,
    Label,
    LinearColorMapper,
    Span,
)
from bokeh.palettes import RdBu11, Turbo256
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Tabs
from ezcharts.plots import BokehPlot
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from sklearn.decomposition import PCA


FEATURE_COUNTS_FILE = "feature_counts.tsv"
EDGE_R_RESULTS_FILE = "edgeR_results.tsv"
CONTRAST_DIRECTORY = re.compile(r"^group_(.+)_vs_(.+)$")
REQUIRED_RESULT_COLUMNS = ("feature_id", "logFC", "logCPM", "PValue", "FDR")
LABEL_COLUMNS = ("gene", "gene_id", "locus_tag", "feature_id")
NONSIGNIFICANT_COLOR = "#B8B8B8"
CONDITION_PALETTE = (
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#F0E442",
    "#000000",
    "#332288",
    "#88CCEE",
    "#44AA99",
    "#117733",
    "#999933",
    "#DDCC77",
    "#CC6677",
    "#882255",
    "#AA4499",
)


@dataclass
class ContrastResult:
    """One target-versus-reference edgeR result."""

    contrast_id: str
    target_group: str
    reference_group: str
    results: pd.DataFrame

    @property
    def label(self) -> str:
        """Human-readable contrast label."""
        return f"{self.target_group} vs {self.reference_group}"


@dataclass
class DifferentialResult:
    """Validated inputs needed by all differential report plots."""

    feature_counts: pd.DataFrame
    sample_metadata: pd.DataFrame
    contrasts: list[ContrastResult]
    condition_colors: dict[str, str]


def safe_name(value: str) -> str:
    """Match the workflow/R safe-name convention used for contrast directories."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", value)


def condition_palette(groups: list[str]) -> dict[str, str]:
    """Create a deterministic condition-to-color mapping."""
    ordered_groups = sorted(set(groups), key=lambda group: (group != "control", group))
    colors = list(CONDITION_PALETTE)
    if len(ordered_groups) > len(colors):
        remaining = len(ordered_groups) - len(colors)
        colors.extend(
            Turbo256[index * 255 // max(remaining - 1, 1)] for index in range(remaining)
        )
    return dict(zip(ordered_groups, colors[: len(ordered_groups)], strict=True))


def _normalise_sample_metadata(samples_df: pd.DataFrame) -> pd.DataFrame:
    """Return sample metadata indexed by unique sample name."""
    if {"Name", "Group"}.issubset(samples_df.columns):
        metadata = samples_df[["Name", "Group"]].rename(
            columns={"Name": "name", "Group": "group"}
        )
    elif {"name", "group"}.issubset(samples_df.columns):
        metadata = samples_df[["name", "group"]].copy()
    else:
        raise ValueError(
            "Sample metadata must contain either Name/Group or name/group columns."
        )

    metadata = metadata.astype(str)
    if metadata["name"].eq("").any() or metadata["group"].eq("").any():
        raise ValueError("Sample names and conditions must not be empty.")
    duplicated = metadata.loc[metadata["name"].duplicated(), "name"].tolist()
    if duplicated:
        raise ValueError(
            "Sample names must be unique for differential analysis: "
            + ", ".join(sorted(set(duplicated)))
        )
    return metadata.set_index("name", drop=False)


def _resolve_contrast_groups(
    directory_name: str,
    groups: list[str],
) -> tuple[str, str]:
    """Resolve a sanitized contrast directory to original condition labels."""
    if not CONTRAST_DIRECTORY.fullmatch(directory_name):
        raise ValueError(f"Invalid edgeR contrast directory name: {directory_name}")

    matches = [
        (target, reference)
        for target in groups
        for reference in groups
        if target != reference
        and directory_name == f"group_{safe_name(target)}_vs_{safe_name(reference)}"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Could not uniquely map edgeR contrast directory '{directory_name}' "
            "to the sample conditions."
        )
    return matches[0]


def _read_contrasts(
    results_dir: Path,
    groups: list[str],
) -> list[ContrastResult]:
    """Read and validate every edgeR contrast table."""
    contrast_dirs = sorted(
        path
        for path in results_dir.iterdir()
        if path.is_dir() and (path / EDGE_R_RESULTS_FILE).is_file()
    )
    if not contrast_dirs:
        raise ValueError(
            f"No contrast directories containing {EDGE_R_RESULTS_FILE} were found "
            f"in {results_dir}."
        )

    contrasts = []
    for contrast_dir in contrast_dirs:
        target_group, reference_group = _resolve_contrast_groups(
            contrast_dir.name,
            groups,
        )
        results = pd.read_csv(contrast_dir / EDGE_R_RESULTS_FILE, sep="\t")
        missing = [
            column for column in REQUIRED_RESULT_COLUMNS if column not in results
        ]
        if missing:
            raise ValueError(
                f"{contrast_dir / EDGE_R_RESULTS_FILE} is missing columns: "
                + ", ".join(missing)
            )

        results["feature_id"] = results["feature_id"].astype(str)
        if (
            results["feature_id"].eq("").any()
            or results["feature_id"].duplicated().any()
        ):
            raise ValueError(
                f"{contrast_dir / EDGE_R_RESULTS_FILE} contains empty or duplicate "
                "feature_id values."
            )
        for column in ("logFC", "logCPM", "PValue", "FDR"):
            results[column] = pd.to_numeric(results[column], errors="coerce")
        finite_values = results[["logFC", "logCPM", "PValue", "FDR"]].dropna()
        if not np.isfinite(finite_values.to_numpy(dtype=float)).all():
            raise ValueError(
                f"{contrast_dir / EDGE_R_RESULTS_FILE} contains infinite statistics."
            )
        probabilities = results[["PValue", "FDR"]].dropna(how="all")
        if probabilities.lt(0).any().any() or probabilities.gt(1).any().any():
            raise ValueError(
                f"{contrast_dir / EDGE_R_RESULTS_FILE} contains probability values "
                "outside [0, 1]."
            )

        contrasts.append(
            ContrastResult(
                contrast_id=contrast_dir.name,
                target_group=target_group,
                reference_group=reference_group,
                results=results,
            )
        )
    return contrasts


def _align_feature_counts(
    feature_counts: pd.DataFrame,
    contrasts: list[ContrastResult],
    sample_names: list[str],
) -> pd.DataFrame:
    """Align CPM rows to feature IDs, including the legacy row-order format."""
    count_columns = [
        column for column in feature_counts.columns if column != "feature_id"
    ]
    missing_samples = sorted(set(sample_names) - set(count_columns))
    extra_columns = sorted(set(count_columns) - set(sample_names))
    if missing_samples or extra_columns:
        details = []
        if missing_samples:
            details.append("missing samples: " + ", ".join(missing_samples))
        if extra_columns:
            details.append("unexpected columns: " + ", ".join(extra_columns))
        raise ValueError(
            "feature_counts.tsv columns do not match report samples ("
            + "; ".join(details)
            + ")."
        )

    first_ids = contrasts[0].results["feature_id"].tolist()

    if "feature_id" in feature_counts:
        feature_counts["feature_id"] = feature_counts["feature_id"].astype(str)
        if (
            feature_counts["feature_id"].eq("").any()
            or feature_counts["feature_id"].duplicated().any()
        ):
            raise ValueError(
                "feature_counts.tsv contains empty or duplicate feature_id values."
            )
        feature_counts = feature_counts.set_index("feature_id")
        for contrast in contrasts:
            if set(feature_counts.index) != set(contrast.results["feature_id"]):
                raise ValueError(
                    "feature_counts.tsv and edgeR_results.tsv contain different "
                    "feature_id values."
                )
        feature_counts = feature_counts.loc[first_ids]
    else:
        for contrast in contrasts[1:]:
            if contrast.results["feature_id"].tolist() != first_ids:
                raise ValueError(
                    "Legacy feature-count alignment requires identical feature "
                    "ordering across every edgeR contrast."
                )
        if len(feature_counts) != len(first_ids):
            raise ValueError(
                "Legacy feature_counts.tsv has no feature_id column and its row "
                "count does not match edgeR_results.tsv."
            )
        feature_counts.index = pd.Index(first_ids, name="feature_id")

    feature_counts = feature_counts[sample_names].apply(
        pd.to_numeric,
        errors="coerce",
    )
    if feature_counts.isna().any().any():
        raise ValueError(
            "feature_counts.tsv contains missing or non-numeric CPM values."
        )
    if not np.isfinite(feature_counts.to_numpy(dtype=float)).all():
        raise ValueError("feature_counts.tsv contains infinite CPM values.")
    if (feature_counts < 0).any().any():
        raise ValueError("feature_counts.tsv contains negative CPM values.")
    return feature_counts


def load_differential_results(
    results_dir: str | Path,
    samples_df: pd.DataFrame,
) -> DifferentialResult:
    """Load edgeR outputs and validate their relationship to report samples."""
    results_dir = Path(results_dir)
    counts_path = results_dir / FEATURE_COUNTS_FILE
    if not counts_path.is_file():
        raise ValueError(f"Missing edgeR feature-count table: {counts_path}")

    metadata = _normalise_sample_metadata(samples_df)
    groups = metadata["group"].drop_duplicates().tolist()
    contrasts = _read_contrasts(results_dir, groups)
    feature_counts = _align_feature_counts(
        pd.read_csv(counts_path, sep="\t"),
        contrasts,
        metadata.index.tolist(),
    )

    return DifferentialResult(
        feature_counts=feature_counts,
        sample_metadata=metadata,
        contrasts=contrasts,
        condition_colors=condition_palette(groups),
    )


def _display_labels(results: pd.DataFrame) -> pd.Series:
    """Choose stable, human-readable labels with duplicate disambiguation."""
    labels = pd.Series(index=results.index, dtype=object)
    for index, row in results.iterrows():
        label = None
        for column in LABEL_COLUMNS:
            if column not in results:
                continue
            value = row[column]
            if pd.notna(value) and str(value).strip().lower() not in {"", "na", "nan"}:
                label = str(value).strip()
                break
        labels.at[index] = label or str(row["feature_id"])

    duplicates = labels.duplicated(keep=False)
    labels.loc[duplicates] = (
        labels.loc[duplicates]
        + " ["
        + results.loc[duplicates, "feature_id"].astype(str)
        + "]"
    )
    return labels


def classify_results(
    contrast: ContrastResult,
    condition_colors: dict[str, str],
    lfc_cutoff: float,
    padj_cutoff: float,
) -> pd.DataFrame:
    """Classify genes by both effect-size and adjusted-p-value cutoffs."""
    results = contrast.results.copy()
    significant = (
        results["logFC"].abs().ge(lfc_cutoff)
        & results["FDR"].le(padj_cutoff)
        & results[["logFC", "FDR"]].notna().all(axis=1)
    )
    positive = significant & results["logFC"].gt(0)
    negative = significant & results["logFC"].lt(0)

    results["direction"] = "Not significant"
    results.loc[positive, "direction"] = f"Up in {contrast.target_group}"
    results.loc[negative, "direction"] = f"Up in {contrast.reference_group}"
    results["plot_color"] = NONSIGNIFICANT_COLOR
    results.loc[positive, "plot_color"] = condition_colors[contrast.target_group]
    results.loc[negative, "plot_color"] = condition_colors[contrast.reference_group]
    results["significant"] = positive | negative
    results["display_label"] = _display_labels(results)

    positive_fdr = results.loc[results["FDR"].gt(0), "FDR"]
    minimum = (
        max(float(positive_fdr.min()), np.finfo(float).tiny)
        if not positive_fdr.empty
        else np.finfo(float).tiny
    )
    results["neg_log10_fdr"] = -np.log10(results["FDR"].clip(lower=minimum))
    return results


def _empty_plot(title: str, message: str, height: int = 360) -> BokehPlot:
    """Create a report-safe placeholder Bokeh plot."""
    plot = BokehPlot(
        title=title,
        x_range=(0, 1),
        y_range=(0, 1),
        height=height,
        sizing_mode="stretch_width",
        tools="",
    )
    plot._fig.xaxis.visible = False
    plot._fig.yaxis.visible = False
    plot._fig.grid.visible = False
    plot._fig.add_layout(Label(x=0.5, y=0.5, text=message, text_align="center"))
    plot.report_height = height
    return plot


def create_pca_plot(data: DifferentialResult) -> BokehPlot:
    """Create a PCA plot from log2(CPM + 1) values for all samples."""
    sample_names = data.sample_metadata.index.tolist()
    transformed = np.log2(data.feature_counts[sample_names].T + 1)
    if min(transformed.shape) < 2:
        return _empty_plot(
            "PCA of log2(CPM + 1)",
            "At least two samples and two retained features are required.",
        )

    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(transformed)
    groups = data.sample_metadata.loc[sample_names, "group"].tolist()
    source = ColumnDataSource(
        {
            "sample": sample_names,
            "group": groups,
            "PC1": coordinates[:, 0],
            "PC2": coordinates[:, 1],
            "color": [data.condition_colors[group] for group in groups],
        }
    )
    plot = BokehPlot(
        title="PCA of log2(CPM + 1)",
        x_axis_label=f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)",
        y_axis_label=f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)",
        height=420,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    hover = HoverTool(
        tooltips=[
            ("Sample", "@sample"),
            ("Condition", "@group"),
            ("PC1", "@PC1{0.000}"),
            ("PC2", "@PC2{0.000}"),
        ]
    )
    plot._fig.add_tools(hover)
    for group in data.condition_colors:
        indices = [index for index, value in enumerate(groups) if value == group]
        if not indices:
            continue
        group_source = ColumnDataSource(
            {key: [source.data[key][index] for index in indices] for key in source.data}
        )
        plot._fig.scatter(
            "PC1",
            "PC2",
            source=group_source,
            size=14,
            fill_color=data.condition_colors[group],
            fill_alpha=0.85,
            line_color="#333333",
            line_alpha=0.45,
            legend_label=group,
        )
    plot._fig.legend.title = "Condition"
    plot._fig.legend.location = "top_right"
    return plot


def _add_classified_points(
    plot: BokehPlot,
    classified: pd.DataFrame,
    x_column: str,
    y_column: str,
    contrast: ContrastResult,
    condition_colors: dict[str, str],
) -> None:
    """Add gene points in deterministic legend order."""
    class_colors = (
        ("Not significant", NONSIGNIFICANT_COLOR),
        (
            f"Up in {contrast.reference_group}",
            condition_colors[contrast.reference_group],
        ),
        (f"Up in {contrast.target_group}", condition_colors[contrast.target_group]),
    )
    for direction, color in class_colors:
        subset = classified.loc[classified["direction"].eq(direction)]
        if subset.empty:
            continue
        plot._fig.scatter(
            x_column,
            y_column,
            source=ColumnDataSource(subset),
            size=6,
            fill_color=color,
            fill_alpha=0.7 if direction != "Not significant" else 0.35,
            line_color=None,
            legend_label=direction,
        )
    plot._fig.legend.location = "top_right"


def create_ma_plot(
    contrast: ContrastResult,
    condition_colors: dict[str, str],
    lfc_cutoff: float,
    padj_cutoff: float,
) -> BokehPlot:
    """Create the edgeR logFC-versus-logCPM (MA) plot."""
    classified = classify_results(
        contrast,
        condition_colors,
        lfc_cutoff,
        padj_cutoff,
    )
    plot = BokehPlot(
        title=f"logFC vs logCPM — {contrast.label}",
        x_axis_label="edgeR logCPM",
        y_axis_label="edgeR log2 fold change",
        height=420,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    plot._fig.add_tools(
        HoverTool(
            tooltips=[
                ("Feature", "@display_label"),
                ("Feature ID", "@feature_id"),
                ("logFC", "@logFC{0.000}"),
                ("logCPM", "@logCPM{0.000}"),
                ("FDR", "@FDR{0.000e}"),
            ]
        )
    )
    _add_classified_points(
        plot,
        classified.dropna(subset=["logCPM", "logFC"]),
        "logCPM",
        "logFC",
        contrast,
        condition_colors,
    )
    for location in (-lfc_cutoff, lfc_cutoff):
        plot._fig.add_layout(
            Span(
                location=location,
                dimension="width",
                line_color="#555555",
                line_dash="dashed",
                line_width=1,
            )
        )
    plot._fig.add_layout(
        Span(location=0, dimension="width", line_color="#888888", line_width=1)
    )
    return plot


def create_volcano_plot(
    contrast: ContrastResult,
    condition_colors: dict[str, str],
    lfc_cutoff: float,
    padj_cutoff: float,
) -> BokehPlot:
    """Create an edgeR volcano plot colored by contrast direction."""
    classified = classify_results(
        contrast,
        condition_colors,
        lfc_cutoff,
        padj_cutoff,
    )
    plot = BokehPlot(
        title=f"Volcano plot — {contrast.label}",
        x_axis_label="edgeR log2 fold change",
        y_axis_label="-log10(edgeR FDR)",
        height=420,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,save,reset",
    )
    plot._fig.add_tools(
        HoverTool(
            tooltips=[
                ("Feature", "@display_label"),
                ("Feature ID", "@feature_id"),
                ("logFC", "@logFC{0.000}"),
                ("logCPM", "@logCPM{0.000}"),
                ("FDR", "@FDR{0.000e}"),
            ]
        )
    )
    _add_classified_points(
        plot,
        classified.dropna(subset=["logFC", "neg_log10_fdr"]),
        "logFC",
        "neg_log10_fdr",
        contrast,
        condition_colors,
    )
    for location in (-lfc_cutoff, lfc_cutoff):
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
        Span(
            location=-np.log10(padj_cutoff),
            dimension="width",
            line_color="#555555",
            line_dash="dashed",
            line_width=1,
        )
    )
    return plot


def select_heatmap_genes(
    contrast: ContrastResult,
    condition_colors: dict[str, str],
    lfc_cutoff: float,
    padj_cutoff: float,
    top_n: int = 20,
) -> pd.DataFrame:
    """Select strongest significant positive and negative features."""
    classified = classify_results(
        contrast,
        condition_colors,
        lfc_cutoff,
        padj_cutoff,
    )
    up = (
        classified.loc[classified["significant"] & classified["logFC"].gt(0)]
        .sort_values(["logFC", "FDR"], ascending=[False, True])
        .head(top_n)
    )
    down = (
        classified.loc[classified["significant"] & classified["logFC"].lt(0)]
        .sort_values(["logFC", "FDR"], ascending=[True, True])
        .head(top_n)
    )
    return pd.concat([down, up], ignore_index=True)


def cluster_heatmap_rows(matrix: pd.DataFrame) -> list[str]:
    """Order rows by average-linkage clustering without drawing a tree."""
    if len(matrix.index) < 2:
        return matrix.index.tolist()
    hierarchy = linkage(
        matrix.to_numpy(dtype=float),
        method="average",
        metric="euclidean",
        optimal_ordering=True,
    )
    return matrix.index.take(leaves_list(hierarchy)).tolist()


def create_heatmap_plot(
    data: DifferentialResult,
    contrast: ContrastResult,
    lfc_cutoff: float,
    padj_cutoff: float,
    top_n: int = 20,
) -> BokehPlot:
    """Create a row-z-scored log2(CPM + 1) heatmap for one contrast."""
    selected = select_heatmap_genes(
        contrast,
        data.condition_colors,
        lfc_cutoff,
        padj_cutoff,
        top_n,
    )
    if selected.empty:
        return _empty_plot(
            f"Top differential genes — {contrast.label}",
            (
                "No genes pass both "
                f"|log2FC| ≥ {lfc_cutoff:g} and FDR ≤ {padj_cutoff:g}."
            ),
            height=440,
        )

    metadata = data.sample_metadata
    sample_order = []
    for group in (contrast.reference_group, contrast.target_group):
        sample_order.extend(
            sorted(metadata.index[metadata["group"].eq(group)].tolist())
        )

    feature_ids = selected["feature_id"].tolist()
    labels = selected.set_index("feature_id")["display_label"].to_dict()
    log_counts = np.log2(data.feature_counts.loc[feature_ids, sample_order] + 1)
    row_means = log_counts.mean(axis=1)
    row_std = log_counts.std(axis=1).replace(0, np.nan)
    z_scores = log_counts.sub(row_means, axis=0).div(row_std, axis=0).fillna(0)
    feature_ids = cluster_heatmap_rows(z_scores)
    z_scores = z_scores.loc[feature_ids]
    melted = (
        z_scores.rename_axis("feature_id")
        .reset_index()
        .melt(id_vars="feature_id", var_name="sample", value_name="z_score")
    )
    melted["display_label"] = melted["feature_id"].map(labels)

    z_limit = max(float(melted["z_score"].abs().max()), 1.0)
    mapper = LinearColorMapper(
        palette=RdBu11,
        low=-z_limit,
        high=z_limit,
    )
    heatmap = BokehPlot(
        title=f"Top differential genes — {contrast.label}",
        x_range=sample_order,
        y_range=list(reversed([labels[feature_id] for feature_id in feature_ids])),
        x_axis_location="above",
        x_axis_label="Sample",
        y_axis_label="Gene",
        height=max(440, 22 * len(feature_ids) + 120),
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    heatmap._fig.grid.grid_line_color = None
    heatmap._fig.xaxis.major_label_orientation = np.pi / 4
    heatmap._fig.add_tools(
        HoverTool(
            tooltips=[
                ("Gene", "@display_label"),
                ("Feature ID", "@feature_id"),
                ("Sample", "@sample"),
                ("Row z-score", "@z_score{0.000}"),
            ]
        )
    )
    heatmap._fig.rect(
        x="sample",
        y="display_label",
        width=0.98,
        height=0.98,
        source=ColumnDataSource(melted),
        fill_color={"field": "z_score", "transform": mapper},
        line_color=None,
    )
    heatmap._fig.add_layout(
        ColorBar(color_mapper=mapper, title="Row z-score"),
        "right",
    )

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
    for group in (contrast.reference_group, contrast.target_group):
        group_samples = [
            sample for sample in sample_order if metadata.at[sample, "group"] == group
        ]
        annotation._fig.rect(
            x=group_samples,
            y=[0.5] * len(group_samples),
            width=0.98,
            height=0.7,
            fill_color=data.condition_colors[group],
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


def add_differential_analysis(
    data: DifferentialResult,
    lfc_cutoff: float,
    padj_cutoff: float,
) -> None:
    """Add overview, contrast, and plot-type tabs to the report section."""
    tabs = Tabs()
    with tabs.add_tab("Overview"):
        EZChart(create_pca_plot(data), "epi2melabs", height="460px")

    with tabs.add_tab("Contrasts"):
        contrast_tabs = Tabs()
        for contrast in data.contrasts:
            with contrast_tabs.add_tab(contrast.label):
                plot_tabs = Tabs()
                with plot_tabs.add_tab("MA plot"):
                    EZChart(
                        create_ma_plot(
                            contrast,
                            data.condition_colors,
                            lfc_cutoff,
                            padj_cutoff,
                        ),
                        "epi2melabs",
                        height="460px",
                    )
                with plot_tabs.add_tab("Volcano"):
                    EZChart(
                        create_volcano_plot(
                            contrast,
                            data.condition_colors,
                            lfc_cutoff,
                            padj_cutoff,
                        ),
                        "epi2melabs",
                        height="460px",
                    )
                with plot_tabs.add_tab("Heatmap"):
                    heatmap = create_heatmap_plot(
                        data,
                        contrast,
                        lfc_cutoff,
                        padj_cutoff,
                    )
                    EZChart(
                        heatmap,
                        "epi2melabs",
                        height=f"{getattr(heatmap, 'report_height', 720)}px",
                    )
