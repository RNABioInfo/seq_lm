"""Tests for edgeR parsing and differential report plots."""

from pathlib import Path

import pytest

pytest.importorskip("ezcharts")
pd = pytest.importorskip("pandas")

from bokeh.models import (  # noqa: E402
    ColumnDataSource,
    Label,
    LinearColorMapper,
    Span,
)
from bokeh.palettes import RdBu11  # noqa: E402

from workflow_glue.qc_report_types import differential_plots as dp  # noqa: E402


def sample_metadata(include_second_target=False):
    """Return report-style sample metadata."""
    rows = [
        {"Name": "control_1", "Group": "control"},
        {"Name": "control_2", "Group": "control"},
        {"Name": "treated_1", "Group": "treated"},
        {"Name": "treated_2", "Group": "treated"},
    ]
    if include_second_target:
        rows.extend(
            [
                {"Name": "rescue_1", "Group": "rescue"},
                {"Name": "rescue_2", "Group": "rescue"},
            ]
        )
    return pd.DataFrame(rows)


def write_edge_r_tree(
    path: Path,
    *,
    legacy_counts: bool = False,
    include_second_target: bool = False,
    second_order_mismatch: bool = False,
):
    """Write compact edgeR outputs for report tests."""
    metadata = sample_metadata(include_second_target)
    samples = metadata["Name"].tolist()
    feature_ids = ["feature_a", "feature_b", "feature_c", "feature_d"]
    counts = pd.DataFrame(
        {
            sample: [1 + index, 20 - index, 5 + index, 8]
            for index, sample in enumerate(samples)
        },
        index=feature_ids,
    )
    if legacy_counts:
        counts.to_csv(path / "feature_counts.tsv", sep="\t", index=False)
    else:
        counts.rename_axis("feature_id").to_csv(
            path / "feature_counts.tsv",
            sep="\t",
        )

    pd.DataFrame(
        {
            "feature_id": feature_ids,
            "average_log_cpm": [2.0, 3.0, 4.0, 5.0],
            "tagwise_dispersion": [0.16, 0.09, 0.04, 0.01],
            "tagwise_bcv": [0.4, 0.3, 0.2, 0.1],
            "trended_dispersion": [0.1225, 0.09, 0.0625, 0.04],
            "trended_bcv": [0.35, 0.3, 0.25, 0.2],
            "common_dispersion": [0.09] * 4,
            "common_bcv": [0.3] * 4,
        }
    ).to_csv(path / "edgeR_bcv_data.tsv", sep="\t", index=False)

    results = pd.DataFrame(
        {
            "feature_id": feature_ids,
            "logFC": [3.0, -2.5, 0.2, 1.5],
            "logCPM": [4.0, 5.0, 3.0, 2.0],
            "PValue": [0.0001, 0.0002, 0.8, 0.01],
            "FDR": [0.001, 0.002, 0.9, 0.2],
            "gene": ["gene_a", "gene_b", "gene_c", "gene_a"],
            "gene_id": ["id_a", "id_b", "id_c", "id_d"],
        }
    )
    contrast = path / "group_treated_vs_control"
    contrast.mkdir()
    results.to_csv(contrast / "edgeR_results.tsv", sep="\t", index=False)

    if include_second_target:
        rescue = path / "group_rescue_vs_control"
        rescue.mkdir()
        rescue_results = results.copy()
        if second_order_mismatch:
            rescue_results = rescue_results.iloc[[1, 0, 2, 3]]
        rescue_results.to_csv(
            rescue / "edgeR_results.tsv",
            sep="\t",
            index=False,
        )
    return metadata


def test_condition_palette_is_deterministic_and_control_first():
    """Condition colors do not depend on sample arrival order."""
    first = dp.condition_palette(["treated", "control", "rescue"])
    second = dp.condition_palette(["rescue", "treated", "control", "treated"])

    assert first == second
    assert list(first) == ["control", "rescue", "treated"]
    assert len(set(first.values())) == 3


def test_load_differential_results_aligns_by_feature_id(tmp_path):
    """Current edgeR counts align explicitly by feature ID and sample name."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)

    assert data.feature_counts.index.tolist() == [
        "feature_a",
        "feature_b",
        "feature_c",
        "feature_d",
    ]
    assert data.feature_counts.columns.tolist() == metadata["Name"].tolist()
    assert data.bcv_data["feature_id"].tolist() == data.feature_counts.index.tolist()
    assert data.contrasts[0].label == "treated vs control"


def test_load_differential_results_supports_valid_legacy_counts(tmp_path):
    """Legacy counts without feature_id use validated edgeR row ordering."""
    metadata = write_edge_r_tree(
        tmp_path,
        legacy_counts=True,
        include_second_target=True,
    )
    data = dp.load_differential_results(tmp_path, metadata)

    assert data.feature_counts.index.name == "feature_id"
    assert (
        data.feature_counts.index.tolist()
        == data.contrasts[0].results["feature_id"].tolist()
    )
    assert len(data.contrasts) == 2


def test_legacy_counts_reject_different_contrast_feature_order(tmp_path):
    """A legacy row-order fallback cannot silently align inconsistent contrasts."""
    metadata = write_edge_r_tree(
        tmp_path,
        legacy_counts=True,
        include_second_target=True,
        second_order_mismatch=True,
    )

    with pytest.raises(ValueError, match="identical feature ordering"):
        dp.load_differential_results(tmp_path, metadata)


def test_load_differential_results_rejects_sample_mismatch(tmp_path):
    """CPM columns must exactly match report sample names."""
    metadata = write_edge_r_tree(tmp_path)
    counts_path = tmp_path / "feature_counts.tsv"
    counts = pd.read_csv(counts_path, sep="\t").drop(columns="treated_2")
    counts.to_csv(counts_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="missing samples: treated_2"):
        dp.load_differential_results(tmp_path, metadata)


def test_classification_uses_both_cutoffs_and_condition_colors(tmp_path):
    """Only genes passing LFC and FDR are colored by contrast direction."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)
    contrast = data.contrasts[0]

    classified = dp.classify_results(
        contrast,
        data.condition_colors,
        lfc_cutoff=1,
        padj_cutoff=0.05,
    ).set_index("feature_id")

    assert classified.at["feature_a", "direction"] == "Up in treated"
    assert classified.at["feature_b", "direction"] == "Up in control"
    assert classified.at["feature_c", "direction"] == "Not significant"
    assert classified.at["feature_d", "direction"] == "Not significant"
    assert classified.at["feature_a", "plot_color"] == data.condition_colors["treated"]
    assert classified.at["feature_b", "plot_color"] == data.condition_colors["control"]
    assert classified.at["feature_c", "plot_color"] == dp.NONSIGNIFICANT_COLOR


def test_heatmap_selection_ranks_effect_size_and_disambiguates_labels(tmp_path):
    """Heatmap genes are significant, direction-balanced, and uniquely labeled."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)

    selected = dp.select_heatmap_genes(
        data.contrasts[0],
        data.condition_colors,
        lfc_cutoff=1,
        padj_cutoff=0.05,
        top_n=20,
    )

    assert selected["feature_id"].tolist() == ["feature_b", "feature_a"]
    assert selected["display_label"].is_unique


def test_pca_uses_log_cpm_and_reports_explained_variance(tmp_path):
    """PCA axes expose explained variance and retain condition metadata."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)
    plot = dp.create_pca_plot(data)

    assert plot._fig.title.text == "PCA of log2(CPM + 1)"
    assert plot._fig.xaxis.axis_label.startswith("PC1 (")
    assert plot._fig.yaxis.axis_label.startswith("PC2 (")
    groups = {
        group
        for source in plot._fig.select({"type": ColumnDataSource})
        for group in source.data.get("group", [])
    }
    assert groups == {"control", "treated"}


def test_bcv_plot_uses_exported_tagwise_trended_and_common_values(tmp_path):
    """BCV overview plot exposes all three edgeR dispersion summaries."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)
    plot = dp.create_bcv_plot(data)

    assert plot._fig.title.text == "edgeR biological coefficient of variation (BCV)"
    assert plot._fig.xaxis.axis_label == "Average log CPM"
    assert plot._fig.yaxis.axis_label == "Biological coefficient of variation"
    assert [item.label["value"] for item in plot._fig.legend[0].items] == [
        "Tagwise",
        "Trended",
        "Common",
    ]
    sources = list(plot._fig.select({"type": ColumnDataSource}))
    assert any(
        list(source.data.get("tagwise_bcv", [])) == [0.4, 0.3, 0.2, 0.1]
        for source in sources
    )


def test_overview_adds_bcv_plot_with_pca(tmp_path, monkeypatch):
    """The differential Overview tab emits both global diagnostic plots."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)
    rendered_titles = []

    def record_chart(plot, *_args, **_kwargs):
        title = getattr(plot._fig, "title", None)
        if title is not None:
            rendered_titles.append(title.text)

    monkeypatch.setattr(dp, "EZChart", record_chart)
    dp.add_differential_analysis(data, lfc_cutoff=1, padj_cutoff=0.05)

    assert rendered_titles[:2] == [
        "PCA of log2(CPM + 1)",
        "edgeR biological coefficient of variation (BCV)",
    ]


def test_ma_and_volcano_include_cutoff_spans(tmp_path):
    """Both contrast plots show their applicable decision thresholds."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)
    contrast = data.contrasts[0]

    ma = dp.create_ma_plot(contrast, data.condition_colors, 1, 0.05)
    volcano = dp.create_volcano_plot(contrast, data.condition_colors, 1, 0.05)

    assert len(list(ma._fig.select({"type": Span}))) == 3
    assert len(list(volcano._fig.select({"type": Span}))) == 3


def test_heatmap_contains_only_compared_conditions(tmp_path):
    """Contrast heatmaps exclude samples from unrelated conditions."""
    metadata = write_edge_r_tree(tmp_path, include_second_target=True)
    data = dp.load_differential_results(tmp_path, metadata)
    contrast = next(item for item in data.contrasts if item.target_group == "treated")

    heatmap = dp.create_heatmap_plot(data, contrast, 1, 0.05)
    samples = {
        sample
        for source in heatmap._fig.select({"type": ColumnDataSource})
        for sample in source.data.get("sample", [])
    }

    assert samples == {"control_1", "control_2", "treated_1", "treated_2"}
    assert not any(sample.startswith("rescue") for sample in samples)

    mapper = next(iter(heatmap._fig.select({"type": LinearColorMapper})))
    assert mapper.palette[0] == RdBu11[0]
    assert mapper.palette[-1] == RdBu11[-1]


def test_heatmap_row_clustering_groups_similar_profiles():
    """Clustering places similarly high and low row profiles together."""
    matrix = pd.DataFrame(
        {
            "sample_1": [2.0, 1.8, -2.0, -1.8],
            "sample_2": [1.0, 0.9, -1.0, -0.9],
            "sample_3": [-2.0, -1.8, 2.0, 1.8],
        },
        index=["high_a", "high_b", "low_a", "low_b"],
    )

    order = dp.cluster_heatmap_rows(matrix)

    positions = {feature: index for index, feature in enumerate(order)}
    assert abs(positions["high_a"] - positions["high_b"]) == 1
    assert abs(positions["low_a"] - positions["low_b"]) == 1


def test_heatmap_placeholder_when_no_gene_is_significant(tmp_path):
    """A valid empty result renders an explanatory placeholder."""
    metadata = write_edge_r_tree(tmp_path)
    data = dp.load_differential_results(tmp_path, metadata)

    heatmap = dp.create_heatmap_plot(data, data.contrasts[0], 10, 0.05)
    labels = list(heatmap._fig.select({"type": Label}))

    assert "No genes pass both" in labels[0].text
