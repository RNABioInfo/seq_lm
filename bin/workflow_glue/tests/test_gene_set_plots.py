"""Tests for fry parsing and gene-set enrichment report plots."""

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ezcharts")
pd = pytest.importorskip("pandas")

from bokeh.models import ColumnDataSource, Div, Label, Select, Span  # noqa: E402

from workflow_glue.qc_report_types import differential_plots as dp  # noqa: E402
from workflow_glue.qc_report_types import gene_set_plots as gp  # noqa: E402
from workflow_glue.qc_report_types import gsva_plots as gsp  # noqa: E402


def write_gene_set_tree(path: Path):
    """Write a compact, report-ready edgeR and fry output tree."""
    metadata = pd.DataFrame(
        [
            {"Name": "control_1", "Group": "control"},
            {"Name": "treated_1", "Group": "treated"},
        ]
    )
    feature_ids = ["feature_a", "feature_b", "feature_c", "feature_d"]
    pd.DataFrame(
        {
            "feature_id": feature_ids,
            "control_1": [2, 40, 8, 7],
            "treated_1": [30, 4, 9, 20],
        }
    ).to_csv(path / "feature_counts.tsv", sep="\t", index=False)

    contrast_dir = path / "group_treated_vs_control"
    contrast_dir.mkdir()
    pd.DataFrame(
        {
            "feature_id": feature_ids,
            "logFC": [2.5, -2.2, -0.2, 1.1],
            "logCPM": [5.0, 5.5, 4.0, 3.0],
            "PValue": [0.0001, 0.0002, 0.5, 0.02],
            "FDR": [0.001, 0.002, 0.8, 0.04],
            "gene": ["gene_a", "gene_b", "gene_c", "gene_d"],
        }
    ).to_csv(contrast_dir / "edgeR_results.tsv", sep="\t", index=False)

    fry = pd.DataFrame(
        {
            "gene_set": ["set_up", "set_down", "set_mixed"],
            "NGenes": [2, 2, 2],
            "Direction": ["Up", "Down", "Down"],
            "PValue": [0.001, 0.002, 0.7],
            "FDR": [0.01, 0.02, 0.8],
            "PValue.Mixed": [0.005, 0.01, 0.0001],
            "FDR.Mixed": [0.02, 0.03, 0.001],
            "description": ["Shared label", "Shared label", "Mixed response"],
            "gmt_members": [3, 2, 2],
            "matched_gmt_members": [3, 2, 2],
            "count_matrix_members": [3, 2, 2],
            "tested_members": [2, 2, 2],
            "tested_gmt_members": [2, 2, 2],
            "count_matrix_coverage": [1.0, 1.0, 1.0],
            "tested_coverage": [2 / 3, 1.0, 1.0],
        }
    )
    fry.to_csv(contrast_dir / "fry_results.tsv", sep="\t", index=False)

    pd.DataFrame(
        {
            "gene_set": [
                "set_up",
                "set_up",
                "set_down",
                "set_down",
                "set_mixed",
                "set_mixed",
            ],
            "description": [
                "Shared label",
                "Shared label",
                "Shared label",
                "Shared label",
                "Mixed response",
                "Mixed response",
            ],
            "n_genes": [2, 2, 2, 2, 2, 2],
            "sample": [
                "control_1",
                "treated_1",
                "control_1",
                "treated_1",
                "control_1",
                "treated_1",
            ],
            "group": [
                "control",
                "treated",
                "control",
                "treated",
                "control",
                "treated",
            ],
            "score": [-0.6, 0.7, 0.5, -0.4, 0.0, 0.1],
        }
    ).to_csv(path / "gsva_scores_long.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "gene_set": ["set_up", "set_down", "set_mixed", "set_small"],
            "description": [
                "Shared label",
                "Shared label",
                "Mixed response",
                "Small set",
            ],
            "resolved_members": [3, 2, 2, 1],
            "retained_members": [2, 2, 2, 1],
            "variable_members": [2, 2, 2, 1],
            "scored_members": [2, 2, 2, 0],
            "status": ["scored", "scored", "scored", "below_min_size"],
        }
    ).to_csv(path / "gsva_gene_set_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "gene_set": ["set_up", "set_down", "set_mixed"],
            "description": ["Shared label", "Shared label", "Mixed response"],
            "n_genes": [2, 2, 2],
            "target_group": ["treated"] * 3,
            "control_group": ["control"] * 3,
            "effect_size": [1.3, -0.9, 0.1],
            "average_score": [0.05, 0.05, 0.05],
            "t_statistic": [5.0, -4.0, 0.2],
            "p_value": [0.001, 0.003, 0.8],
            "adjusted_p_value": [0.003, 0.0045, 0.8],
            "log_odds": [2.0, 1.5, -4.0],
        }
    ).to_csv(contrast_dir / "gsva_limma_results.tsv", sep="\t", index=False)

    pd.DataFrame(
        {
            "gene_set": [
                "set_up",
                "set_up",
                "set_up",
                "set_up",
                "set_down",
                "set_down",
                "set_mixed",
                "set_mixed",
            ],
            "feature_id": [
                "feature_a",
                "feature_d",
                "feature_unretained",
                "feature_a",
                "feature_b",
                "feature_c",
                "feature_a",
                "feature_c",
            ],
        }
    ).to_csv(path / "gene_set_resolution.tsv", sep="\t", index=False)
    return metadata


def load_tree(path: Path):
    """Load the compact differential and fry fixture."""
    metadata = write_gene_set_tree(path)
    differential = dp.load_differential_results(path, metadata)
    analyses = gp.load_gene_set_results(path, differential)
    return differential, analyses[0]


def load_gsva_tree(path: Path):
    """Load the compact differential and GSVA fixture."""
    metadata = write_gene_set_tree(path)
    differential = dp.load_differential_results(path, metadata)
    return differential, gsp.load_gsva_results(path, differential)


def test_load_gsva_results_aligns_scores_coverage_and_limma(tmp_path):
    """GSVA report inputs retain score and contrast ordering."""
    differential, gsva = load_gsva_tree(tmp_path)

    assert gsva.gene_set_order == ["set_up", "set_down", "set_mixed"]
    assert gsva.scores_long["sample"].tolist() == [
        "control_1",
        "treated_1",
    ] * 3
    assert gsva.coverage["status"].tolist() == [
        "scored",
        "scored",
        "scored",
        "below_min_size",
    ]
    assert gsva.contrasts[0].label == "treated vs control"
    assert gsva.contrasts[0].results["effect_size"].tolist() == [1.3, -0.9, 0.1]
    assert set(differential.condition_colors) == {"control", "treated"}


def test_gsva_plots_use_score_differences_and_adjusted_p_values(tmp_path):
    """Limma visuals use GSVA score units rather than fold-change labels."""
    differential, gsva = load_gsva_tree(tmp_path)
    analysis = gsva.contrasts[0]
    prepared = gsp.prepare_limma_volcano(
        analysis,
        differential.condition_colors,
        0.05,
    ).set_index("gene_set")

    assert prepared.at["set_up", "significant"]
    assert prepared.at["set_down", "significant"]
    assert not prepared.at["set_mixed", "significant"]
    assert prepared.at["set_up", "plot_color"] == differential.condition_colors["treated"]
    assert prepared.at["set_down", "plot_color"] == differential.condition_colors["control"]

    volcano = gsp.create_limma_volcano(
        analysis,
        differential.condition_colors,
        0.05,
    )
    assert "GSVA score difference" in volcano._fig.xaxis[0].axis_label
    assert "adjusted p-value" in volcano._fig.yaxis[0].axis_label


def test_gsva_raw_and_differential_heatmaps_render(tmp_path):
    """Raw-score, distribution, and differential views contain data marks."""
    differential, gsva = load_gsva_tree(tmp_path)

    raw_heatmap = gsp.create_score_heatmap(gsva, differential.condition_colors)
    distribution = gsp.create_score_distribution(
        gsva,
        "set_up",
        differential.condition_colors,
    )
    differential_heatmap = gsp.create_limma_heatmap(
        gsva,
        gsva.contrasts[0],
        differential.condition_colors,
        0.05,
    )
    summary = gsp.create_multi_contrast_dot_plot(gsva, 0.05)

    assert raw_heatmap.report_height > 0
    assert len(distribution._fig.renderers) > 0
    assert differential_heatmap.report_height > 0
    assert summary._fig.y_range.factors == [
        "Mixed response",
        "Shared label [set_down]",
        "Shared label [set_up]",
    ]


@pytest.mark.parametrize(
    ("relative_path", "column", "value", "message"),
    [
        ("gsva_scores_long.tsv", "score", np.inf, "nonfinite score"),
        (
            "group_treated_vs_control/gsva_limma_results.tsv",
            "adjusted_p_value",
            1.2,
            "invalid adjusted_p_value",
        ),
    ],
)
def test_load_gsva_results_rejects_invalid_numeric_values(
    tmp_path,
    relative_path,
    column,
    value,
    message,
):
    """Malformed score and limma values fail before report rendering."""
    metadata = write_gene_set_tree(tmp_path)
    path = tmp_path / relative_path
    table = pd.read_csv(path, sep="\t")
    table.loc[0, column] = value
    table.to_csv(path, sep="\t", index=False)
    differential = dp.load_differential_results(tmp_path, metadata)

    with pytest.raises(ValueError, match=message):
        gsp.load_gsva_results(tmp_path, differential)


def test_load_gene_set_results_resolves_retained_unique_members(tmp_path):
    """Resolution rows are deduplicated and intersected with retained features."""
    _, analysis = load_tree(tmp_path)

    assert analysis.members["set_up"] == ("feature_a", "feature_d")
    assert analysis.members["set_down"] == ("feature_b", "feature_c")
    assert analysis.fry_results["display_label"].is_unique
    assert analysis.fry_results.loc[0, "display_label"].endswith("[set_up]")


def test_load_gene_set_results_rejects_membership_count_mismatch(tmp_path):
    """The fry NGenes field must agree with retained resolved feature IDs."""
    metadata = write_gene_set_tree(tmp_path)
    fry_path = tmp_path / "group_treated_vs_control/fry_results.tsv"
    fry = pd.read_csv(fry_path, sep="\t")
    fry.loc[fry["gene_set"].eq("set_up"), "NGenes"] = 3
    fry.to_csv(fry_path, sep="\t", index=False)
    differential = dp.load_differential_results(tmp_path, metadata)

    with pytest.raises(ValueError, match="reports NGenes=3, but 2"):
        gp.load_gene_set_results(tmp_path, differential)


def test_load_gene_set_results_requires_resolution_table(tmp_path):
    """A missing membership mapping fails with the exact required path."""
    metadata = write_gene_set_tree(tmp_path)
    (tmp_path / "gene_set_resolution.tsv").unlink()
    differential = dp.load_differential_results(tmp_path, metadata)

    with pytest.raises(ValueError, match="Missing edgeR gene-set resolution"):
        gp.load_gene_set_results(tmp_path, differential)


def test_load_gene_set_results_rejects_duplicate_gene_sets(tmp_path):
    """Each contrast must contain one fry row per gene-set identifier."""
    metadata = write_gene_set_tree(tmp_path)
    fry_path = tmp_path / "group_treated_vs_control/fry_results.tsv"
    fry = pd.read_csv(fry_path, sep="\t")
    pd.concat([fry, fry.iloc[[0]]], ignore_index=True).to_csv(
        fry_path,
        sep="\t",
        index=False,
    )
    differential = dp.load_differential_results(tmp_path, metadata)

    with pytest.raises(ValueError, match="duplicate gene_set values: set_up"):
        gp.load_gene_set_results(tmp_path, differential)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("Direction", "Sideways", "invalid Direction"),
        ("FDR", 1.2, "invalid FDR"),
        ("NGenes", 1, "invalid NGenes"),
    ],
)
def test_load_gene_set_results_rejects_malformed_fry(
    tmp_path,
    column,
    value,
    message,
):
    """Invalid categorical, probability, and membership values fail clearly."""
    metadata = write_gene_set_tree(tmp_path)
    fry_path = tmp_path / "group_treated_vs_control/fry_results.tsv"
    fry = pd.read_csv(fry_path, sep="\t")
    fry.loc[0, column] = value
    fry.to_csv(fry_path, sep="\t", index=False)
    differential = dp.load_differential_results(tmp_path, metadata)

    with pytest.raises(ValueError, match=message):
        gp.load_gene_set_results(tmp_path, differential)


def test_signed_significance_uses_directional_fdr_and_condition_colors(tmp_path):
    """Only directional significance controls sign and condition coloring."""
    differential, analysis = load_tree(tmp_path)
    prepared = gp.prepare_signed_gene_sets(
        analysis,
        differential.condition_colors,
        padj_cutoff=0.05,
    ).set_index("gene_set")

    assert prepared.at["set_up", "signed_significance"] == pytest.approx(2.0)
    assert prepared.at["set_down", "signed_significance"] == pytest.approx(
        np.log10(0.02)
    )
    assert prepared.at["set_mixed", "signed_significance"] < 0
    assert (
        prepared.at["set_up", "plot_color"]
        == differential.condition_colors["treated"]
    )
    assert (
        prepared.at["set_down", "plot_color"]
        == differential.condition_colors["control"]
    )
    assert prepared.at["set_mixed", "plot_color"] == gp.NONSIGNIFICANT_COLOR

    plot = gp.create_signed_significance_plot(
        analysis,
        differential.condition_colors,
        0.05,
    )
    assert len(list(plot._fig.select({"type": Span}))) == 3
    assert plot._fig.y_range.factors == ["set_mixed", "set_down", "set_up"]


def test_zero_fdr_is_finite_for_plotting(tmp_path):
    """A reported zero is clamped for plotting without changing the raw FDR."""
    differential, analysis = load_tree(tmp_path)
    analysis.fry_results.loc[
        analysis.fry_results["gene_set"].eq("set_up"),
        "FDR",
    ] = 0
    prepared = gp.prepare_signed_gene_sets(
        analysis,
        differential.condition_colors,
        0.05,
    )

    assert np.isfinite(prepared["signed_significance"]).all()
    assert analysis.fry_results.loc[0, "FDR"] == 0


def test_signed_plot_limits_summary_to_top_30(tmp_path):
    """Large collections state that the signed summary is truncated."""
    differential, analysis = load_tree(tmp_path)
    template = analysis.fry_results.iloc[0]
    rows = []
    members = {}
    for index in range(31):
        row = template.copy()
        row["gene_set"] = f"set_{index:02d}"
        row["display_label"] = f"Set {index:02d}"
        row["FDR"] = (index + 1) / 100
        rows.append(row)
        members[row["gene_set"]] = ("feature_a", "feature_d")
    analysis.fry_results = pd.DataFrame(rows)
    analysis.members = members

    prepared = gp.prepare_signed_gene_sets(
        analysis,
        differential.condition_colors,
        0.05,
    )
    plot = gp.create_signed_significance_plot(
        analysis,
        differential.condition_colors,
        0.05,
    )

    assert len(prepared) == 30
    assert "top 30 of 31" in plot._fig.title.text


def test_empty_fry_table_renders_informative_placeholders(tmp_path):
    """A valid empty fry table does not prevent report generation."""
    metadata = write_gene_set_tree(tmp_path)
    fry_path = tmp_path / "group_treated_vs_control/fry_results.tsv"
    fry = pd.read_csv(fry_path, sep="\t").iloc[0:0]
    fry.to_csv(fry_path, sep="\t", index=False)
    differential = dp.load_differential_results(tmp_path, metadata)
    analysis = gp.load_gene_set_results(tmp_path, differential)[0]

    signed = gp.create_signed_significance_plot(
        analysis,
        differential.condition_colors,
        0.05,
    )
    barcode = gp.create_barcode_plot(analysis, differential.condition_colors)
    messages = [
        label.text
        for plot in (signed, barcode)
        for label in plot._fig.select({"type": Label})
    ]

    assert messages == [
        "No gene sets were tested for this contrast.",
        "No gene sets were tested for this contrast.",
    ]


def test_tricube_moving_average_matches_limma_fixture():
    """The Python worm smoother matches limma::tricubeMovingAverage."""
    observed = gp.tricube_moving_average(
        [0, 1, 0, 0, 1, 0, 0, 0, 1, 0],
        span=0.45,
    )
    expected = np.array(
        [
            0.426199123201002117,
            0.374337240186720766,
            0.317927849979762756,
            0.317927849979762756,
            0.364144300040474489,
            0.290698552298447488,
            0.054458595362630535,
            0.290698552298447488,
            0.374337240186720766,
            0.426199123201002117,
        ]
    )

    np.testing.assert_allclose(observed, expected, rtol=1e-14, atol=1e-14)


def test_barcode_renders_one_set_without_a_bokeh_selector(tmp_path):
    """The report-level dropdown replaces the visually inconsistent widget."""
    differential, analysis = load_tree(tmp_path)
    barcode = gp.create_barcode_plot(analysis, differential.condition_colors)

    selectors = list(barcode._fig.select({"type": Select}))
    assert selectors == []

    sources = list(barcode._fig.select({"type": ColumnDataSource}))
    member_sources = [
        source for source in sources if "feature_id" in source.data
    ]
    assert len(member_sources) == 1
    assert member_sources[0].data["feature_id"] == ["feature_d", "feature_a"]
    assert member_sources[0].data["plot_color"] == [
        differential.condition_colors["treated"],
        differential.condition_colors["treated"],
    ]
    details = [
        div.text
        for div in barcode._fig.select({"type": Div})
        if "Gene-set details" in div.text
    ]
    assert len(details) == 1
    assert "<details" in details[0]
    assert "Shared label [set_up]" in details[0]


def test_gene_set_selector_labels_use_gmt_identifiers(tmp_path):
    """Dropdown labels are the exact first-column GMT identifiers."""
    _, analysis = load_tree(tmp_path)
    long_description = (
        "A very long carbon-stress response description used in report details"
    )
    analysis.fry_results.loc[
        analysis.fry_results["gene_set"].eq("set_mixed"),
        "description",
    ] = long_description
    analysis.fry_results["display_label"] = gp._gene_set_display_labels(
        analysis.fry_results
    )

    labels = gp._gene_set_selector_labels(analysis.fry_results)

    assert labels == {
        "set_up": "set_up",
        "set_down": "set_down",
        "set_mixed": "set_mixed",
    }
    row = analysis.fry_results.loc[
        analysis.fry_results["gene_set"].eq("set_mixed")
    ].iloc[0]
    assert long_description in gp._gene_set_summary_html(row)
