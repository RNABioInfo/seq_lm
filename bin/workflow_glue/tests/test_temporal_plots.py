"""Tests for descriptive temporal gene-set report plots."""

import json

import numpy as np
import pytest

pytest.importorskip("ezcharts")
pd = pytest.importorskip("pandas")

from bokeh.models import HoverTool, Legend  # noqa: E402
from ezcharts.components.reports import labs  # noqa: E402
from ezcharts.layout.snippets import Tabs  # noqa: E402

from workflow_glue.report_compat import labs_report  # noqa: E402
from workflow_glue.qc_report_types import differential_plots as dp  # noqa: E402
from workflow_glue.qc_report_types import temporal_plots as tp  # noqa: E402


def make_differential(samples):
    """Build the metadata portion needed by temporal validation."""
    metadata = pd.DataFrame(samples)[["Name", "Group"]].rename(
        columns={"Name": "name", "Group": "group"}
    )
    metadata = metadata.set_index("name", drop=False)
    return dp.DifferentialResult(
        feature_counts=pd.DataFrame(),
        bcv_data=pd.DataFrame(),
        mds_data=pd.DataFrame(),
        sample_metadata=metadata,
        contrasts=[],
        condition_colors=dp.condition_palette(metadata["group"].tolist()),
    )


def make_temporal_result():
    """Build an irregular three-time-point trajectory with a singleton."""
    samples = ["s0a", "s0b", "s15a", "s15b", "s45a"]
    groups = ["t0", "t0", "t15", "t15", "t45"]
    times = [0, 0, 15, 15, 45]
    metadata = pd.DataFrame(
        {
            "sample": samples,
            "group": groups,
            "time_minutes": times,
            "sample_order": np.arange(len(samples)),
        }
    ).set_index("sample", drop=False)
    scores = pd.DataFrame(
        {
            "gene_set": ["set_a"] * len(samples),
            "sample": samples,
            "group": groups,
            "score": [-1.0, 1.0, 2.0, 4.0, 6.0],
            "time_minutes": times,
            "sample_order": np.arange(len(samples)),
        }
    )
    log_cpm = pd.DataFrame(
        {
            "s0a": [1.0, 8.0],
            "s0b": [3.0, 6.0],
            "s15a": [5.0, 4.0],
            "s15b": [7.0, 2.0],
            "s45a": [9.0, 1.0],
        },
        index=["feature_a", "feature_b"],
    )
    return tp.TemporalResult(
        metadata=metadata,
        scores_long=scores,
        log_cpm=log_cpm,
        members={"set_a": ("feature_a", "feature_b")},
        feature_labels={"feature_a": "gene_a", "feature_b": "gene_b"},
        gene_set_order=["set_a"],
    )


def test_temporal_metadata_sorts_signed_integer_minutes():
    """Elapsed minutes are numeric and may be irregular or negative."""
    samples = pd.DataFrame(
        [
            {"Name": "late", "Group": "late_group", "Time (min)": "45"},
            {"Name": "early", "Group": "early_group", "Time (min)": "-5"},
            {"Name": "middle", "Group": "middle_group", "Time (min)": "15"},
        ]
    )
    differential = make_differential(samples.to_dict(orient="records"))

    observed = tp._temporal_metadata(samples, differential)

    assert observed["sample"].tolist() == ["early", "middle", "late"]
    assert observed["time_minutes"].tolist() == [-5, 15, 45]


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        (
            [
                {"Name": "a", "Group": "same", "Time (min)": 0},
                {"Name": "b", "Group": "same", "Time (min)": 15},
            ],
            "every group to use one elapsed minute",
        ),
        (
            [
                {"Name": "a", "Group": "first", "Time (min)": 0},
                {"Name": "b", "Group": "second", "Time (min)": 0},
            ],
            "every elapsed minute to identify one group",
        ),
        (
            [
                {"Name": "a", "Group": "same", "Time (min)": 0},
                {"Name": "b", "Group": "same", "Time (min)": 0},
            ],
            "at least two distinct elapsed minutes",
        ),
    ],
)
def test_temporal_metadata_rejects_invalid_single_trajectory(rows, message):
    """The first temporal version represents exactly one ordered trajectory."""
    samples = pd.DataFrame(rows)
    differential = make_differential(rows)

    with pytest.raises(ValueError, match=message):
        tp._temporal_metadata(samples, differential)


def test_temporal_metadata_rejects_noninteger_minutes():
    """Fractional or labelled time values are not silently reinterpreted."""
    rows = [
        {"Name": "a", "Group": "first", "Time (min)": "0"},
        {"Name": "b", "Group": "second", "Time (min)": "1.5"},
    ]
    samples = pd.DataFrame(rows)

    with pytest.raises(ValueError, match="signed integer order values"):
        tp._temporal_metadata(samples, make_differential(rows))


def test_temporal_score_summary_uses_sample_sd_and_stable_jitter():
    """Score means, sample SDs, replicate counts, and jitter are deterministic."""
    data = make_temporal_result()
    colors = {"t0": "#111111", "t15": "#222222", "t45": "#333333"}

    summary, raw = tp.prepare_temporal_scores(data, "set_a", colors)
    repeated = tp.prepare_temporal_scores(data, "set_a", colors)[1]

    assert summary["time_minutes"].tolist() == [0, 15, 45]
    assert summary["mean_score"].tolist() == pytest.approx([0.0, 3.0, 6.0])
    assert summary["sd_score"].iloc[:2].tolist() == pytest.approx(
        [np.sqrt(2), np.sqrt(2)]
    )
    assert pd.isna(summary["sd_score"].iloc[2])
    assert summary["n"].tolist() == [2, 2, 1]
    assert raw["plot_time"].tolist() == repeated["plot_time"].tolist()
    assert raw.loc[raw["time_minutes"].eq(45), "plot_time"].iloc[0] == 45


def test_temporal_gene_summary_and_plots_include_sd_and_hover():
    """Gene lines use time-point mean log CPM, SD, labels, and interactive legends."""
    data = make_temporal_result()
    summary = tp.prepare_temporal_gene_expression(data, "set_a")
    feature_a = summary.loc[summary["feature_id"].eq("feature_a")]

    assert feature_a["time_minutes"].tolist() == [0, 15, 45]
    assert feature_a["mean_log_cpm"].tolist() == pytest.approx([2.0, 6.0, 9.0])
    assert feature_a["sd_log_cpm"].iloc[:2].tolist() == pytest.approx(
        [np.sqrt(2), np.sqrt(2)]
    )
    assert feature_a["sd_label"].iloc[2] == "unavailable (n=1)"

    score_plot = tp.create_temporal_score_plot(
        data,
        "set_a",
        {"t0": "#111111", "t15": "#222222", "t45": "#333333"},
    )
    gene_plot = tp.create_temporal_gene_plot(data, "set_a")

    assert score_plot._fig.title.text == "GSVA score over time — set_a"
    assert gene_plot._fig.title.text == "Gene expression over time — set_a"
    assert gene_plot._fig.yaxis[0].axis_label == "log2(TMM-normalized CPM + 1)"
    assert any("Sample" in dict(tool.tooltips) for tool in score_plot._fig.select({"type": HoverTool}))
    legends = list(gene_plot._fig.select({"type": Legend}))
    assert len(legends) == 1
    assert legends[0].click_policy == "hide"
    assert [item.label["value"] for item in legends[0].items] == ["gene_a", "gene_b"]


def test_large_gene_set_omits_oversized_legend():
    """Large temporal line plots remain hover-driven instead of adding huge legends."""
    data = make_temporal_result()
    members = tuple(f"feature_{index}" for index in range(21))
    data.members["set_a"] = members
    data.log_cpm = pd.DataFrame(
        np.arange(21 * 5, dtype=float).reshape(21, 5),
        index=members,
        columns=data.metadata["sample"],
    )
    data.feature_labels = {feature_id: feature_id for feature_id in members}

    plot = tp.create_temporal_gene_plot(data, "set_a")

    assert list(plot._fig.select({"type": Legend})) == []


def test_temporal_analysis_writes_one_shared_gene_set_dropdown(tmp_path):
    """A report tab serializes both figures beneath one gene-set selector."""
    versions = tmp_path / "versions"
    versions.mkdir()
    (versions / "versions.txt").write_text("qc_report,workflow\n")
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"timeline_analysis": True}))
    report = labs_report(
        labs,
        "Temporal test report",
        "qc_report",
        params,
        versions,
        "workflow",
    )
    with report.add_section("Analysis", "Analysis"):
        primary_tabs = Tabs()
        with primary_tabs.add_tab("Temporal Analysis"):
            tp.add_temporal_analysis(
                make_temporal_result(),
                {"t0": "#111111", "t15": "#222222", "t45": "#333333"},
            )

    output = tmp_path / "temporal_report.html"
    report.write(output)
    html = output.read_text()

    assert html.count(">Gene set<") == 1
    assert "GSVA score over time" in html
    assert "Gene expression over time" in html
    assert "Descriptive temporal view" in html
