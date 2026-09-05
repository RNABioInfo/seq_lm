"""Test ICA snapshot validation and report plot preparation."""

import json
import math

import numpy as np
from scipy import stats
from types import SimpleNamespace

import pytest

pytest.importorskip("ezcharts")
pd = pytest.importorskip("pandas")

from workflow_glue.qc_report_types.imodulon_plots import (  # noqa: E402
    create_activity_heatmap,
    create_component_distribution,
    create_multi_contrast_effects,
    create_time_effects,
    create_time_heatmap,
    create_timecourse,
    create_volcano,
    load_imodulon_results,
    imodulon_timecourse_enabled,
    validate_imodulon_timecourse,
)
from workflow_glue import qc_report  # noqa: E402


def write_ica_snapshot(root, *, status="ready", groups=None, times=None, values=None):
    """Write a compact, internally consistent ICA report fixture."""
    root.mkdir()
    groups = groups or ["control", "control", "early", "early", "late", "late"]
    times = times or [-5, -5, 10, 10, 40, 40]
    aliases = [f"rep{i % 2}" for i in range(len(groups))]
    ids = [f"sample-{i}" for i in range(len(groups))]
    samples = pd.DataFrame(
        {
            "sample_id": ids,
            "alias": aliases,
            "group": groups,
            "order": times,
            "source_batch_index": [7] * len(groups),
            "assigned_abundance": [20000 if status == "ready" else 5000] * len(groups),
            "ready": [status == "ready"] * len(groups),
        }
    )
    samples.to_csv(root / "sample_metadata.tsv", sep="\t", index=False)
    model = {
        "components": ["iM-A", "NA"],
        "diagnostics": {
            "rank": 2,
            "singular_values": [3.0, 1.0],
            "rank_tolerance": 1e-15,
            "condition_number": 3.0,
            "shared_gene_count": 2,
            "model_gene_count": 2,
            "gene_coverage": 1.0,
        },
    }
    provenance = {
        "schema_version": 1,
        "batch_index": 7,
        "report_sequence": 3,
        "analysis_index": 11,
        "model": model,
        "settings": {
            "log_base": 2,
            "pseudocount": 1,
            "min_read_count": 10000,
            "padj_cutoff": 0.05,
        },
        "control_sample_ids": ids[:2],
    }
    (root / "provenance.json").write_text(json.dumps(provenance))
    if status == "deferred":
        (root / "status.json").write_text(
            json.dumps(
                {
                    "status": "deferred",
                    "statistical_availability": "unavailable",
                    "sample_ids": ids,
                }
            )
        )
        return root
    (root / "status.json").write_text(
        json.dumps(
            {
                "status": "ready",
                "statistical_availability": "complete",
                "tested_count": 4,
                "untested_count": 0,
            }
        )
    )
    values = values or {
        "iM-A": [0.0, 0.2, 1.0, 1.2, -0.5, -0.3],
        "NA": [0.1, -0.1, 0.4, 0.5, 1.4, 1.5],
    }
    activity_rows = []
    for component, component_values in values.items():
        for sid, alias, group, time, value in zip(
            ids, aliases, groups, times, component_values[: len(ids)], strict=True
        ):
            activity_rows.append(
                {
                    "component_id": component,
                    "sample_id": sid,
                    "alias": alias,
                    "group": group,
                    "order": time,
                    "activity": value,
                }
            )
    activities = pd.DataFrame(activity_rows)
    activities.to_csv(root / "activities_long.tsv", sep="\t", index=False)
    activities.pivot(
        index="component_id", columns="sample_id", values="activity"
    ).to_csv(root / "activities.tsv", sep="\t")
    summary = (
        activities.groupby(["component_id", "group"], sort=False)
        .activity.agg(["mean", "std", "count"])
        .rename(columns={"std": "sd", "count": "n"})
        .reset_index()
    )
    summary.to_csv(root / "activity_summary.tsv", sep="\t", index=False)
    rows = []
    targets = list(dict.fromkeys(group for group in groups if group != "control"))
    for component in model["components"]:
        for target in targets:
            x = activities.loc[
                activities.component_id.eq(component) & activities.group.eq(target),
                "activity",
            ]
            y = activities.loc[
                activities.component_id.eq(component) & activities.group.eq("control"),
                "activity",
            ]
            effect = x.mean() - y.mean()
            row = {
                "component_id": component,
                "target_group": target,
                "control_group": "control",
                "activity_difference": effect,
                "target_mean": x.mean(),
                "control_mean": y.mean(),
                "target_sd": x.std(),
                "control_sd": y.std(),
                "target_n": len(x),
                "control_n": len(y),
                "status": "insufficient_replicates",
            }
            if min(len(x), len(y)) >= 2:
                se = math.sqrt(x.var() / len(x) + y.var() / len(y))
                row["status"] = "zero_variance" if se == 0 else "tested"
                if se > 0:
                    test = stats.ttest_ind(x, y, equal_var=False)
                    ci = test.confidence_interval()
                    row.update(
                        standard_error=se,
                        degrees_of_freedom=test.df,
                        t_statistic=test.statistic,
                        ci_lower=ci.low,
                        ci_upper=ci.high,
                        p_value=test.pvalue,
                    )
            rows.append(row)
    statistic_columns = [
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
    ]
    differential = pd.DataFrame(rows, columns=statistic_columns)
    for target in targets:
        contrast = differential.target_group.eq(target)
        adjusted = stats.false_discovery_control(
            differential.loc[contrast, "p_value"].fillna(1).to_numpy(dtype=float)
        )
        differential.loc[contrast, "adjusted_p_value"] = adjusted
    untested = differential.status.ne("tested")
    differential.loc[untested, "adjusted_p_value"] = np.nan
    differential["significant"] = differential.adjusted_p_value.le(0.05).astype(object)
    differential.loc[untested, "significant"] = None
    differential.to_csv(root / "differential_activity.tsv", sep="\t", index=False)
    tested = int((~untested).sum())
    availability = (
        "complete"
        if tested and tested == len(differential)
        else "partial"
        if tested
        else "unavailable"
    )
    (root / "status.json").write_text(
        json.dumps(
            {
                "status": "ready",
                "statistical_availability": availability,
                "tested_count": tested,
                "untested_count": int(untested.sum()),
            }
        )
    )
    pd.DataFrame(
        {
            "component_id": model["components"],
            "gene_coverage": [1, 1],
            "retained_squared_weight_fraction": [1, 0.8],
        }
    ).to_csv(root / "component_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "sample_id": ids,
            "residual_sum_squares": [0.1] * len(ids),
            "centered_sum_squares": [1] * len(ids),
            "residual_rmse": [0.2] * len(ids),
            "normalized_residual": [0.1] * len(ids),
        }
    ).to_csv(root / "projection_qc.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "model_gene_id": ["g1", "g2"],
            "gene_id": ["a", "b"],
            "transcript_id": ["t1", "t2"],
            "method": ["annotation", "annotation"],
            "status": ["mapped", "mapped"],
        }
    ).to_csv(root / "gene_mapping.tsv", sep="\t", index=False)
    return root


def test_ready_snapshot_and_plots_support_literal_na_and_irregular_time(tmp_path):
    data = load_imodulon_results(write_ica_snapshot(tmp_path / "ica"), 7, 3, 11)
    assert data.ready
    assert data.components == ("iM-A", "NA")
    validate_imodulon_timecourse(data)
    assert create_activity_heatmap(
        data, list(data.components)
    )._fig.y_range.factors == ["NA", "iM-A"]
    assert (
        create_component_distribution(data, "iM-A")._fig.yaxis.axis_label
        == "Control-centered activity"
    )
    assert (
        create_volcano(data.differential, "all", data.cutoff)._fig.xaxis.axis_label
        == "Activity difference (target − control)"
    )
    assert create_multi_contrast_effects(data)._fig.x_range.factors == [
        "early vs control",
        "late vs control",
    ]
    assert create_timecourse(data, "iM-A")._fig.xaxis.axis_label == "Elapsed time (min)"
    assert (
        create_time_effects(data, "iM-A")._fig.xaxis.axis_label
        == "Target elapsed time (min)"
    )
    heatmap = create_time_heatmap(data)
    assert not hasattr(heatmap._fig.x_range, "factors")
    assert heatmap._fig.xaxis.axis_label == "Elapsed time (min)"


def test_deferred_snapshot_loads_without_stale_activity_files(tmp_path):
    data = load_imodulon_results(
        write_ica_snapshot(tmp_path / "ica", status="deferred"), 7, 3, 11
    )
    assert not data.ready
    assert data.activities is None


def test_single_comparison_uses_one_direct_contrast(tmp_path):
    data = load_imodulon_results(
        write_ica_snapshot(
            tmp_path / "ica",
            groups=["control", "control", "treated", "treated"],
            times=[0, 0, 25, 25],
        )
    )
    contrasts = list(data.differential.groupby(["target_group", "control_group"]))
    assert len(contrasts) == 1
    assert contrasts[0][0] == ("treated", "control")
    assert (
        create_volcano(
            contrasts[0][1], "treated vs control", data.cutoff
        )._fig.title.text
        == "Differential activity — treated vs control"
    )
    validate_imodulon_timecourse(data)


@pytest.mark.parametrize(
    ("groups", "times", "message"),
    [
        (
            ["control", "control", "early", "early"],
            [0, 0, 5, 6],
            "one elapsed minute per group",
        ),
        (
            ["control", "control", "early", "late"],
            [0, 0, 5, 5],
            "one group per elapsed minute",
        ),
        (["control", "control"], [0, 0], "at least two elapsed minutes"),
        (
            ["control", "control", "early", "early", "late", "late"],
            [None, None, 5, 5, 20, 20],
            "order value for every sample",
        ),
    ],
)
def test_invalid_timecourse_metadata_is_rejected(tmp_path, groups, times, message):
    data = load_imodulon_results(
        write_ica_snapshot(tmp_path / "ica", groups=groups, times=times)
    )
    with pytest.raises(ValueError, match=message):
        validate_imodulon_timecourse(data)


def test_snapshot_identity_mismatch_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="analysis_index"):
        load_imodulon_results(write_ica_snapshot(tmp_path / "ica"), 7, 3, 12)


def test_qc_report_adds_active_ica_timecourse_main_tab(tmp_path, monkeypatch):
    snapshot = write_ica_snapshot(tmp_path / "ica")
    samples_frame = pd.DataFrame(
        {
            "Name": ["rep0", "rep1", "rep0", "rep1", "rep0", "rep1"],
            "Group": ["control", "control", "early", "early", "late", "late"],
            "Time (min)": [-5, -5, 10, 10, 40, 40],
        }
    )
    monkeypatch.setattr(
        qc_report,
        "load_qc_samples",
        lambda _path: SimpleNamespace(samples_df=samples_frame, sample_results=[]),
    )
    monkeypatch.setattr(
        qc_report, "add_sample_read_fate_sankeys", lambda _samples: None
    )
    monkeypatch.setattr(qc_report, "add_sample_hists", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(qc_report, "add_sample_2d_kdes", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        qc_report, "create_nanoplot_metrics_table", lambda _samples: pd.DataFrame()
    )
    samples = tmp_path / "samples.tsv"
    samples.write_text("unused\n")
    versions = tmp_path / "versions"
    versions.mkdir()
    (versions / "versions.txt").write_text("qc_report,workflow\n")
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"ica_analysis": True}))
    report = tmp_path / "report.html"
    args = qc_report.argparser().parse_args(
        [
            str(report),
            "--samples",
            str(samples),
            "--versions",
            str(versions),
            "--params",
            str(params),
            "--latest-batch",
            "7",
            "--refresh-seconds",
            "0",
            "--imodulon-results",
            str(snapshot),
            "--imodulon-batch",
            "7",
            "--imodulon-sequence",
            "3",
            "--imodulon-analysis-index",
            "11",
        ]
    )
    monkeypatch.chdir(tmp_path)
    qc_report.main(args)
    html = report.read_text()
    assert ">iModulon Analysis</button>" in html
    assert html.find("Quality Control") < html.find("iModulon Analysis")
    assert "Time course" in html
    assert "Differential activity" in html
    assert "Elapsed time (min)" in html
    assert "Component signs and scales belong to the supplied model" in html
    assert "Search components" in html
    assert '"iModulon Analysis"' in html


def test_ica_timecourse_is_inferred_from_complete_order_metadata(tmp_path):
    temporal = load_imodulon_results(write_ica_snapshot(tmp_path / "temporal"))
    comparison = load_imodulon_results(
        write_ica_snapshot(tmp_path / "comparison", times=[None] * 6)
    )
    assert imodulon_timecourse_enabled(temporal)
    assert not imodulon_timecourse_enabled(comparison)
    comparison.samples.loc[0, "order"] = 0
    with pytest.raises(ValueError, match="order for every sample"):
        imodulon_timecourse_enabled(comparison)


def test_removed_ica_timecourse_cli_option_is_rejected(capsys):
    with pytest.raises(SystemExit):
        qc_report.argparser().parse_args(
            [
                "report.html",
                "--samples",
                "samples.tsv",
                "--versions",
                "versions",
                "--params",
                "params.json",
                "--imodulon-timecourse",
            ]
        )
    assert "unrecognized arguments: --imodulon-timecourse" in capsys.readouterr().err


@pytest.mark.parametrize(
    "problem",
    [
        "missing_activity",
        "missing_contrast",
        "wrong_time",
        "garbage_statistic",
        "missing_tested_p",
        "invalid_probability",
        "wrong_significance",
        "wrong_count",
        "wrong_summary",
        "wrong_availability",
        "incomplete_controls",
        "wrong_readiness",
    ],
)
def test_corrupt_snapshot_is_rejected(tmp_path, problem):
    root = write_ica_snapshot(tmp_path / "ica")
    if problem in {"wrong_availability", "incomplete_controls"}:
        path = root / (
            "status.json" if problem == "wrong_availability" else "provenance.json"
        )
        contents = json.loads(path.read_text())
        if problem == "wrong_availability":
            contents["tested_count"] = 99
        else:
            contents["control_sample_ids"] = contents["control_sample_ids"][:1]
        path.write_text(json.dumps(contents))
    else:
        filename = (
            "activities_long.tsv"
            if problem in {"missing_activity", "wrong_time"}
            else "activity_summary.tsv"
            if problem == "wrong_summary"
            else "sample_metadata.tsv"
            if problem == "wrong_readiness"
            else "differential_activity.tsv"
        )
        frame = pd.read_csv(root / filename, sep="\t", dtype=str, keep_default_na=False)
        if problem == "missing_activity":
            frame = frame.iloc[1:]
        elif problem == "missing_contrast":
            frame = frame.loc[frame.target_group.ne("early")]
        else:
            column, value = {
                "wrong_time": ("order", "999"),
                "garbage_statistic": ("adjusted_p_value", "garbage"),
                "missing_tested_p": ("p_value", ""),
                "invalid_probability": ("adjusted_p_value", "2"),
                "wrong_significance": ("significant", ""),
                "wrong_count": ("target_n", "9"),
                "wrong_summary": ("mean", "99"),
                "wrong_readiness": ("ready", "false"),
            }[problem]
            frame.loc[0, column] = value
        frame.to_csv(root / filename, sep="\t", index=False)
    with pytest.raises(ValueError, match="ICA|invalid or missing numeric"):
        load_imodulon_results(root)


@pytest.mark.parametrize("singleton", [True, False])
def test_unavailable_inference_preserves_descriptive_values(tmp_path, singleton):
    kwargs = (
        {"groups": ["control", "control", "treated"], "times": [0, 0, 10]}
        if singleton
        else {"values": {"iM-A": [0, 0, 1, 1, 2, 2], "NA": [0, 0, -1, -1, -2, -2]}}
    )
    data = load_imodulon_results(write_ica_snapshot(tmp_path / "ica", **kwargs))
    assert data.differential.activity_difference.notna().all()
    assert data.differential.p_value.isna().all()
    assert data.differential.significant.eq("unavailable").all()
    assert not data.differential.is_significant.any()
    assert data.status["statistical_availability"] == "unavailable"
    create_timecourse(data, "NA")
    create_time_effects(data, "NA")


def test_volcano_only_clips_zero_and_retains_actual_probability(tmp_path):
    from bokeh.models import Span

    data = load_imodulon_results(write_ica_snapshot(tmp_path / "ica"))
    rows = data.differential.copy()
    rows["adjusted_p_value"] = [0, 0.2, 1e-310, 0.4]
    fig = create_volcano(rows, "contrast", 0.05)._fig
    source = next(
        renderer.data_source
        for renderer in fig.renderers
        if hasattr(renderer, "data_source")
    )
    assert source.data["adjusted_p_value"][0] == 0
    y = source.data["minus_log10_adjusted_p"]
    assert np.isfinite(y).all()
    assert y[0] > y[2] > y[1]
    np.testing.assert_allclose(y[1:], -np.log10([0.2, 1e-310, 0.4]))
    assert any(
        np.isclose(span.location, -np.log10(0.05)) for span in fig.select(type=Span)
    )


def test_heatmap_uses_sample_ids_even_with_colliding_display_labels(tmp_path):
    data = load_imodulon_results(write_ica_snapshot(tmp_path / "ica"))
    data.samples.loc[2, ["group", "alias"]] = ["a / b", "c"]
    data.samples.loc[3, ["group", "alias"]] = ["a", "b / c"]
    fig = create_activity_heatmap(data, list(data.components))._fig
    assert fig.x_range.factors == data.samples.sample_id.tolist()
    assert len(set(fig.x_range.factors)) == len(data.samples)
    assert (
        fig.xaxis[0].major_label_overrides["sample-2"]
        == fig.xaxis[0].major_label_overrides["sample-3"]
    )


def test_loader_accepts_actual_analysis_cli_outputs(tmp_path):
    import runpy
    import subprocess
    import sys
    from pathlib import Path

    bin_dir = Path(__file__).resolve().parents[2]
    fixture = runpy.run_path(str(bin_dir / "imodulon_analysis/tests/smoke.py"))[
        "fixture"
    ]
    fixture(tmp_path)
    cli = [sys.executable, str(bin_dir / "imodulon-analysis")]
    subprocess.run(
        cli
        + [
            "prepare",
            "--matrix",
            str(tmp_path / "matrix.csv"),
            "--annotation",
            str(tmp_path / "annotation.gtf"),
            "--output",
            str(tmp_path / "prepared"),
        ],
        check=True,
    )
    for threshold, name in ((0, "ready"), (10000, "deferred")):
        subprocess.run(
            cli
            + [
                "analyze",
                "--prepared",
                str(tmp_path / "prepared"),
                "--manifest",
                str(tmp_path / "manifest.tsv"),
                "--counts-dir",
                str(tmp_path),
                "--min-reads",
                str(threshold),
                "--output",
                str(tmp_path / name),
            ],
            check=True,
        )
        data = load_imodulon_results(tmp_path / name, 0, 0, 0)
        assert data.ready == (name == "ready")
