import argparse
import importlib.machinery
import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).parents[2] / "stability_de_analysis"


def _load_module():
    loader = importlib.machinery.SourceFileLoader("stability_de_analysis", str(SCRIPT))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules[loader.name] = module
    loader.exec_module(module)
    return module


stability = _load_module()


def _write_snapshot(root: Path, shifts=None, calls=None):
    shifts = shifts or {}
    calls = calls or {}
    root.mkdir()
    pd.DataFrame(
        {"sample": ["c1", "c2", "a1", "a2", "b1", "b2"], "group": ["control", "control", "A", "A", "B", "B"]}
    ).to_csv(root / "sample_metadata.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_id": ["f1", "f2", "f3"], "c1": [1, 2, 3]}).to_csv(
        root / "feature_counts.tsv", sep="\t", index=False
    )
    for group in ("A", "B"):
        contrast = root / f"group_{group}_vs_control"
        contrast.mkdir()
        called = set(calls.get(group, {"f1"}))
        pd.DataFrame(
            {
                "feature_id": ["f1", "f2", "f3"],
                "logFC": [2.0 + shifts.get(group, 0.0), 0.1, -0.1],
                "FDR": [0.01 if "f1" in called else 0.5, 0.5, 0.5],
            }
        ).to_csv(contrast / "edgeR_results.tsv", sep="\t", index=False)


def _args(**overrides):
    values = dict(
        max_feature_diff_fraction=0.05,
        max_median_abs_lfc_delta=0.1,
        min_jaccard_similarity=0.9,
        max_call_churn_fraction=0.1,
        max_lost_call_fraction=0.1,
        max_fdr=0.05,
        min_abs_lfc=1.0,
        min_de_calls_for_fraction_metrics=20,
        max_small_set_call_changes=2,
    )
    values.update(overrides)
    return argparse.Namespace(**values)


def test_assesses_every_contrast_and_accepts_threshold_boundary(tmp_path):
    previous = tmp_path / "previous"
    current = tmp_path / "current"
    _write_snapshot(previous)
    _write_snapshot(current, shifts={"A": 0.1})
    rows = stability.assess_results(previous, current, _args())
    assert [row.contrast_id for row in rows] == ["group_A_vs_control", "group_B_vs_control"]
    assert all(row.stable for row in rows)
    assert rows[0].median_abs_lfc_delta == pytest.approx(0.0)


def test_large_call_set_uses_fraction_metrics():
    features = [f"f{index}" for index in range(20)]
    previous = pd.DataFrame({"feature_id": features, "logFC": [2.0] * 20, "FDR": [0.01] * 20})
    current = previous.copy()
    current.loc[0:2, "FDR"] = 0.5
    de_results = stability.DEResults(current, previous)
    assert not stability.de_calls_are_stable(0.9, 0.1, 0.1, 0.05, 1.0, 20, 2, de_results)


def test_small_call_set_uses_absolute_change_limit():
    previous = pd.DataFrame({"feature_id": ["a", "b", "c"], "logFC": [2.0] * 3, "FDR": [0.01] * 3})
    current = previous.copy()
    current.loc[0:1, "FDR"] = 0.5
    de_results = stability.DEResults(current, previous)
    assert stability.de_calls_are_stable(0.99, 0.01, 0.01, 0.05, 1.0, 20, 2, de_results)
    current.loc[2, "FDR"] = 0.5
    assert not stability.de_calls_are_stable(
        0.99, 0.01, 0.01, 0.05, 1.0, 20, 2, stability.DEResults(current, previous)
    )


def test_rejects_nonfinite_statistics(tmp_path):
    previous = tmp_path / "previous"
    current = tmp_path / "current"
    _write_snapshot(previous)
    _write_snapshot(current)
    path = current / "group_A_vs_control" / "edgeR_results.tsv"
    table = pd.read_csv(path, sep="\t")
    table.loc[0, "logFC"] = float("inf")
    table.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="non-finite"):
        stability.assess_results(previous, current, _args())


def test_cli_writes_structured_multi_contrast_tsv(tmp_path):
    previous = tmp_path / "previous"
    current = tmp_path / "current"
    output = tmp_path / "assessment.tsv"
    _write_snapshot(previous)
    _write_snapshot(current)
    command = [
        sys.executable, str(SCRIPT), "--previous-results", str(previous), "--current-results", str(current), "--output", str(output),
        "--max-feature-diff-fraction", "0.05", "--max-median-abs-lfc-delta", "0.1",
        "--min-jaccard-similarity", "0.9", "--max-call-churn-fraction", "0.1",
        "--max-lost-call-fraction", "0.1", "--max-fdr", "0.05", "--min-abs-lfc", "1",
        "--min-de-calls-for-fraction-metrics", "20", "--max-small-set-call-changes", "2",
    ]
    subprocess.run(command, check=True)
    result = pd.read_csv(output, sep="\t")
    assert result["contrast_id"].tolist() == ["group_A_vs_control", "group_B_vs_control"]
    assert result["stable"].tolist() == [True, True]
