import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from imodulon_analysis.analysis import (
    analyze,
    annotation_targets,
    bh_adjust,
    differential_activity,
    prepare,
)
from smoke import fixture


def setup_model(tmp_path):
    fixture(tmp_path)
    prepared = tmp_path / "prepared"
    prepare(tmp_path / "matrix.csv", tmp_path / "annotation.gtf", None, 1, prepared)
    return prepared


def run(prepared, root, name="results", **kwargs):
    output = root / name
    analyze(prepared, root / "manifest.tsv", root, output, min_reads=0, **kwargs)
    return output


def test_end_to_end_denominator_aggregation_baseline_and_sign(tmp_path):
    prepared = setup_model(tmp_path)
    out = run(prepared, tmp_path)
    a = pd.read_csv(out / "activities.tsv", sep="\t", index_col=0)
    logs = np.log2(np.array([[25, 30, 60, 65], [40, 35, 15, 10]]) * 10000 + 1)
    centered = logs - logs[:, :2].mean(axis=1)[:, None]
    expected = centered * np.array([[1], [-0.5]])
    np.testing.assert_allclose(a.values, expected)
    np.testing.assert_allclose(a.values[:, :2].mean(axis=1), 0, atol=1e-14)
    metadata = pd.read_csv(out / "sample_metadata.tsv", sep="\t")
    assert metadata.sample_id.nunique() == 4  # aliases repeat across groups
    assert set(metadata.assigned_abundance) == {100}
    assert json.loads((out / "status.json").read_text())["status"] == "ready"
    assert (
        pd.read_csv(out / "projection_qc.tsv", sep="\t").residual_sum_squares.max()
        < 1e-20
    )
    assert {p.name for p in out.iterdir()} == {
        "activities.tsv",
        "activities_long.tsv",
        "activity_summary.tsv",
        "differential_activity.tsv",
        "sample_metadata.tsv",
        "gene_mapping.tsv",
        "component_coverage.tsv",
        "projection_qc.tsv",
        "centered_expression.tsv",
        "reference_expression.tsv",
        "provenance.json",
        "status.json",
    }


def test_exact_recovery_and_shuffled_matrix(tmp_path):
    prepared = setup_model(tmp_path)
    m = np.array([[1.0, 2.0], [3.0, -4.0]])
    with np.load(prepared / "basis.npz") as b:
        np.testing.assert_allclose(b["inverse"] @ b["weights"], np.eye(2))
    out = run(prepared, tmp_path)
    (tmp_path / "matrix.csv").write_text(",positive,negative\ng2,0,-2\ng1,1,0\n")
    second = tmp_path / "prepared2"
    prepare(tmp_path / "matrix.csv", tmp_path / "annotation.gtf", None, 1, second)
    out2 = run(second, tmp_path, "results2")
    np.testing.assert_allclose(
        pd.read_csv(out / "activities.tsv", sep="\t", index_col=0),
        pd.read_csv(out2 / "activities.tsv", sep="\t", index_col=0),
    )
    (tmp_path / "matrix.csv").write_text(",positive,negative\ng1,1,2\ng2,3,-4\n")
    third = tmp_path / "prepared3"
    prepare(tmp_path / "matrix.csv", tmp_path / "annotation.gtf", None, 1, third)
    known = np.array([[2.0, -4, 3], [-1.0, 0.1, 8]])
    with np.load(third / "basis.npz") as b:
        np.testing.assert_allclose(b["inverse"] @ (m @ known), known, atol=1e-13)


def test_live_baseline_recomputed_and_deferred(tmp_path):
    prepared = setup_model(tmp_path)
    a = run(prepared, tmp_path)
    analyze(
        prepared,
        tmp_path / "manifest.tsv",
        tmp_path,
        tmp_path / "deferred",
        min_reads=101,
    )
    assert (
        json.loads((tmp_path / "deferred/status.json").read_text())["status"]
        == "deferred"
    )
    assert not (tmp_path / "deferred/activities.tsv").exists()
    p = tmp_path / "q0.quant"
    p.write_text(p.read_text().replace("\t20\n", "\t40\n"))
    b = run(
        prepared,
        tmp_path,
        "changed",
        batch_index=3,
        report_sequence=1,
        analysis_index=8,
    )
    aa = pd.read_csv(a / "activities.tsv", sep="\t", index_col=0)
    bb = pd.read_csv(b / "activities.tsv", sep="\t", index_col=0)
    assert not np.allclose(aa.values[:, 2:], bb.values[:, 2:])
    assert json.loads((b / "provenance.json").read_text())["analysis_index"] == 8


@pytest.mark.parametrize(
    "content,message",
    [
        (",a,a\ng1,1,2\ng2,3,4\n", "duplicate"),
        (",a\ng1,1\ng1,2\n", "duplicate"),
        (",a\ng1,nan\ng2,1\n", "finite"),
        (",a,b\ng1,1,2\ng2,2,4\n", "rank deficient"),
        (",a\ng1,1\nmissing,2\n", "coverage"),
    ],
)
def test_invalid_matrices(tmp_path, content, message):
    fixture(tmp_path)
    (tmp_path / "matrix.csv").write_text(content)
    with pytest.raises(ValueError, match=message):
        prepare(
            tmp_path / "matrix.csv",
            tmp_path / "annotation.gtf",
            None,
            1,
            tmp_path / "p",
        )


def test_partial_coverage_and_explicit_map(tmp_path):
    fixture(tmp_path)
    (tmp_path / "matrix.csv").write_text(",a\nmodel1,2\nunmapped,3\n")
    (tmp_path / "map.tsv").write_text("gene_id\tmodel_gene_id\ng1\tmodel1\n")
    prepare(
        tmp_path / "matrix.csv",
        tmp_path / "annotation.gtf",
        tmp_path / "map.tsv",
        0.5,
        tmp_path / "p",
    )
    coverage = pd.read_csv(tmp_path / "p/component_coverage.tsv", sep="\t")
    assert coverage.gene_coverage.iloc[0] == 0.5
    assert coverage.retained_squared_weight_fraction.iloc[0] == pytest.approx(4 / 13)
    # Explicit maps do not acquire automatic matches for g2.
    (tmp_path / "matrix.csv").write_text(",a\nmodel1,2\ng2,3\n")
    with pytest.raises(ValueError, match="coverage"):
        prepare(
            tmp_path / "matrix.csv",
            tmp_path / "annotation.gtf",
            tmp_path / "map.tsv",
            1,
            tmp_path / "p2",
        )


@pytest.mark.parametrize("rows", [["g1\tm1", "g1\tm2"], ["g1\tm1", "g2\tm1"]])
def test_map_collisions(tmp_path, rows):
    fixture(tmp_path)
    (tmp_path / "map.tsv").write_text(
        "gene_id\tmodel_gene_id\n" + "\n".join(rows) + "\n"
    )
    with pytest.raises(ValueError, match="duplicate"):
        prepare(
            tmp_path / "matrix.csv",
            tmp_path / "annotation.gtf",
            tmp_path / "map.tsv",
            1,
            tmp_path / "p",
        )


def test_gff_parentage_and_aliases(tmp_path):
    (tmp_path / "a.gff").write_text(
        "c\ts\tgene\t1\t100\t.\t+\t.\tID=gene:g1;locus_tag=LOC1\n"
        "c\ts\tmRNA\t1\t100\t.\t+\t.\tID=transcript:t1;Parent=gene:g1\n"
        "c\ts\texon\t1\t100\t.\t+\t.\tParent=transcript:t1\n"
    )
    aliases, tx, _ = annotation_targets(tmp_path / "a.gff")
    assert aliases["LOC1"] == {"gene:g1"}
    assert tx == {"transcript:t1": "gene:g1"}
    (tmp_path / "m.csv").write_text(",component\nLOC1,1\n")
    prepare(tmp_path / "m.csv", tmp_path / "a.gff", None, 1, tmp_path / "p")
    assert (
        json.loads((tmp_path / "p/model.json").read_text())["transcript_aliases"]["t1"]
        == "transcript:t1"
    )


def test_ambiguous_parent(tmp_path):
    fixture(tmp_path)
    p = tmp_path / "annotation.gtf"
    p.write_text(
        p.read_text() + p.read_text().splitlines()[0].replace('"g1"', '"g2"') + "\n"
    )
    with pytest.raises(ValueError, match="parentage"):
        annotation_targets(p)


@pytest.mark.parametrize(
    "mutation,message",
    [
        (lambda s: s.replace("t1b\t100\t5\n", ""), "missing expected"),
        (lambda s: s.replace("t1b\t100\t5", "t1\t100\t5"), "duplicate"),
        (lambda s: s.replace("t1b\t100\t5", "t1b\t100\t-1"), "nonnegative"),
        (lambda s: s.replace("t1b\t100\t5", "t1b\t100\tinf"), "finite"),
    ],
)
def test_bad_quantifications(tmp_path, mutation, message):
    p = setup_model(tmp_path)
    q = tmp_path / "q0.quant"
    q.write_text(mutation(q.read_text()))
    with pytest.raises(ValueError, match=message):
        run(p, tmp_path)


def test_zero_target_is_measured_and_zero_total_defers(tmp_path):
    p = setup_model(tmp_path)
    q = tmp_path / "q0.quant"
    q.write_text(q.read_text().replace("\t20\n", "\t0\n"))
    run(p, tmp_path)
    q.write_text("tname\tnum_reads\nt1\t0\nt1b\t0\nt2\t0\n")
    out = run(p, tmp_path, "zero")
    assert json.loads((out / "status.json").read_text())["status"] == "deferred"


def test_welch_reference_and_translation_invariance():
    raw = np.array(
        [[1.0, 2.0, 5.0, 2.0, 3.0, 8.0], [7, 7, 7, 7, 7, 7], [1, 1, 1, 2, 3, 4]]
    )
    meta = pd.DataFrame({"group": ["control"] * 3 + ["treated"] * 3})
    centered = raw - raw[:, :3].mean(axis=1)[:, None]
    result = differential_activity(
        raw, centered, meta, ["a", "b", "c"], "control", 0.05
    )
    reference = stats.ttest_ind(raw[0, 3:], raw[0, :3], equal_var=False)
    assert result.iloc[0].p_value == pytest.approx(reference.pvalue)
    assert result.iloc[0].t_statistic == pytest.approx(reference.statistic)
    assert result.iloc[0].ci_lower == pytest.approx(reference.confidence_interval().low)
    assert result.iloc[0].ci_upper == pytest.approx(
        reference.confidence_interval().high
    )
    assert result.iloc[1].status == "zero_variance"
    assert result.iloc[2].status == "tested"
    shift = differential_activity(
        centered, centered, meta, ["a", "b", "c"], "control", 0.05
    )
    np.testing.assert_allclose(result.p_value, shift.p_value, equal_nan=True)
    np.testing.assert_allclose(
        result.adjusted_p_value, shift.adjusted_p_value, equal_nan=True
    )
    np.testing.assert_allclose(bh_adjust([0.01, 1, 0.03]), [0.03, 1, 0.045])
    assert result.iloc[1].significant is None


def test_small_group_and_control_only(tmp_path):
    p = setup_model(tmp_path)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("\n".join(manifest.read_text().splitlines()[:-1]) + "\n")
    out = run(p, tmp_path)
    result = pd.read_csv(out / "differential_activity.tsv", sep="\t")
    assert set(result.status) == {"insufficient_replicates"}
    manifest.write_text("\n".join(manifest.read_text().splitlines()[:-1]) + "\n")
    out = run(p, tmp_path, "control_only")
    assert pd.read_csv(out / "differential_activity.tsv", sep="\t").empty


def test_ambiguous_control_labels(tmp_path):
    p = setup_model(tmp_path)
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(manifest.read_text().replace("control", "Control", 1))
    with pytest.raises(ValueError, match="unambiguous"):
        run(p, tmp_path)


def test_gff_transcript_id_and_id_are_aliases(tmp_path):
    (tmp_path / "a.gff").write_text(
        "c\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        "c\ts\tmRNA\t1\t100\t.\t+\t.\tID=rna1;transcript_id=t1;Parent=g1\n"
    )
    (tmp_path / "m.csv").write_text(",a\ng1,1\n")
    prepare(tmp_path / "m.csv", tmp_path / "a.gff", None, 1, tmp_path / "p")
    model = json.loads((tmp_path / "p/model.json").read_text())
    assert model["transcript_map"] == {"t1": "g1"}
    assert model["transcript_aliases"]["rna1"] == "t1"


def test_config_schema_defaults():
    import re

    root = Path(__file__).resolve().parents[3]
    schema = json.loads((root / "nextflow_schema.json").read_text())
    properties = schema["definitions"]["imodulon_analysis_options"]["properties"]
    config = (root / "nextflow.config").read_text()
    assert len(properties) == 8
    assert "ica_timecourse" not in properties
    for name, definition in properties.items():
        match = re.search(r"^\s*" + name + r" = (.+)$", config, re.M)
        assert match, name
        if definition["type"] == "boolean":
            assert match.group(1) == str(definition["default"]).lower()
        elif "default" in definition:
            assert float(match.group(1)) == definition["default"]
        else:
            assert match.group(1) == "null"
    assert {"$ref": "#/definitions/imodulon_analysis_options"} in schema["allOf"]


def test_noninteger_constant_groups_are_untestable():
    raw = np.array([[0.1, 0.1, 0.1, 0.2, 0.2, 0.2]])
    metadata = pd.DataFrame({"group": ["control"] * 3 + ["treated"] * 3})
    result = differential_activity(
        raw, raw - 0.1, metadata, ["constant"], "control", 0.05
    )
    assert result.iloc[0].status == "zero_variance"
    assert result.iloc[0].target_sd == 0
    assert result.iloc[0].control_sd == 0
    assert result.iloc[0].p_value is None
