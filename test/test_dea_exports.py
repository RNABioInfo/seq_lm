"""Native edgeR exports: independent pair fits and unfiltered uncertainty rows."""

import csv
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "bin/de_analysis/edgeR_analysis.R"


def run_fixture(root, third=False, zero=False, flat=False, min_reads=10000,
                stopping_gmt=None, annotation=None, zero_features=(), error=None,
                warning=None):
    import random

    random.seed(18)
    rows = []
    for group in ["T0", "T1"] + (["T4"] if third else []):
        for j in range(3):
            name = f"{group}_{j}"
            file = root / (name + ".quant")
            with file.open("w") as f:
                f.write("tname\tnum_reads\n")
                for i in range(300):
                    count = (
                        0
                        if zero and group == "T4"
                        else random.randint(50, 110)
                        * (3 if group == "T1" and i < 30 else 1)
                    )
                    if flat:
                        count = 1
                    if i in zero_features:
                        count = 0
                    f.write(f"g{i}\t{count}\n")
            rows.append({"name": name, "group": group, "count_file": str(file)})
    manifest = root / "manifest.tsv"
    with manifest.open("w") as f:
        w = csv.DictWriter(
            f, fieldnames=["name", "group", "count_file"], delimiter="\t"
        )
        w.writeheader()
        w.writerows(rows)
    definitions = "contrast_id\ttarget_group\treference_group\npair\tT1\tT0\n" + (
        "other\tT4\tT0\n" if third else ""
    )
    (root / "contrasts.tsv").write_text(definitions)
    extra_args = []
    if stopping_gmt is not None:
        path = root / "stopping genes.gmt"
        path.write_text(stopping_gmt)
        extra_args += ["--dea_stopping_genes", str(path)]
    if annotation is not None:
        path = root / "annotation.gtf"
        path.write_text(annotation)
        extra_args += ["--annotation", str(path)]
    result = subprocess.run(
        [
            "Rscript",
            str(SCRIPT),
            "--quant_manifest",
            str(manifest),
            "--contrasts",
            str(root / "contrasts.tsv"),
            "--min_read_count",
            str(min_reads),
            "--output_dir",
            str(root / "out"),
        ] + extra_args,
        check=False,
        text=True,
        capture_output=True,
    )
    if error:
        assert result.returncode != 0
        assert error in result.stdout + result.stderr
    else:
        assert result.returncode == 0, result.stdout + result.stderr
    if warning:
        assert warning in result.stderr
    return root / "out"


def test_pair_is_independent_of_other_group_and_exports_intervals(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    first = run_fixture(a)
    second = run_fixture(b, third=True)
    assert (first / "pair/edgeR_results.tsv").read_text() == (
        second / "pair/edgeR_results.tsv"
    ).read_text()
    assert (first / "pair/precision.tsv").read_text() == (
        second / "pair/precision.tsv"
    ).read_text()
    assert (second / "pair/precision_model.rds").stat().st_size > 1000
    with (second / "pair/precision_samples.tsv").open() as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    assert len(samples) == 6 and {s["group"] for s in samples} == {"T0", "T1"}
    with (second / "pair/precision.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 300
    assert set(rows[0]) == {
        "feature_id", "R", "ci_lower", "ci_upper", "unshrunk_logFC"
    }


def test_zero_depth_other_target_is_deferred_without_blocking_pair(tmp_path):
    result = run_fixture(tmp_path, third=True, zero=True)
    assert (result / "pair/precision.tsv").exists()
    assert (
        result / "other/analysis_deferred.txt"
    ).read_text().strip() == "insufficient_read_depth"
    assert not (result / "other/precision.tsv").exists()


def test_unavailable_cohort_is_not_passed_to_enrichment(tmp_path):
    result = run_fixture(tmp_path, third=True, zero=True)
    with (result / "sample_metadata.tsv").open() as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    assert {r["group"] for r in samples} == {"T0", "T1"}
    with (result / "contrasts.tsv").open() as f:
        fitted = list(csv.DictReader(f, delimiter="\t"))
    assert [r["contrast_id"] for r in fitted] == ["pair"]
    with (result / "contrast_status.tsv").open() as f:
        statuses = list(csv.DictReader(f, delimiter="\t"))
    assert len(statuses) == 2 and statuses[1]["reason"] == "insufficient_read_depth"


def test_read_depth_without_retained_genes_defers_the_entire_snapshot(tmp_path):
    result = run_fixture(tmp_path, flat=True, min_reads=100)
    assert (result / "analysis_deferred.txt").exists()
    assert not (result / "pair/precision.tsv").exists()
    assert (
        result / "pair/analysis_deferred.txt"
    ).read_text().strip() == "no_retained_features"


def test_stopping_union_reuses_annotation_mapping_without_changing_fit(tmp_path):
    annotation = (
        'chr\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "shared"; transcript_id "g1";\n'
        'chr\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id "shared"; transcript_id "g2";\n'
    )
    default = tmp_path / "default"
    selected = tmp_path / "selected"
    default.mkdir()
    selected.mkdir()
    first = run_fixture(default, third=True, annotation=annotation, zero_features=(299,))
    second = run_fixture(
        selected, third=True, annotation=annotation, zero_features=(299,),
        stopping_gmt="A\tfirst\tshared\tg1\tunknown\nB\tsecond\tg2\tg299\tshared\n",
        warning="1 --dea_stopping_genes member(s) did not resolve: unknown",
    )
    assert (second / "dea_stopping_features.tsv").read_text() == "feature_id\ng1\ng2\ng299\n"
    with (second / "dea_stopping_resolution.tsv").open() as f:
        resolution = list(csv.DictReader(f, delimiter="\t"))
    assert {r["feature_id"] for r in resolution if r["gmt_member"] == "shared"} == {"g1", "g2"}
    assert [r["feature_id"] for r in resolution if r["gmt_member"] == "unknown"] == ["NA"]
    for pair in ("pair", "other"):
        for name in ("edgeR_results.tsv", "precision.tsv", "precision_samples.tsv", "precision_methods.tsv"):
            assert (first / pair / name).read_bytes() == (second / pair / name).read_bytes()
        with (second / pair / "precision.tsv").open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        assert len(rows) == 300
        assert next(r for r in rows if r["feature_id"] == "g299")["R"] == "NA"


@pytest.mark.parametrize("gmt,error", [
    ("", "must contain at least one GMT member"),
    ("A\tempty\n", "must contain at least one GMT member"),
    ("invalid\n", "Invalid GMT record"),
    ("A\tx\tg1\nA\ty\tg2\n", "Duplicate gene-set name"),
    ("A\tx\tunknown\n", "No --dea_stopping_genes members resolved"),
])
def test_invalid_stopping_gmt_is_rejected(tmp_path, gmt, error):
    run_fixture(tmp_path, stopping_gmt=gmt, error=error)


def test_deferred_fit_does_not_require_selection_exports(tmp_path):
    output = run_fixture(tmp_path, flat=True, min_reads=100, stopping_gmt="A\tx\tg1\n")
    assert (output / "analysis_deferred.txt").exists()
    assert not (output / "dea_stopping_features.tsv").exists()
