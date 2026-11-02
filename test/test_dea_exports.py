"""Native edgeR exports: independent pair fits and unfiltered uncertainty rows."""

import csv
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "bin/de_analysis/edgeR_analysis.R"


def run_fixture(root, third=False, zero=False, flat=False, min_reads=10000):
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
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
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
