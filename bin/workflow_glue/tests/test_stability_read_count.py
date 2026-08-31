"""Tests for the DEA read-count stability command."""

from pathlib import Path
import subprocess
import sys

SCRIPT = Path(__file__).parents[2] / "stability_read_count"


def _write_count(path: Path, reads: int, feature_id: str = "tx1") -> None:
    path.write_text(
        f"tname\tlen\tnum_reads\n{feature_id}\t100\t{reads}\n"
    )


def _run_check(
    tmp_path: Path,
    reads: list[int],
    feature_ids: list[str] | None = None,
) -> str:
    counts_dir = tmp_path / "counts"
    counts_dir.mkdir()
    rows = []
    for index, read_count in enumerate(reads):
        count_file = f"sample_{index}.quant"
        _write_count(
            counts_dir / count_file,
            read_count,
            feature_ids[index] if feature_ids else "tx1",
        )
        rows.append(
            {
                "name": f"sample_{index}",
                "group": "control" if index < 2 else "treated",
                "count_file": count_file,
            }
        )
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "name\tgroup\tcount_file\n"
        + "".join(
            f"{row['name']}\t{row['group']}\t{row['count_file']}\n"
            for row in rows
        )
    )

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--min_read_count",
            "10000",
            "--min_replicate_sample_count",
            "2",
            "--metadata",
            str(metadata),
            "--counts_dir",
            str(counts_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def test_stability_read_count_requires_enough_deep_samples_per_group(tmp_path):
    assert (
        _run_check(tmp_path, [10000, 12000, 11000, 9000])
        == "insufficient_read_depth"
    )


def test_stability_read_count_accepts_all_groups_at_threshold(tmp_path):
    assert _run_check(tmp_path, [10000, 12000, 11000, 13000]) == "ready"


def test_stability_read_count_rejects_no_shared_feature_ids(tmp_path):
    """Disjoint quantification identifiers defer edgeR instead of crashing it."""
    assert (
        _run_check(
            tmp_path,
            [10000, 12000, 11000, 13000],
            ["tx1", "tx2", "tx3", "tx4"],
        )
        == "no_matching_feature_ids"
    )


def test_stability_read_count_accepts_at_least_one_shared_feature_id(tmp_path):
    assert (
        _run_check(
            tmp_path,
            [10000, 12000, 11000, 13000],
            ["tx1", "tx1", "tx1", "tx1"],
        )
        == "ready"
    )
