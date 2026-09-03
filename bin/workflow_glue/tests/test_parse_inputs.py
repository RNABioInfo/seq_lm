"""Tests for QC report input parsing."""

import gzip

import pytest

pd = pytest.importorskip("pandas")

from workflow_glue.qc_report_types import parse_inputs  # noqa: E402


def write_nanoplot(path, rows):
    """Write a gzipped NanoPlot-data style TSV."""
    with gzip.open(path, "wt") as handle:
        handle.write("\t".join(parse_inputs.NANOPLOT_COLUMNS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(str(row[column]) for column in parse_inputs.NANOPLOT_COLUMNS)
                + "\n"
            )


def test_aggregate_nanoplot_accepts_legacy_headerless_empty_chunk(tmp_path):
    """A terminal chunk with no usable reads does not break aggregation."""
    populated = tmp_path / "nanoplot_data_chunk_0.tsv.gz"
    empty = tmp_path / "nanoplot_data_chunk_1.tsv.gz"
    row = {
        "readIDs": "read_1",
        "quals": 10,
        "aligned_quals": 9,
        "lengths": 100,
        "aligned_lengths": 90,
        "mapQ": 30,
        "percentIdentity": 98.0,
    }
    write_nanoplot(populated, [row])
    with gzip.open(empty, "wt"):
        pass

    result = parse_inputs.aggregate_nanoplot([populated, empty])

    assert result is not None
    assert result.to_dict(orient="records") == [row]


def test_aggregate_nanoplot_all_empty_chunks_keep_expected_schema(tmp_path):
    """Samples with only empty QC chunks retain a report-compatible schema."""
    empty = tmp_path / "nanoplot_data_chunk_0.tsv.gz"
    with gzip.open(empty, "wt"):
        pass

    result = parse_inputs.aggregate_nanoplot([empty])

    assert result is not None
    assert result.empty
    assert result.columns.tolist() == parse_inputs.NANOPLOT_COLUMNS


def test_aggregate_nanoplot_rejects_nonempty_malformed_chunk(tmp_path):
    """The empty-file fallback does not hide malformed non-empty output."""
    malformed = tmp_path / "nanoplot_data_chunk_0.tsv.gz"
    with gzip.open(malformed, "wt") as handle:
        handle.write("wrong\tcolumns\n1\t2\n")

    with pytest.raises(ValueError, match="NanoPlot data table is missing columns"):
        parse_inputs.aggregate_nanoplot([malformed])


def test_load_qc_samples_preserves_optional_elapsed_minutes(tmp_path, monkeypatch):
    """The report exposes order metadata without requiring it in legacy TSVs."""
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "name\tgroup\torder\tchunks_seen\tlatest_batch_index\tqc_dir\n"
        "rep_1\tcontrol\t-5\t1\t0\tqc\n"
    )
    monkeypatch.setattr(
        parse_inputs,
        "retrieve_flagstat_paths",
        lambda *_args: [tmp_path / "flagstat.tsv"],
    )
    monkeypatch.setattr(
        parse_inputs,
        "retrieve_nanoplot_paths",
        lambda *_args: [tmp_path / "nanoplot.tsv.gz"],
    )
    monkeypatch.setattr(
        parse_inputs,
        "aggregate_flagstat",
        lambda _paths: parse_inputs.FlagstatResult(1, 1, 0, 0, 1, 1),
    )
    monkeypatch.setattr(
        parse_inputs,
        "aggregate_nanoplot",
        lambda _paths: pd.DataFrame(columns=parse_inputs.NANOPLOT_COLUMNS),
    )

    result = parse_inputs.load_qc_samples(samples)

    assert result.samples_df["Time (min)"].tolist() == ["-5"]
