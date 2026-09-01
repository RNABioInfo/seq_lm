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
