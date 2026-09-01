"""Test check_sample_sheet.py."""

import os
import csv
from pathlib import Path

import pytest
from workflow_glue import check_sample_sheet


# define a list of error messages to be tested.
ERROR_MESSAGES = [
    ("sample_sheet_1.csv", "Sample sheet requires at least 2 'control' group samples"),
    ("sample_sheet_2.csv", "Parsing error: Sample aliases must be unique within each group"),
]


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return os.path.join(
        request.config.getoption("--test_data"),
        "workflow_glue",
        "check_sample_sheet")


@pytest.mark.parametrize("sample_sheet_name,error_msg", ERROR_MESSAGES)
def test_check_sample_sheet(
        capsys, test_data, sample_sheet_name, error_msg, tmp_path):
    """Test the sample sheets."""
    expected_error_message = error_msg
    source = Path(test_data) / sample_sheet_name
    sample_sheet_path = tmp_path / sample_sheet_name
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
        fieldnames = list(rows[0])
    for index, row in enumerate(rows):
        bam_dir = tmp_path / f"bams-{index}"
        bam_dir.mkdir()
        row["bam_dir"] = str(bam_dir)
    with sample_sheet_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    args = check_sample_sheet.argparser().parse_args([sample_sheet_path])
    try:
        check_sample_sheet.main(args)
    except SystemExit:
        pass
    out, _ = capsys.readouterr()
    assert out.startswith(expected_error_message)


@pytest.mark.parametrize(
    "value,expected",
    [
        (None, True),
        ("", True),
        ("  ", True),
        ("true", True),
        ("TrUe", True),
        ("false", False),
        ("FALSE", False),
    ],
)
def test_parse_is_live(value, expected):
    """The optional flag is strict but case-insensitive and defaults to live."""
    assert check_sample_sheet.parse_is_live(value) is expected


def test_parse_is_live_rejects_invalid_value():
    """Ambiguous truthy and falsey spellings are rejected."""
    with pytest.raises(ValueError, match="Invalid is_live value 'yes'"):
        check_sample_sheet.parse_is_live("yes")


@pytest.mark.parametrize("header,row", [
    ("alias,group,bam_dir", "rep_1,control,/tmp/rep_1"),
    ("alias,group,bam_dir,is_live", "rep_1,control,/tmp/rep_1,"),
])
def test_sample_sheet_without_explicit_is_live_is_valid(tmp_path, header, row):
    """Legacy and blank-valued samplesheets retain all-live semantics."""
    sample_sheet = Path(tmp_path) / "sample_sheet.csv"
    first = tmp_path / "rep_1"
    second = tmp_path / "rep_2"
    first.mkdir()
    second.mkdir()
    row = row.replace("/tmp/rep_1", str(first))
    sample_sheet.write_text(
        f"{header}\n"
        f"{row}\n"
        f"rep_2,control,{second}"
        f"{',' if 'is_live' in header else ''}\n"
    )
    args = check_sample_sheet.argparser().parse_args([str(sample_sheet)])
    check_sample_sheet.main(args)


def test_sample_sheet_rejects_invalid_is_live(capsys, tmp_path):
    """The command reports the row-level boolean parsing error."""
    sample_sheet = Path(tmp_path) / "sample_sheet.csv"
    first = tmp_path / "rep_1"
    second = tmp_path / "rep_2"
    first.mkdir()
    second.mkdir()
    sample_sheet.write_text(
        "alias,group,bam_dir,is_live\n"
        f"rep_1,control,{first},true\n"
        f"rep_2,control,{second},yes\n"
    )
    args = check_sample_sheet.argparser().parse_args([str(sample_sheet)])
    with pytest.raises(SystemExit):
        check_sample_sheet.main(args)
    out, _ = capsys.readouterr()
    assert "Invalid is_live value 'yes'" in out
