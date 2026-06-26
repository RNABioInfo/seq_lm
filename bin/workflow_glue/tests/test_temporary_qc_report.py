"""Test TEMPORARY temporary_qc_report.py."""
import json

import pytest

pytest.importorskip("ezcharts")
pytest.importorskip("pandas")

from workflow_glue import temporary_qc_report  # noqa: E402


def test_temporary_qc_report_writes_html(tmp_path):
    """TEMPORARY: report lists the current QC result samples."""
    samples = tmp_path / "temporary_qc_report_samples.tsv"
    samples.write_text(
        "sample_id\talias\tgroup\ttype\tchunks_seen\tlatest_batch_index\n"
        "sample_1\tcontrol_1\tcontrol\tCONTROL\t2\t2\n"
        "sample_3\tcondition_1\tcondition\tCONDITION\t1\t1\n"
    )

    versions = tmp_path / "versions"
    versions.mkdir()
    (versions / "versions.txt").write_text("temporary_qc_report,temporary\n")

    params = tmp_path / "params.json"
    params.write_text(json.dumps({"temporary_qc_report": True}))

    report = tmp_path / "temporary_qc_report.html"
    args = temporary_qc_report.argparser().parse_args([
        str(report),
        "--samples", str(samples),
        "--versions", str(versions),
        "--params", str(params),
        "--latest-batch", "2",
        "--refresh-seconds", "0",
    ])

    temporary_qc_report.main(args)

    html = report.read_text()
    assert "TEMPORARY QC REPORT" in html
    assert "sample_1" in html
    assert "sample_3" in html
