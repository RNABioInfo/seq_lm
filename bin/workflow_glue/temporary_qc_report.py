"""TEMPORARY: create a live QC placeholder report."""
from pathlib import Path

from ezcharts.components.reports import labs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd

from .report_compat import labs_report
from .util import get_named_logger, wf_parser  # noqa: ABS101


TEMPORARY_NOTICE = (
    "TEMPORARY QC REPORT - remove this placeholder when the permanent live "
    "QC report is implemented."
)
EXPECTED_COLUMNS = [
    "sample_id",
    "alias",
    "group",
    "type",
    "chunks_seen",
    "latest_batch_index",
]
DISPLAY_COLUMNS = {
    "sample_id": "Sample ID",
    "alias": "Alias",
    "group": "Group",
    "type": "Type",
    "chunks_seen": "QC chunks seen",
    "latest_batch_index": "Latest batch index",
}


def load_temporary_qc_samples(samples_path):
    """Load the temporary QC report sample table."""
    samples = pd.read_csv(samples_path, sep="\t", dtype=str).fillna("")
    missing_columns = [column for column in EXPECTED_COLUMNS if column not in samples]
    if missing_columns:
        raise ValueError(
            "Temporary QC report samples table is missing columns: "
            + ", ".join(missing_columns)
        )
    return samples[EXPECTED_COLUMNS].rename(columns=DISPLAY_COLUMNS)


def add_temporary_auto_refresh(report_path, refresh_seconds):
    """TEMPORARY: refresh the open report while live chunks are still arriving."""
    if refresh_seconds <= 0:
        return

    path = Path(report_path)
    refresh_tag = f'<meta http-equiv="refresh" content="{refresh_seconds}">'
    html = path.read_text()
    if refresh_tag in html:
        return
    if "</head>" in html:
        html = html.replace("</head>", f"    {refresh_tag}\n</head>", 1)
    else:
        html = refresh_tag + "\n" + html
    path.write_text(html)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("TemporaryQC")
    samples = load_temporary_qc_samples(args.samples)

    report = labs_report(
        labs,
        "TEMPORARY seq_lm QC report - remove later",
        "temporary_qc_report",
        args.params,
        args.versions,
        "temporary",
    )

    with report.add_section("TEMPORARY notice", "Temporary notice"):
        DataTable.from_pandas(pd.DataFrame([{
            "Temporary notice": TEMPORARY_NOTICE,
            "Latest batch index": args.latest_batch,
        }]))

    with report.add_section(
        "TEMPORARY current QC result samples",
        "Temporary QC samples",
    ):
        DataTable.from_pandas(samples)

    report.write(args.report)
    add_temporary_auto_refresh(args.report, args.refresh_seconds)
    logger.info(f"Temporary QC report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("temporary_qc_report")
    parser.add_argument("report", help="Temporary QC report output HTML file")
    parser.add_argument(
        "--samples",
        required=True,
        help="TSV containing the current temporary QC report sample rows.",
    )
    parser.add_argument(
        "--versions",
        required=True,
        help="Directory containing CSVs containing name,version.",
    )
    parser.add_argument(
        "--params",
        required=True,
        help="JSON file containing the workflow parameter key/values.",
    )
    parser.add_argument(
        "--latest-batch",
        default="unknown",
        help="Latest live QC batch index represented in the temporary report.",
    )
    parser.add_argument(
        "--refresh-seconds",
        default=15,
        type=int,
        help="TEMPORARY browser auto-refresh interval; use 0 to disable.",
    )
    return parser
