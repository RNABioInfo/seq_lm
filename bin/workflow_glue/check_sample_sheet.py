"""Check if a sample sheet is valid."""
import sys

from seq_run_manager.utils.sample_sheet import (
    SampleSheetError,
    parse_is_live,
    parse_sample_sheet,
)

from .util import get_named_logger, wf_parser  # noqa: ABS101

__all__ = ["argparser", "main", "parse_is_live"]


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheet")

    try:
        samples = parse_sample_sheet(args.sample_sheet)
    except SampleSheetError as error:
        sys.stdout.write(f"Parsing error: {error}")
        raise SystemExit from error

    control_count = sum(
        sample.group.lower() == args.control_group.lower() for sample in samples
    )
    if control_count < args.min_control_samples:
        sys.stdout.write(
            f"Sample sheet requires at least {args.min_control_samples} "
            f"'{args.control_group}' group samples"
        )
        raise SystemExit

    logger.info(f"Checked sample sheet {args.sample_sheet}.")


def argparser():
    """Run the entry point argument parser."""
    parser = wf_parser("check_sample_sheet")
    parser.add_argument("sample_sheet", help="Sample sheet to check")
    parser.add_argument(
        "--control_group",
        default="control",
        help="Group name used to identify control samples.",
    )
    parser.add_argument(
        "--min_control_samples",
        default=2,
        type=int,
        help="Minimum number of control-group rows required in the sample sheet.",
    )
    return parser
