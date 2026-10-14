"""Check if a sample sheet is valid."""
import argparse
import csv
import sys

from seq_run_manager.utils.sample_sheet import (
    SampleSheetError,
    parse_is_live,
    parse_sample_sheet,
)

from .gene_contrasts import read_gene_contrasts, resolve_ica_groups
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
    if not args.gene_contrast_file and (args.differential_expression or not args.ica_analysis) and control_count < args.min_control_samples:
        sys.stdout.write(
            f"Sample sheet requires at least {args.min_control_samples} "
            f"'{args.control_group}' group samples"
        )
        raise SystemExit

    try:
        groups = {sample.group for sample in samples}
        if args.gene_contrast_file:
            read_gene_contrasts(args.gene_contrast_file, groups)
        if args.ica_analysis:
            resolve_ica_groups(groups, args.ica_control_group, args.ica_condition_groups)
            if sum(sample.group == args.ica_control_group.strip() for sample in samples) < 2:
                raise ValueError("ICA requires at least two samples in its control group")
    except (ValueError, OSError, csv.Error) as error:
        sys.stdout.write(f"Parsing error: {error}")
        raise SystemExit(1) from error

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
    parser.add_argument("--gene_contrast_file", help="Headerless CSV of target,reference gene contrasts")
    parser.add_argument("--differential_expression", action=argparse.BooleanOptionalAction, default=True,
                        help="Validate DEA controls; disable for an ICA-only cohort.")
    parser.add_argument("--ica_analysis", action="store_true")
    parser.add_argument("--ica_control_group", default="control")
    parser.add_argument("--ica_condition_groups")
    return parser
