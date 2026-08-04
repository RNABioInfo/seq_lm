"""Check if a sample sheet is valid."""
import csv
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def parse_is_live(value):
    """Parse the optional per-sample live-analysis flag."""
    normalized = (value or "").strip().lower()
    if normalized in {"", "true"}:
        return True
    if normalized == "false":
        return False
    raise ValueError(
        f"Invalid is_live value '{value}'. Expected true, false, or a blank value."
    )


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheet")

    group_names = []
    control_count = 0

    try:
        with open(args.sample_sheet, "r") as f:
            csv_reader = csv.DictReader(f)
            removed_fields = {"id", "type"}.intersection(csv_reader.fieldnames or [])
            if removed_fields:
                sys.stdout.write(
                    "Sample sheet contains removed fields: "
                    + ", ".join(sorted(removed_fields))
                )
                sys.exit()
            n_row = 0
            for row in csv_reader:
                n_row += 1
                if n_row == 1:
                    n_cols = len(row)
                else:
                    # check we got the same number of fields
                    if len(row) != n_cols:
                        raise ValueError(
                            f"Unexpected number of cells in row number {n_row}."
                        )
                try:
                    alias = row["alias"].strip()
                except KeyError:
                    sys.stdout.write("'alias' column missing")
                    sys.exit()
                try:
                    group = row["group"].strip()
                except KeyError:
                    sys.stdout.write("'group' column missing")
                    sys.exit()
                if not alias:
                    sys.stdout.write("empty value in 'alias' column")
                    sys.exit()
                if not group:
                    sys.stdout.write("empty value in 'group' column")
                    sys.exit()
                parse_is_live(row.get("is_live"))
                group_names.append((group, alias))
                if group.lower() == args.control_group.lower():
                    control_count += 1
    except Exception as e:
        sys.stdout.write(f"Parsing error: {e}")
        sys.exit()

    if len(group_names) > len(set(group_names)):
        sys.stdout.write("values in 'alias' column not unique within group")
        sys.exit()

    if control_count < args.min_control_samples:
        sys.stdout.write(
            f"Sample sheet requires at least {args.min_control_samples} "
            f"'{args.control_group}' group samples"
        )
        sys.exit()

    logger.info(f"Checked sample sheet {args.sample_sheet}.")


def argparser():
    """Argument parser for entrypoint."""
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
