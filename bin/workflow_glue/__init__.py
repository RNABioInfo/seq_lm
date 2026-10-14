"""Workflow Python code."""
import argparse
import importlib

from .util import _log_level, get_main_logger  # noqa: ABS101


__version__ = "0.0.1"
_package_name = "workflow_glue"


def get_components():
    """List the supported workflow commands."""
    return ["check_sample_sheet", "qc_report", "transcript_biotypes"]


def cli():
    """Run workflow entry points."""
    parser = argparse.ArgumentParser(
        "wf-glue",
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(__version__)
    )

    subparsers = parser.add_subparsers(
        title="subcommands",
        description="valid commands",
        help="additional help",
        dest="command",
    )
    subparsers.required = True

    # Register the supported commands.
    components = [f"{_package_name}.{comp}" for comp in get_components()]
    for module in components:
        mod = importlib.import_module(module)
        p = subparsers.add_parser(module.split(".")[-1], parents=[mod.argparser()])
        p.set_defaults(func=mod.main)

    logger = get_main_logger(_package_name)
    args = parser.parse_args()

    logger.info("Starting entrypoint.")
    args.func(args)
