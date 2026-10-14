"""A dummy test."""

import argparse

from workflow_glue import qc_report


def test():
    """Just showing that we can import using the workflow-glue."""
    assert isinstance(qc_report.argparser(), argparse.ArgumentParser)
