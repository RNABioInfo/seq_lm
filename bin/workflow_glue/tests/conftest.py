#!/usr/bin/env python
"""Pytest argument definitions."""

from pathlib import Path


def pytest_addoption(parser):
    """Define command line arguments for pytest."""
    parser.addoption(
        "--test_data",
        action="store",
        default=str(Path(__file__).resolve().parents[3] / "test_data")
    )
