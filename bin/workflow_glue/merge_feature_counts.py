#!/usr/bin/env python
"""Merge feature counts from new bam files into a single table."""

import pandas as pd
from .util import get_named_logger, wf_parser


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("Merge feature counts.")
    parser.add_argument("-n", "--new_counts", type=str)
    parser.add_argument("-a", "--all_counts", type=str)
    parser.add_argument("-o", "--output", type=str)

    return parser


def main(args):
    newCountsFile = args.new_counts
    allCountsFile = args.all_counts
    outputCountsFile = args.output

    new_counts = pd.read_csv(newCountsFile, sep="\t", index_col=0, header=0)
    new_counts.rename(columns={new_counts.columns[-1]: "counts"}, inplace=True)

    all_counts = pd.read_csv(allCountsFile, sep="\t", index_col=0, header=0)
    all_counts.rename(columns={all_counts.columns[-1]: "counts"}, inplace=True)

    all_counts.counts = all_counts.counts + new_counts.counts

    all_counts.to_csv(outputCountsFile, sep="\t", index=True)
