import pandas as pd
import sys

from pathlib import Path
from functools import reduce
import re

from .result_types import FlagstatResult, QCResult, SampleQCResult
from .constants import (
    DISPLAY_COLUMNS,
    EXPECTED_COLUMNS,
    NANOPLOT_COLUMNS,
    OPTIONAL_COLUMNS,
)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", value)


def retrieve_flagstat_paths(qc_dir: str, group: str, name: str) -> list[Path]:
    return list(
        Path(qc_dir, "flagstat", safe_name(group), safe_name(name)).glob("*.tsv")
    )


def retrieve_nanoplot_paths(qc_dir: str, group: str, name: str) -> list[Path]:
    return list(
        Path(qc_dir, "nanoplot", safe_name(group), safe_name(name)).glob("*.tsv.gz")
    )


def aggregate_flagstat(flagstat_paths: list[Path]) -> FlagstatResult | None:
    flagstat_results = [FlagstatResult.from_tsv(path) for path in flagstat_paths]
    if not flagstat_results:
        return None

    return reduce(lambda prev, cur: prev + cur, flagstat_results)  # type: ignore


def aggregate_nanoplot(nanoplot_paths: list[Path]) -> pd.DataFrame | None:
    if not nanoplot_paths:
        return None

    nanoplot_results = []
    for path in nanoplot_paths:
        try:
            nanoplot_results.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            # QC versions before the producer wrote an explicit header created
            # a valid, headerless gzip stream when a BAM had no usable reads.
            # Preserve that chunk as an empty table so cumulative live reports
            # can still include its flagstat counts and later QC chunks.
            nanoplot_results.append(pd.DataFrame(columns=NANOPLOT_COLUMNS))
    missing_columns = sorted(
        {
            column
            for result in nanoplot_results
            for column in NANOPLOT_COLUMNS
            if column not in result
        }
    )
    if missing_columns:
        raise ValueError(
            "NanoPlot data table is missing columns: " + ", ".join(missing_columns)
        )

    return pd.concat(nanoplot_results, ignore_index=True)


def load_qc_samples(samples_path) -> QCResult:
    """Load the QC report sample table."""
    samples = pd.read_csv(samples_path, sep="\t", dtype=str).fillna("")
    missing_columns = [column for column in EXPECTED_COLUMNS if column not in samples]
    if missing_columns:
        raise ValueError(
            "QC report samples table is missing columns: " + ", ".join(missing_columns)
        )
    columns = EXPECTED_COLUMNS + [
        column for column in OPTIONAL_COLUMNS if column in samples
    ]
    samples = samples[columns].rename(columns=DISPLAY_COLUMNS)

    qc_results: list[SampleQCResult] = []

    for _, row in samples.iterrows():
        name = row["Name"]
        group = row["Group"]
        qc_dir = row["QC input directory"]

        flagstat_paths = retrieve_flagstat_paths(qc_dir, group, name)
        nanoplot_paths = retrieve_nanoplot_paths(qc_dir, group, name)
        if flagstat_res := aggregate_flagstat(flagstat_paths):
            nanoplot_res = aggregate_nanoplot(nanoplot_paths)
            if nanoplot_res is None:
                sys.exit(f"Could not aggregate NanoPlot results for: {group}/{name}")

            qc_results.append(
                SampleQCResult(
                    name,
                    group,
                    flagstat_res,
                    nanoplot_res,
                )
            )
        else:
            sys.exit(f"Could not aggregate flagstat results for: {group}/{name}")

    return QCResult(qc_results, samples)
