import pandas as pd

from .constants import NANOPLOT_METRICS
from .result_types import SampleQCResult

def _format_integer_metric(value: float | int | None) -> str:
    if value is None or pd.isna(value):
        return ""
    return f"{value:.0f}"


def _format_float_metric(value: float | int | None) -> str:
    if value is None or pd.isna(value):
        return ""
    return f"{value:.2f}"


def create_nanoplot_metrics_table(sample_results: list[SampleQCResult]) -> pd.DataFrame:
    """Create a cross-sample NanoPlot metrics table."""
    rows = [{"Metric": metric} for metric in NANOPLOT_METRICS]

    for sample in sample_results:
        nanoplot = sample.nanoplot.copy()
        lengths = pd.to_numeric(nanoplot["lengths"], errors="coerce")
        aligned_lengths = pd.to_numeric(nanoplot["aligned_lengths"], errors="coerce")
        quals = pd.to_numeric(nanoplot["quals"], errors="coerce")
        mapq = pd.to_numeric(nanoplot["mapQ"], errors="coerce")
        metrics = {
            "Number of mapped reads": _format_integer_metric(len(nanoplot)),
            "Number of bases": _format_integer_metric(lengths.sum(min_count=1)),
            "Number of aligned bases": _format_integer_metric(
                aligned_lengths.sum(min_count=1)
            ),
            "Median read length": _format_float_metric(lengths.median()),
            "Mean read length": _format_float_metric(lengths.mean()),
            "Median read quals": _format_float_metric(quals.median()),
            "Mean read quals": _format_float_metric(quals.mean()),
            "Median MapQ": _format_float_metric(mapq.median()),
            "Mean MapQ": _format_float_metric(mapq.mean()),
        }
        for row in rows:
            row[sample.label] = metrics[row["Metric"]]

    return pd.DataFrame(rows)
