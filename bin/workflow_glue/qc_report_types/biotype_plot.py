"""Transcript-biotype composition plots for the integrated QC report."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool, NumeralTickFormatter, Title
from ezcharts.components.ezchart import EZChart
from ezcharts.plots import BokehPlot

from .result_types import SampleQCResult
from ..transcript_biotypes import BIOTYPE_ORDER


REQUIRED_COLUMNS = ("name", "group", "biotype", "num_reads", "fraction")
BIOTYPE_COLORS = {
    "Protein-coding": "#0072B2",
    "rRNA": "#D55E00",
    "tRNA": "#009E73",
    "lncRNA": "#CC79A7",
    "Other ncRNA": "#56B4E9",
    "Pseudogene": "#E69F00",
    "Other": "#999999",
    "Unknown": "#4D4D4D",
}
CAPTION = "Unknown denotes targets without one unambiguous annotation biotype."


def load_transcript_biotypes(path: Path | str) -> pd.DataFrame:
    """Load and validate the fixed-category transcript composition table."""
    data = pd.read_csv(path, sep="\t", dtype={"name": str, "group": str})
    missing = [column for column in REQUIRED_COLUMNS if column not in data]
    if missing:
        raise ValueError(
            "Transcript-biotype table is missing columns: " + ", ".join(missing)
        )
    data = data[list(REQUIRED_COLUMNS)].copy()
    if data.empty:
        raise ValueError("Transcript-biotype table contains no rows.")
    if data[["name", "group", "biotype"]].isna().any().any():
        raise ValueError("Transcript-biotype table contains empty identifiers.")
    unsupported = sorted(set(data["biotype"]).difference(BIOTYPE_ORDER))
    if unsupported:
        raise ValueError(
            "Transcript-biotype table contains unsupported classes: "
            + ", ".join(unsupported)
        )
    if data.duplicated(["group", "name", "biotype"]).any():
        raise ValueError("Transcript-biotype table contains duplicate sample classes.")

    for column in ("num_reads", "fraction"):
        data[column] = pd.to_numeric(data[column], errors="coerce")
        if not np.isfinite(data[column]).all():
            raise ValueError(
                f"Transcript-biotype table contains non-finite {column} values."
            )
    if (data["num_reads"] < 0).any():
        raise ValueError("Transcript-biotype table contains negative num_reads values.")
    if ((data["fraction"] < 0) | (data["fraction"] > 1)).any():
        raise ValueError("Transcript-biotype fractions must be between 0 and 1.")

    for (group, name), sample in data.groupby(["group", "name"], sort=False):
        observed = set(sample["biotype"])
        if observed != set(BIOTYPE_ORDER):
            raise ValueError(
                f"Transcript-biotype table has incomplete classes for {group}/{name}."
            )
        total_reads = float(sample["num_reads"].sum())
        total_fraction = float(sample["fraction"].sum())
        expected = 1.0 if total_reads > 0 else 0.0
        if not np.isclose(total_fraction, expected, atol=1e-8):
            raise ValueError(
                f"Transcript-biotype fractions for {group}/{name} sum to "
                f"{total_fraction}, expected {expected}."
            )
    data["sample"] = data["group"] + "/" + data["name"]
    return data


def create_transcript_biotype_plot(
    data: pd.DataFrame,
    sample_order: list[str],
) -> BokehPlot:
    """Create a 100% stacked horizontal bar chart across samples."""
    observed_samples = set(data["sample"])
    if observed_samples != set(sample_order):
        missing = sorted(set(sample_order).difference(observed_samples))
        unexpected = sorted(observed_samples.difference(sample_order))
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if unexpected:
            details.append("unexpected " + ", ".join(unexpected))
        raise ValueError(
            "Transcript-biotype samples do not match QC samples: " + "; ".join(details)
        )

    plot_height = max(380, 42 * len(sample_order) + 190)
    plot = BokehPlot(
        title="Transcript biotype composition",
        x_axis_label="Fraction of Oarfish-assigned abundance",
        y_axis_label="Sample",
        x_range=(0, 1),
        y_range=list(reversed(sample_order)),
        height=plot_height,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,save,reset",
    )
    cumulative = {sample: 0.0 for sample in sample_order}
    for biotype in BIOTYPE_ORDER:
        selected = data.loc[data["biotype"].eq(biotype)].set_index("sample")
        rows = []
        for sample in sample_order:
            row = selected.loc[sample]
            left = cumulative[sample]
            right = left + float(row["fraction"])
            rows.append(
                {
                    "sample": sample,
                    "group": row["group"],
                    "name": row["name"],
                    "biotype": biotype,
                    "num_reads": float(row["num_reads"]),
                    "fraction": float(row["fraction"]),
                    "left": left,
                    "right": right,
                }
            )
            cumulative[sample] = right
        source = ColumnDataSource(pd.DataFrame(rows))
        bars = plot._fig.hbar(
            y="sample",
            left="left",
            right="right",
            height=0.72,
            fill_color=BIOTYPE_COLORS[biotype],
            fill_alpha=0.9,
            line_color="white",
            line_width=0.5,
            legend_label=biotype,
            source=source,
        )
        plot._fig.add_tools(
            HoverTool(
                renderers=[bars],
                tooltips=[
                    ("Sample", "@sample"),
                    ("Biotype", "@biotype"),
                    ("Estimated reads", "@num_reads{0,0.000}"),
                    ("Fraction", "@fraction{0.00%}"),
                ],
            )
        )

    zero_samples = [
        sample
        for sample in sample_order
        if float(data.loc[data["sample"].eq(sample), "num_reads"].sum()) == 0
    ]
    if zero_samples:
        plot._fig.text(
            x=[0.5] * len(zero_samples),
            y=zero_samples,
            text=["No assigned reads"] * len(zero_samples),
            text_align="center",
            text_baseline="middle",
            text_color="#4D4D4D",
        )

    plot._fig.xaxis.formatter = NumeralTickFormatter(format="0%")
    plot._fig.grid.grid_line_alpha = 0.15
    plot._fig.legend.title = "Transcript biotype"
    plot._fig.legend.orientation = "horizontal"
    plot._fig.legend.location = "top_left"
    plot._fig.legend.click_policy = "hide"
    plot._fig.add_layout(
        Title(
            text=CAPTION,
            text_font_size="10pt",
            text_color="#4B5563",
            align="left",
        ),
        "below",
    )
    plot.report_height = plot_height
    return plot


def add_transcript_biotype_composition(
    path: Path | str,
    sample_results: list[SampleQCResult],
) -> None:
    """Add cross-sample transcript composition to the current report panel."""
    data = load_transcript_biotypes(path)
    sample_order = [sample.label for sample in sample_results]
    EZChart(
        create_transcript_biotype_plot(data, sample_order),
        "epi2melabs",
        height=f"{max(380, 42 * len(sample_order) + 190)}px",
    )
