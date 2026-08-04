import ezcharts as ezc
import pandas as pd

from ezcharts.layout.snippets import Tabs
from ezcharts.components.ezchart import EZChart

from .result_types import SampleQCResult


def create_hist_plot(
    data: pd.DataFrame,
    x_column: str,
    title: str,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
    bins: int = 100,
):
    """Create a one-dimensional histogram from one numeric column."""
    x_axis_label = x_axis_label or x_column
    y_axis_label = y_axis_label or "Count"

    if x_column in data.columns:
        plot_data = pd.to_numeric(data[x_column], errors="coerce").dropna()
    else:
        plot_data = pd.Series(dtype=float)

    plot = ezc.histplot(data=plot_data, bins=bins)
    plot.title = {"text": title}
    plot._fig.xaxis.axis_label = x_axis_label
    plot._fig.yaxis.axis_label = y_axis_label
    return plot


def add_sample_hists(
    sample_results: list[SampleQCResult],
    x_column: str,
    title: str,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
    height: str = "360px",
) -> None:
    """Add per-sample histograms to the current report section."""
    tabs = Tabs()
    with tabs.add_dropdown_menu("Sample", change_header=True):  # type: ignore
        for sample in sample_results:
            with tabs.add_dropdown_tab(sample.label):  # type: ignore
                EZChart(
                    create_hist_plot(
                        sample.nanoplot,
                        x_column=x_column,
                        title=title,
                        x_axis_label=x_axis_label,
                        y_axis_label=y_axis_label,
                    ),
                    "epi2melabs",
                    height=height,
                )
