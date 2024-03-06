from bokeh.plotting import figure
from bokeh.model import Model

from deseq_analysis import DeseqStats


def create_ma_plot(title: str, stat_results: DeseqStats) -> Model:
    colors = ["red" if x < 0.05 else "black" for x in stat_results.results_df["padj"]]

    ma_plot = figure(
        title=f"MA Plot {title}",
        x_axis_label="log2 Fold Change",
        y_axis_label="log2 Counts",
        sizing_mode="stretch_width",
    )
    ma_plot.circle(
        stat_results.results_df["log2FoldChange"],
        stat_results.results_df["baseMean"],
        line_color=None,
        fill_color=colors,
        fill_alpha=0.6,
        size=6,
    )

    return ma_plot
