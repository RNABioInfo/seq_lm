from bokeh import model

from bokeh.models import ColumnDataSource, Whisker
from bokeh.transform import factor_cmap
from bokeh.plotting import figure
from bokeh.palettes import Pastel1

import pandas as pd


def create_library_size_plot(
    counts_df: pd.DataFrame, title: str, yAxisLabel: str
) -> model.Model:
    df = counts_df.T.melt(var_name="columns")
    df = df[["columns", "value"]].rename(columns={"columns": "kind"})

    kinds = df.kind.unique()

    # compute quantiles
    qs = df.groupby("kind").value.quantile([0.25, 0.5, 0.75])
    qs = qs.unstack().reset_index()
    qs.columns = ["kind", "q1", "q2", "q3"]
    df = pd.merge(df, qs, on="kind", how="left")

    # compute IQR outlier bounds
    iqr = df.q3 - df.q1
    df["upper"] = df.q3 + 1.5 * iqr
    df["lower"] = df.q1 - 1.5 * iqr

    source = ColumnDataSource(df)

    p = figure(
        x_range=kinds,
        toolbar_location="right",
        title=title,
        background_fill_color="#eaefef",
        y_axis_label=yAxisLabel,
        x_axis_label="Sample [Run_Replicate]",
    )

    # outlier range
    whisker = Whisker(base="kind", upper="upper", lower="lower", source=source)
    whisker.upper_head.size = whisker.lower_head.size = 20
    p.add_layout(whisker)

    # quantile boxes
    cmap = factor_cmap("kind", Pastel1[len(kinds)], kinds)  # type: ignore
    p.vbar("kind", 0.7, "q2", "q3", source=source, color=cmap, line_color="black")
    p.vbar("kind", 0.7, "q1", "q2", source=source, color=cmap, line_color="black")

    # outliers
    outliers = df[~df.value.between(df.lower, df.upper)]
    p.scatter("kind", "value", source=outliers, size=6, color="black", alpha=0.3)

    p.xgrid.grid_line_color = None
    p.axis.major_label_text_font_size = "14px"
    p.axis.axis_label_text_font_size = "12px"

    return p
