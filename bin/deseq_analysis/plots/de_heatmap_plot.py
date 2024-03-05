from bokeh import model

from bokeh.models import ColumnDataSource, ColorBar, LinearColorMapper
from bokeh.transform import factor_cmap, transform
from bokeh.plotting import figure, gridplot
from bokeh.models import HoverTool
import pandas as pd

import colorcet as cc


def create_de_heatmap_plot(title: str, norm_counts_df: pd.DataFrame) -> model.Model:
    heatmap_df = z_score(norm_counts_df)

    column_names: list[str] = []
    for col in heatmap_df.columns:
        elements = col.split("_")
        name = f"run_{elements[0]}_replicate_{elements[1]}"
        column_names.append(name)

    heatmap_df.columns = column_names

    gene_ids = heatmap_df.index.tolist()
    sample_ids = heatmap_df.columns.tolist()

    heatmap_df = pd.melt(heatmap_df.reset_index(), id_vars=["Geneid"])

    TOOLS = "save,ypan,reset,ywheel_zoom"
    hover = HoverTool()
    hover.tooltips = [
        ("Gene ID", "@Geneid"),
        ("Sample", "@variable"),
        ("z-score", "@value"),
    ]
    heatmap_plot = figure(
        tools=TOOLS,
        tooltips=hover.tooltips,
        title=f"Differential gene expression {name}",
        y_axis_label="Gene ID",
        sizing_mode="stretch_width",
        x_range=sample_ids,
        y_range=gene_ids,
        x_axis_location="above",
        width=400,
        min_border_bottom=5,
    )
    heatmap_plot.grid.grid_line_color = None
    heatmap_plot.axis.axis_line_color = None
    heatmap_plot.axis.major_tick_line_color = None
    heatmap_plot.axis.major_label_text_font_size = "7px"
    heatmap_plot.axis.major_label_standoff = 0
    heatmap_plot.xaxis.visible = False

    mapper = LinearColorMapper(
        palette=cc.b_diverging_bwr_20_95_c54,
        low=heatmap_df.value.min(),
        high=heatmap_df.value.max(),
    )
    colors = ["#1f77b4", "#1f77b4", "#2ca02c", "#2ca02c"]
    # Create a data source with your samples and their types
    source = ColumnDataSource(data=dict(samples=sample_ids))
    # Create a new figure for the annotation bar
    annotation_bar: model.Model = figure(
        tools="",
        x_axis_label="Samples",
        y_axis_label="Run",
        y_axis_location="left",
        x_range=heatmap_plot.x_range,
        y_range=(0, 1),
        frame_height=15,
        outline_line_color="white",
        background_fill_color="white",
    )
    # Add rectangles to the figure, colored based on the sample type
    annotation_bar.rect(
        source=source,
        x="samples",
        y=0.5,
        width=0.98,
        height=20,
        height_units="screen",
        fill_color=factor_cmap("samples", palette=colors, factors=sample_ids),
    )

    annotation_bar.xaxis.major_tick_line_color = None
    annotation_bar.xaxis.minor_tick_line_color = None
    annotation_bar.yaxis.major_tick_line_color = None
    annotation_bar.yaxis.minor_tick_line_color = None
    annotation_bar.yaxis.visible = True
    annotation_bar.yaxis.major_label_text_font_size = "0px"
    annotation_bar.yaxis.major_label_text_alpha = 0.0
    annotation_bar.yaxis.axis_label_orientation = "horizontal"
    annotation_bar.grid.grid_line_color = None

    fill_color = transform("value", mapper)  # type: ignore
    heatmap_plot.rect(
        x="variable",
        y="Geneid",
        width=0.98,
        height=1,
        source=heatmap_df,
        fill_color=fill_color,
        line_color=None,
    )
    color_bar = ColorBar(
        color_mapper=mapper,
        location=(0, 0),
        orientation="vertical",
        title="z-score",
    )
    heatmap_plot.add_layout(color_bar, "right")
    p = gridplot(
        [[heatmap_plot], [annotation_bar]],  # type: ignore
        merge_tools=True,
        toolbar_location="right",
    )

    return p


def z_score(df) -> pd.DataFrame:
    # copy the dataframe
    df_std = df.copy()
    # apply the z-score method
    for column in df_std.columns:
        df_std[column] = (df_std[column] - df_std[column].mean()) / df_std[column].std()

    return df_std
