from bokeh import model

from bokeh.models import ColumnDataSource, ColorBar, LinearColorMapper
from bokeh.transform import factor_cmap, transform
from bokeh.plotting import figure, gridplot
from bokeh.models import HoverTool
from bokeh.palettes import Turbo256
import pandas as pd

import colorcet as cc


def create_de_heatmap_plot(
    title: str,
    norm_counts_df: pd.DataFrame,
    metadata: pd.DataFrame,
) -> model.Model:
    heatmap_df = z_score(norm_counts_df)

    gene_ids = heatmap_df.index.tolist()
    sample_ids = heatmap_df.columns.tolist()

    # Oarfish names its identifier column ``tname``, while featureCounts uses
    # ``Geneid``. Normalize the index name here so the plot supports both.
    heatmap_df = pd.melt(
        heatmap_df.rename_axis("feature_id").reset_index(),
        id_vars=["feature_id"],
        var_name="sample",
        value_name="z_score",
    )

    TOOLS = "save,ypan,reset,ywheel_zoom"
    hover = HoverTool()
    hover.tooltips = [
        ("Feature ID", "@feature_id"),
        ("Sample", "@sample"),
        ("z-score", "@z_score"),
    ]
    heatmap_plot = figure(
        tools=TOOLS,
        tooltips=hover.tooltips,
        title=f"Differential gene expression {title}",
        y_axis_label="Feature ID",
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
        low=heatmap_df.z_score.min(),
        high=heatmap_df.z_score.max(),
    )
    group_names = tuple(metadata["group"].unique())
    group_palette = [
        Turbo256[index * 255 // max(len(group_names) - 1, 1)]
        for index in range(len(group_names))
    ]
    # Create a data source with your samples and their types
    source = ColumnDataSource(
        data=dict(
            samples=sample_ids,
            groups=[metadata.loc[sample_id, "group"] for sample_id in sample_ids],
        )
    )
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
        fill_color=factor_cmap("groups", palette=group_palette, factors=group_names),
    )

    annotation_bar.xaxis.major_tick_line_color = None
    annotation_bar.xaxis.minor_tick_line_color = None
    annotation_bar.yaxis.major_tick_line_color = None
    annotation_bar.yaxis.minor_tick_line_color = None
    annotation_bar.yaxis.visible = True
    annotation_bar.yaxis.major_label_text_font_size = "0px"
    annotation_bar.yaxis.major_label_text_alpha = 0.0
    annotation_bar.grid.grid_line_color = None

    fill_color = transform("z_score", mapper)  # type: ignore
    heatmap_plot.rect(
        x="sample",
        y="feature_id",
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


def z_score(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize each feature across samples for heatmap display."""
    feature_means = df.mean(axis=1)
    feature_std = df.std(axis=1).replace(0, float("nan"))
    return df.sub(feature_means, axis=0).div(feature_std, axis=0).fillna(0.0)
