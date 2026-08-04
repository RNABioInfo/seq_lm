from bokeh.plotting import figure
from bokeh.model import Model
from bokeh.palettes import Turbo256
from bokeh.transform import factor_cmap

from sklearn.decomposition import PCA

import pandas as pd

import numpy as np


def create_pca_plot(title: str, pca_df: pd.DataFrame, metadata: pd.DataFrame) -> Model:
    log2_norm_data = np.log2(pca_df.T + 1)

    pca = PCA(n_components=2)
    pca_df = pd.DataFrame(pca.fit_transform(log2_norm_data), columns=["PC1", "PC2"])
    group_names: tuple[str, ...] = tuple(metadata["group"].unique())
    group_palette = [
        Turbo256[index * 255 // max(len(group_names) - 1, 1)]
        for index in range(len(group_names))
    ]
    pca_df["group"] = metadata["group"].tolist()
    pca_df["sample"] = metadata.index.tolist()

    pca_plot = figure(
        title=f"PCA Plot {title}",
        x_axis_label="PC1",
        y_axis_label="PC2",
        sizing_mode="stretch_width",
    )
    pca_plot.scatter(
        "PC1",
        "PC2",
        source=pca_df,
        line_color="black",
        line_alpha=0.3,
        fill_color=factor_cmap("group", group_palette, group_names),
        fill_alpha=0.8,
        legend_field="group",
        size=18,
    )

    return pca_plot
