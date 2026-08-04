import pandas as pd
import numpy as np
import json

from ezcharts.plots import BokehPlot
from ezcharts.layout.snippets import Tabs
from ezcharts.components.ezchart import EZChart
from scipy.stats import gaussian_kde  # type: ignore
from bokeh.palettes import Blues9
from bokeh.embed import json_item

from .result_types import SampleQCResult

MAX_KDE_POINTS = 20_000
KDE_SAMPLE_SEED = 0
KDE_GRID_SIZE = 100


def _empty_kde_plot(title: str, x_axis_label: str, y_axis_label: str) -> BokehPlot:
    """Create an empty 2D KDE plot for sparse live QC batches."""
    plot = BokehPlot(
        title=f"{title} (insufficient data)",
        x_axis_label=x_axis_label,
        y_axis_label=y_axis_label,
        sizing_mode="stretch_width",
    )
    plot._fig.background_fill_color = "#fafafa"
    plot._fig.grid.grid_line_alpha = 0.25
    plot._fig.line([], [], alpha=0)
    return plot


def _sample_kde_values(
    values: np.ndarray,
    max_points: int = MAX_KDE_POINTS,
) -> np.ndarray:
    """Deterministically cap the observations used to fit a KDE."""
    if len(values) <= max_points:
        return values

    rng = np.random.default_rng(KDE_SAMPLE_SEED)
    indices = rng.choice(len(values), size=max_points, replace=False)
    return values[indices]


def _kde(x_values: np.ndarray, y_values: np.ndarray, grid_size: int):
    """Evaluate a bivariate Gaussian KDE over a regular grid."""
    x_min, x_max = float(x_values.min()), float(x_values.max())
    y_min, y_max = float(y_values.min()), float(y_values.max())

    x_padding = max((x_max - x_min) * 0.05, 1.0)
    y_padding = max((y_max - y_min) * 0.05, 0.5)
    x_min -= x_padding
    x_max += x_padding
    y_min -= y_padding
    y_max += y_padding

    x_grid, y_grid = np.mgrid[
        x_min : x_max : complex(grid_size),
        y_min : y_max : complex(grid_size),
    ]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    values = np.vstack([x_values, y_values])
    kernel = gaussian_kde(values)
    density = np.reshape(kernel(positions).T, x_grid.shape)

    return x_grid, y_grid, density


def create_2d_kde_plot(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
    grid_size: int = KDE_GRID_SIZE,
) -> BokehPlot:
    """Create a Bokeh contour-style 2D KDE plot from two numeric columns."""
    x_axis_label = x_axis_label or x_column
    y_axis_label = y_axis_label or y_column

    if x_column not in data.columns or y_column not in data.columns:
        return _empty_kde_plot(title, x_axis_label, y_axis_label)

    plot_data = data[[x_column, y_column]].copy()
    plot_data[x_column] = pd.to_numeric(plot_data[x_column], errors="coerce")
    plot_data[y_column] = pd.to_numeric(plot_data[y_column], errors="coerce")
    plot_data = plot_data.dropna(subset=[x_column, y_column])

    values = plot_data[[x_column, y_column]].to_numpy(dtype=float)
    values = _sample_kde_values(values)
    if len(values) <= 2 or np.linalg.matrix_rank(values - values.mean(axis=0)) < 2:
        return _empty_kde_plot(title, x_axis_label, y_axis_label)

    x_values = values[:, 0]
    y_values = values[:, 1]

    try:
        x_grid, y_grid, density = _kde(x_values, y_values, grid_size)
    except (ValueError, np.linalg.LinAlgError):
        return _empty_kde_plot(title, x_axis_label, y_axis_label)

    density_min = float(np.nanmin(density))
    density_max = float(np.nanmax(density))
    if density_min == density_max:
        return _empty_kde_plot(title, x_axis_label, y_axis_label)

    plot = BokehPlot(
        title=title,
        x_axis_label=x_axis_label,
        y_axis_label=y_axis_label,
        sizing_mode="stretch_width",
    )
    plot._fig.background_fill_color = "#fafafa"
    plot._fig.grid.level = "overlay"
    plot._fig.grid.grid_line_color = "black"
    plot._fig.grid.grid_line_alpha = 0.05

    palette = Blues9[::-1]
    levels = np.linspace(density_min, density_max, len(palette) + 1)
    plot._fig.contour(
        x_grid,
        y_grid,
        density,
        levels[1:],
        fill_color=palette,
        line_color=palette,
    )
    return plot

def create_2d_kde_html(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
) -> str:
    """Return generic 2D KDE plot JSON for tests/probes."""
    plot = create_2d_kde_plot(
        data,
        x_column=x_column,
        y_column=y_column,
        title=title,
        x_axis_label=x_axis_label,
        y_axis_label=y_axis_label,
    )
    return json.dumps(json_item(plot._fig))


def add_sample_2d_kdes(
    sample_results: list[SampleQCResult],
    x_column: str,
    y_column: str,
    title: str,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
    height: str = "360px",
) -> None:
    """Add per-sample 2D KDE plots to the current report section."""
    tabs = Tabs()
    with tabs.add_dropdown_menu("Sample", change_header=True): # type: ignore
        for sample in sample_results:
            with tabs.add_dropdown_tab(sample.label):  # type: ignore
                EZChart(
                    create_2d_kde_plot(
                        sample.nanoplot,
                        x_column=x_column,
                        y_column=y_column,
                        title=title,
                        x_axis_label=x_axis_label,
                        y_axis_label=y_axis_label,
                    ),
                    "epi2melabs",
                    height=height,
                )
