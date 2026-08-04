import json

from ezcharts.plots import Plot
from ezcharts.layout.snippets import Tabs
from ezcharts.components.ezchart import EZChart

from .result_types import FlagstatResult, SampleQCResult

class ReadFateSankeyPlot(Plot):
    """ECharts Sankey plot without axis finalization."""

    def finalise(self):
        """Skip axis defaults because Sankey charts do not define axes."""
        return None


def create_read_fate_sankey_plot(flagstat_res: FlagstatResult) -> ReadFateSankeyPlot:
    """Create a native ECharts Sankey showing read alignment fate."""
    primary_unmapped = max(flagstat_res.primary_reads - flagstat_res.primary_mapped, 0)

    nodes = [
        {"name": "Total entries", "itemStyle": {"color": "#334155"}},
        {"name": "Primary entries", "itemStyle": {"color": "#2563eb"}},
        {"name": "Secondary entries", "itemStyle": {"color": "#878787"}},
        {"name": "Primary mapped", "itemStyle": {"color": "#16a34a"}},
        {"name": "Primary unmapped", "itemStyle": {"color": "#dc2626"}},
    ]
    links = [
        ("Total entries", "Primary entries", flagstat_res.primary_reads),
        ("Total entries", "Secondary entries", flagstat_res.secondary_reads),
        ("Primary entries", "Primary mapped", flagstat_res.primary_mapped),
        ("Primary entries", "Primary unmapped", primary_unmapped),
    ]
    plot = ReadFateSankeyPlot(
        tooltip={"trigger": "item", "triggerOn": "mousemove"},
        series=[
            {
                "type": "sankey",
                "data": nodes,
                "links": [
                    {"source": source, "target": target, "value": value}
                    for source, target, value in links
                    if value > 0
                ],
                "nodeWidth": 18,
                "nodeGap": 18,
                "draggable": False,
                "label": {"fontSize": 12},
                "lineStyle": {
                    "color": "source",
                    "opacity": 0.28,
                    "curveness": 0.5,
                },
                "emphasis": {"focus": "adjacency"},
            }
        ],
    )
    plot.toolbox = {"show": False}
    return plot


def add_sample_read_fate_sankeys(sample_results: list[SampleQCResult]) -> None:
    """Add one native ECharts Sankey per sample behind an ezCharts dropdown."""
    tabs = Tabs()
    with tabs.add_dropdown_menu("Sample", change_header=True): # type: ignore
        for sample in sample_results:
            with tabs.add_dropdown_tab(sample.label):  # type: ignore
                EZChart(
                    create_read_fate_sankey_plot(sample.flagstat),
                    "epi2melabs",
                    height="300px",
                )


def create_read_fate_sankey_html(flagstat_res: FlagstatResult) -> str:
    """Return the native ECharts Sankey options as JSON for tests/probes."""
    return create_read_fate_sankey_plot(flagstat_res).to_json()


def create_sample_read_fate_sankey_html(sample_results: list[SampleQCResult]) -> str:
    """Return sample labels and Sankey options as JSON for tests/probes."""
    sample_payload = []
    for sample in sample_results:
        sample_payload.append(
            {
                "sample": sample.label,
                "plot": create_read_fate_sankey_plot(sample.flagstat).dict(
                    exclude_unset=True
                ),
            }
        )
    return json.dumps(sample_payload, separators=(",", ":"))
