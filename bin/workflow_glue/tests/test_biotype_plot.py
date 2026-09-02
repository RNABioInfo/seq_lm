"""Tests for the transcript-biotype QC report plot."""

from pathlib import Path

import pytest

pytest.importorskip("ezcharts")
pd = pytest.importorskip("pandas")

from bokeh.models import ColumnDataSource, HoverTool, Title  # noqa: E402

from workflow_glue.qc_report_types import biotype_plot as bp  # noqa: E402
from workflow_glue.transcript_biotypes import BIOTYPE_ORDER  # noqa: E402


def composition_rows(samples=("control/A", "treated/B")):
    rows = []
    for sample_index, sample in enumerate(samples):
        group, name = sample.split("/", 1)
        for index, biotype in enumerate(BIOTYPE_ORDER):
            fraction = 0.0
            reads = 0.0
            if index == sample_index:
                fraction = 0.75
                reads = 75.5
            elif biotype == "Unknown":
                fraction = 0.25
                reads = 25.0
            rows.append(
                {
                    "name": name,
                    "group": group,
                    "biotype": biotype,
                    "num_reads": reads,
                    "fraction": fraction,
                }
            )
    return pd.DataFrame(rows)


def test_load_transcript_biotypes_validates_and_labels_samples(tmp_path):
    path = Path(tmp_path, "biotypes.tsv")
    composition_rows().to_csv(path, sep="\t", index=False)

    loaded = bp.load_transcript_biotypes(path)

    assert loaded["sample"].drop_duplicates().tolist() == ["control/A", "treated/B"]
    assert loaded.groupby("sample")["fraction"].sum().tolist() == [1.0, 1.0]


def test_plot_has_fixed_stack_order_hover_caption_and_sample_order():
    data = composition_rows()
    data["sample"] = data["group"] + "/" + data["name"]

    plot = bp.create_transcript_biotype_plot(data, ["treated/B", "control/A"])

    assert plot._fig.y_range.factors == ["control/A", "treated/B"]
    sources = [
        renderer.data_source
        for renderer in plot._fig.renderers
        if hasattr(renderer, "data_source")
        and "biotype" in renderer.data_source.data
    ]
    assert [source.data["biotype"][0] for source in sources] == list(BIOTYPE_ORDER)
    colors = [
        renderer.glyph.fill_color
        for renderer in plot._fig.renderers
        if hasattr(renderer, "data_source")
        and "biotype" in renderer.data_source.data
    ]
    assert colors == [bp.BIOTYPE_COLORS[biotype] for biotype in BIOTYPE_ORDER]
    hover_tools = list(plot._fig.select({"type": HoverTool}))
    assert len(hover_tools) == len(BIOTYPE_ORDER)
    tooltip_labels = {label for label, _value in hover_tools[0].tooltips}
    assert tooltip_labels == {"Sample", "Biotype", "Estimated reads", "Fraction"}
    titles = [title.text for title in plot._fig.select({"type": Title})]
    assert bp.CAPTION in titles
    assert plot._fig.xaxis[0].formatter.format == "0%"


def test_zero_read_sample_is_annotated():
    data = composition_rows(samples=("control/zero",))
    data["num_reads"] = 0.0
    data["fraction"] = 0.0
    data["sample"] = data["group"] + "/" + data["name"]

    plot = bp.create_transcript_biotype_plot(data, ["control/zero"])

    texts = [
        text
        for source in plot._fig.select({"type": ColumnDataSource})
        for text in source.data.get("text", [])
    ]
    assert "No assigned reads" in texts


def test_incomplete_categories_are_rejected(tmp_path):
    path = Path(tmp_path, "biotypes.tsv")
    composition_rows().iloc[:-1].to_csv(path, sep="\t", index=False)

    with pytest.raises(ValueError, match="incomplete classes"):
        bp.load_transcript_biotypes(path)
