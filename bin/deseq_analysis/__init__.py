#!/usr/bin/env python
"""Perform differential expression analysis using DESeq2."""

import re
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from bokeh.io import output_file, save
from bokeh.plotting import figure, row
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from .plots.de_heatmap_plot import create_de_heatmap_plot
from .plots.library_size_plot import create_library_size_plot
from .plots.ma_plot import create_ma_plot
from .plots.pca_plot import create_pca_plot
from .util.metadata import get_metadata, read_quantification_counts

__version__ = "0.0.1"
_package_name = "deseq_analysis"


def argparser() -> ArgumentParser:
    """Create the command-line parser."""
    parser = ArgumentParser(description="Differential expression analysis.")
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-m",
        "--quant-manifest",
        help="TSV with name, group, and count_file columns.",
    )
    inputs.add_argument(
        "-q",
        "--quant-files",
        nargs="+",
        help="Legacy run_N_replicate_N featureCounts files.",
    )
    parser.add_argument("-o", "--output-dir", default=".")
    parser.add_argument("-t", "--threads", type=int, default=1)
    return parser


def _safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", value)


def _placeholder_plot(title: str, message: str):
    plot = figure(title=title, width=800, height=250, toolbar_location=None)
    plot.text(x=[0.5], y=[0.5], text=[message], text_align="center")
    plot.xaxis.visible = False
    plot.yaxis.visible = False
    plot.grid.visible = False
    return plot


def _load_count_matrix(metadata: pd.DataFrame) -> pd.DataFrame:
    count_series: list[pd.Series] = []
    for sample in metadata.itertuples(index=False):
        counts = read_quantification_counts(sample.count_file)
        count_series.append(counts.rename(sample.sample))

    counts_df = pd.concat(count_series, axis=1, join="outer").fillna(0)
    counts_df = counts_df.apply(pd.to_numeric, errors="coerce").fillna(0)
    if (counts_df < 0).any().any():
        raise ValueError("Quantification files contain negative counts.")

    # Oarfish emits EM-estimated read counts, which can be fractional. PyDESeq2
    # requires a non-negative integer count matrix.
    counts_df = counts_df.round().astype(int)
    counts_df = counts_df[counts_df.sum(axis=1) > 20]
    if counts_df.empty:
        raise ValueError("No quantified targets remain after the total-count > 20 filter.")

    return counts_df.T


def main() -> None:
    args = argparser().parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_df = get_metadata(args.quant_files, args.quant_manifest)
    if "control" not in set(metadata_df["group"]):
        raise ValueError("Differential expression requires a control group.")
    comparison_groups = sorted(set(metadata_df["group"]) - {"control"})
    if not comparison_groups:
        raise ValueError("Differential expression requires at least one non-control group.")

    counts_df = _load_count_matrix(metadata_df)
    metadata_df = metadata_df.set_index("sample")
    metadata_df["group"] = pd.Categorical(
        metadata_df["group"],
        categories=["control", *comparison_groups],
    )

    inference = DefaultInference(n_cpus=args.threads)
    data_set = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design="~group",
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    data_set.deseq2()

    output_file(str(output_dir / "library_sizes.html"), title="Library Sizes")
    non_norm_plot = create_library_size_plot(
        counts_df,
        title="Non-normalized library counts",
        y_axis_label="Counts",
    )
    norm_data = pd.DataFrame(
        data_set.layers["normed_counts"],
        index=counts_df.index,
        columns=counts_df.columns,
    )
    norm_plot = create_library_size_plot(
        norm_data,
        title="Normalized library counts",
        y_axis_label="Counts",
    )
    save(row(non_norm_plot, norm_plot, sizing_mode="stretch_width"))

    for comparison_group in comparison_groups:
        comparison_title = f"{comparison_group} vs control"
        comparison_dir = output_dir / (
            f"group_{_safe_label(comparison_group)}_vs_control"
        )
        comparison_dir.mkdir(parents=True, exist_ok=True)

        stat_res = DeseqStats(
            data_set,
            contrast=["group", comparison_group, "control"],
            inference=inference,
            quiet=True,
        )
        stat_res.summary()
        stat_res.lfc_shrink(f"group[T.{comparison_group}]")
        stat_res.results_df.to_csv(
            comparison_dir / "deseq2_results.tsv",
            sep="\t",
            index=True,
            header=True,
            float_format="%.10f",
        )

        output_file(str(comparison_dir / "ma_plot.html"))
        save(create_ma_plot(comparison_title, stat_res))

        comparison_metadata = metadata_df[
            metadata_df["group"].isin([comparison_group, "control"])
        ]
        comparison_counts = norm_data.loc[comparison_metadata.index].T

        output_file(str(comparison_dir / "pca_plot.html"))
        if min(comparison_counts.shape) >= 2:
            pca_plot = create_pca_plot(
                comparison_title,
                comparison_counts,
                comparison_metadata,
            )
        else:
            pca_plot = _placeholder_plot(
                f"PCA Plot {comparison_title}",
                "Insufficient samples or quantified targets for PCA.",
            )
        save(pca_plot)

        significant_ids = stat_res.results_df.index[
            stat_res.results_df["padj"].notna()
            & (stat_res.results_df["padj"] < 0.05)
        ]
        significant_counts = comparison_counts.loc[
            comparison_counts.index.intersection(significant_ids)
        ]
        output_file(str(comparison_dir / "heatmap.html"))
        if significant_counts.empty:
            heatmap_plot = _placeholder_plot(
                f"Differential gene expression {comparison_title}",
                "No targets pass adjusted p-value < 0.05 at this live checkpoint.",
            )
        else:
            heatmap_plot = create_de_heatmap_plot(
                comparison_title,
                significant_counts,
                comparison_metadata,
            )
        save(heatmap_plot)
