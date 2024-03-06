#!/usr/bin/env python
"""Perform differential expression analysis using DESeq2."""
import os
from argparse import ArgumentParser

import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from bokeh.plotting import row
from bokeh.io import save, output_file

from .plots.library_size_plot import create_library_size_plot
from .plots.de_heatmap_plot import create_de_heatmap_plot
from .plots.ma_plot import create_ma_plot
from .plots.pca_plot import create_pca_plot
from .util.metadata import get_metadata

__version__ = "0.0.1"
_package_name = "deseq_analysis"


def argparser():
    """Argument parser for entrypoint."""
    parser = ArgumentParser("Differential expression analysis.")
    parser.add_argument("-q", "--quant_files", type=str, nargs="+", default=[])
    parser.add_argument("-t", "--threads", type=int, default=1)

    return parser


def main():
    args = argparser().parse_args()

    quantFiles: list[str] = args.quant_files
    quantFiles.sort()
    threads: int = args.threads

    metadata_df = get_metadata(quantFiles)

    # Read in the count files and merge them into a single data frame by the "Name" column
    count_files = metadata_df["count_file"].tolist()
    sample_names = metadata_df["sample"].tolist()
    metadata_df.set_index("sample", inplace=True)

    dfs = []

    # Read and merge the count files
    for count_file, sample_name in zip(count_files, sample_names):
        df: pd.DataFrame = pd.read_csv(count_file, sep="\t", header=0, index_col=0)
        df = df.iloc[:, [-1]]
        df.columns = [sample_name]
        dfs.append(df)

    counts_df: pd.DataFrame = dfs[0].join(dfs[1:])

    counts_df = counts_df.apply(pd.to_numeric, errors="coerce")
    counts_df = counts_df[counts_df.sum(axis=1) > 20]
    counts_df = counts_df.T

    # Run DESeq2 normalization
    inference = DefaultInference(n_cpus=threads)

    data_set = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design_factors=["run"],
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )

    data_set.deseq2()

    output_file("library_sizes.html", title="Library Sizes")
    # Plot non normalized library counts
    non_norm_plot = create_library_size_plot(
        counts_df, title="Non-normalized library counts", yAxisLabel="Counts"
    )

    norm_data: pd.DataFrame = counts_df.div(data_set.obsm["size_factors"], axis=0)  # type: ignore
    norm_plot = create_library_size_plot(
        norm_data, title="Normalized library counts", yAxisLabel="Counts"
    )

    library_size_plots = row(non_norm_plot, norm_plot, sizing_mode="stretch_width")  # type: ignore
    save(library_size_plots)

    # Run DESeq2 analysis for each run compared to the first run as the ground truth
    runs = sorted(metadata_df["run"].unique())
    ground_truth = runs.pop(0)

    for run in runs:
        comparison_name: str = f"run_{run}_vs_{ground_truth}"

        out_path = comparison_name
        if not os.path.exists(out_path):
            os.mkdir(out_path)

        contrast = ["run", run, ground_truth]
        stat_res = DeseqStats(
            data_set, contrast=contrast, inference=inference, quiet=True
        )
        stat_res.summary()
        stat_res.lfc_shrink(comparison_name)

        output_file(f"{out_path}/ma_plot.html")
        ma_plot = create_ma_plot(comparison_name, stat_res)
        save(ma_plot)  # type: ignore

        comparison_metadata = metadata_df[
            metadata_df["run"].str.contains(f"{run}|{ground_truth}")
        ]
        norm_filtered_counts_df = norm_data.loc[comparison_metadata.index].T

        norm_filtered_counts_df = norm_filtered_counts_df.loc[
            stat_res.results_df["padj"] < 0.05
        ]

        output_file(f"{out_path}/pca_plot.html")
        pca_plot = create_pca_plot(
            comparison_name, norm_filtered_counts_df, comparison_metadata
        )
        save(pca_plot)  # type: ignore

        output_file(f"{out_path}/heatmap.html")
        heatmap_plot = create_de_heatmap_plot(comparison_name, norm_filtered_counts_df)
        save(heatmap_plot)  # type: ignore

        deseq_path = f"{out_path}/deseq2_results_{run}_vs_{ground_truth}.tsv"

        with open(deseq_path, "w") as f:
            stat_res.results_df.to_csv(
                f, sep="\t", index=True, header=True, float_format="%.10f"
            )
