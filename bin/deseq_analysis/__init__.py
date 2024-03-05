#!/usr/bin/env python
"""Perform differential expression analysis using DESeq2."""
import os
from argparse import ArgumentParser
import numpy as np
import pandas as pd

from math import pi

from pathlib import Path
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from bokeh.plotting import figure, row
from bokeh import model
from bokeh.transform import factor_cmap
from bokeh.palettes import Pastel1
from bokeh.io import save, output_file

from .plots.library_size_plot import createLibrarySizePlot
from .plots.de_heatmap_plot import create_de_heatmap_plot

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

    replicate_names: list[str] = []
    run_names: list[str] = []
    sample_names: list[str] = []
    count_files: list[str] = []

    # Expects quant files being prepared by featureCounts with merge_feature_counts.py
    for quant_file in quantFiles:
        path = Path(quant_file)
        file_name: str = path.stem
        name_components: list[str] = file_name.split("_")

        run_name = name_components[1]
        replicate_name = name_components[3]
        sample_name = f"{run_name}_{replicate_name}"

        replicate_names.append(replicate_name)
        run_names.append(run_name)
        sample_names.append(sample_name)
        count_files.append(quant_file)

    metadata: dict[str, list[str]] = {
        "replicate": replicate_names,
        "run": run_names,
        "sample": sample_names,
        "count_file": count_files,
    }

    metadata_df = pd.DataFrame(metadata)

    # Read in the count files and merge them into a single data frame by the "Name" column
    count_files = metadata_df["count_file"].tolist()
    sample_names = metadata_df["sample"].tolist()

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

    metadata_df.set_index("sample", inplace=True)

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
    nonNormPlot = createLibrarySizePlot(
        counts_df, title="Non-normalized library counts", yAxisLabel="Counts"
    )
    norm_data: pd.DataFrame = counts_df.div(data_set.obsm["size_factors"], axis=0)  # type: ignore

    normPlot = createLibrarySizePlot(
        norm_data, title="Normalized library counts", yAxisLabel="Counts"
    )

    librarySizePlots = row(nonNormPlot, normPlot, sizing_mode="stretch_width")  # type: ignore
    save(librarySizePlots)

    # Run DESeq2 analysis for each run compared to the first run as the ground truth
    runs_df = metadata_df["run"].unique()
    runs_df.sort()

    runs: list[str] = runs_df.tolist()
    ground_truth = runs.pop(0)

    for run in runs:
        name = f"run_{run}_vs_{ground_truth}"
        # Create dir for comparison#
        out_path = name
        if not os.path.exists(out_path):
            os.mkdir(out_path)

        contrast = ["run", run, ground_truth]
        stat_res = DeseqStats(
            data_set, contrast=contrast, inference=inference, quiet=True
        )
        stat_res.summary()
        stat_res.lfc_shrink(name)

        # Plot MA plot
        output_file(f"{out_path}/ma_plot.html")
        colors = ["red" if x < 0.05 else "black" for x in stat_res.results_df["padj"]]
        ma_plot = figure(
            title=f"MA Plot {name}",
            x_axis_label="log2 Fold Change",
            y_axis_label="log2 Counts",
            sizing_mode="stretch_width",
        )
        ma_plot.circle(
            stat_res.results_df["log2FoldChange"],
            stat_res.results_df["baseMean"],
            line_color=None,
            fill_color=colors,
            fill_alpha=0.6,
            size=6,
        )
        save(ma_plot)  # type: ignore

        # Plot PCA for the run vs ground truth
        # Filter rows by sample names for run and ground truth
        comparison_metadata = metadata_df[
            metadata_df["run"].str.contains(f"{run}|{ground_truth}")
        ]
        norm_filtered_counts_df = norm_data.loc[comparison_metadata.index].T

        # Filter out the genes with padj > 0.05
        norm_filtered_counts_df = norm_filtered_counts_df.loc[
            stat_res.results_df["padj"] < 0.05
        ]

        log2_norm_data = np.log2(norm_filtered_counts_df.T + 1)

        pca = PCA(n_components=2)
        pca_df = pd.DataFrame(pca.fit_transform(log2_norm_data), columns=["PC1", "PC2"])

        pca_df["run"] = comparison_metadata["run"].tolist()
        pca_df["sample"] = pca_df.index.tolist()

        output_file(f"{out_path}/pca_plot.html")
        pca_plot = figure(
            title=f"PCA Plot {name}",
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
            fill_color=factor_cmap("run", Pastel1[3], runs),
            fill_alpha=0.8,
            legend_field="run",
            size=18,
        )
        save(pca_plot)  # type: ignore

        output_file(f"{out_path}/heatmap.html")
        p = create_de_heatmap_plot(
            f"Differential gene expression {name}", norm_filtered_counts_df
        )
        save(p)  # type: ignore

        deseq_path = f"{out_path}/deseq2_results_{run}_vs_{ground_truth}.tsv"

        with open(deseq_path, "w") as f:
            stat_res.results_df.to_csv(
                f, sep="\t", index=True, header=True, float_format="%.10f"
            )
