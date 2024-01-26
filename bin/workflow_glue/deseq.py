#!/usr/bin/env python
"""Perform differential expression analysis using DESeq2."""
import os

os.environ["MPLCONFIGDIR"] = "/tmp"

from .util import get_named_logger, wf_parser
import numpy as np
import pandas as pd

from pathlib import Path
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

from bokeh.models import ColumnDataSource, Whisker
from bokeh.plotting import figure, show, row
from bokeh import model
from bokeh.transform import factor_cmap
from bokeh.palettes import Pastel1


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("Differential expression analysis.")
    parser.add_argument("-q", "--quant_files", type=str, nargs="+", default=[])
    parser.add_argument("-t", "--threads", type=int, default=1)

    return parser


def createLibrarySizePlot(
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


def main(args):
    quantFiles: list[str] = args.quant_files
    threads: int = args.threads

    replicate_names: list[str] = []
    run_names: list[str] = []
    sample_names: list[str] = []
    count_files: list[str] = []

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
        df: pd.DataFrame = pd.read_csv(
            count_file, sep="\t", skiprows=1, header=1, index_col=0
        )
        df = df.iloc[:, [-1]]
        df.columns = [sample_name]
        dfs.append(df)

    counts_df: pd.DataFrame = dfs[0].join(dfs[1:])

    print(counts_df)
    # set all columns data type to int
    counts_df = counts_df.astype(int)
    counts_df = counts_df[counts_df.sum(axis=1) > 50]
    counts_df = counts_df.T

    metadata_df.set_index("sample", inplace=True)

    # Run DESeq2 analysis
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

    # Plot non normalized library counts
    nonNormPlot = createLibrarySizePlot(
        counts_df, title="Non-normalized library counts", yAxisLabel="Counts"
    )

    norm_data: pd.DataFrame = counts_df.div(data_set.obsm["size_factors"], axis=0)  # type: ignore

    normPlot = createLibrarySizePlot(
        norm_data, title="Normalized library counts", yAxisLabel="Counts"
    )

    # show(row(nonNormPlot, normPlot, sizing_mode="stretch_width"))  # type: ignore

    # Run DESeq2 analysis for each run compared to the first run as the ground truth
    runs_df = metadata_df["run"].unique()
    runs_df.sort()

    runs: list[str] = runs_df.tolist()
    ground_truth = runs.pop(0)

    for run in runs:
        # Create dir for comparison#
        out_path = f"run_{run}_vs_{ground_truth}"
        if not os.path.exists(out_path):
            os.mkdir(out_path)

        contrast = ["run", run, ground_truth]
        stat_res = DeseqStats(
            data_set, contrast=contrast, inference=inference, quiet=True
        )
        stat_res.summary()
        stat_res.lfc_shrink(f"run_{run}_vs_{ground_truth}")

        # FIXME: This is not working currently
        # stat_res.results_df.sort_values("padj", inplace=True)
        # ma_path = f"{out_path}/run_{run}_vs_{ground_truth}_MA.png"
        # stat_res.plot_MA(save_path=ma_path)
        # plt.clf()

        # Plot PCA for the run vs ground truth
        # Filter rows by sample names for run and ground truth
        comparison_metadata = metadata_df[
            metadata_df["run"].str.contains(f"{run}|{ground_truth}")
        ]
        filtered_counts_df = counts_df.loc[comparison_metadata.index].T

        # Filter out the genes with padj > 0.05
        filtered_count_df = filtered_counts_df.loc[stat_res.results_df["padj"] < 0.05]

        log2_norm_data = np.log2(filtered_count_df.T + 1)

        pca = PCA(n_components=2)
        pca_df = pd.DataFrame(pca.fit_transform(log2_norm_data), columns=["PC1", "PC2"])

        pca_df["run"] = comparison_metadata["run"].tolist()
        pca_df["sample"] = pca_df.index.tolist()

        # FIXME: This is not working currently
        # pca_plot = sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="run")

        # pca_path = f"{out_path}/run_{run}_vs_{ground_truth}_PCA.png"
        # plt.savefig(pca_path)
        # plt.clf()

        deseq_path = f"{out_path}/deseq2_results_{run}_vs_{ground_truth}.tsv"

        with open(deseq_path, "w") as f:
            stat_res.results_df.to_csv(
                f, sep="\t", index=True, header=True, float_format="%.10f"
            )


if __name__ == "__main__":
    main(argparser().parse_args())
