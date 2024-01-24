# Differential enrichment analysis using DESeq2

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

parser = argparse.ArgumentParser(description="Experiment metadata file.")
parser.add_argument("--metadata_file", type=str, help="Path to the metadata file")

args = parser.parse_args()
metadata_file_path = args.metadata_file

# The metadata should have following columns:
# 1. Sample name "sample"
# 2. Condition / Run number "run"
# 3. Replicate number "replicate"
# 4. Path to the count file "count_file"
metadata: pd.DataFrame = pd.read_csv(
    metadata_file_path, sep="\t", index_col=None, header=0
)

# Read in the count files and merge them into a single data frame by the "Name" column
count_files = metadata["count_file"].tolist()
sample_names = metadata["sample"].tolist()

dfs = []

# Read and merge the count files
for count_file, sample_name in zip(count_files, sample_names):
    df = pd.read_csv(count_file, sep="\t", index_col=0, header=None)
    df = df.iloc[:-5]
    df.columns = [sample_name]
    dfs.append(df)


counts_df: pd.DataFrame = dfs[0].join(dfs[1:])
counts_df = counts_df[counts_df.sum(axis=1) > 50]
counts_df = counts_df.T

metadata.set_index("sample", inplace=True)

# Run DESeq2 analysis
inference = DefaultInference(n_cpus=8)

data_set = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors=["run"],
    refit_cooks=True,
    inference=inference,
    quiet=True,
)
data_set.deseq2()

# Plot raw and normalized data library sizes to validate normalization
raw_data = counts_df.T.melt(var_name="columns")
raw_plot = sns.boxplot(data=raw_data, x="columns", y="value")
plt.yscale("log")
plt.savefig("raw_data.png")
plt.clf()


norm_data: pd.DataFrame = counts_df.div(data_set.obsm["size_factors"], axis=0) # type: ignore

norm_plot_data = norm_data.copy()
norm_plot_data["run"] = metadata["run"].tolist()
norm_plot_data["sample"] = norm_plot_data.index.tolist()
norm_data_melted = norm_plot_data.melt(id_vars=["sample", "run"], var_name="columns")

norm_plot = sns.boxplot(data=norm_data_melted, x="sample", y="value", hue="run")
plt.yscale("log")
plt.savefig("norm_data.png")
plt.clf()


# Run DESeq2 analysis for each run compared to the first run as the ground truth
runs_df = metadata["run"].unique()
runs_df.sort()

runs: list[str] = runs_df.tolist()
ground_truth = runs.pop(0)

for run in runs:
    # Create dir for comparison#
    out_path = f"run_{run}_vs_{ground_truth}"
    os.mkdir(out_path)

    contrast = ["run", run, ground_truth]
    stat_res = DeseqStats(data_set, contrast=contrast, inference=inference, quiet=True)
    stat_res.summary()
    stat_res.lfc_shrink(f"run_{run}_vs_{ground_truth}")
    # stat_res.results_df.sort_values("padj", inplace=True)
    ma_path = f"{out_path}/run_{run}_vs_{ground_truth}_MA.png"
    stat_res.plot_MA(save_path=ma_path)
    plt.clf()

    # Plot PCA for the run vs ground truth
    # Filter rows by sample names for run and ground truth
    comparison_metadata = metadata[
        metadata["run"].str.contains(f"{run}|{ground_truth}")
    ]
    filtered_counts_df = counts_df.loc[comparison_metadata.index].T

    # Filter out the genes with padj > 0.05
    filtered_count_df = filtered_counts_df.loc[stat_res.results_df["padj"] < 0.05]

    log2_norm_data = np.log2(filtered_count_df.T + 1)

    pca = PCA(n_components=2)
    pca_df = pd.DataFrame(pca.fit_transform(log2_norm_data), columns=["PC1", "PC2"])

    pca_df["run"] = comparison_metadata["run"].tolist()
    pca_df["sample"] = pca_df.index.tolist()

    pca_plot = sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="run")

    pca_path = f"{out_path}/run_{run}_vs_{ground_truth}_PCA.png"
    plt.savefig(pca_path)
    plt.clf()

    deseq_path = f"{out_path}/deseq2_results_{run}_vs_{ground_truth}.tsv"

    with open(deseq_path, "w") as f:
        stat_res.results_df.to_csv(
            f, sep="\t", index=True, header=True, float_format="%.10f"
        )
