from pathlib import Path

import pandas as pd


MANIFEST_COLUMNS = {"name", "group", "count_file"}


def get_metadata(
    quant_files: list[str] | None = None,
    quant_manifest: str | None = None,
) -> pd.DataFrame:
    """Build DESeq2 sample metadata from a live manifest or legacy filenames."""
    if quant_manifest is not None:
        return _metadata_from_manifest(Path(quant_manifest))

    return _metadata_from_legacy_quant_files(quant_files or [])


def _metadata_from_manifest(manifest_path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(manifest_path, sep="\t", dtype=str)
    missing_columns = MANIFEST_COLUMNS - set(metadata.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Quantification manifest is missing required columns: {missing}")
    if metadata.empty:
        raise ValueError("Quantification manifest contains no samples.")

    metadata = metadata[["name", "group", "count_file"]].copy()
    for column in ("name", "group", "count_file"):
        metadata[column] = metadata[column].fillna("").str.strip()
        if (metadata[column] == "").any():
            raise ValueError(f"Quantification manifest contains an empty {column} value.")

    metadata["group"] = metadata["group"].map(
        lambda group: "control" if group.casefold() == "control" else group
    )
    metadata["sample"] = metadata["group"] + "/" + metadata["name"]
    duplicate_samples = metadata.loc[metadata["sample"].duplicated(), "sample"].tolist()
    if duplicate_samples:
        duplicates = ", ".join(sorted(set(duplicate_samples)))
        raise ValueError(f"Quantification manifest contains duplicate samples: {duplicates}")

    manifest_dir = manifest_path.resolve().parent
    metadata["count_file"] = metadata["count_file"].map(
        lambda value: str(
            (manifest_dir / value).resolve() if not Path(value).is_absolute() else Path(value)
        )
    )
    return metadata


def _metadata_from_legacy_quant_files(quant_files: list[str]) -> pd.DataFrame:
    if not quant_files:
        raise ValueError("Provide --quant-manifest or at least one --quant-file.")

    replicate_names: list[str] = []
    run_names: list[str] = []
    sample_names: list[str] = []
    count_files: list[str] = []

    for quant_file in sorted(quant_files):
        path = Path(quant_file)
        name_components = path.stem.split("_")
        if len(name_components) < 4:
            raise ValueError(
                f"Legacy quantification filename does not match run_N_replicate_N: {path.name}"
            )

        run_name = name_components[1]
        replicate_name = name_components[3]
        sample_name = f"{run_name}_{replicate_name}"

        replicate_names.append(replicate_name)
        run_names.append(run_name)
        sample_names.append(sample_name)
        count_files.append(str(path))

    metadata = pd.DataFrame(
        {
            "replicate": replicate_names,
            "run": run_names,
            "sample": sample_names,
            "count_file": count_files,
        }
    )
    baseline_run = sorted(metadata["run"].unique())[0]
    metadata["group"] = metadata["run"].map(
        lambda run: "control" if run == baseline_run else run
    )
    return metadata


def read_quantification_counts(count_file: str) -> pd.Series:
    """Read Oarfish ``.quant`` or legacy featureCounts output as one count vector."""
    counts = pd.read_csv(count_file, sep="\t", header=0)

    if {"tname", "num_reads"}.issubset(counts.columns):
        return counts.set_index("tname")["num_reads"]

    if counts.shape[1] < 2:
        raise ValueError(f"Quantification file has no count column: {count_file}")

    return counts.set_index(counts.columns[0])[counts.columns[-1]]
