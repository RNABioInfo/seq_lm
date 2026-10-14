from pathlib import Path

import pandas as pd


MANIFEST_COLUMNS = {"name", "group", "count_file"}


def get_metadata(quant_manifest: str) -> pd.DataFrame:
    """Build DESeq2 sample metadata from the quantification manifest."""
    manifest_path = Path(quant_manifest)
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


def read_quantification_counts(count_file: str) -> pd.Series:
    """Read an Oarfish ``.quant`` count vector."""
    counts = pd.read_csv(count_file, sep="\t", header=0)
    if not {"tname", "num_reads"}.issubset(counts.columns):
        raise ValueError(f"Quantification file requires tname and num_reads: {count_file}")
    return counts.set_index("tname")["num_reads"]
