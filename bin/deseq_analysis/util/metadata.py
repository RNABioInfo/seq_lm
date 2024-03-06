import pandas as pd

from pathlib import Path


def get_metadata(quant_files: list[str]) -> pd.DataFrame:
    replicate_names: list[str] = []
    run_names: list[str] = []
    sample_names: list[str] = []
    count_files: list[str] = []

    for quant_file in quant_files:
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

    return metadata_df
