"""Evaluate a completed workflow snapshot against its committed R state."""

import argparse
import json
from pathlib import Path

from .engine import FORECAST_METHOD, atomic_json, evaluate, read_tsv, write_tsv


def load_precision(root, definitions, selected=False):
    """Load complete exports, restricting only the monitoring universe when enabled."""
    precision = {}
    selected_ids = None
    for definition in definitions:
        folder = root / definition["contrast_id"]
        if not (folder / "precision.tsv").exists():
            continue
        samples = read_tsv(folder / "precision_samples.tsv")
        rows = read_tsv(folder / "precision.tsv")
        if selected:
            if selected_ids is None:
                selection = root / "dea_stopping_features.tsv"
                if not selection.is_file():
                    raise ValueError("Selected DEA stopping requires dea_stopping_features.tsv")
                selected_ids = [r.get("feature_id", "") for r in read_tsv(selection)]
                if not selected_ids or any(not x for x in selected_ids):
                    raise ValueError("DEA stopping feature selection is empty or invalid")
                if len(selected_ids) != len(set(selected_ids)):
                    raise ValueError("Duplicate DEA stopping feature IDs")
            selected_set = set(selected_ids)
            missing = selected_set - {r["feature_id"] for r in rows}
            if missing:
                raise ValueError(f"Selected DEA stopping features absent from precision export: {sorted(missing)}")
            rows = [r for r in rows if r["feature_id"] in selected_set]
        precision[definition["contrast_id"]] = {
            "rows": rows,
            "samples": samples,
            "depths": {x["key"]: float(x["assigned_reads"]) for x in samples},
            "universe": selected_ids if selected else [r["feature_id"] for r in rows],
        }
    return precision


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--state", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    payload = json.loads(a.input.read_text())
    root = Path(payload.pop("results"))
    payload["precision"] = load_precision(
        root, payload["contrasts"],
        selected=payload["identities"].get("dea_stopping_genes") is not None,
    )
    previous = json.loads(a.state.read_text()) if a.state.exists() else None
    state, contrasts, samples, rates = evaluate(payload, previous)
    a.output.mkdir(parents=True, exist_ok=True)
    # State is staged beside the audit. The caller commits it after action auditing.
    atomic_json(a.output / "state.json", state)
    atomic_json(
        a.output / "audit.json",
        {"forecast_method": FORECAST_METHOD, "contrasts": contrasts,
         "rates": rates, "config": state["config"]},
    )
    write_tsv(a.output / "contrast_monitoring.tsv", contrasts)
    write_tsv(a.output / "sample_monitoring.tsv", samples)
    write_tsv(a.output / "read_rates.tsv", rates)
    write_tsv(a.output / "forecast_history.tsv", state["audit_history"])
    write_tsv(
        a.output / "read_rate_history.tsv",
        state["rate_history"],
        [
            "contrast_id",
            "elapsed_seconds",
            "sample",
            "predicted_increment",
            "previous_horizon_increment",
            "eligible_rate_intervals",
            "last_interval_reads_per_minute",
            "productive_reads_per_minute",
            "interrupted",
            "recovery_intervals",
            "rate_coverage_seconds",
            "cadence_adequate",
            "last_interval_seconds",
            "data_available_seconds",
            "recovery_seconds",
        ],
    )
    history = []
    for cid, s in state["contrasts"].items():
        history.extend(dict(contrast_id=cid, **x) for x in s["points"])
    write_tsv(
        a.output / "R_history.tsv", history, ["contrast_id", "time", "depth", "R"]
    )


if __name__ == "__main__":
    main()
