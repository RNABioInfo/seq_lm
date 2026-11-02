"""Evaluate a completed workflow snapshot against its committed R state."""

import argparse
import json
from pathlib import Path

from .engine import FORECAST_METHOD, atomic_json, evaluate, read_tsv, write_tsv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--state", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    payload = json.loads(a.input.read_text())
    root = Path(payload.pop("results"))
    precision = {}
    for d in payload["contrasts"]:
        folder = root / d["contrast_id"]
        if not (folder / "precision.tsv").exists():
            continue
        samples = read_tsv(folder / "precision_samples.tsv")
        rows = read_tsv(folder / "precision.tsv")
        precision[d["contrast_id"]] = {
            "rows": rows,
            "samples": samples,
            "depths": {x["key"]: float(x["assigned_reads"]) for x in samples},
            "universe": [r["feature_id"] for r in rows],
        }
    payload["precision"] = precision
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
