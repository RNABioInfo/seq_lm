#!/usr/bin/env python3
"""Evaluate saved edgeR profile exports; never invoke Nextflow or sequencing tools.

Data availability is taken from the saved checkpoint clock. Real delivery and
processing latency are unavailable retrospectively and are tested separately.
"""

import argparse
import csv
import gzip
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from dea_termination.engine import atomic_json, clean, evaluate, write_tsv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--horizon-minutes", type=int, default=30)
    p.add_argument("--sample-every-minutes", type=int, default=10)
    a = p.parse_args()
    if a.horizon_minutes < 1 or a.sample_every_minutes < 1:
        p.error("horizon and observation spacing must be positive")
    root = a.reference
    a.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((root / "input_manifest.json").read_text())
    with (root / manifest["files"]["universe"]).open() as f:
        universe = [r["feature_id"] for r in csv.DictReader(f, delimiter="\t")]
    with (root / manifest["files"]["acquisitions"]).open() as f:
        ends = {
            r["run_id"]: float(r["acquisition_duration_seconds"])
            for r in csv.DictReader(f, delimiter="\t")
        }
    definitions = {c["id"]: c for c in manifest["contrasts"]}
    diagnostics = {}
    for cid in definitions:
        with np.load(root / f"{cid}_elbow_diagnostics.npz", allow_pickle=False) as z:
            diagnostics[cid] = {
                "elapsed": z["elapsed"],
                "totals": z["assigned_totals"],
                "R": z["R"],
                "FDR": z["FDR"],
                "FC": z["logFC"],
            }
    cuts = (0.005, 0.01, 0.02)
    states = {}
    records = []
    summary = []
    started = time.monotonic()
    with gzip.open(root / manifest["files"]["observations"], "rt") as f:
        grouped = itertools.groupby(
            csv.DictReader(f, delimiter="\t"),
            lambda r: (r["session"], int(r["batch_index"])),
        )
        previous = None
        for (cid, batch), group in grouped:
            saved = list(group)
            t = float(saved[0]["elapsed_seconds"])
            if (
                a.sample_every_minutes != 10
                and t % (a.sample_every_minutes * 60) > 1e-6
            ):
                continue
            definition = definitions[cid]
            diag = diagnostics[cid]
            i = int(np.flatnonzero(np.abs(diag["elapsed"] - t) < 0.001)[0])
            totals = diag["totals"][i]
            if previous != cid:
                print(f"Evaluating {cid} at all three cutoffs", flush=True)
                previous = cid
            precision_rows = {
                r["feature_id"]: dict(
                    feature_id=r["feature_id"],
                    R=r["R_profile"],
                )
                for r in saved
            }
            rows = [
                precision_rows.get(
                    g,
                    dict(
                        feature_id=g,
                        R="NA",
                    ),
                )
                for g in universe
            ]
            samples = []
            contrast_samples = []
            for j, s in enumerate(definition["samples"]):
                group_name = "T0" if j < 3 else cid
                key = f"{group_name}/{s['id']}"
                active = j >= 3 and t < ends[s["acquisition"]]
                samples.append(
                    {
                        "key": key,
                        "name": s["id"],
                        "group": group_name,
                        "bam_dir": "/unused",
                        "protocol_run_id": s.get("acquisition", ""),
                        "live": active,
                        "declared_live": j >= 3,
                    }
                )
                contrast_samples.append({"key": key, "group": group_name})
            base = {
                "samples": samples,
                "contrasts": [
                    {"contrast_id": cid, "target_group": cid, "reference_group": "T0"}
                ],
                "identities": {"evaluation": "saved-profile-exports"},
                "started_at": 0.0,
                "observed_at": t,
                "completed_at": t,
                "behavior": "log",
                "precision": {
                    cid: {
                        "rows": rows,
                        "samples": contrast_samples,
                        "depths": {s["key"]: float(d) for s, d in zip(samples, totals)},
                        "universe": universe,
                    }
                },
            }
            for cutoff in cuts:
                state, metrics, _, _ = evaluate(
                    base
                    | {
                        "config": {
                            "gain_cutoff": cutoff,
                            "horizon_seconds": a.horizon_minutes * 60,
                        }
                    },
                    states.get((cid, cutoff)),
                )
                states[cid, cutoff] = state
                records.append(
                    dict(
                        contrast=cid,
                        batch_index=batch,
                        gain_cutoff=cutoff,
                        **metrics[0],
                    )
                )
    for cid, definition in definitions.items():
        diag = diagnostics[cid]
        times = diag["elapsed"]
        final_calls = (diag["FDR"][-1] <= 0.05) & (np.abs(diag["FC"][-1]) >= 1)
        for cutoff in cuts:
            candidates = [
                r
                for r in records
                if r["contrast"] == cid and r["gain_cutoff"] == cutoff and r["eligible"]
            ]
            stop = candidates[0] if candidates else None
            actual = None
            recall = None
            if stop:
                i = int(
                    np.flatnonzero(np.abs(times - stop["elapsed_seconds"]) < 0.001)[0]
                )
                state = states[cid, cutoff]
                old = next(
                    (
                        r
                        for r in state["audit_history"]
                        if r["elapsed_seconds"] == stop["elapsed_seconds"]
                    ),
                    None,
                )
                actual = old.get("actual_gain_next_horizon") if old else None
                calls = (diag["FDR"][i] <= 0.05) & (np.abs(diag["FC"][i]) >= 1)
                recall = (
                    float(np.sum(calls & final_calls) / np.sum(final_calls))
                    if np.any(final_calls)
                    else None
                )
            summary.append(
                {
                    "contrast": cid,
                    "gain_cutoff": cutoff,
                    "stop_seconds": stop["elapsed_seconds"] if stop else None,
                    "reason": stop["reason"] if stop else "no_precision_stop_observed",
                    "predicted_conservative_gain": stop["predicted_conservative_gain"]
                    if stop
                    else None,
                    "actual_next_horizon_gain": actual,
                    "endpoint_DE_call_recall": recall,
                    "eligible_checkpoints": len(
                        [
                            r
                            for r in records
                            if r["contrast"] == cid and r["gain_cutoff"] == cutoff and r["eligible"]
                        ]
                    ),
                }
            )
    write_tsv(a.output / "checkpoint_decisions.tsv", records)
    write_tsv(a.output / "cutoff_comparison.tsv", summary)
    atomic_json(
        a.output / "methods.json",
        {
            "method": "Curve-only R forecast with causal error calibration; interruption guards enabled; no yield fallback",
            "cutoffs": cuts,
            "horizon_minutes": a.horizon_minutes,
            "sample_every_minutes": a.sample_every_minutes,
            "source": str(root),
            "elapsed_seconds": time.monotonic() - started,
            "limitations": [
                "Saved acquisition/checkpoint times approximate data availability; live latency is not reconstructed.",
                "Final-depth calls are a reference, not biological truth.",
                "No automatic low-yield stop is used for T4.",
            ],
        },
    )
    table = "".join(
        "<tr>"
        + "".join(
            f"<td>{v}</td>"
            for v in [
                r["contrast"],
                f"{r['gain_cutoff'] * 100:g}%",
                "No precision stop"
                if r["stop_seconds"] is None
                else f"{r['stop_seconds'] / 3600:.2f} h",
                "—"
                if r["actual_next_horizon_gain"] is None
                else f"{r['actual_next_horizon_gain'] * 100:.3f}%",
                "—"
                if r["endpoint_DE_call_recall"] is None
                else f"{r['endpoint_DE_call_recall'] * 100:.1f}%",
            ]
        )
        + "</tr>"
        for r in summary
    )
    html = (
        """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>seqLM R stopping validation</title><style>body{font:17px/1.6 system-ui;max-width:1100px;margin:32px auto;padding:16px;color:#172b42}table{border-collapse:collapse;width:100%}td,th{padding:8px;border-bottom:1px solid #ccd}th{text-align:left}svg{max-width:100%;height:auto}</style><h1>R-based stopping without pore-check timestamps</h1><p>The evaluated workflow uses only completed checkpoints and read increments. A throughput interruption blocks precision stopping until recovery. There is no low-yield stopping fallback.</p><p>The forecast covers the configured horizon. It describes the custom aggregate uncertainty score, not the percentage of genes that are correct. The full-depth endpoint is a comparison reference, not ground truth.</p><table><tr><th>Contrast</th><th>Cutoff</th><th>First stop</th><th>Actual next horizon improvement</th><th>Final DE calls retained</th></tr>"""
        + f"<caption>Forecast horizon: {a.horizon_minutes} minutes; observations every {a.sample_every_minutes} minutes.</caption>"
        + table
        + "</table><p>T4 has no eligible edgeR uncertainty models. It cannot receive a precision-based stopping decision. Saved checkpoints cannot establish live file-delivery latency; interruption and processing-delay behavior are tested separately.</p>"
    )
    # Scientific figures are standalone SVGs and embedded for offline viewing.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for cid in definitions:
        rows = [r for r in records if r["contrast"] == cid and r["gain_cutoff"] == 0.01]
        if not rows:
            continue
        fig, ax = plt.subplots(figsize=(10, 3.5))
        ax.plot(
            [r["elapsed_seconds"] / 3600 for r in rows],
            [r["R_information"] for r in rows],
            label="Observed aggregate R",
        )
        for cutoff, color in zip(cuts, ["#0072B2", "#D55E00", "#009E73"]):
            stop = next(
                r
                for r in summary
                if r["contrast"] == cid and r["gain_cutoff"] == cutoff
            )
            if stop["stop_seconds"] is not None:
                ax.axvline(
                    stop["stop_seconds"] / 3600,
                    color=color,
                    label=f"{cutoff * 100:g}% stop",
                )
        ax.set(
            xlabel="Saved sequencing checkpoint time (hours)",
            ylabel="Aggregate R",
            title=f"{cid} versus fixed full-depth T0",
        )
        ax.set_ylim(bottom=0)
        ax.legend()
        fig.tight_layout()
        path = a.output / f"{cid}_R.svg"
        fig.savefig(path)
        plt.close(fig)
        svg = path.read_text()
        html += svg[svg.index("<svg") :]
    (a.output / "report.html").write_text(html + "</html>")
    print(json.dumps(clean(summary), indent=2), flush=True)


if __name__ == "__main__":
    main()
