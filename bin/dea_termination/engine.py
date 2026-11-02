"""DEA sequencing run termination engine
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import tempfile
from itertools import pairwise
from pathlib import Path

from .forecast import curve_summary, information_radius, quantile

STATE_SCHEMA_VERSION = 3
FORECAST_METHOD = "depth_curve"

DEFAULTS = {
    "gain_cutoff": 0.01,
    "information_radius_scale": 1.0,
    "max_processing_lag_seconds": 600,
    "min_points": 6,
    "earliest_seconds": 3600,
    "confirmations": 3,
    "horizon_seconds": 1800,
    "interruption_fraction": 0.25,
    "recovery_fraction": 0.5,
    "recovery_intervals": 2,
    "rate_quantile": 0.75,
    "minimum_rate_intervals": 3,
    "max_fit_discrepancy": 0.1,
    "elbow_min_strength": 0.1,
    "elbow_min_observed_drop": 0.05,
    "post_elbow_points": 3,
    "max_raw_increase_30min": 0.03,
    "minimum_forecast_margin": 0.002,
    "calibration_points": 6,
    "calibration_quantile": 0.9,
}


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, float):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix="." + path.name)
    try:
        with os.fdopen(fd, "w") as out:
            json.dump(clean(value), out, allow_nan=False, sort_keys=True)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def finite(x):
    return x is not None and math.isfinite(float(x))


def number(x):
    return float(x) if x not in ("", None, "NA", "NaN") else math.nan


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, columns=None):
    columns = columns or list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            r = clean(row)
            w.writerow(
                {
                    k: (
                        "true"
                        if v is True
                        else "false"
                        if v is False
                        else ""
                        if v is None
                        else v
                    )
                    for k, v in r.items()
                }
            )


def validate_config(config):
    unknown = config.keys() - DEFAULTS.keys()
    if unknown:
        raise ValueError(f"unknown R settings: {sorted(unknown)}")
    c = DEFAULTS | config
    for k, v in c.items():
        if (
            isinstance(v, bool)
            or not isinstance(v, (int, float))
            or not math.isfinite(v)
        ):
            raise ValueError(f"invalid R setting {k}")
    for k in (
        "horizon_seconds",
        "max_processing_lag_seconds",
        "min_points",
        "earliest_seconds",
        "confirmations",
        "recovery_intervals",
    ):
        if c[k] < 1 or int(c[k]) != c[k]:
            raise ValueError(f"{k} must be a positive integer")
    if c["min_points"] < 6 or c["earliest_seconds"] < 3600:
        raise ValueError("R monitoring requires at least six points and one hour")
    for k in ("gain_cutoff", "interruption_fraction", "recovery_fraction"):
        if not 0 < c[k] < 1:
            raise ValueError(f"{k} must be between zero and one")
    if c["information_radius_scale"] <= 0:
        raise ValueError("information_radius_scale must be positive")
    return c


def weighted_quantile(pairs, q):
    pairs = sorted((v, w) for v, w in pairs if w > 0)
    target = sum(w for v, w in pairs) * q
    total = 0.0
    for value, weight in pairs:
        total += weight
        if total >= target:
            return value
    return math.nan


def rate_state(history, t, state, c):
    """Integrate actual arrival intervals. Missing output is never a zero bin."""
    H = c["horizon_seconds"]
    intervals = []

    for left, right in pairwise(history):
        dt = right["time"] - left["time"]
        delta = right["depth"] - left["depth"]
        if dt <= 0 or delta < -1e-6:
            raise ValueError(
                "invalid read-arrival interval or decreasing assigned-read counter"
            )
        intervals.append((left["time"], right["time"], max(0.0, delta) / dt))

    recent = [
        (rate, min(end, t) - max(start, t - 3 * H))
        for start, end, rate in intervals
        if end > t - 3 * H and start < t and end - start <= H
    ]

    prod = (
        weighted_quantile(recent, c["rate_quantile"])
        if len(recent) >= c["minimum_rate_intervals"]
        else math.nan
    )

    window = [
        (rate, min(end, t) - max(start, t - H))
        for start, end, rate in intervals
        if end > t - H and start < t and end - start <= H
    ]

    coverage = sum(dt for rate, dt in window)
    covered = abs(coverage - H) <= 1e-6
    last = sum(rate * dt for rate, dt in window) if covered else math.nan
    reference = state.get("reference")
    interrupted = state.get("interrupted", False)

    # Also detect a long low-yield interval on return, without treating silence as data.
    newest = intervals[-1] if intervals else None
    changed = newest is not None and newest[1] != state.get("last_interval_time")

    if changed:
        start, end, rate = newest
        if (
            finite(reference)
            and reference > 0
            and (
                (finite(last) and last / H < reference * c["interruption_fraction"])
                or (end - start >= H and rate < reference * c["interruption_fraction"])
            )
        ):
            interrupted = True

        if interrupted:
            if (
                end - start <= H
                and finite(reference)
                and rate >= reference * c["recovery_fraction"]
            ):
                state["recovery"] = state.get("recovery", 0) + 1
                state["recovery_seconds"] = (
                    state.get("recovery_seconds", 0) + end - start
                )
            else:
                state["recovery"] = 0
                state["recovery_seconds"] = 0
            if (
                state["recovery"] >= c["recovery_intervals"]
                and state["recovery_seconds"] >= H / 3
            ):
                interrupted = False
                state["recovery"] = 0
                state["recovery_seconds"] = 0
        else:
            state["recovery"] = 0
            state["recovery_seconds"] = 0

        state["last_interval_time"] = end

    if not interrupted and finite(prod) and prod > 0:
        reference = prod

    state.update(reference=reference, interrupted=interrupted)

    return {
        "predicted_increment": H * max(last / H, prod)
        if finite(last) and finite(prod)
        else math.nan,
        "previous_horizon_increment": last,
        "eligible_rate_intervals": len(recent),
        "rate_coverage_seconds": coverage,
        "cadence_adequate": covered and bool(history) and t - history[-1]["time"] <= H,
        "data_age_seconds": t - history[-1]["time"] if history else math.nan,
        "last_interval_seconds": newest[1] - newest[0] if newest else math.nan,
        "last_interval_reads_per_minute": newest[2] * 60 if newest else math.nan,
        "productive_reads_per_minute": reference * 60
        if finite(reference)
        else math.nan,
        "interrupted": interrupted,
        "recovery_intervals": state.get("recovery", 0),
        "recovery_seconds": state.get("recovery_seconds", 0),
    }


def interpolate_R(points, when, horizon):
    """Causal log-R interpolation; the caller supplies only already observed points."""
    for p in points:
        if p["time"] == when and finite(p["R"]) and p["R"] > 0:
            return p["R"]

    for left, right in pairwise(points):
        dt = right["time"] - left["time"]
        if (
            left["time"] < when < right["time"]
            and 0 < dt <= horizon
            and min(left["R"], right["R"]) > 0
        ):
            f = (when - left["time"]) / dt
            return math.exp((1 - f) * math.log(left["R"]) + f * math.log(right["R"]))

    return math.nan


def assess_contrast(
    state,
    definition,
    rows,
    contrast_samples,
    depths,
    live,
    t,
    completed,
    c,
    universe,
    sequenced=None,
    arrivals=None,
):
    cid = definition["contrast_id"]
    s = state.setdefault(
        cid, {"points": [], "samples": {}, "errors": [], "forecasts": [], "streak": 0}
    )
    H = c["horizon_seconds"]
    live_keys = sorted(x["key"] for x in contrast_samples if live.get(x["key"], False))

    if s.get("live_keys", live_keys) != live_keys:
        s["streak"] = 0
        s.pop("streak_started", None)

    s["live_keys"] = live_keys
    # Reject an obsolete snapshot before mutating any member history.
    for sample in contrast_samples:
        key = sample["key"]
        old_version = s["samples"].get(key, {}).get("version")
        version = (arrivals or {}).get(key, {}).get("version", t)

        if old_version is not None and version < old_version:
            previous = dict(s.get("latest_metric", {}))
            previous.update(
                contrast_id=cid,
                reason="stale_snapshot",
                eligible=False,
                informative=False,
            )
            return clean(previous), []

    points = s["points"]
    baseline = not points
    rates = []
    stale = False

    for sample in contrast_samples:
        key = sample["key"]
        d = depths[key]

        if not finite(d) or d < 0:
            raise ValueError("assigned-read depth must be finite and nonnegative")

        arrival = (arrivals or {}).get(key, {"time": t, "version": t})
        at = arrival["time"]
        version = arrival["version"]

        if not finite(at) or at < 0 or at > t or not finite(version) or version < 0:
            raise ValueError("invalid per-sample arrival timestamp or version")

        ss = s["samples"].setdefault(key, {"history": []})
        hist = ss["history"]
        old_version = ss.get("version")

        if old_version is not None and version < old_version:
            stale = True
            continue

        if old_version is not None and version == old_version:
            if d != hist[-1]["depth"] or at != hist[-1]["time"]:
                raise ValueError(
                    "sample version reused with different timestamp or depth"
                )
        else:
            if arrival.get("baseline") and hist:
                if d < hist[-1]["depth"] - 1e-6:
                    raise ValueError(
                        "decreasing assigned-read counter at resume baseline"
                    )
                hist.clear()
                s["streak"] = 0
                s.pop("streak_started", None)

            if hist and (at <= hist[-1]["time"] or d < hist[-1]["depth"] - 1e-6):
                raise ValueError(
                    "non-increasing observation time or decreasing assigned-read counter"
                )

            hist.append({"time": at, "depth": d})
            ss["version"] = version

            # Keep the boundary point for all rolling windows.
            while len(hist) > 2 and hist[1]["time"] < t - 3 * H:
                hist.pop(0)

        if live.get(key, False):
            rates.append(
                dict(
                    sample=key, **rate_state(hist, at, ss, c), data_available_seconds=at
                )
            )
        else:
            rates.append(
                {
                    "sample": key,
                    "predicted_increment": 0.0,
                    "interrupted": False,
                    "cadence_adequate": True,
                }
            )

    if stale or (
        s.get("last_evaluated_time") is not None and t < s["last_evaluated_time"]
    ):
        previous = dict(s.get("latest_metric", {}))
        previous.update(
            contrast_id=cid, reason="stale_snapshot", eligible=False, informative=False
        )
        return clean(previous), []

    previous_depths = s.get("evaluated_depths", {})
    fresh = all(depths[k] > previous_depths.get(k, depths[k]) + 1e-6 for k in live_keys)
    universe = sorted(universe)
    ids = [r["feature_id"] for r in rows]

    if len(ids) != len(set(ids)):
        raise ValueError("duplicate precision feature IDs")

    if s.get("universe") not in (None, universe):
        raise ValueError("R feature universe changed")

    s["universe"] = universe
    if set(ids) != set(universe):
        raise ValueError("precision export must include the fixed unfiltered universe")

    radius = [number(r["R"]) for r in rows]
    R, I = information_radius(radius, c["information_radius_scale"])
    D = sum(
        depths[x["key"]]
        for x in contrast_samples
        if (sequenced or live).get(x["key"], False)
    )

    finite_r = finite(R) and R > 0
    informative = bool((fresh or baseline) and finite_r and live_keys)

    # No duplicate R observation when a model is rerun for an unrelated arrival.
    if points and t <= points[-1]["time"]:
        informative = False
    curve_points = points + [{"time": t, "depth": D, "R": R}] if informative else points
    extra = sum(x["predicted_increment"] for x in rates)
    curve = curve_summary(
        [x["depth"] for x in curve_points],
        [x["R"] for x in curve_points],
        [x["time"] for x in curve_points],
        extra,
        c,
    )

    # Settle forecasts only when their deadlines are actually bracketed.
    matured = []
    for old in s["forecasts"]:
        if old.get("settled") or old["deadline"] > t:
            continue
        projected = interpolate_R(curve_points, old["deadline"], H)
        if finite(projected):
            gain = 1 - projected / old["R"]
            old.update(settled=True, actual_gain=gain)
            if finite(old["curve_gain"]):
                s["errors"].append(gain - old["curve_gain"])
            matured.append({"time": old["time"], "actual_gain_next_horizon": gain})
        elif curve_points and curve_points[-1]["time"] > old["deadline"]:
            old["settled"] = True  # An oversized gap cannot later become evidence.

    s["errors"] = s["errors"][-c["calibration_points"] :]
    errors = [max(0.0, x) for x in s["errors"]]
    margin = (
        max(c["minimum_forecast_margin"], quantile(errors, c["calibration_quantile"]))
        if errors
        else c["minimum_forecast_margin"]
    )
    conservative = (
        curve["curve_gain"] + margin if finite(curve["curve_gain"]) else math.nan
    )
    past_R = interpolate_R(curve_points, t - H, H)
    observed = 1 - R / past_R if finite(past_R) and finite_r else math.nan
    interrupted = any(x["interrupted"] for x in rates)
    cadence = all(x["cadence_adequate"] for x in rates) and all(
        completed - s["samples"][k]["history"][-1]["time"] <= H for k in live_keys
    )
    ages = [completed - s["samples"][k]["history"][-1]["time"] for k in live_keys]
    # Completion lag measures pipeline time, independently from the file cadence.
    lag = completed - t
    reason = (
        "no_live_samples"
        if not live_keys
        else "processing_lag"
        if lag > c["max_processing_lag_seconds"]
        else "interrupted"
        if interrupted
        else "insufficient_cadence"
        if not cadence
        else "waiting_for_fresh_data"
        if not fresh
        else "unavailable_precision"
        if not finite_r
        else "warming_up"
        if t < c["earliest_seconds"] or len(curve_points) < c["min_points"]
        else "unusable_curve"
        if not curve["curve_usable"]
        else "elbow_not_confirmed"
        if not curve["elbow_pass"]
        else "raw_R_unstable"
        if not finite(observed) or -observed > c["max_raw_increase_30min"]
        else "continue_observed_gain"
        if observed >= c["gain_cutoff"]
        else "continue_expected_gain"
        if not finite(conservative) or conservative >= c["gain_cutoff"]
        else "marginal_precision"
    )

    passed = reason == "marginal_precision" and informative
    if passed:
        if not s["streak"]:
            s["streak_started"] = t
        s["streak"] += 1
    elif reason != "waiting_for_fresh_data":
        s["streak"] = 0
        s.pop("streak_started", None)

    span = t - s["streak_started"] if s.get("streak_started") is not None else 0
    if informative:
        points.append({"time": t, "depth": D, "R": R})
        s["evaluated_depths"] = dict(depths)
        s["last_evaluated_time"] = t
        s["forecasts"].append(
            {"time": t, "deadline": t + H, "curve_gain": curve["curve_gain"], "R": R}
        )

    s["forecasts"] = [x for x in s["forecasts"] if not x.get("settled")]
    result = dict(
        contrast_id=cid,
        target_group=definition["target_group"],
        reference_group=definition["reference_group"],
        elapsed_seconds=t,
        analysis_completed_seconds=completed,
        processing_lag_seconds=lag,
        maximum_sample_data_age_seconds=max(ages, default=0),
        forecast_horizon_seconds=H,
        assigned_monitored_reads=D,
        R_information=R,
        information_score=I,
        fraction_R_le_05=sum(x <= 0.5 for x in radius) / len(radius),
        fraction_R_le_1=sum(x <= 1 for x in radius) / len(radius),
        finite_R_fraction=sum(math.isfinite(x) for x in radius) / len(radius),
        informative=informative,
        observed_gain_previous_horizon=observed,
        forecast_method=FORECAST_METHOD,
        forecast_margin=margin,
        predicted_conservative_gain=conservative,
        reason=reason,
        cadence_adequate=cadence,
        consecutive_confirmations=s["streak"],
        confirmation_span_seconds=span,
        eligible=(passed and s["streak"] >= c["confirmations"] and span >= 2 * H / 3)
        or (
            reason == "waiting_for_fresh_data"
            and all(depths[k] == previous_depths.get(k) for k in live_keys)
            and s.get("latest_metric", {}).get("eligible", False)
        ),
        **curve,
    )

    s["latest_metric"] = result
    s["matured"] = matured

    return clean(result), clean(rates)


def fingerprint(payload, c):
    fields = {k: payload[k] for k in ("contrasts", "identities")}
    fields["samples"] = [
        {k: v for k, v in s.items() if k not in ("live", "available_at", "version")}
        for s in payload["samples"]
    ]

    return hashlib.sha256(
        json.dumps(
            clean(dict(forecast_method=FORECAST_METHOD, config=c, **fields)),
            sort_keys=True,
        ).encode()
    ).hexdigest()


def evaluate(payload, state=None):
    c = validate_config(payload.get("config", {}))
    identity = fingerprint(payload, c)

    if state is not None and (
        state.get("schema_version") != STATE_SCHEMA_VERSION
        or state.get("forecast_method") != FORECAST_METHOD
    ):
        raise ValueError(
            "unsupported DEA state format or forecast method; use a fresh output directory"
        )

    if state is None or state.get("fingerprint") != identity:
        state = {
            "schema_version": STATE_SCHEMA_VERSION,
            "forecast_method": FORECAST_METHOD,
            "fingerprint": identity,
            "contrasts": {},
            "started_at": payload["started_at"],
        }

    if state["started_at"] != payload["started_at"]:
        raise ValueError("monitoring clock changed")
    t = payload["observed_at"] - state["started_at"]
    completed = payload["completed_at"] - state["started_at"]
    if t < 0 or completed < t or not math.isfinite(t + completed):
        raise ValueError("invalid monitoring clock")

    samples = payload["samples"]
    live = {x["key"]: x["live"] for x in samples}
    metrics = []
    rate_rows = []
    sequenced = {x["key"]: x.get("declared_live", x["live"]) for x in samples}

    for definition in payload["contrasts"]:
        cid = definition["contrast_id"]
        data = payload["precision"].get(cid)

        if data:
            result, rates = assess_contrast(
                state["contrasts"],
                definition,
                data["rows"],
                data["samples"],
                data["depths"],
                live,
                max(
                    (
                        payload.get("arrivals", {})
                        .get(x["key"], {})
                        .get("available_at", payload["observed_at"])
                        - state["started_at"]
                        for x in data["samples"]
                        if live.get(x["key"], False)
                    ),
                    default=t,
                ),
                completed,
                c,
                data["universe"],
                sequenced,
                {
                    k: {
                        "time": v["available_at"] - state["started_at"],
                        "version": v["version"],
                        "baseline": v.get("baseline", False),
                    }
                    for k, v in payload.get("arrivals", {}).items()
                },
            )

            metrics.append(result)
            rate_rows.extend(
                dict(contrast_id=cid, elapsed_seconds=result["elapsed_seconds"], **r)
                for r in rates
            )

        else:
            if cid in state["contrasts"]:
                state["contrasts"][cid]["streak"] = 0
            metrics.append(
                {
                    "contrast_id": cid,
                    "elapsed_seconds": t,
                    "reason": "analysis_deferred",
                    "eligible": False,
                    "consecutive_confirmations": 0,
                }
            )

    by = {r["contrast_id"]: r for r in metrics}
    sample_rows = []

    for sample in samples:
        required = [
            x["contrast_id"]
            for x in payload["contrasts"]
            if sample["group"] in (x["target_group"], x["reference_group"])
        ]

        eligible = bool(
            sample["live"] and required and all(by[x]["eligible"] for x in required)
        )

        sample_rows.append(
            {
                "group": sample["group"],
                "sample": sample["name"],
                "bam_dir": sample["bam_dir"],
                "protocol_run_id": sample.get("protocol_run_id") or "",
                "effectively_live": sample["live"],
                "method": "dea",
                "metric": "R_information",
                "required_contrasts": ",".join(required),
                "consecutive_confirmations": min(
                    [by[x]["consecutive_confirmations"] for x in required], default=0
                ),
                "eligible": eligible,
                "newly_eligible": eligible
                and not state.get("sample_eligible", {}).get(sample["key"], False),
                "reason": ",".join(sorted({by[x]["reason"] for x in required})),
                "behavior": payload["behavior"],
                "action_result": "none",
                "processing_lag_seconds": max(
                    (by[x].get("processing_lag_seconds", 0) for x in required),
                    default=0,
                ),
                "result_checked_at": payload["completed_at"],
            }
        )

    state["sample_eligible"] = {
        x["key"]: row["eligible"] for x, row in zip(samples, sample_rows)
    }

    history = state.setdefault("audit_history", [])
    for row in metrics:
        for matured in (
            state["contrasts"].get(row["contrast_id"], {}).pop("matured", [])
        ):
            for old in reversed(history):
                if (
                    old["contrast_id"] == row["contrast_id"]
                    and old["elapsed_seconds"] == matured["time"]
                    and old.get("informative")
                ):
                    old["actual_gain_next_horizon"] = matured[
                        "actual_gain_next_horizon"
                    ]
                    break
                
    history.extend(metrics)
    state.setdefault("rate_history", []).extend(rate_rows)
    state["config"] = c
    return clean(state), metrics, sample_rows, rate_rows
