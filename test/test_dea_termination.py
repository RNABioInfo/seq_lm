"""Causal forecasts, interruption protection, current-format state and cohorts."""

import copy
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "bin"))
from dea_termination.engine import (
    DEFAULTS,
    atomic_json,
    evaluate,
    interpolate_R,
    rate_state,
    validate_config,
    weighted_quantile,
)
from dea_termination.forecast import curve_summary, fit_curve, information_radius


def payload(t=0, R=1.0, depth=1000.0, extra=None):
    samples = [
        {
            "key": f"{g}/{i}",
            "name": str(i),
            "group": g,
            "bam_dir": f"/tmp/{g}_{i}",
            "live": g == "T1",
            "protocol_run_id": "shared",
        }
        for g in ("T0", "T1")
        for i in range(3)
    ]
    contrast_samples = [{"key": x["key"], "group": x["group"]} for x in samples]
    rows = [
        dict(feature_id=f"g{i}", R=R)
        for i in range(10)
    ]
    return {
        "config": extra or {},
        "samples": samples,
        "contrasts": [
            {"contrast_id": "T1", "target_group": "T1", "reference_group": "T0"}
        ],
        "identities": {"reference": "fixed"},
        "started_at": 100000.0,
        "observed_at": 100000.0 + t,
        "completed_at": 100000.0 + t,
        "behavior": "log",
        "precision": {
            "T1": {
                "rows": rows,
                "samples": contrast_samples,
                "depths": {
                    x["key"]: 100000 if x["group"] == "T0" else depth for x in samples
                },
                "universe": [r["feature_id"] for r in rows],
            }
        },
    }


def test_aggregate_scale_and_missing_denominator():
    r, i = information_radius([1, 1])
    assert r == 1 and i == 0.5
    r, i = information_radius([0, np.inf])
    assert r == 1 and i == 0.5
    assert information_radius([np.nan, np.inf])[1] == 0
    assert information_radius([1], 0.5)[1] == pytest.approx(0.2)
    with pytest.raises(ValueError):
        information_radius([-1])


def test_nonnegative_curve_recovers_known_model():
    depths = np.arange(1, 30) * 10000
    a, b = fit_curve(depths, np.sqrt(0.25 + 20000 / depths))
    assert a == pytest.approx(0.25)
    assert b == pytest.approx(0.02)
    c = curve_summary(
        np.arange(1, 30) * 10000,
        np.sqrt(0.25 + 20000 / (np.arange(1, 30) * 10000)),
        np.arange(1, 30) * 600,
        10000,
        DEFAULTS,
    )
    assert c["curve_usable"] and c["curve_gain"] > 0
    assert not curve_summary([100] * 6, [1] * 6, list(range(6)), 0, DEFAULTS)[
        "curve_usable"
    ]


def test_long_pause_does_not_erase_reference_and_recovery_requires_two_bins():
    state = {}
    history = []
    depth = 0
    for i in range(20):
        depth += 100 if i < 6 or i >= 16 else 0
        history.append({"time": i * 600, "depth": depth})
        rate = rate_state(history[-12:], i * 600, state, DEFAULTS)
        if i in (9, 15):
            assert rate["interrupted"] and rate[
                "productive_reads_per_minute"
            ] == pytest.approx(10)
        if i == 16:
            assert rate["interrupted"] and rate["recovery_intervals"] == 1
        if i == 17:
            assert not rate["interrupted"]


@pytest.mark.parametrize("spacing", [6, 12, 24])
def test_pore_check_spacing_is_not_an_input(spacing):
    state = {}
    history = []
    depth = 0
    for i in range(40):
        depth += 0 if i and i % spacing == 0 else 100
        history.append({"time": i * 600, "depth": depth})
        rate = rate_state(history[-12:], i * 600, state, DEFAULTS)
        if i >= 3:
            assert rate["predicted_increment"] >= 200


def test_zero_yield_cannot_create_stop_or_informative_observations():
    state = None
    for i in range(15):
        state, metrics, samples, _ = evaluate(payload(i * 600), state)
        assert metrics[0]["informative"] == (i == 0)
        assert not samples[-1]["eligible"]
    assert len(state["contrasts"]["T1"]["points"]) == 1


def test_first_hour_warmup_and_lag_cannot_stop():
    state = None
    for i in range(7):
        p = payload(i * 600, R=np.sqrt(0.25 + 1 / (i + 1)), depth=1000 * (i + 1))
        state, metrics, _, _ = evaluate(p, state)
        if i < 6:
            assert (
                metrics[0]["reason"] in ("warming_up", "insufficient_cadence")
                or metrics[0]["reason"] == "waiting_for_fresh_data"
            )
    p = payload(4200, depth=8000)
    p["completed_at"] += 601
    state, metrics, _, _ = evaluate(p, state)
    assert metrics[0]["reason"] == "processing_lag"
    assert not metrics[0]["eligible"]


def test_gap_resets_confirmation_and_is_not_a_zero_bin():
    state, *_ = evaluate(payload())
    state, *_ = evaluate(payload(600, depth=2000), state)
    state["contrasts"]["T1"]["streak"] = 2
    state, metrics, _, _ = evaluate(payload(3000, depth=3000), state)
    assert metrics[0]["reason"] == "insufficient_cadence"
    assert metrics[0]["consecutive_confirmations"] == 0
    h = [{"time": 0, "depth": 0}, {"time": 2400, "depth": 100}]
    assert rate_state(h, 2400, {}, DEFAULTS)["eligible_rate_intervals"] == 0


def test_config_and_identity_changes_reset_state_and_legacy_is_rejected():
    state, *_ = evaluate(payload())
    state["contrasts"]["T1"]["streak"] = 2
    p = payload(600, depth=2000, extra={"gain_cutoff": 0.02})
    state, metrics, _, _ = evaluate(p, state)
    assert metrics[0]["consecutive_confirmations"] == 0
    with pytest.raises(ValueError, match="unsupported"):
        evaluate(payload(), {"schema_version": 2})
    _, repeated, *_ = evaluate(
        payload(600, depth=2000, extra={"gain_cutoff": 0.02}), state
    )
    assert not repeated[0]["informative"]
    with pytest.raises(ValueError, match="version reused"):
        evaluate(payload(600, depth=2100, extra={"gain_cutoff": 0.02}), state)


def test_feature_universe_drift_and_counter_regression_are_errors():
    state, *_ = evaluate(payload())
    p = payload(600, depth=500)
    with pytest.raises(ValueError, match="decreasing"):
        evaluate(p, copy.deepcopy(state))
    p = payload(600, depth=2000)
    p["precision"]["T1"]["universe"] = ["other"]
    with pytest.raises(ValueError, match="universe"):
        evaluate(p, copy.deepcopy(state))


def test_independent_contrast_does_not_require_unrelated_precision():
    p = payload()
    p["contrasts"].append(
        {"contrast_id": "T4", "target_group": "T4", "reference_group": "T0"}
    )
    _state, metrics, samples, _ = evaluate(p)
    assert metrics[1]["reason"] == "analysis_deferred"
    # T1 sample eligibility depends only on its own contrast.
    assert samples[-1]["required_contrasts"] == "T1"
    assert (
        samples[0]["required_contrasts"] == "T1,T4"
        and not samples[0]["effectively_live"]
    )


def test_state_roundtrip_is_strict_json_and_atomic(tmp_path):
    state, *_ = evaluate(payload())
    path = tmp_path / "state.json"
    atomic_json(path, state)
    reread = json.loads(path.read_text())
    assert reread == state
    assert not list(tmp_path.glob(".state*"))
    _state2, metrics, *_ = evaluate(payload(600, depth=2000), reread)
    assert metrics[0]["elapsed_seconds"] == 600


@pytest.mark.parametrize(
    "config",
    [
        {"min_points": 5},
        {"earliest_seconds": 3599},
        {"gain_cutoff": 1},
        {"gain_cutoff": 0},
        {"information_radius_scale": 0},
        {"checkpoint_seconds": 300},
        {"horizon_seconds": 0},
    ],
)
def test_invalid_config(config):
    with pytest.raises(ValueError):
        validate_config(config)


def test_threshold_equality_and_three_confirmations(monkeypatch):
    from dea_termination import engine

    predicted = 0.008

    def curve(*args):
        return {
            "curve_gain": predicted,
            "curve_usable": True,
            "elbow_pass": True,
            "finite_points": 10,
        }

    monkeypatch.setattr(engine, "curve_summary", curve)
    state = None
    for i in range(7):
        state, metrics, _, _ = evaluate(payload(i * 600, depth=10000 * (i + 1)), state)
    assert metrics[0]["predicted_conservative_gain"] == 0.01
    assert (
        not metrics[0]["eligible"] and metrics[0]["reason"] == "continue_expected_gain"
    )
    predicted = 0.007
    for i in range(7, 10):
        state, metrics, _, _ = evaluate(payload(i * 600, depth=10000 * (i + 1)), state)
        assert metrics[0]["eligible"] == (i == 9)


def test_natural_completion_changes_only_affected_live_cohort():
    p = payload()
    p["samples"][3]["declared_live"] = True
    state, *_ = evaluate(p)
    p = payload(600, depth=2000)
    p["samples"][3].update(live=False, declared_live=True)
    state, metrics, _, _ = evaluate(p, state)
    assert len(state["contrasts"]["T1"]["points"]) == 2
    assert metrics[0]["assigned_monitored_reads"] == 6000
    assert state["contrasts"]["T1"]["live_keys"] == ["T1/1", "T1/2"]
    assert metrics[0]["consecutive_confirmations"] == 0


def test_feature_order_is_not_a_change_in_the_feature_universe():
    state, *_ = evaluate(payload())
    p = payload(600, depth=2000)
    p["precision"]["T1"]["rows"].reverse()
    p["precision"]["T1"]["universe"].reverse()
    _, metrics, _, _ = evaluate(p, state)
    assert metrics[0]["R_information"] == 1


def stamp(p, times=None, versions=None):
    """Give every sample explicit input provenance, independent of analysis time."""
    t = p["observed_at"]
    p["arrivals"] = {
        x["key"]: {
            "available_at": (times or {}).get(x["key"], t),
            "version": (versions or {}).get(x["key"], t),
        }
        for x in p["samples"]
    }
    return p


def test_duration_weighted_quantile_and_clipped_windows():
    assert weighted_quantile([(100, 1), (2, 99)], 0.75) == 2
    h = [{"time": t, "depth": 2 * t} for t in (0, 500, 1100, 1700, 2300)]
    rate = rate_state(h, 2300, {}, DEFAULTS)
    assert rate["previous_horizon_increment"] == pytest.approx(3600)
    assert rate["predicted_increment"] == pytest.approx(3600)
    assert rate["cadence_adequate"]
    assert rate["last_interval_reads_per_minute"] == 120


@pytest.mark.parametrize(
    "step,horizon,adequate",
    [(60, 1800, True), (600, 1800, True), (3600, 1800, False), (3600, 3600, True)],
)
def test_output_cadence_and_configurable_horizon(step, horizon, adequate):
    c = validate_config({"horizon_seconds": horizon})
    h = [
        {"time": i * step, "depth": i * step * 2}
        for i in range(max(7, 3 * horizon // step + 1))
    ]
    result = rate_state(h, h[-1]["time"], {}, c)
    assert result["cadence_adequate"] == adequate
    assert math.isfinite(result["predicted_increment"]) == adequate


def test_log_R_interpolation_does_not_extrapolate_or_bridge_long_gaps():
    points = [{"time": 0, "R": 4}, {"time": 1200, "R": 1}]
    assert interpolate_R(points, 600, 1800) == pytest.approx(2)
    assert math.isnan(interpolate_R(points, -1, 1800))
    assert math.isnan(interpolate_R(points, 1201, 1800))
    assert math.isnan(interpolate_R(points, 600, 600))


def test_stale_snapshot_cannot_rewind_sample_history():
    state, *_ = evaluate(stamp(payload()))
    state, *_ = evaluate(stamp(payload(600, depth=2000)), state)
    before = copy.deepcopy(state["contrasts"]["T1"]["samples"])
    state, metrics, *_ = evaluate(stamp(payload(300, depth=1500)), state)
    assert metrics[0]["reason"] == "stale_snapshot"
    assert state["contrasts"]["T1"]["samples"] == before


def test_partial_member_updates_accumulate_until_every_live_member_is_fresh():
    p = stamp(payload())
    state, *_ = evaluate(p)
    versions = {x["key"]: 100000 for x in p["samples"]}
    times = dict(versions)
    for i, key in enumerate(("T1/0", "T1/1", "T1/2"), 1):
        p = payload(i * 100, depth=1000)
        for advanced in ("T1/0", "T1/1", "T1/2")[:i]:
            p["precision"]["T1"]["depths"][advanced] = 2000
        versions[key] = i
        times[key] = 100000 + i * 100
        # Versions must be monotonic independently of timestamp units.
        versions[key] = 100000 + i
        state, metrics, *_ = evaluate(stamp(p, times, versions), state)
        assert metrics[0]["informative"] == (i == 3)
    assert len(state["contrasts"]["T1"]["points"]) == 2


def test_hourly_output_can_qualify_only_with_hourly_horizon(monkeypatch):
    from dea_termination import engine

    monkeypatch.setattr(
        engine,
        "curve_summary",
        lambda *a: {"curve_gain": 0.001, "curve_usable": True, "elbow_pass": True},
    )
    for horizon in (1800, 3600):
        state = None
        for i in range(10):
            state, metrics, *_ = evaluate(
                payload(
                    i * 3600, depth=(i + 1) * 10000, extra={"horizon_seconds": horizon}
                ),
                state,
            )
        assert metrics[0]["eligible"] == (horizon == 3600)


def test_burst_delivery_cannot_shortcut_confirmation_span(monkeypatch):
    from dea_termination import engine

    monkeypatch.setattr(
        engine,
        "curve_summary",
        lambda *a: {"curve_gain": 0.001, "curve_usable": True, "elbow_pass": True},
    )
    state = None
    for t in range(0, 3601, 600):
        state, metrics, *_ = evaluate(payload(t, depth=10000 + t), state)
    assert metrics[0]["consecutive_confirmations"] == 1
    for t in (3601, 3602, 3603):
        state, metrics, *_ = evaluate(payload(t, depth=10000 + t), state)
        assert not metrics[0]["eligible"]
    state, metrics, *_ = evaluate(payload(4800, depth=14800), state)
    assert metrics[0]["eligible"]


def test_forecasts_mature_once_after_deadline_is_bracketed():
    state = None
    for t in (0, 500, 1100, 1700, 2300, 2900, 3500):
        state, _metrics, *_ = evaluate(
            payload(t, R=math.sqrt(0.25 + 1000 / (t + 1000)), depth=t + 1000), state
        )
    matured = [r for r in state["audit_history"] if "actual_gain_next_horizon" in r]
    assert matured
    first = matured[0]
    assert first["elapsed_seconds"] + 1800 <= 3500
    before = copy.deepcopy(state["contrasts"]["T1"]["errors"])
    state, *_ = evaluate(
        payload(3500, R=math.sqrt(0.25 + 1000 / 4500), depth=4500), state
    )
    assert state["contrasts"]["T1"]["errors"] == before


def test_curve_only_margin_uses_completed_curve_errors(monkeypatch):
    from dea_termination import engine

    monkeypatch.setattr(
        engine,
        "curve_summary",
        lambda *args: {
            "curve_gain": 0.001, "curve_usable": True, "elbow_pass": True,
        },
    )
    state = None
    for t in (0, 600, 1200):
        state, metrics, *_ = evaluate(
            payload(t, R=math.exp(-t / 10000), depth=1000 + t), state
        )
        assert metrics[0]["forecast_margin"] == 0.002
    state, metrics, *_ = evaluate(
        payload(1800, R=math.exp(-0.18), depth=2800), state
    )
    actual_gain = 1 - math.exp(-0.18)
    assert metrics[0]["forecast_margin"] == pytest.approx(actual_gain - 0.001)
    assert metrics[0]["predicted_conservative_gain"] == pytest.approx(actual_gain)
    assert metrics[0]["forecast_method"] == "depth_curve"
    assert "model_gain" not in metrics[0]
    assert state["schema_version"] == 3
    assert state["forecast_method"] == "depth_curve"
    assert all("curve_gain" in f and "base" not in f
               for f in state["contrasts"]["T1"]["forecasts"])


def test_curve_forecast_needs_only_radii_and_assigned_depths():
    state = None
    for i in range(12):
        depth = (i + 1) * 10000
        state, metrics, *_ = evaluate(
            payload(i * 600, R=math.sqrt(0.25 + 20000 / depth), depth=depth), state
        )
    metric = metrics[0]
    assert math.isfinite(metric["curve_gain"])
    assert metric["predicted_conservative_gain"] == pytest.approx(
        metric["curve_gain"] + metric["forecast_margin"]
    )


def test_resume_baseline_does_not_invent_a_zero_production_interval():
    state = None
    for t in range(0, 3601, 600):
        state, *_ = evaluate(stamp(payload(t, depth=1000 + t)), state)
    p = stamp(payload(10000, depth=4600))
    for arrival in p["arrivals"].values():
        arrival["baseline"] = True
    state, metrics, _rates, _ = evaluate(p, state)
    hist = state["contrasts"]["T1"]["samples"]["T1/0"]["history"]
    assert hist == [{"time": 10000, "depth": 4600}]
    assert not metrics[0]["eligible"]
    assert metrics[0]["consecutive_confirmations"] == 0


def test_old_timed_state_and_removed_settings_are_rejected():
    with pytest.raises(ValueError, match="unsupported"):
        evaluate(payload(), {"schema_version": 1})
    for name in ("checkpoint_seconds", "recovery_bins"):
        with pytest.raises(ValueError, match="unknown"):
            validate_config({name: 600})


def test_resume_baseline_rejects_counter_regression():
    state, *_ = evaluate(stamp(payload(depth=2000)))
    p = stamp(payload(10000, depth=1000))
    for arrival in p["arrivals"].values():
        arrival["baseline"] = True
    with pytest.raises(ValueError, match="decreasing"):
        evaluate(p, state)
