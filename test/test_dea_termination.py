"""Causal forecasts, interruption protection, current-format state and cohorts."""

import copy
import json
import math
import subprocess
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
    write_tsv,
)
from dea_termination.cli import load_precision
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


def export_precision(root, p, selection=None):
    """Write the same full-feature artifacts consumed by the monitoring CLI."""
    root.mkdir(parents=True, exist_ok=True)
    for cid, data in p["precision"].items():
        folder = root / cid
        folder.mkdir(exist_ok=True)
        write_tsv(folder / "precision.tsv", data["rows"])
        write_tsv(folder / "precision_samples.tsv", [
            row | {"assigned_reads": data["depths"][row["key"]]}
            for row in data["samples"]
        ])
    if selection is not None:
        write_tsv(root / "dea_stopping_features.tsv", [
            {"feature_id": feature} for feature in selection
        ], ["feature_id"])


def test_selected_score_ignores_unselected_radii_and_keeps_library_depths(tmp_path):
    p = payload()
    p["identities"]["dea_stopping_genes"] = "gmt-hash"
    p["precision"]["T1"]["rows"][0]["R"] = 0.5
    p["precision"]["T1"]["rows"][1]["R"] = 2
    export_precision(tmp_path, p, ["g0", "g1"])
    full = load_precision(tmp_path, p["contrasts"])
    p["precision"] = load_precision(tmp_path, p["contrasts"], selected=True)
    _, metrics, _, _ = evaluate(p)
    expected, info = information_radius([0.5, 2])
    assert metrics[0]["R_information"] == pytest.approx(expected)
    assert metrics[0]["information_score"] == pytest.approx(info)
    assert metrics[0]["uncertainty_scope"] == "selected_features"
    assert metrics[0]["uncertainty_feature_count"] == 2
    assert metrics[0]["assigned_monitored_reads"] == 3000
    assert p["precision"]["T1"]["depths"] == full["T1"]["depths"]
    full["T1"]["rows"][2]["R"] = 1000
    p["precision"] = full
    export_precision(tmp_path, p, ["g1", "g0"])
    p["precision"] = load_precision(tmp_path, p["contrasts"], selected=True)
    _, changed, _, _ = evaluate(p)
    assert changed[0]["R_information"] == metrics[0]["R_information"]
    p["precision"] = load_precision(tmp_path, p["contrasts"])
    del p["identities"]["dea_stopping_genes"]
    _, default, _, _ = evaluate(p)
    assert default[0]["uncertainty_feature_count"] == 10
    assert default[0]["uncertainty_scope"] == "all_features"
    assert default[0]["R_information"] != expected


@pytest.mark.parametrize("radii", [(1, "NA"), (1, "Inf"), ("NA", "Inf")])
def test_unestimable_selected_features_stay_in_denominator(tmp_path, radii):
    p = payload()
    p["identities"]["dea_stopping_genes"] = "gmt-hash"
    for row, radius in zip(p["precision"]["T1"]["rows"], radii):
        row["R"] = radius
    export_precision(tmp_path, p, ["g0", "g1"])
    p["precision"] = load_precision(tmp_path, p["contrasts"], selected=True)
    _, metrics, samples, _ = evaluate(p)
    expected = 0.25 if radii[0] == 1 else 0
    assert metrics[0]["information_score"] == expected
    assert metrics[0]["finite_R_fraction"] == (0.5 if expected else 0)
    assert metrics[0]["uncertainty_feature_count"] == 2
    if not expected:
        assert metrics[0]["R_information"] is None
        assert not metrics[0]["eligible"]
        assert not any(row["eligible"] for row in samples)


@pytest.mark.parametrize("selection,error", [
    (None, "requires dea_stopping_features.tsv"),
    ([], "selection is empty or invalid"),
    (["g0", "g0"], "Duplicate DEA stopping feature IDs"),
    (["absent"], "absent from precision export"),
])
def test_selection_export_is_required_for_fitted_snapshots(tmp_path, selection, error):
    p = payload()
    export_precision(tmp_path, p, selection)
    with pytest.raises(ValueError, match=error):
        load_precision(tmp_path, p["contrasts"], selected=True)


def test_selection_change_resets_state_and_unexpected_drift_fails(tmp_path):
    p = payload()
    p["identities"]["dea_stopping_genes"] = "first-content-hash"
    export_precision(tmp_path, p, ["g0", "g1"])
    p["precision"] = load_precision(tmp_path, p["contrasts"], selected=True)
    state, *_ = evaluate(p)
    state["contrasts"]["T1"]["streak"] = 2
    later = payload(600, depth=2000)
    later["identities"]["dea_stopping_genes"] = "first-content-hash"
    export_precision(tmp_path, later, ["g2", "g3"])
    later["precision"] = load_precision(tmp_path, later["contrasts"], selected=True)
    with pytest.raises(ValueError, match="R feature universe changed"):
        evaluate(later, copy.deepcopy(state))
    later["identities"]["dea_stopping_genes"] = "second-content-hash"
    reset, metrics, *_ = evaluate(later, state)
    assert reset["fingerprint"] != state["fingerprint"]
    assert metrics[0]["consecutive_confirmations"] == 0
    assert len(reset["contrasts"]["T1"]["points"]) == 1


def test_selected_cli_audits_all_contrasts_and_handles_deferral(tmp_path):
    p = payload()
    p["identities"]["dea_stopping_genes"] = "gmt-content-hash"
    other = {"contrast_id": "other", "target_group": "T1", "reference_group": "T0"}
    p["contrasts"].append(other)
    p["precision"]["other"] = copy.deepcopy(p["precision"]["T1"])
    root = tmp_path / "fit"
    export_precision(root, p, ["g0", "g1"])
    p.pop("precision")
    p["results"] = str(root)
    request = tmp_path / "request.json"
    request.write_text(json.dumps(p))
    script = Path(__file__).parents[1] / "bin/dea-termination"
    for deferred in (False, True):
        if deferred:
            for cid in ("T1", "other"):
                (root / cid / "precision.tsv").unlink()
            (root / "dea_stopping_features.tsv").unlink()
        output = tmp_path / ("deferred" if deferred else "selected")
        result = subprocess.run([
            sys.executable, str(script), "--input", str(request),
            "--state", str(tmp_path / "no-previous-state.json"), "--output", str(output),
        ], capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        audit = json.loads((output / "audit.json").read_text())
        assert len(audit["contrasts"]) == 2
        for row in audit["contrasts"]:
            assert row["uncertainty_scope"] == "selected_features"
            assert row["uncertainty_feature_count"] == (None if deferred else 2)
            if deferred:
                assert row["reason"] == "analysis_deferred" and not row["eligible"]


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
