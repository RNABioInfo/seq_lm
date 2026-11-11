"""ICA adapter regressions for the shared precision stopping engine."""
import csv
import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from precision_termination.engine import evaluate
from precision_termination.ica import component_radius, load_ica_precision


def test_ica_radius_uses_larger_side_and_accounts_for_log_base_and_weight():
    interval = {"status": "tested", "activity_difference": "2", "ci_lower": "0", "ci_upper": "3"}
    assert component_radius(interval, 0.25, 16) == 2
    assert component_radius({"status": "zero_variance"}, 0.25, 16) is None
    with pytest.raises(ValueError, match="Invalid tested ICA"):
        component_radius(interval | {"ci_upper": "1"}, 0.25, 16)


def test_unestimable_selected_ica_component_stays_in_denominator(tmp_path):
    definitions = [{"contrast_id": "pair", "target_group": "case", "reference_group": "control"}]
    samples = [
        {"key": f"{group}\trep{i}", "group": group, "name": f"rep{i}", "live": group == "case", "bam_dir": ""}
        for group in ("control", "case") for i in range(1, 3 if group == "control" else 4)
    ]
    provenance = {
        "schema_version": 2,
        "model": {"components": ["x", "y"], "stopping_components": ["x"],
                  "component_max_abs_weights": {"x": 1, "y": 1}, "hashes": {}},
        "settings": {"control_group": "control", "condition_groups": ["case", "control"],
                     "log_base": 2, "min_read_count": 10000, "min_replicate_sample_count": 2},
        "control_sample_ids": ["fixed-1", "fixed-2"],
        "software": {},
    }
    (tmp_path / "provenance.json").write_text(json.dumps(provenance))
    (tmp_path / "status.json").write_text('{"status":"ready"}')
    with (tmp_path / "sample_metadata.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["alias", "group", "assigned_abundance"])
        for sample in samples:
            depth = 1200 if sample["name"] == "rep3" else 10000
            writer.writerow([sample["name"], sample["group"], depth])
    (tmp_path / "differential_activity.tsv").write_text(
        "component_id\ttarget_group\tcontrol_group\tstatus\tactivity_difference\tci_lower\tci_upper\n"
        "x\tcase\tcontrol\tzero_variance\t0\t\t\n"
        "y\tcase\tcontrol\ttested\t1\t0\t2\n"
    )

    precision, identities = load_ica_precision(tmp_path, definitions, samples)
    assert precision["pair"]["universe"] == ["x"]
    assert precision["pair"]["depths"]["case\trep3"] == 1200
    payload = {
        "method": "ica", "samples": samples, "contrasts": definitions, "precision": precision,
        "identities": identities, "started_at": 0, "observed_at": 600, "completed_at": 600,
        "behavior": "log",
    }
    _, metrics, rows, _ = evaluate(payload)
    assert metrics[0]["information_score"] == 0
    assert metrics[0]["uncertainty_feature_count"] == 1
    assert metrics[0]["uncertainty_scope"] == "selected_features"
    assert all(row["method"] == "ica" and not row["eligible"] for row in rows)

    del provenance["model"]["stopping_components"]
    (tmp_path / "provenance.json").write_text(json.dumps(provenance))
    precision, identities = load_ica_precision(tmp_path, definitions, samples)
    _, metrics, _, _ = evaluate(payload | {"precision": precision, "identities": identities})
    assert metrics[0]["information_score"] == 0.25
    assert metrics[0]["R_information"] == pytest.approx(math.sqrt(3))
