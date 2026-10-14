"""ICA audit validation without requiring HTML plotting dependencies."""

from types import SimpleNamespace

import pandas as pd
import pytest

from workflow_glue.qc_report_types.ica_monitoring import load_ica_monitoring_results


def fixture_audit(tmp_path):
    data = SimpleNamespace(
        samples=pd.DataFrame({"group": ["control", "treated"], "alias": ["rep1", "rep1"]}),
        provenance={"batch_index": 7, "analysis_index": 11},
    )
    rows = pd.DataFrame({
        "group": ["control", "treated"], "sample": ["rep1", "rep1"],
        "batch_index": 7, "analysis_index": 11, "method": "ica", "metric": "ica_monitoring",
        "behavior": "log", "required_components": 84, "passing_components": [84, 83],
        "blocking_components": ["", "stress"], "consecutive_stable_batches": [3, 0],
        "eligible": ["true", "false"], "reason": ["all_components_converged", "component_change_exceeds_tolerance"],
        "action_result": ["waiting_for_run_stability", "none"],
    })
    path = tmp_path / "sample_stability.tsv"
    rows.to_csv(path, sep="\t", index=False)
    return data, rows, path


def test_repeated_aliases_and_component_blockers(tmp_path):
    data, _, path = fixture_audit(tmp_path)
    result = load_ica_monitoring_results(path, data, "log")
    assert result["Sample"].tolist() == ["rep1", "rep1"]
    assert result["Passing components"].tolist() == ["84", "83"]
    assert result["Blocking components"].tolist() == ["", "stress"]


@pytest.mark.parametrize("field,value,message", [
    ("analysis_index", 12, "analysis_index does not match"),
    ("batch_index", 8, "batch_index does not match"),
    ("method", "dea", "Unexpected ICA monitoring metric"),
    ("metric", "unknown", "Unexpected ICA monitoring metric"),
    ("behavior", "terminate", "behavior or eligibility is invalid"),
    ("eligible", "yes", "behavior or eligibility is invalid"),
    ("group", "other", "sample identities do not match"),
])
def test_audit_must_match_the_report(tmp_path, field, value, message):
    data, rows, path = fixture_audit(tmp_path)
    rows[field] = value
    rows.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match=message):
        load_ica_monitoring_results(path, data, "log")


def test_duplicate_and_missing_audit_rows(tmp_path):
    data, rows, path = fixture_audit(tmp_path)
    for invalid in (rows.iloc[[0]], pd.concat([rows, rows.iloc[[0]]])):
        invalid.to_csv(path, sep="\t", index=False)
        with pytest.raises(ValueError, match="sample identities do not match"):
            load_ica_monitoring_results(path, data, "log")
    rows.drop(columns=["passing_components"]).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing required columns"):
        load_ica_monitoring_results(path, data, "log")


def test_full_audit_omits_unselected_groups_from_ica_tab(tmp_path):
    data, rows, path = fixture_audit(tmp_path)
    extra = rows.iloc[[1]].copy()
    extra["group"] = "QC only"
    extra["eligible"] = "false"
    audit = pd.concat([rows, extra], ignore_index=True)
    report_samples = pd.DataFrame({"Name": audit["sample"], "Group": audit["group"]})
    audit.to_csv(path, sep="\t", index=False)
    table = load_ica_monitoring_results(path, data, "log", report_samples)
    assert table["Group"].tolist() == ["control", "treated"]
    audit.loc[2, "eligible"] = "true"
    audit.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="outside the ICA cohort"):
        load_ica_monitoring_results(path, data, "log", report_samples)
