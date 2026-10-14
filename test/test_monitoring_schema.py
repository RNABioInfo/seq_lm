"""Schema-level validation of selectable sequencing stopping methods."""

import json
from pathlib import Path

import pytest
from jsonschema import Draft7Validator


SCHEMA = json.loads((Path(__file__).resolve().parents[1] / "nextflow_schema.json").read_text())
VALIDATOR = Draft7Validator(SCHEMA)
BASE = {
    "sample_sheet": "samples.csv",
    "monitoring_behavior": "log", "stopping_method": "ica",
    "ica_analysis": True, "ica_matrix": "weights.csv",
    "differential_expression": False, "gene_set_enrichment": False,
    "reference_genome": "genome.fa", "reference_annotation": "annotation.gtf",
    "live_analysis": True, "num_stable_batches": 3,
    "ica_max_fold_change": 1.05, "ica_min_read_growth_fraction": 0.1,
}


def test_ica_stopping_does_not_require_dea():
    VALIDATOR.validate(BASE)
    VALIDATOR.validate(BASE | {"monitoring_behavior": "terminate"})
    VALIDATOR.validate(BASE | {"ica_max_fold_change": 1})
    VALIDATOR.validate(BASE | {"stopping_method": "dea", "differential_expression": True})
    Draft7Validator.check_schema(SCHEMA)


@pytest.mark.parametrize("values", [
    {"stopping_method": "unknown"},
    {"stopping_method": "dea"},
    {"ica_analysis": False},
    {"ica_matrix": ""},
    {"ica_max_fold_change": 0.95},
    {"ica_max_fold_change": "1.05"},
    {"ica_min_read_growth_fraction": 0},
    {"ica_min_read_growth_fraction": -0.1},
    {"num_stable_batches": 0},
    {"num_stable_batches": 1.5},
    {"monitoring_behavior": "terminate", "live_analysis": False},
])
def test_invalid_stopping_parameters_are_rejected(values):
    assert list(VALIDATOR.iter_errors(BASE | values))


def test_ica_model_required_only_when_monitoring_enabled():
    without_model = {key: value for key, value in BASE.items() if key != "ica_matrix"}
    assert list(VALIDATOR.iter_errors(without_model))
    VALIDATOR.validate(without_model | {"monitoring_behavior": "disabled", "ica_analysis": False})


def test_minknow_completion_options():
    VALIDATOR.validate(BASE | {"detect_minknow_completion": True, "minknow_completion_settle_seconds": 30})
    VALIDATOR.validate(BASE | {"detect_minknow_completion": False, "minknow_completion_settle_seconds": 1})
    for value in [0, -1, 1.5, "30", None, True]:
        assert list(VALIDATOR.iter_errors(BASE | {"minknow_completion_settle_seconds": value}))
    assert list(VALIDATOR.iter_errors(BASE | {"detect_minknow_completion": "true"}))
