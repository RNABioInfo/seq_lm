"""Validate ICA monitoring audits independently of plotting dependencies."""

from __future__ import annotations

from typing import TYPE_CHECKING
import pandas as pd

if TYPE_CHECKING:
    from .imodulon_plots import IModulonResult


def load_ica_monitoring_results(path, data: IModulonResult, behavior: str, report_samples=None):
    """Validate the matching sample audit, including aliases repeated across groups."""
    audit = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"analysis_index", "batch_index", "group", "sample", "metric", "method",
                "required_components", "passing_components", "blocking_components",
                "consecutive_stable_batches", "eligible", "reason", "behavior", "action_result"}
    if not required.issubset(audit.columns):
        raise ValueError("ICA monitoring audit is missing required columns")
    identities = list(zip(audit.group, audit["sample"]))
    expected = list(zip(data.samples.group, data.samples.alias))
    if report_samples is not None:
        qc_pairs = list(zip(report_samples["Group"], report_samples["Name"]))
        if len(set(identities)) != len(identities) or set(identities) != set(qc_pairs):
            raise ValueError("ICA monitoring audit does not match QC sample identities")
        selected = pd.Series([pair in set(expected) for pair in identities], index=audit.index)
        if audit.loc[~selected, "eligible"].ne("false").any():
            raise ValueError("Samples outside the ICA cohort cannot be eligible")
        audit = audit.loc[selected].copy()
        identities = list(zip(audit.group, audit["sample"]))
    if len(set(identities)) != len(identities) or set(identities) != set(expected):
        raise ValueError("ICA monitoring audit sample identities do not match")
    for field in ("analysis_index", "batch_index"):
        if set(audit[field]) != {str(data.provenance[field])}:
            raise ValueError(f"ICA monitoring audit {field} does not match")
    if set(audit.metric) != {"ica_monitoring"} or set(audit.method) != {"ica"}:
        raise ValueError("Unexpected ICA monitoring metric")
    if set(audit.behavior) != {behavior} or not set(audit.eligible).issubset({"true", "false"}):
        raise ValueError("ICA monitoring audit behavior or eligibility is invalid")
    columns = ["sample", "group", "required_components", "passing_components", "blocking_components",
               "consecutive_stable_batches", "eligible", "reason", "action_result"]
    return audit[columns].rename(columns={
        "sample": "Sample", "group": "Group", "required_components": "Required components",
        "passing_components": "Passing components", "blocking_components": "Blocking components",
        "consecutive_stable_batches": "Informative stable updates", "eligible": "Eligible",
        "reason": "Reason", "action_result": "Action",
    })

