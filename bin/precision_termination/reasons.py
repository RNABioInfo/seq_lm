"""Public stopping reasons shared by audits and reports."""

_DETAILED_REASONS = {
    "precision_met": ("marginal_precision",),
    "improvement_continues": ("continue_observed_gain", "continue_expected_gain"),
    "stability_not_confirmed": ("unusable_curve", "elbow_not_confirmed", "raw_R_unstable"),
    "waiting_for_reads": (
        "warming_up", "insufficient_cadence", "waiting_for_fresh_data", "interrupted",
    ),
    "results_not_ready": (
        "analysis_deferred", "unavailable_precision", "processing_lag", "stale_snapshot",
        "newer_input_or_stale_result",
    ),
    "not_monitored": ("no_live_samples", "terminal_snapshot", "fixed_reference", "fixed", "not_required"),
}
_PUBLIC_REASONS = {
    detail: code for code, details in _DETAILED_REASONS.items() for detail in details
}


def audit_reason(reason: str) -> str:
    """Group detailed reasons and remove duplicates across a sample's contrasts."""
    reasons = (value.strip() for value in reason.split(","))
    codes = {_PUBLIC_REASONS.get(value, value) for value in reasons if value}
    return ",".join(sorted(codes)) or "not_monitored"
