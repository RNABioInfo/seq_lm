"""Adapt ICA contrast intervals to the shared precision forecast."""
import json
import math
from collections import defaultdict

from workflow_glue.readiness import read_depth_readiness

from .engine import read_tsv


def component_radius(row, weight, log_base):
    if row["status"] != "tested":
        return None

    effect = float(row["activity_difference"])
    lower = float(row["ci_lower"])
    upper = float(row["ci_upper"])
    if not all(map(math.isfinite, (effect, lower, upper))) or not lower <= effect <= upper:
        raise ValueError("Invalid tested ICA activity interval")

    radius = max(effect - lower, upper - effect) * weight * math.log2(log_base)
    return radius if math.isfinite(radius) else None


def load_ica_precision(root, definitions, samples):
    provenance = json.loads((root / "provenance.json").read_text())
    if provenance.get("schema_version") != 2:
        raise ValueError("Unsupported ICA provenance schema version")
    model = provenance["model"]
    settings = provenance["settings"]
    control = settings["control_group"]
    components = model["components"]
    selected = model.get("stopping_components")
    universe = selected if selected is not None else components
    weights = model["component_max_abs_weights"]
    log_base = float(settings["log_base"])

    if not math.isfinite(log_base) or log_base <= 1:
        raise ValueError("Invalid ICA log base")
    if not universe or len(set(universe)) != len(universe) or not set(universe) <= set(components):
        raise ValueError("Invalid ICA stopping component universe")
    if any(not math.isfinite(weights[k]) or weights[k] <= 0 for k in universe):
        raise ValueError("Invalid ICA stopping component weights")

    expected = {s["key"] for s in samples if s["group"] in settings["condition_groups"]}
    members = {}
    for row in read_tsv(root / "sample_metadata.tsv"):
        key = f'{row["group"]}\t{row["alias"]}'
        if key in members:
            raise ValueError("Duplicate ICA monitoring sample")
        depth = float(row["assigned_abundance"])
        if not math.isfinite(depth) or depth < 0:
            raise ValueError("Invalid ICA assigned abundance")
        members[key] = {
            "key": key,
            "name": row["alias"],
            "group": row["group"],
            "depth": depth,
        }
    if members.keys() != expected:
        raise ValueError("ICA monitoring samples do not match the selected cohort")

    groups = {sample["group"] for sample in members.values()}
    _, ready_groups = read_depth_readiness(
        ((sample["group"], sample["depth"]) for sample in members.values()),
        settings["min_read_count"],
        settings["min_replicate_sample_count"],
    )
    status = json.loads((root / "status.json").read_text())["status"]
    if status not in ("ready", "deferred") or (status == "ready") != (ready_groups == groups):
        raise ValueError("Invalid ICA snapshot status")
    table = root / ("differential_activity.tsv" if status == "ready" else "monitoring_differential_activity.tsv")
    rows_by_group = defaultdict(dict)
    if table.exists():
        for row in read_tsv(table):
            group, component = row["target_group"], row["component_id"]
            if row["control_group"] != control or component in rows_by_group[group]:
                raise ValueError("Unexpected or duplicate ICA monitoring interval")
            rows_by_group[group][component] = row

    precision = {}
    for definition in definitions:
        target = definition["target_group"]
        if definition["reference_group"] != control:
            raise ValueError("ICA monitoring control does not match the snapshot")

        pair_groups = {target, control}
        if not pair_groups <= groups:
            raise ValueError("ICA monitoring groups do not match the snapshot")
        if not pair_groups <= ready_groups:
            continue

        pair = [sample for sample in members.values() if sample["group"] in pair_groups]
        intervals = rows_by_group[target]
        if intervals.keys() != set(components):
            raise ValueError(f"Missing or unexpected ICA intervals for {target}")

        precision[definition["contrast_id"]] = {
            "rows": [
                {"feature_id": key, "R": component_radius(intervals[key], weights[key], log_base)}
                for key in universe
            ],
            "samples": pair,
            "depths": {sample["key"]: sample["depth"] for sample in pair},
            "universe": universe,
        }

    identities = {
        "model": model["hashes"],
        "settings": settings,
        "software": provenance["software"],
        "controls": sorted(provenance["control_sample_ids"]),
        "ica_stopping_components": selected,
    }
    return precision, identities
