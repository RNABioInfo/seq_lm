"""Validate explicit gene contrasts and independent ICA group selections."""

import csv
import re
from pathlib import Path


def read_gene_contrasts(path, groups):
    """Read headerless, ordered target-minus-reference CSV pairs."""
    pairs = []
    seen = set()
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        for row_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            field = r'(?:"(?:[^"]|"")*"|[^",\r\n]*)'
            if not re.fullmatch(rf"\s*{field}\s*,\s*{field}\s*", line):
                raise ValueError(f"Contrast row {row_number} requires exactly two CSV fields")
            # Permit surrounding spaces even after a quoted field closes.
            line = re.sub(r' *("(?:[^"]|"")*") *', r"\1", line)
            row = next(csv.reader([line], strict=True, skipinitialspace=True))
            pair = tuple(value.strip() for value in row)
            if len(pair) != 2 or not all(pair):
                raise ValueError(f"Contrast row {row_number} requires exactly two nonempty groups")
            if any(any(char in group for char in "\t\r\n") for group in pair):
                raise ValueError(f"Contrast row {row_number} contains tabs or newlines")
            if not set(pair).issubset(groups):
                raise ValueError(f"Contrast row {row_number} references unknown groups: {pair}")
            if pair[0] == pair[1] or pair in seen:
                raise ValueError(f"Contrast row {row_number} is a self-comparison or duplicate: {pair}")
            seen.add(pair)
            pairs.append(pair)
    if not pairs:
        raise ValueError("At least one gene contrast is required")
    return pairs


def resolve_ica_groups(groups, control_group, condition_groups):
    """Include the configured reference automatically in an explicit condition list."""
    control_group = control_group.strip()
    if not control_group or control_group not in groups:
        raise ValueError(f"Unknown ICA control group: {control_group!r}")
    if condition_groups is None or not condition_groups.strip():
        return set(groups)
    conditions = [group.strip() for group in condition_groups.split(",")]
    if (not all(conditions) or len(set(conditions)) != len(conditions)
            or not set(conditions).issubset(groups)):
        raise ValueError("ICA condition groups must be distinct, nonempty sample-sheet group names")
    return set(conditions) | {control_group}
