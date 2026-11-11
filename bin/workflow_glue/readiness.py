"""Shared library-depth precondition for DEA and ICA."""

from collections import Counter
from collections.abc import Iterable


def read_depth_readiness(
    libraries: Iterable[tuple[str, float]],
    min_read_count: int,
    min_replicate_sample_count: int,
) -> tuple[list[bool], set[str]]:
    """Return per-library threshold flags and the groups eligible for analysis."""
    if (
        type(min_read_count) is not int
        or min_read_count < 0
        or type(min_replicate_sample_count) is not int
        or min_replicate_sample_count < 1
    ):
        raise ValueError("Readiness requires a nonnegative integer read count and a positive integer replicate count.")

    passing_replicates: Counter[str] = Counter()
    groups_with_empty_libraries = set()
    passing_libraries = []

    for group, depth in libraries:
        positive = depth > 0
        passes = positive and depth >= min_read_count
        passing_libraries.append(passes)
        passing_replicates[group] += passes
        if not positive:
            groups_with_empty_libraries.add(group)

    ready_groups = {
        group
        for group, count in passing_replicates.items()
        if count >= min_replicate_sample_count and group not in groups_with_empty_libraries
    }
    return passing_libraries, ready_groups
