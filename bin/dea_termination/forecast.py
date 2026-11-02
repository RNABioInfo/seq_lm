"""Dependency-free implementation of the evaluated aggregate and depth curve."""

import math


def quantile(values, q):
    values = sorted(values)
    if not values:
        return math.nan
    position = (len(values) - 1) * q
    lo = int(position)
    hi = min(lo + 1, len(values) - 1)
    return values[lo] + (position - lo) * (values[hi] - values[lo])


def information_radius(r, scale=1.0):
    r = list(r)
    if any(math.isfinite(x) and x < 0 for x in r):
        raise ValueError("negative R")
    credits = [
        1 / (1 + (x / scale) ** 2)
        if math.isfinite(x) and abs(x / scale) < 1e154
        else 0.0
        for x in r
    ]
    info = math.fsum(credits) / len(r) if r else 0.0
    return (scale * math.sqrt(1 / info - 1) if info > 0 else math.inf), info


def fit_curve(depths, radii):
    pairs = [
        (1e6 / d, r * r)
        for d, r in zip(depths, radii)
        if math.isfinite(d + r) and d > 0 and r > 0
    ]
    if not pairs:
        return math.nan, math.nan
    n = len(pairs)
    sx = math.fsum(x for x, y in pairs)
    sxx = math.fsum(x * x for x, y in pairs)
    sy = math.fsum(y for x, y in pairs)
    sxy = math.fsum(x * y for x, y in pairs)
    det = n * sxx - sx * sx
    candidates = [(sy / n, 0.0), (0.0, sxy / sxx)]
    if det > max(n * sxx, 1) * 1e-12:
        b = (n * sxy - sx * sy) / det
        a = (sy - b * sx) / n
        if a >= 0 and b >= 0:
            candidates.append((a, b))
    return min(
        candidates,
        key=lambda ab: math.fsum((y - ab[0] - ab[1] * x) ** 2 for x, y in pairs),
    )


def curve_summary(depths, radii, times, extra, config):
    valid = [
        (d, r, t)
        for d, r, t in zip(depths, radii, times)
        if math.isfinite(d + r) and d > 0 and r > 0
    ]
    result = {
        "finite_points": len(valid),
        "curve_a": math.nan,
        "curve_b": math.nan,
        "fit_R": math.nan,
        "fit_discrepancy": math.nan,
        "curve_gain": math.nan,
        "elbow_depth": math.nan,
        "elbow_seconds": math.nan,
        "elbow_strength": 0.0,
        "post_elbow_points": 0,
        "observed_relative_drop": math.nan,
        "tail_slope_per_million_reads": math.nan,
        "tail_to_initial_slope_ratio": math.nan,
        "curve_usable": False,
        "elbow_pass": False,
    }

    if not valid:
        return result

    d, r, ts = map(list, zip(*valid))
    a, b = fit_curve(d, r)
    current = depths[-1]
    observed = radii[-1]

    if not math.isfinite(a + b + current + observed) or current <= 0 or observed <= 0:
        return result

    fitted = math.sqrt(a + b * 1e6 / current)
    disc = abs(fitted - observed) / observed
    gain = (
        1 - math.sqrt(a + b * 1e6 / (current + extra)) / fitted
        if math.isfinite(extra) and extra >= 0
        else math.nan
    )
    result.update(
        curve_a=a,
        curve_b=b,
        fit_R=fitted,
        fit_discrepancy=disc,
        curve_gain=gain,
        observed_relative_drop=1 - observed / r[0],
        curve_usable=len(r) >= config["min_points"]
        and max(d) - min(d) > 1e-6
        and disc <= config["max_fit_discrepancy"]
        and math.isfinite(gain),
    )

    if len(r) < 2 or b <= 0 or max(d) - min(d) <= 1e-6:
        return result
    first = d[0]
    end = d[-1]
    y0 = math.sqrt(a + b * 1e6 / first)
    yend = math.sqrt(a + b * 1e6 / end)

    if y0 <= yend:
        return result

    grid = [first + (end - first) * k / 1000 for k in range(1001)]
    chord = [
        1
        - (x - first) / (end - first)
        - (math.sqrt(a + b * 1e6 / x) - yend) / (y0 - yend)
        for x in grid
    ]

    k = max(range(1001), key=chord.__getitem__)
    kd = grid[k]
    kt = ts[-1]

    for j in range(1, len(d)):
        if d[j] > d[j - 1] and d[j - 1] <= kd <= d[j]:
            kt = ts[j - 1] + (kd - d[j - 1]) / (d[j] - d[j - 1]) * (ts[j] - ts[j - 1])
            break

    slope = lambda dd: -(b * 1e6) / (2 * dd * dd * math.sqrt(a + b * 1e6 / dd))
    posts = sum(x > kd for x in d)
    
    result.update(
        elbow_depth=kd,
        elbow_seconds=kt,
        elbow_strength=chord[k],
        post_elbow_points=posts,
        tail_slope_per_million_reads=slope(current) * 1e6,
        tail_to_initial_slope_ratio=slope(current) / slope(first),
        elbow_pass=chord[k] >= config["elbow_min_strength"]
        and posts >= config["post_elbow_points"]
        and result["observed_relative_drop"] >= config["elbow_min_observed_drop"],
    )
    return result
