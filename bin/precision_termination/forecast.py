"""Dependency-free implementation of the evaluated aggregate and depth curve."""

import math


FIT_HALF_LIFE_SECONDS = 2 * 60 * 60
HUBER_RELATIVE_LIMIT = 0.04
ROBUST_FIT_ITERATIONS = 8


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


def nonnegative_fit(rows, weights):
    """Fit two nonnegative coefficients to a response of one."""
    weighted = [(x, y, w * x, w * y) for (x, y), w in zip(rows, weights, strict=True)]
    s00 = math.fsum(x * wx for x, _, wx, _ in weighted)
    s01 = math.fsum(x * wy for x, _, _, wy in weighted)
    s11 = math.fsum(y * wy for _, y, _, wy in weighted)
    s0 = math.fsum(wx for _, _, wx, _ in weighted)
    s1 = math.fsum(wy for _, _, _, wy in weighted)

    if s00 == 0:
        return math.nan, math.nan

    candidates = [(s0 / s00, 0.0)]
    if s11 > 0:
        candidates.append((0.0, s1 / s11))
    det = s00 * s11 - s01 * s01
    if det > s00 * s11 * 1e-12:
        a = (s0 * s11 - s1 * s01) / det
        b = (s1 * s00 - s0 * s01) / det
        if a >= 0 and b >= 0:
            candidates.append((a, b))

    return min(
        candidates,
        key=lambda ab: math.fsum(w * (1 - ab[0] * x - ab[1] * y) ** 2
                                 for (x, y), w in zip(rows, weights, strict=True)),
    )


def fit_curve(depths, radii, times):
    """Fit the raw scores with recency weights and a robust relative error loss."""
    points = [(d, r, t) for d, r, t in zip(depths, radii, times, strict=True)
              if math.isfinite(d + r + t) and d > 0 and r > 0]
    if not points:
        return math.nan, math.nan

    min_depth = min(d for d, _, _ in points)
    min_radius = min(r for _, r, _ in points)
    radius_squared = min_radius * min_radius
    if not 0 < radius_squared < math.inf:
        return math.nan, math.nan

    # Divide each residual by observed R². Scale the coefficients so the
    # two design columns are bounded by one and the response is one.
    inverse_radii = [(min_radius / r) ** 2 for _, r, _ in points]
    rows = [(x, x * min_depth / d) for x, (d, _, _) in zip(inverse_radii, points, strict=True)]
    latest = max(t for _, _, t in points)
    recency = [2 ** (-(latest - t) / FIT_HALF_LIFE_SECONDS) for _, _, t in points]

    a, b = nonnegative_fit(rows, recency)
    for _ in range(ROBUST_FIT_ITERATIONS):
        # Huber loss downweights large relative misses without removing points.
        weights = [w / max(1.0, abs(1 - a * x - b * y) / HUBER_RELATIVE_LIMIT)
                   for (x, y), w in zip(rows, recency, strict=True)]
        a, b = nonnegative_fit(rows, weights)

    return a * radius_squared, b * radius_squared * (min_depth / 1e6)


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
    a, b = fit_curve(d, r, ts)
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
