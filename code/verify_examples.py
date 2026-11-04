"""Deterministic numerical sanity checks for the Gaussian hull candidate.

These finite evaluations do not prove the general theorems.  They use exact
one-dimensional normal probabilities and deterministic root finding, not Monte
Carlo.  Run with Python 3, numpy and scipy.
"""

from __future__ import annotations

import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq, minimize_scalar
from scipy.special import log_ndtr, logsumexp, ndtr


def logdiffexp(big: float, small: float) -> float:
    return big + math.log(-math.expm1(small - big))


def log_interval_mass(x: float, left: float, right: float) -> float:
    """log P{N(x,1) in (left,right)}, stable in the tail."""
    lo, hi = left - x, right - x
    if lo > 0:
        return logdiffexp(float(log_ndtr(-lo)), float(log_ndtr(-hi)))
    return logdiffexp(float(log_ndtr(hi)), float(log_ndtr(lo)))


def log_q(x: float) -> float:
    return float(logsumexp([
        log_interval_mass(x, -2.0, -1.0),
        log_interval_mass(x, 1.0, 2.0),
    ]))


def tail_cutoff(log_t: float) -> float:
    # All chosen t <= .01. On [-2,2], q >= Phi(2)-Phi(1) > .13;
    # on [2,infinity), q decreases strictly. Thus {q(X)<=t}={|X|>=r}.
    return float(brentq(
        lambda x: log_q(x) - log_t,
        2.0,
        5.0 + math.sqrt(-2.0 * log_t),
        xtol=1e-13,
    ))


def main() -> None:
    max_result = minimize_scalar(lambda x: -math.exp(log_q(x)), bounds=(0, 3), method="bounded")
    rows = []
    for theta in [0.0, 0.5, 0.99, 2.0, 3.0]:
        signed_hull_distance = 2.0 - abs(theta)
        for log_t in [math.log(1e-2), math.log(1e-5), math.log(1e-20), math.log(1e-100), math.log(1e-300)]:
            cutoff = tail_cutoff(log_t)
            log_f = float(logsumexp([log_ndtr(-cutoff - theta), log_ndtr(theta - cutoff)]))
            rows.append({
                "theta": theta,
                "log_t": log_t,
                "t": math.exp(log_t),
                "cutoff": cutoff,
                "log_F": log_f,
                "F_over_t": math.exp(log_f - log_t),
                "scaled_log_ratio": (log_f - log_t) / math.sqrt(-2.0 * log_t),
                "predicted_limit": -signed_hull_distance,
                "strict_false_confidence": log_f < log_t,
                "root_log_residual": abs(log_q(cutoff) - log_t),
            })
    # Halfspace at its boundary must be exactly uniform: posterior Q_X(H)=Phi(-Z).
    threshold_values = [0.01, 0.05, 0.5, 0.95, 0.99]
    from scipy.special import ndtri
    halfspace_residuals = [abs(float(ndtr(ndtri(t))) - t) for t in threshold_values]
    # Severe exterior-interval example has a deterministic high-confidence floor.
    exterior_floors = [{"a": a, "posterior_infimum": float(2 * ndtr(-a))} for a in [.01, .1, 1.0]]
    result = {
        "scope": "finite deterministic sanity checks only; no formal proof verification",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "hypothesis": "(-2,-1) union (1,2)",
        "sampling": "N(theta,1)",
        "posterior": "N(x,1)",
        "essential_convex_hull": "[-2,2]",
        "numerical_maximum_posterior_mass": -float(max_result.fun),
        "maximiser_positive": float(max_result.x),
        "tail_checks": rows,
        "halfspace_max_uniform_residual": max(halfspace_residuals),
        "exterior_interval_checks": exterior_floors,
        "max_root_log_residual": max(r["root_log_residual"] for r in rows),
    }
    out = Path(__file__).with_name("gaussian_false_confidence_checks.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(out),
        "rows": len(rows),
        "max_root_log_residual": result["max_root_log_residual"],
        "posterior_maximum": result["numerical_maximum_posterior_mass"],
        "at_t_1e_minus_300": [r for r in rows if r["log_t"] < -600],
        "halfspace_max_uniform_residual": result["halfspace_max_uniform_residual"],
    }, indent=2))


if __name__ == "__main__":
    main()
