"""Finite deterministic checks for the p>1 extension, not a proof.

For H=(-2,-1) union (1,2), posterior masses and sampling tails use the
incomplete-gamma formula for the exactly specified exponential-power density.
No data download or Monte Carlo is used.
"""
from __future__ import annotations

import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gammaincc, gammaln, logsumexp


def log_sf_positive(x: float, p: float) -> float:
    if not (x >= 0):
        raise SystemExit("check failed: x >= 0")
    value = float(gammaincc(1 / p, x**p / p))
    if value > 1e-250:
        return math.log(0.5 * value)
    # Exact rescaling of the upper incomplete-gamma integral prevents
    # underflow. The remaining integral is of order one and is evaluated
    # by deterministic quadrature; this remains a floating-point check.
    shape, lower = 1 / p, x**p / p
    scaled, _ = quad(lambda s: math.exp(-s) * (1 + s / lower)**(shape - 1),
                     0, math.inf, epsabs=1e-13, epsrel=1e-13)
    return -lower + (shape - 1) * math.log(lower) + math.log(scaled) - float(gammaln(shape)) - math.log(2)


def logdiff(big: float, small: float) -> float:
    return big + math.log(-math.expm1(small - big))


def log_q(x: float, p: float) -> float:
    # Only x >= 2 is used, so both interval masses are differences of
    # small lower-tail probabilities, evaluated by reflected upper tails.
    if not (x >= 2):
        raise SystemExit("check failed: x >= 2")
    return float(logsumexp([
        logdiff(log_sf_positive(x + 1, p), log_sf_positive(x + 2, p)),
        logdiff(log_sf_positive(x - 2, p), log_sf_positive(x - 1, p)),
    ]))


def main() -> None:
    rows = []
    for p in [1.5, 3.0]:
        compact_lower_bound = math.exp(logdiff(log_sf_positive(1, p), log_sf_positive(2, p)))
        for log_t in [-10.0, -40.0, -100.0, -200.0]:
            if not (math.exp(log_t) < compact_lower_bound):
                raise SystemExit("check failed: math.exp(log_t) < compact_lower_bound")
            # q(x) >= compact_lower_bound on [-2,2] by unimodality of
            # each unit interval mass, and decreases strictly on x>=2.
            cutoff = float(brentq(lambda x: log_q(x, p) - log_t,
                                 2, (-p * log_t)**(1 / p) + 3,
                                 xtol=1e-13))
            for theta in [0.0, 2.0, 3.0]:
                if not (cutoff > abs(theta)):
                    raise SystemExit("check failed: cutoff > abs(theta)")
                log_f = float(logsumexp([
                    log_sf_positive(cutoff + theta, p),
                    log_sf_positive(cutoff - theta, p),
                ]))
                rows.append({
                    "p": p,
                    "theta": theta,
                    "log_t": log_t,
                    "cutoff": cutoff,
                    "log_F": log_f,
                    "scaled_log_ratio": (log_f - log_t) / (-p * log_t)**((p - 1) / p),
                    "predicted_limit": abs(theta) - 2,
                    "root_log_residual": abs(log_q(cutoff, p) - log_t),
                })
    result = {
        "scope": "24 finite deterministic sanity checks, not formal or interval proof verification",
        "density": "c_p exp(-abs(z)^p/p)",
        "hypothesis": "(-2,-1) union (1,2)",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "tail_checks": rows,
        "max_root_log_residual": max(row["root_log_residual"] for row in rows),
        "laplace_coefficient": math.sinh(math.log(2.5)) - math.sinh(math.log(2)),
        "laplace_coefficient_exact": "3/10 (proved analytically in manuscript)",
    }
    out = Path(__file__).with_name("exponential_power_checks.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(out), "rows": len(rows),
                      "max_root_log_residual": result["max_root_log_residual"],
                      "at_log_t_minus_200": [r for r in rows if r["log_t"] == -200]}, indent=2))


if __name__ == "__main__":
    main()
