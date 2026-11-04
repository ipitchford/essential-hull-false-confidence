# SPDX-License-Identifier: MIT
"""Compare regenerated check JSON with the shipped records, allowing last-digit float differences."""
import json, math, sys
def cmp(x, y, path):
    if isinstance(x, dict):
        if set(x) != set(y): sys.exit(f"keys differ at {path}")
        for k in x:
            if k == "environment": continue
            cmp(x[k], y[k], path + "/" + k)
    elif isinstance(x, list):
        if len(x) != len(y): sys.exit(f"length differs at {path}")
        for i, (u, v) in enumerate(zip(x, y)): cmp(u, v, f"{path}[{i}]")
    elif isinstance(x, float) and isinstance(y, (int, float)):
        if not math.isclose(x, y, rel_tol=1e-7, abs_tol=1e-9): sys.exit(f"value differs at {path}: {x} vs {y}")
    elif x != y:
        sys.exit(f"value differs at {path}: {x!r} vs {y!r}")
for n in ["gaussian_false_confidence_checks", "exponential_power_checks", "extension_checks"]:
    a, b = (json.load(open(f"{d}/{n}.json")) for d in sys.argv[1:3])
    cmp(a, b, n)
print("JSON agrees with the shipped records")
