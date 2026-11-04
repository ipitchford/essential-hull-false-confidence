"""Deterministic checks of the extension results. Floating-point sanity checks, not proofs.
Every check raises SystemExit on failure, so results do not depend on python -O."""
import json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gammaln, ndtr
from scipy.stats import chi2, ncx2

def need(c, m):
    if not c: raise SystemExit("check failed: " + m)
out = {}

# 1. Two-interval Gaussian example: last threshold violating F(t)>=t at theta=0.
Phi = ndtr
q = lambda x: Phi(2-x)-Phi(1-x)+Phi(-1-x)-Phi(-2-x)
xs = np.linspace(0, 12, 240001); qs = q(xs); i = int(np.argmax(qs)); xm = xs[i]
def F(t):
    # q increases on [0,xm], decreases after; |X| folded normal
    p = 0.0
    if q(0.0) <= t:
        p += 2*(Phi(brentq(lambda x: q(x)-t, 0, xm))-0.5)
    hi = brentq(lambda x: q(x)-t, xm, 40)
    return p + 2*(1-Phi(hi))
need(F(1e-3) < 1e-3, "low-threshold violation")
need(F(0.35) > 0.35, "high end conservative")
tc = brentq(lambda t: F(t)-t, 1e-3, 0.35)
need(abs(q(xm)-0.3893500291) < 1e-8, "max mass")
out["two_interval_crossover_t"] = tc

# 2. Laplace bounded failure: F(t)=t/K when K=sinh B - sinh d.
B, d = 1.5, 0.1
K = math.sinh(B)-math.sinh(d); need(K > 1, "K>1")
qL = lambda x: 0.5*quad(lambda u: math.exp(-abs(x-u)), d, B)[0]+0.5*quad(lambda u: math.exp(-abs(x-u)), -B, -d)[0]
x = 30.0; need(abs(qL(x)/(K*math.exp(-x))-1) < 1e-8, "tail constant")
t = K*math.exp(-x)
need(abs(math.exp(-math.log(K/t))-t/K) < 1e-20, "F formula")
out["laplace_failure"] = {"B": B, "delta": d, "K": K, "F_over_t": 1/K}

# 3. Tail criterion: log ratio log S(R)-log g(R-b).
def logS_ratio(psi, R, b, dim):
    # S(R)/g(R-b) with common normaliser and sphere area removed from S to compare growth
    I = quad(lambda s: (R+s)**(dim-1)*math.exp(psi(R)-psi(R+s)), 0, 200, limit=400)[0]
    return math.log(I)-psi(R)+psi(R-b)
pow_ = lambda p: (lambda r: r**p/p)
xl = lambda r: (1+r)*math.log(1+r)
tab = {}
for name, psi, dim, b, trend in [("p=1.5,d=2", pow_(1.5), 2, .5, -1), ("p=1,d=3", pow_(1.0), 3, .5, 1),
                                  ("xlogx,d=1", xl, 1, .5, -1), ("xlogx,d=3", xl, 3, .5, 1)]:
    v = [logS_ratio(psi, R, b, dim) for R in (50., 500., 5000.)]
    need((v[2]-v[0])*trend > 0, name)
    tab[name] = v
out["tail_log_ratio_R50_500_5000"] = tab

# 4. Conjunction disc. C = disc radius b in sigma units, true miss m. Theorem bound and exact values.
rows = []
for sr, mr in [(1,0),(2,0),(2,.5),(10,0),(10,.5),(10,.9),(20,.5)]:
    b = 1/sr; m = mr*b
    Pc = lambda r: ncx2.cdf(b*b, 2, r*r) if r > 0 else chi2.cdf(b*b, 2)
    atruth = float(Pc(m)); row = dict(sigma_over_R=sr, miss_over_R=mr, Pc_at_truth=atruth, P_report_at_most={})
    for s in (.1, .01, .001):
        if s >= Pc(0): pr = 1.0
        else:
            r = brentq(lambda r: Pc(r)-s, 1e-9, 60)
            pr = float(ncx2.sf(r*r, 2, m*m) if m > 0 else chi2.sf(r*r, 2))
        if atruth <= s: need(pr >= .5, "finite-level theorem")
        row["P_report_at_most"][str(s)] = pr
    rows.append(row)
out["conjunction"] = rows

Path(__file__).with_name("extension_checks.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=0)[:1800]); print("ALL PASSED")
