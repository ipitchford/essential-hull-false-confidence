"""Standalone calibration figure from normal CDFs and deterministic roots.

The plotted hypothesis is H=(-2,-1) union (1,2), theta=0, unit Gaussian
sampling/posterior variance. q is even and, for x>0, its derivative has the
sign of exp(3/2)-2cosh(x), giving exactly two monotone branches. This permits
CDF evaluation without simulation or discretisation of the sampling law.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/false_confidence_matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_examples import log_q


PEAK_X = math.acosh(math.exp(1.5) / 2)
PEAK_Q = math.exp(log_q(PEAK_X))
Q_ZERO = math.exp(log_q(0))


def cdf_q(t: float) -> float:
    if t <= 0:
        return 0.0
    if t >= PEAK_Q:
        return 1.0
    log_t = math.log(t)
    outer = brentq(lambda x: log_q(x) - log_t, PEAK_X,
                   5 + math.sqrt(-2 * log_t), xtol=1e-13)
    value = float(2 * ndtr(-outer))
    if t > Q_ZERO:
        inner = brentq(lambda x: log_q(x) - log_t, 0, PEAK_X, xtol=1e-13)
        value += float(2 * ndtr(inner) - 1)
    return value


def main() -> None:
    aubergine, ochre, teal = "#65334D", "#B07D26", "#247C7B"
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.labelcolor": "#252525", "text.color": "#252525",
        "axes.titleweight": "semibold", "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6), gridspec_kw={"width_ratios": [1.18, 1]})
    t_full = np.unique(np.r_[np.linspace(0, 1, 600),
                            PEAK_Q - np.geomspace(1e-9, PEAK_Q - Q_ZERO, 150),
                            Q_ZERO, PEAK_Q])
    f_full = np.array([cdf_q(float(t)) for t in t_full])
    ax = axes[0]
    ax.plot(t_full, f_full, color=aubergine, lw=2.5, label=r"$F(t)=P_0\{q_H(X)\leq t\}$")
    ax.plot([0, 1], [0, 1], color=ochre, lw=1.7, ls="--", label="Uniform reference")
    ax.axvline(PEAK_Q, color=teal, lw=1, ls=":", alpha=.8)
    ax.scatter([PEAK_Q], [1], color=teal, s=27, zorder=4)
    ax.annotate(f"Maximum posterior mass\n{PEAK_Q:.5f}", xy=(PEAK_Q, 1),
                xytext=(.53, .76), color=teal, fontsize=10,
                arrowprops={"arrowstyle": "-", "color": teal, "lw": 1},
                ha="left", va="center")
    ax.set(xlim=(0, 1), ylim=(0, 1.035), xlabel="Posterior threshold  t",
           ylabel="Sampling cumulative probability  F(t)", title="Full calibration curve")
    ax.legend(loc="lower right", frameon=False, fontsize=9.5)
    ax.grid(alpha=.15)

    ax = axes[1]
    t_tail = np.geomspace(1e-6, .05, 220)
    f_tail = np.array([cdf_q(float(t)) for t in t_tail])
    ax.loglog(t_tail, f_tail, color=aubergine, lw=2.5)
    ax.loglog(t_tail, t_tail, color=ochre, lw=1.7, ls="--")
    value = cdf_q(.01)
    ax.scatter([.01], [value], color=teal, s=30, zorder=4)
    ax.annotate(r"$F(0.01)\approx1.64\times10^{-5}$" + "\n" + r"$F(0.01)<0.01$",
                xy=(.01, value), xytext=(8e-6, 5e-4), color=teal,
                arrowprops={"arrowstyle": "-", "color": teal, "lw": 1},
                ha="left", va="center", fontsize=10.5)
    ax.set(xlim=(1e-6, .05), xlabel="Posterior threshold  t  (log scale)",
           ylabel="F(t)  (log scale)", title="Lower-tail violation")
    ax.grid(alpha=.15, which="major")
    fig.suptitle("One fixed false hypothesis can cross the uniform calibration line", fontsize=13, y=.985)
    fig.text(.5, .012, r"$H=(-2,-1)\cup(1,2)$;  $X\sim N(0,1)$;  $Q_x=N(x,1)$.  Normal-CDF roots; no simulation.",
             ha="center", fontsize=9.5, color="#555555")
    fig.tight_layout(rect=(0, .055, 1, .94), w_pad=2.4)
    out = Path(__file__).parent / "figures"
    out.mkdir(exist_ok=True)
    fig.savefig(out / "calibration_cdf.png", dpi=220)
    fig.savefig(out / "calibration_cdf.pdf")
    plt.close(fig)
    print(f"Figure written; exact maximiser x={PEAK_X:.12g}, numerical qmax={PEAK_Q:.12g}, F(.01)={value:.12g}")


if __name__ == "__main__":
    main()
