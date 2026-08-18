"""Session brief N, part 4 -- the intra-fibre curl constant.

n3_curl.py Sec. 5 finds that the supremum of |curl A| over triples with a COMMON
sigma is about 0.1438, far below the global (log 2)/2 = 0.34657.  The candidate
closed form is

    (1/2) log(4/3) = 0.1438410362258904 ... ,

i.e. an intra-fibre cycle asymmetry ratio Omega = 4/3 against the global 2.

This script refines the search: a long multi-restart differential evolution
followed by a deep pattern search, at k = 2..5 lines, plus a targeted local
polish of the best point found, and reports the agreement with (1/2)log(4/3).
It also reports the extremal configuration, so the family can be read off.

    python research/birkhoff/n4_fibre_curl.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "realizability"))

import i_cone as T                                        # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

LOG2 = math.log(2.0)
TARGET = 0.5 * math.log(4.0 / 3.0)


def parts(a, b):
    bs = T.candidates(a, b)
    v = T.phi(a, b, bs)
    e0 = math.log(b.R) - math.log(a.R)
    e1 = math.log(b.Lam) - math.log(a.Lam)
    hi = max(float(v.max()), e0, e1)
    lo = min(float(v.min()), e0, e1)
    return 0.5 * ((hi - max(e0, e1)) - (min(e0, e1) - lo))


def curl_of(fs):
    return sum(parts(fs[i], fs[j]) for i, j in ((0, 1), (1, 2), (2, 0)))


def mk(z, k):
    z = np.asarray(z, float).reshape(3, k, 2)
    out = []
    for row in z:
        c = np.maximum(row[:, 0], 0.0)
        x = np.maximum(row[:, 1], 0.0)
        if c.max() <= 1e-9:
            c = c + 1e-3
        if x.max() <= 1e-9:
            x = x + 1e-3
        t = T.Trop(c, x)
        out.append(T.Trop(t.c, t.x * (t.R / t.Lam)))       # push sigma to 0
    return out


def neg(z, k):
    try:
        return -abs(curl_of(mk(z, k)))
    except ValueError:
        return 1e3


def main():
    print("=== the supremum of |curl A| on a single sigma-fibre ===")
    print(f"  candidate closed form (1/2) log(4/3) = {TARGET:.15f}")
    print(f"  {'k':>3} {'best |curl| found':>22} {'- (1/2)log(4/3)':>18} "
          f"{'/(log2/2)':>11}")
    best_all, best_z, best_k = -1.0, None, None
    for k in (2, 3, 4, 5):
        bounds = [(0.0, 12.0)] * (3 * k * 2)
        best = math.inf
        bz = None
        for t in range(16):
            z, f = differential_evolution(neg, bounds, args=(k,), seed=4021 + 733 * t + k,
                                          maxiter=700, popsize=18, F=(0.3, 1.2), CR=0.9)
            for step in (1.0, 0.3, 0.08, 0.02, 5e-3, 1e-3, 2e-4, 4e-5, 8e-6,
                         1.6e-6, 3e-7, 6e-8, 1e-8):
                z, f = pattern_search(neg, z, args=(k,), step=step, min_step=1e-14,
                                      maxiter=120000, bounds=bounds)
            if f < best:
                best, bz = f, z
        print(f"  {k:>3} {-best:22.12f} {-best - TARGET:18.3e} "
              f"{-best/(LOG2/2):11.7f}", flush=True)
        if -best > best_all:
            best_all, best_z, best_k = -best, bz, k
    print(f"\n  best overall: {best_all:.12f} at k = {best_k}; "
          f"(1/2)log(4/3) - best = {TARGET - best_all:.3e}")
    fs = mk(best_z, best_k)
    print("  extremal configuration (lines (intercept, slope), sigma = 0 each):")
    for i, t in enumerate(fs):
        print(f"    Phi_{i+1}: c = {np.array2string(t.c, precision=6)}  "
              f"x = {np.array2string(t.x, precision=6)}")
    print("  pairwise D:", [f"{parts(fs[i], fs[j]):.9f}"
                            for i, j in ((0, 1), (1, 2), (2, 0))])


if __name__ == "__main__":
    main()
