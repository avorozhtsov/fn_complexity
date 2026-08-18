"""Extremal search over signatures with free multiplicities.

The sharp witnesses of the structure theorem need very large multiplicities
(the thin signature (r,1,...,1) with r -> infinity), so a search over a fixed
number of *atoms* cannot find them.  Here a signature is parametrised by k
distinct log-values x_i >= 0 and k log-multiplicities c_i >= 0, which reaches
r = e^{sum c_i} up to astronomically large values at constant cost.

Quantities climbed:
    |D|      the non-exact part of the flow A;   proved  < (log 2)/2
    eps      = P + Q = d - |sigma_b - sigma_a|;  proved  <= 2 log(1+e^{-Delta})
    curl     |A(a,b)+A(b,c)+A(c,a)| over a triangle

    python research/realizability/extremes.py
"""
from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

HERE = Path(__file__).resolve().parent
GRID = C.make_grid(-30.0, 60.0, 0.004)
XMAX, CMAX = 6.0, 90.0        # log-value and log-multiplicity ranges


def sig_from(v):
    """v = [x_1..x_k, c_1..c_k] -> Sig with values e^{x_i}, multiplicities e^{c_i}."""
    k = len(v) // 2
    xs = np.clip(v[:k], 0.0, XMAX)
    cs = np.clip(v[k:], 0.0, CMAX)
    if xs.max() <= 1e-9:
        xs = xs.copy()
        xs[int(np.argmax(v[:k]))] = 1e-3
    return C.Sig.from_logs(xs, np.exp(cs))


def pair_parts(v, k):
    a = sig_from(v[:2 * k])
    b = sig_from(v[2 * k:])
    return C.parts(a, b, GRID), a, b


def neg_absD(v, k):
    try:
        p, _, _ = pair_parts(v, k)
    except ValueError:
        return 1e3
    return -abs(p["D"])


def neg_eps(v, k):
    try:
        p, _, _ = pair_parts(v, k)
    except ValueError:
        return 1e3
    return -p["eps"]


def neg_curl(v, k):
    try:
        ss = [sig_from(v[2 * k * i:2 * k * (i + 1)]) for i in range(3)]
    except ValueError:
        return 1e3
    tot = 0.0
    for a, b in zip(ss, ss[1:] + ss[:1]):
        tot += C.parts(a, b, GRID)["A"]
    return -abs(tot)


def climb(fun, dim, k, seed, maxiter=250, popsize=14):
    bounds = [(0.0, XMAX)] * k + [(0.0, CMAX)] * k
    bounds = bounds * (dim // (2 * k))
    x, f = differential_evolution(fun, bounds, args=(k,), seed=seed,
                                  maxiter=maxiter, popsize=popsize,
                                  F=(0.3, 1.2), CR=0.9)
    for step in (1.0, 0.2, 0.03, 3e-3, 3e-4):
        x, f = pattern_search(fun, x, args=(k,), step=step, min_step=1e-9,
                              maxiter=30000, bounds=bounds)
    return x, -f


def describe(s):
    return "{" + ", ".join(f"e^{x:.4g}^(e^{math.log(m):.4g})"
                           for x, m in zip(s.xs, s.mults)) + "}"


def main():
    rows = []
    LOG2 = C.LOG2

    print("=== sup |D| over pairs, k distinct values with free multiplicity ===")
    print(f"  {'k':>2} {'max |D|':>12} {'/(log2/2)':>11} {'Delta':>9} "
          f"{'bound .5log(1+e^-D)':>20}")
    for k in (1, 2, 3):
        t0 = time.time()
        x, v = climb(neg_absD, 4 * k, k, seed=20 + k)
        p, a, b = pair_parts(x, k)
        D = abs(p["dsigma"])
        bd = 0.5 * math.log1p(math.exp(-D))
        print(f"  {k:>2} {v:12.8f} {v/(LOG2/2):11.4f} {D:9.5f} {bd:20.8f}"
              f"   {time.time()-t0:.0f}s")
        print(f"       a = {describe(a)}   r_a = {a.r:.4g}")
        print(f"       b = {describe(b)}   r_b = {b.r:.4g}")
        rows.append(["max|D|", k, v, LOG2 / 2, D, bd])

    print("\n=== sup eps = P + Q ===")
    print(f"  {'k':>2} {'max eps':>12} {'/log2':>9} {'Delta':>9} "
          f"{'bound 2log(1+e^-D)':>20}")
    for k in (1, 2, 3):
        t0 = time.time()
        x, v = climb(neg_eps, 4 * k, k, seed=40 + k)
        p, a, b = pair_parts(x, k)
        D = abs(p["dsigma"])
        bd = 2 * math.log1p(math.exp(-D))
        print(f"  {k:>2} {v:12.8f} {v/LOG2:9.4f} {D:9.5f} {bd:20.8f}"
              f"   P={p['P']:.6f} Q={p['Q']:.6f}   {time.time()-t0:.0f}s")
        rows.append(["max eps", k, v, LOG2, D, bd])

    print("\n=== sup |curl A| over triangles ===")
    print(f"  {'k':>2} {'max |curl|':>12} {'/(3log2/2)':>12} {'mean|A|':>10}")
    best = None
    for k in (1, 2, 3):
        t0 = time.time()
        x, v = climb(neg_curl, 6 * k, k, seed=60 + k, maxiter=200)
        ss = [sig_from(x[2 * k * i:2 * k * (i + 1)]) for i in range(3)]
        As = [C.parts(a, b, GRID)["A"] for a, b in zip(ss, ss[1:] + ss[:1])]
        cyc = all(t > 0 for t in As) or all(t < 0 for t in As)
        print(f"  {k:>2} {v:12.8f} {v/(3*LOG2/2):12.5f} "
              f"{sum(abs(t) for t in As)/3:10.7f}  directed cycle={cyc}"
              f"   {time.time()-t0:.0f}s")
        for s in ss:
            print(f"       {describe(s)}   r = {s.r:.6g}, sigma = {s.sigma:.5f}")
        rows.append(["max curl", k, v, 3 * LOG2 / 2, sum(abs(t) for t in As) / 3, cyc])
        if best is None or v > best[0]:
            best = (v, ss)

    with (HERE / "extremes.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["quantity", "k", "value", "bound", "extra1", "extra2"])
        for row in rows:
            wr.writerow(row)
    print(f"\nwrote {HERE/'extremes.csv'}")


if __name__ == "__main__":
    main()
