"""G3 -- how much curl is possible?

The structure theorem gives  A = d(psi) + D  with  |D| < (log 2)/2, sharp.
Since curl(d psi) = 0 on every cycle, curl A = curl D, so:

  * |curl A| <= k * (log 2)/2 around any k-cycle -- and on a *directed*
    (i.e. tournament) cycle curl A = sum |A|, so the mean strength of
    preference around any cycle is at most (log 2)/2 = 0.34657;
  * equivalently the geometric mean of C(b->a)/C(a->b) around any
    preference cycle is at most 2.

This script (i) reproduces the seed Hodge split, (ii) hill-climbs the
triangle curl, (iii) hill-climbs the global Hodge residual ||curl||/||A|| at
several n, (iv) measures the slack in the no-arbitrage inequality and in the
sharpened bound  |curl A| <= sum S - (max sigma - min sigma).

    python research/realizability/g3_curl.py
"""
from __future__ import annotations

import csv
import itertools
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import realize as R  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

HERE = Path(__file__).resolve().parent


def _A(x, n, r):
    return R._A_matrix(R._sigs_from(x, n, r))


def neg_triangle_curl(x, n, r):
    try:
        A = _A(x, n, r)
    except ValueError:
        return 1e3
    return -abs(A[0, 1] + A[1, 2] + A[2, 0])


def neg_hodge_ratio(x, n, r):
    try:
        A = _A(x, n, r)
    except ValueError:
        return 1e3
    na = np.linalg.norm(A)
    if na < 1e-9:
        return 0.0
    psi = -A.mean(axis=1)
    G = psi[None, :] - psi[:, None]
    return -float(np.linalg.norm(A - G) / na)


def neg_absD(x):
    """|D| for a single pair, each signature with 3 atoms."""
    try:
        a, b = R._sigs_from(x, 2, 3)
    except ValueError:
        return 1e3
    return -abs(C.parts(a, b, C.make_grid(-25.0, 25.0, 0.004))["D"])


def neg_eps(x):
    try:
        a, b = R._sigs_from(x, 2, 3)
    except ValueError:
        return 1e3
    return -C.parts(a, b, C.make_grid(-25.0, 25.0, 0.004))["eps"]


def climb(fun, dim, seed=0, maxiter=300, popsize=20, xmax=5.0, args=()):
    bounds = [(0.0, xmax)] * dim
    x, f = differential_evolution(fun, bounds, args=args, seed=seed,
                                  maxiter=maxiter, popsize=popsize,
                                  F=(0.3, 1.2), CR=0.9)
    for step in (0.3, 0.05, 5e-3, 5e-4):
        x, f = pattern_search(fun, x, args=args, step=step, min_step=1e-9,
                              maxiter=20000, bounds=bounds)
    return x, -f


def main():
    rows = []

    print("=== 1. seed reproduction: Hodge split of a random integer pool ===")
    import gpools
    pool = gpools.integer_pool(298, seed=11)
    for n in (8, 16, 24):
        rng = np.random.default_rng(n)
        sel = [pool[i] for i in rng.choice(len(pool), n, replace=False)]
        A = R._A_matrix(sel, C.make_grid(-25.0, 25.0, 0.004))
        grad, curl = C.hodge(A)
        cyc = len(C.three_cycles(A))
        print(f"  n={n:3d}  |grad|/|A|={grad:.4f}  |curl|/|A|={curl:.4f}  "
              f"3-cycles={cyc} of {math.comb(n,3)}")
        rows.append(["random pool", n, grad, curl, cyc])

    print("\n=== 2. extremal |D| and eps = P+Q for a single pair (3 atoms) ===")
    t0 = time.time()
    x, v = climb(neg_absD, 6, seed=5, maxiter=250)
    a, b = R._sigs_from(x, 2, 3)
    print(f"  max |D| found = {v:.9f}   bound (log2)/2 = {C.LOG2/2:.9f}"
          f"   ({100*v/(C.LOG2/2):.1f}% of the bound)   {time.time()-t0:.0f}s")
    print(f"    a = {[f'{math.exp(t):.6g}' for t in a.xs]}"
          f"  b = {[f'{math.exp(t):.6g}' for t in b.xs]}")
    rows.append(["max |D| (3 atoms)", 2, v, C.LOG2 / 2, ""])
    x, v = climb(neg_eps, 6, seed=6, maxiter=250)
    print(f"  max eps found = {v:.9f}   crude bound 2 log 2 = {2*C.LOG2:.9f}"
          f"   Lipschitz bound 1.1252")
    rows.append(["max eps (3 atoms)", 2, v, 2 * C.LOG2, ""])

    print("\n=== 3. the largest triangle curl ===")
    print(f"  {'r':>3} {'max |curl A|':>14} {'/ (3 log2/2)':>14} {'sum|A|/3':>12}")
    best_tri = None
    for r in (2, 3, 4, 5, 6):
        t0 = time.time()
        x, v = climb(neg_triangle_curl, 3 * r, seed=7 + r, maxiter=250,
                     args=(3, r))
        A = _A(x, 3, r)
        s = abs(A[0, 1]) + abs(A[1, 2]) + abs(A[2, 0])
        print(f"  {r:>3} {v:14.9f} {v/(3*C.LOG2/2):14.4f} {s/3:12.9f}"
              f"   cycle={abs(v - s) < 1e-9}  {time.time()-t0:.0f}s")
        rows.append([f"max triangle curl r={r}", 3, v, 3 * C.LOG2 / 2, s])
        if best_tri is None or v > best_tri[1]:
            best_tri = (x, v, r)

    print("\n=== 4. the largest Hodge residual ||curl||/||A|| ===")
    print(f"  {'n':>3} {'r':>3} {'max ratio':>12}  {'3-cycles':>9}")
    for n in (3, 4, 5, 6, 8, 12, 16, 24):
        r = 4
        t0 = time.time()
        x, v = climb(neg_hodge_ratio, n * r, seed=100 + n, maxiter=200,
                     args=(n, r))
        A = _A(x, n, r)
        cyc = len(C.three_cycles(A))
        print(f"  {n:>3} {r:>3} {v:12.6f}  {cyc:>5}/{math.comb(n,3):<5}"
              f"  {time.time()-t0:.0f}s")
        rows.append([f"max hodge ratio", n, v, "", cyc])

    print("\n=== 5. no arbitrage and the sharpened bound, on the known cycle ===")
    CYCLE = [C.Sig.of(t) for t in ((6, 3, 3), (7, 2, 1), (6, 5, 1))]
    tot_s = tot_a = 0.0
    sig_vals = [s.sigma for s in CYCLE]
    for a, b in zip(CYCLE, CYCLE[1:] + CYCLE[:1]):
        p = C.parts(a, b)
        tot_s += p["d"] / 2
        tot_a += p["A"]
    spread = max(sig_vals) - min(sig_vals)
    print(f"  |curl A|           = {abs(tot_a):.9f}")
    print(f"  sum S              = {tot_s:.9f}      (brief D(d), slack x"
          f"{tot_s/abs(tot_a):.1f})")
    print(f"  sum S - spread(sigma) = {tot_s - spread:.9f}  (sharpened, slack x"
          f"{(tot_s-spread)/abs(tot_a):.1f})")
    print(f"  3 (log2)/2         = {3*C.LOG2/2:.9f}      (universal, slack x"
          f"{(3*C.LOG2/2)/abs(tot_a):.1f})")
    rows.append(["known cycle", 3, abs(tot_a), tot_s, tot_s - spread])

    with (HERE / "g3_curl.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["what", "n", "value", "bound", "extra"])
        for row in rows:
            wr.writerow(row)
    print(f"\nwrote {HERE/'g3_curl.csv'}")


if __name__ == "__main__":
    main()
