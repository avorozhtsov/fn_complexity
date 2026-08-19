"""G2 -- which metrics does the exchange metric realise?

Two parts.

(a) The oscillation count.  Brief G asks how many sign changes phi = u_b - u_a
    can have, because that caps how many independent directions the sup-norm
    can use.  Measured here directly, together with the convexity of U in
    s = log beta on which the count depends.

(b) Realisation.  For a target metric delta, minimise the distortion
    max_ij (d_ij/delta_ij) / min_ij (d_ij/delta_ij) over signature families:
    this is scale-free, so it answers "realisable up to a scale factor".
    Targets: the uniform metric K_n, the cycles C_4 and C_5, the complete
    bipartite graph metrics K_{2,3} and K_{3,3}, the Petersen graph, and a
    3-regular expander -- i.e. exactly the metrics of T1.2/T1.3.

    python research/realizability/g2_metrics.py
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
import gpools  # noqa: E402
import realize as R  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

HERE = Path(__file__).resolve().parent


# --------------------------------------------------------------------------
# (a) oscillation count
# --------------------------------------------------------------------------

def count_extrema(a, b, grid, thresh):
    """Interior local extrema of phi = U_b - U_a with amplitude >= thresh."""
    phi = b.U(grid) - a.U(grid)
    dp = np.diff(phi)
    sgn = np.sign(dp)
    sgn = sgn[sgn != 0]
    if len(sgn) == 0:
        return 0
    turns = np.flatnonzero(np.diff(sgn) != 0)
    if len(turns) == 0:
        return 0
    # keep only turning points whose excursion from the neighbouring turning
    # points exceeds thresh
    idx = [0] + [int(t) + 1 for t in turns] + [len(phi) - 1]
    vals = phi[idx]
    keep = 0
    for k in range(1, len(vals) - 1):
        amp = min(abs(vals[k] - vals[k - 1]), abs(vals[k] - vals[k + 1]))
        if amp >= thresh:
            keep += 1
    return keep


def convexity_of_U(sigs, grid):
    """max violation of U'' >= 0 in s (relative to the grid's own error)."""
    h = grid[1] - grid[0]
    worst = 0.0
    for s in sigs:
        U = s.U(grid)
        d2 = (U[2:] - 2 * U[1:-1] + U[:-2]) / (h * h)
        worst = max(worst, float(-d2.min()))
    return worst


def part_a():
    grid = C.make_grid(-25.0, 25.0, 0.002)
    print("=== (a) shape of U and the oscillation count ===")
    rows = []
    for name, pool in (
            ("random r<=7, values<=40", gpools.integer_pool(60, seed=11)),
            ("random r<=12, values<=400", gpools.integer_pool(60, seed=12, rmax=12, vmax=400)),
            ("random r<=20, values<=10^5", gpools.integer_pool(60, seed=13, rmax=20, vmax=10 ** 5)),
    ):
        cv = convexity_of_U(pool, grid)
        counts = {}
        rmax_seen = 0
        for a, b in itertools.combinations(pool, 2):
            k = count_extrema(a, b, grid, 1e-9)
            counts[k] = counts.get(k, 0) + 1
            rmax_seen = max(rmax_seen, len(a.xs) + len(b.xs))
        print(f"  {name}")
        print(f"    max(-U'') on the grid = {cv:+.3e}  (grid second-difference "
              f"noise ~ 1e-6)")
        print(f"    interior local extrema of phi, over "
              f"{sum(counts.values())} pairs: "
              + ", ".join(f"{k}:{v}" for k, v in sorted(counts.items())))
        rows.append(["oscillation", name, cv, max(counts), sum(counts.values())])
    return rows


# --------------------------------------------------------------------------
# (b) realisation of prescribed metrics
# --------------------------------------------------------------------------

def graph_metric(n, edges):
    INF = 10 ** 6
    D = np.full((n, n), INF)
    np.fill_diagonal(D, 0)
    for i, j in edges:
        D[i, j] = D[j, i] = 1
    for k in range(n):
        D = np.minimum(D, D[:, k][:, None] + D[k][None, :])
    return D.astype(float)


def targets():
    t = {}
    for n in (4, 6, 8):
        M = np.ones((n, n))
        np.fill_diagonal(M, 0)
        t[f"uniform K_{n}"] = M
    t["C_4"] = graph_metric(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
    t["C_5"] = graph_metric(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
    t["K_2,3"] = graph_metric(5, [(i, 2 + j) for i in range(2) for j in range(3)])
    t["K_3,3"] = graph_metric(6, [(i, 3 + j) for i in range(3) for j in range(3)])
    t["Petersen"] = graph_metric(10, [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
        (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
        (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)])
    # 3-regular expander on 8 nodes: the cube plus a perfect matching (Mobius-Kantor
    # style) -- the Wagner graph V8
    t["Wagner V8"] = graph_metric(8, [(i, (i + 1) % 8) for i in range(8)]
                                  + [(i, i + 4) for i in range(4)])
    return t


def distortion(x, n, r, delta, iu):
    try:
        sigs = R._sigs_from(x, n, r)
    except ValueError:
        return 1e6
    grid = R.GRID
    tab = np.vstack([s.U(grid) for s in sigs])
    eR = np.array([math.log(s.R) for s in sigs])
    eL = np.array([math.log(s.Lam) for s in sigs])
    d = np.zeros((n, n))
    for i in range(n):
        diff = tab - tab[i]
        hi = np.maximum(diff.max(axis=1), np.maximum(eR - eR[i], eL - eL[i]))
        lo = np.minimum(diff.min(axis=1), np.minimum(eR - eR[i], eL - eL[i]))
        d[i] = hi - lo
    ratios = d[iu] / delta[iu]
    if ratios.min() <= 1e-12:
        return 1e6
    return float(ratios.max() / ratios.min())


def realise_metric(delta, r, seed=0, maxiter=300, restarts=2, xmax=5.0):
    n = delta.shape[0]
    iu = np.triu_indices(n, 1)
    bounds = [(0.0, xmax)] * (n * r)
    best_x, best = None, math.inf
    for k in range(restarts):
        x, f = differential_evolution(distortion, bounds, args=(n, r, delta, iu),
                                      seed=seed + 613 * k, maxiter=maxiter,
                                      popsize=18, F=(0.3, 1.2), CR=0.9)
        for step in (0.3, 0.05, 5e-3, 5e-4):
            x, f = pattern_search(distortion, x, args=(n, r, delta, iu),
                                  step=step, min_step=1e-9, maxiter=20000,
                                  bounds=bounds)
        if f < best:
            best_x, best = x, f
    return best_x, best


def part_b(rmax=6):
    print("\n=== (b) realisation of prescribed metrics, up to scale ===")
    print(f"  {'target':<14} {'n':>3} {'r':>3} {'best distortion':>16} "
          f"{'scale (d/delta)':>16}   time")
    rows = []
    for name, delta in targets().items():
        n = delta.shape[0]
        best = (math.inf, None, None)
        for r in (3, 4, rmax):
            if r > rmax:
                continue
            t0 = time.time()
            x, f = realise_metric(delta, r, seed=17 + n, maxiter=250 if n <= 6 else 150)
            sc = ""
            if f < best[0]:
                best = (f, x, r)
            sigs = R._sigs_from(x, n, r)
            iu = np.triu_indices(n, 1)
            dm, _ = C.matrices(sigs, R.GRID)
            sc = f"{np.median(dm[iu] / delta[iu]):.4f}"
            print(f"  {name:<14} {n:>3} {r:>3} {f:16.6f} {sc:>16}   "
                  f"{time.time()-t0:.0f}s")
            rows.append(["metric", name, n, r, f, sc])
    return rows


def main():
    rows = part_a()
    rows += part_b()
    with (HERE / "g2_metrics.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind", "name", "a", "b", "value", "extra"])
        for row in rows:
            wr.writerow(row)
    print(f"\nwrote {HERE/'g2_metrics.csv'}")


if __name__ == "__main__":
    main()
