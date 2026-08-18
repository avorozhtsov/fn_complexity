"""Session brief I -- metric realisation searched inside the tropical cone.

Because  inf over signatures  =  inf over the cone C  (every F_a lies in C, and
every element of C is a locally uniform limit of (1/K)F_{a^(K)} -- see
i_cone.py and i_validate_cone.py), the realisability question for a target
metric can be asked directly in C, where

  * d(a,b) is an exact finite maximum (no grid, no Lipschitz bracket), and
  * the parameters are (intercept, slope) pairs, which is a far better
    conditioned search than atoms of a signature.

    python research/realizability/i_search.py [target ...]
"""
from __future__ import annotations

import itertools
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import i_cone as T  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402


def graph_metric(n, edges):
    INF = 10 ** 6
    D = np.full((n, n), INF, dtype=float)
    np.fill_diagonal(D, 0)
    for i, j in edges:
        D[i, j] = D[j, i] = 1
    for k in range(n):
        D = np.minimum(D, D[:, k][:, None] + D[k][None, :])
    return D


def targets():
    t = {}
    for n in (4, 5, 6, 8):
        M = np.ones((n, n))
        np.fill_diagonal(M, 0)
        t[f"K_{n}"] = M
    t["C_4"] = graph_metric(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
    t["C_5"] = graph_metric(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
    t["K_2,3"] = graph_metric(5, [(i, 2 + j) for i in range(2) for j in range(3)])
    t["K_3,3"] = graph_metric(6, [(i, 3 + j) for i in range(3) for j in range(3)])
    t["path P_3"] = graph_metric(3, [(0, 1), (1, 2)])
    t["path P_4"] = graph_metric(4, [(0, 1), (1, 2), (2, 3)])
    t["star K_1,3"] = graph_metric(4, [(0, 1), (0, 2), (0, 3)])
    return t


def unpack(z, n, k):
    """z -> n tropical functions with k lines each.  c = z[..,0], x = z[..,1]."""
    z = np.asarray(z, float).reshape(n, k, 2)
    out = []
    for row in z:
        c = np.maximum(row[:, 0], 0.0)
        x = np.maximum(row[:, 1], 0.0)
        if c.max() <= 1e-9:
            c = c + 1e-3
        if x.max() <= 1e-9:
            x = x + 1e-3
        out.append(T.Trop(c, x))
    return out


def distortion(z, n, k, delta, iu):
    try:
        fs = unpack(z, n, k)
    except ValueError:
        return 1e6
    D = T.dmatrix(fs)
    ratios = D[iu] / delta[iu]
    if ratios.min() <= 1e-13:
        return 1e6
    return float(ratios.max() / ratios.min())


def search(delta, k=4, seed=0, maxiter=400, restarts=3, xmax=6.0):
    n = delta.shape[0]
    iu = np.triu_indices(n, 1)
    bounds = [(0.0, xmax)] * (n * k * 2)
    best_z, best = None, math.inf
    for t in range(restarts):
        z, f = differential_evolution(distortion, bounds,
                                      args=(n, k, delta, iu), seed=seed + 4021 * t,
                                      maxiter=maxiter, popsize=16, F=(0.3, 1.2),
                                      CR=0.9)
        for step in (0.4, 0.1, 0.02, 4e-3, 8e-4, 1.6e-4, 3e-5, 6e-6, 1e-6):
            z, f = pattern_search(distortion, z, args=(n, k, delta, iu),
                                  step=step, min_step=1e-12, maxiter=60000,
                                  bounds=bounds)
        if f < best:
            best_z, best = z, f
    return best_z, best


def main(names=None):
    tg = targets()
    if names:
        tg = {nm: tg[nm] for nm in names}
    print(f"  {'target':<12} {'n':>3} {'lines k':>8} {'distortion':>14} "
          f"{'scale':>10}   time")
    for name, delta in tg.items():
        n = delta.shape[0]
        iu = np.triu_indices(n, 1)
        for k in (2, 3, 4, 6):
            t0 = time.time()
            z, f = search(delta, k=k, seed=101 + 7 * n + k,
                          maxiter=350 if n <= 6 else 200,
                          restarts=3 if n <= 6 else 2)
            fs = unpack(z, n, k)
            D = T.dmatrix(fs)
            sc = float(np.median(D[iu] / delta[iu]))
            print(f"  {name:<12} {n:>3} {k:>8} {f:14.8f} {sc:10.4f}   "
                  f"{time.time()-t0:.0f}s")
            np.save(Path(__file__).resolve().parent /
                    f"i_best_{name.replace(',', '').replace(' ', '_')}_{k}.npy", z)


if __name__ == "__main__":
    main(sys.argv[1:] or None)
