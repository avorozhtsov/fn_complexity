"""Session brief I -- the other targets of FINDINGS Sec. 4.3, re-examined.

The C_4 argument generalises into a SUFFICIENT test for any target metric
delta on n points:

  (1) find an isometric copy of delta in the oscillation norm on m+1 nodes,
      i.e. a matrix Y (n x (m+1)) with  osc_k (Y[b,k] - Y[a,k]) = delta_ab;
  (2) put the nodes at theta_1 < ... < theta_{m+1} and let S_a be constant on
      each cell.  By Corollary 3.1 of OBSTRUCTION.md the only constraints are

        u_{a,k} = Y[a,k+1] - Y[a,k] + (kappa_{k+1} - kappa_k)  in (0, L_k),
        S_{a,k} = theta_k + delta_{L_k}^{-1}(u_{a,k})   nonincreasing in k,

      with kappa a free common profile and L_k = theta_{k+1} - theta_k free.
      A feasible point realises delta EXACTLY in the cone, hence (Theorem 1)
      with distortion 1 + O(1/log r) by signatures.

For the uniform metric K_n step (1) is free: Y = identity, since
osc(e_b - e_a) = 2 for all a != b.  For the graph metrics an oscillation
embedding is found numerically first.

    python research/realizability/i_targets.py
"""
from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import i_pattern as P  # noqa: E402
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


TARGETS = {
    "K_5": (1.0 - np.eye(5), None),
    "K_6": (1.0 - np.eye(6), None),
    "K_8": (1.0 - np.eye(8), None),
    "C_5": (graph_metric(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]), None),
    "K_2,3": (graph_metric(5, [(i, 2 + j) for i in range(2)
                               for j in range(3)]), None),
    "K_3,3": (graph_metric(6, [(i, 3 + j) for i in range(3)
                               for j in range(3)]), None),
    "Petersen": (graph_metric(10, [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
        (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
        (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]), None),
    "Wagner V8": (graph_metric(8, [(i, (i + 1) % 8) for i in range(8)]
                               + [(i, i + 4) for i in range(4)]), None),
}


# ---------------------------------------------------------------------------
# (1) an isometric copy in the oscillation norm
# ---------------------------------------------------------------------------

def osc_embed(delta, m, seed=0, restarts=6, maxiter=600):
    """Y (n x m) with osc_k (Y[b,k]-Y[a,k]) = delta_ab, up to a common scale."""
    n = delta.shape[0]
    iu = np.triu_indices(n, 1)

    def err(z):
        Y = z.reshape(n, m)
        Y = Y - Y[:, [0]]
        D = np.empty(len(iu[0]))
        for t, (i, j) in enumerate(zip(*iu)):
            v = Y[j] - Y[i]
            D[t] = v.max() - v.min()
        sc = (D * delta[iu]).sum() / (delta[iu] ** 2).sum()
        return float(np.abs(D - sc * delta[iu]).max() / max(sc, 1e-9))

    b = [(-1.5, 1.5)] * (n * m)
    best_z, best = None, math.inf
    for t in range(restarts):
        z, f = differential_evolution(err, b, seed=seed + 1237 * t,
                                      maxiter=maxiter, popsize=10,
                                      F=(0.3, 1.2), CR=0.9)
        for st in (0.3, 0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5,
                   3e-6, 1e-6, 3e-7, 1e-7):
            z, f = pattern_search(err, z, step=st, min_step=1e-13,
                                  maxiter=200000, bounds=b)
        if f < best:
            best_z, best = z, f
    Y = best_z.reshape(n, m)
    return Y - Y[:, [0]], best


# ---------------------------------------------------------------------------
# (2) feasibility of the node system for a given Y
# ---------------------------------------------------------------------------

def feasibility(z, Y):
    n, mm = Y.shape
    m = mm - 1                       # number of cells between the nodes
    L = np.exp(np.clip(z[:m], -12.0, 6.0))
    rho = z[m:2 * m]
    s = math.exp(np.clip(z[2 * m], -14.0, 3.0))
    theta = np.concatenate([[0.0], np.cumsum(L)])
    viol = 0.0
    S = np.zeros((n, m))
    for k in range(m):
        u = rho[k] + s * (Y[:, k + 1] - Y[:, k])
        for a in range(n):
            if u[a] <= 1e-12:
                viol += (1e-12 - u[a]) + 1e-9
            elif u[a] >= L[k] - 1e-12:
                viol += (u[a] - L[k] + 1e-12) + 1e-9
            else:
                S[a, k] = theta[k] + P.delta_inv(u[a], L[k])
    if viol > 0:
        return 1e3 + viol
    for a in range(n):
        for k in range(m - 1):
            if S[a, k] < S[a, k + 1]:
                viol += S[a, k + 1] - S[a, k]
    if viol > 0:
        return 1e3 + viol
    return -s


def test(Y, seed=0, restarts=4, maxiter=350):
    m = Y.shape[1] - 1
    b = [(-3.0, 3.0)] * m + [(-3.0, 3.0)] * m + [(-10.0, 1.0)]
    best = math.inf
    for t in range(restarts):
        z, f = differential_evolution(feasibility, b, args=(Y,),
                                      seed=seed + 733 * t, maxiter=maxiter,
                                      popsize=14, F=(0.3, 1.2), CR=0.9)
        for st in (0.4, 0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5):
            z, f = pattern_search(feasibility, z, args=(Y,), step=st,
                                  min_step=1e-12, maxiter=100000, bounds=b)
        best = min(best, f)
    return best


def antichain(n, m, k):
    """y_a = indicator of the a-th k-subset of {1..m}; osc(y_b - y_a) = 2."""
    sets = list(itertools.combinations(range(m), k))
    if len(sets) < n:
        return None
    Y = np.zeros((n, m))
    for a in range(n):
        for j in sets[a]:
            Y[a, j] = 1.0
    return Y


def main(quick=False):
    print("=== self-check: the C_4 pattern through the same routine ===")
    f = test(P.T0[:, [0, 1, 3, 2]], seed=5, restarts=2, maxiter=200)
    print("  C_4: " + (f"REALISABLE, s = {-f:.8f}" if f < 0
                       else "infeasible (BUG: it must be feasible)"), flush=True)

    print("\n=== uniform K_n, identity ansatz (n nodes, Y = I) ===")
    for n in ((4, 5, 6) if quick else (4, 5, 6, 8)):
        f = test(np.eye(n), seed=31 * n, restarts=2, maxiter=200)
        print(f"  K_{n}: " + (f"REALISABLE, s = {-f:.8f}" if f < 0
                             else "this ansatz infeasible"), flush=True)

    print("\n=== uniform K_n, antichain ansatz (k-subsets of m nodes) ===")
    combos = ((5, 4, 2), (6, 4, 2)) if quick else (
        (5, 4, 2), (6, 4, 2), (8, 5, 2), (6, 5, 2), (8, 6, 3))
    for n, m, k in combos:
        Y = antichain(n, m, k)
        if Y is None:
            continue
        f = test(Y, seed=101 + 7 * n + m, restarts=2, maxiter=200)
        print(f"  K_{n}  (m = {m}, k = {k}): "
              + (f"REALISABLE, s = {-f:.8f}" if f < 0
                 else "this ansatz infeasible"), flush=True)

    if quick:
        return
    print("\n=== graph metrics: an oscillation embedding, then the node test ===")
    for name, (delta, _) in TARGETS.items():
        n = delta.shape[0]
        if np.allclose(delta, 1.0 - np.eye(n)) or n > 6:
            continue
        best = (math.inf, None)
        for m in (n, n + 1):
            Y, e = osc_embed(delta, m, seed=17 * n + m, restarts=3, maxiter=350)
            if e < best[0]:
                best = (e, Y)
        emb, Y = best
        if emb > 1e-7:
            print(f"  {name:<8}: no exact oscillation embedding found "
                  f"(residual {emb:.2e})", flush=True)
            continue
        f = test(Y, seed=31 * n, restarts=2, maxiter=200)
        print(f"  {name:<8}: " + (f"REALISABLE, s = {-f:.8f}" if f < 0
                                 else "this ansatz infeasible"), flush=True)


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)
