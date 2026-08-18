"""Session brief I -- the exact finite model of the exchange metric.

THE SCALE FUNCTION.  For any signature put  F = log Z,  theta = log beta, and

    T(beta) = ( F(beta) - beta F'(beta) ) / F'(beta)      (intercept / slope of
                                                           the tangent to F)
    S(theta) = log T(e^theta).

Then (proved in OBSTRUCTION.md Sec. 2)

    U'(theta) = sigmoid( theta - S(theta) ),     S nonincreasing,

so  sign (U_b - U_a)'(theta) = sign ( S_a(theta) - S_b(theta) ),  and U is
determined by S up to an additive constant.  Conversely EVERY nonincreasing S
arises (Sec. 2.3), so the exchange metric on n points is exactly

    d(a,b) = osc over theta of  int ( sigmoid(t - S_b) - sigmoid(t - S_a) ) dt

over n nonincreasing functions S_1, ..., S_n.

FINITE MODEL.  Take S_a constant on each of m+1 cells cut by common nodes
theta_1 < ... < theta_m.  Then every U_b - U_a is monotone on each cell, so its
oscillation is a maximum over the finite set {-inf, theta_1, ..., theta_m, +inf}
and the model is EXACT for that many cells; letting m grow exhausts the cone.

    python research/realizability/i_nodes.py [--m M] [--restarts N] [--target T]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from optimizers import differential_evolution, pattern_search  # noqa: E402


def sp(t):
    return np.logaddexp(0.0, t)


def y_table(theta, S):
    """Values of y_a at (-inf, theta_1..theta_m, +inf-slope), given S.

    theta : (m,)          nodes
    S     : (n, m+1)      S_a on cell k = (theta_k, theta_{k+1}), k = 0..m
    returns (n, m+2) with column 0 = y(-inf), 1..m = y(theta_k),
    column m+1 = lim (y(theta) - theta).
    """
    n = S.shape[0]
    m = len(theta)
    Y = np.zeros((n, m + 2))
    Y[:, 1] = 0.0
    for k in range(1, m):
        Y[:, k + 1] = Y[:, k] + sp(theta[k] - S[:, k]) - sp(theta[k - 1] - S[:, k])
    Y[:, 0] = Y[:, 1] - sp(theta[0] - S[:, 0])
    Y[:, m + 1] = Y[:, m] - sp(theta[m - 1] - S[:, m]) - S[:, m]
    return Y


def dmat(theta, S):
    Y = y_table(theta, S)
    n = Y.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        diff = Y - Y[i]
        D[i] = diff.max(axis=1) - diff.min(axis=1)
    return D


# ---------------------------------------------------------------------------

def unpack(z, n, m):
    theta = np.empty(m)
    theta[0] = z[0]
    for k in range(1, m):
        theta[k] = theta[k - 1] + math.exp(np.clip(z[k], -20.0, 4.0))
    S = np.empty((n, m + 1))
    off = m
    for a in range(n):
        S[a, 0] = z[off]
        for k in range(1, m + 1):
            S[a, k] = S[a, k - 1] - math.exp(np.clip(z[off + k], -20.0, 4.0))
        off += m + 1
    return theta, S


def distortion(z, n, m, delta, iu):
    theta, S = unpack(z, n, m)
    D = dmat(theta, S)
    r = D[iu] / delta[iu]
    if r.min() <= 1e-14:
        return 1e6
    return float(r.max() / r.min())


def bounds_for(n, m):
    b = [(-4.0, 4.0)] + [(-6.0, 2.0)] * (m - 1)
    for _ in range(n):
        b += [(-4.0, 6.0)] + [(-8.0, 2.0)] * m
    return b


def run(delta, m, seed, maxiter=700, restarts=8, pop=14, init=None):
    n = delta.shape[0]
    iu = np.triu_indices(n, 1)
    b = bounds_for(n, m)
    best_z, best = None, math.inf
    for t in range(restarts):
        z, f = differential_evolution(distortion, b, args=(n, m, delta, iu),
                                      seed=seed + 3319 * t, maxiter=maxiter,
                                      popsize=pop, F=(0.3, 1.2), CR=0.9,
                                      init=init if t == 0 else None)
        for step in (0.6, 0.2, 0.06, 0.02, 6e-3, 2e-3, 6e-4, 2e-4, 6e-5, 2e-5,
                     6e-6, 2e-6, 6e-7, 2e-7, 6e-8, 2e-8):
            z, f = pattern_search(distortion, z, args=(n, m, delta, iu),
                                  step=step, min_step=1e-15, maxiter=400000,
                                  bounds=b)
        if f < best:
            best_z, best = z, f
    return best_z, best


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
    "C_4": graph_metric(4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "K_4": 1.0 - np.eye(4),
    "K_5": 1.0 - np.eye(5),
    "K_6": 1.0 - np.eye(6),
    "C_5": graph_metric(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
    "K_2,3": graph_metric(5, [(i, 2 + j) for i in range(2) for j in range(3)]),
    "P_4": graph_metric(4, [(0, 1), (1, 2), (2, 3)]),
    "K_1,3": graph_metric(4, [(0, 1), (0, 2), (0, 3)]),
}


def report(z, n, m, delta, name):
    iu = np.triu_indices(n, 1)
    theta, S = unpack(z, n, m)
    D = dmat(theta, S)
    r = D[iu] / delta[iu]
    print(f"\n  --- {name}, m = {m} nodes, distortion = {r.max()/r.min():.10f} ---")
    print("  nodes theta: " + ", ".join(f"{t:+.5f}" for t in theta))
    for a in range(n):
        print(f"   S_{a} = " + ", ".join(f"{v:+.5f}" for v in S[a]))
    print("  d matrix:")
    for i in range(n):
        print("   " + "  ".join(f"{D[i, j]:10.7f}" for j in range(n)))
    print("  ratios: " + "  ".join(f"{v:.7f}" for v in r))
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms", type=int, nargs="*", default=[3, 4, 5, 6, 8])
    ap.add_argument("--restarts", type=int, default=8)
    ap.add_argument("--target", default="C_4")
    args = ap.parse_args()
    delta = TARGETS[args.target]
    n = delta.shape[0]
    best = (math.inf, None, None)
    for m in args.ms:
        z, f = run(delta, m, seed=2207 + 53 * m, restarts=args.restarts)
        print(f"  {args.target}  m={m}  distortion = {f:.10f}", flush=True)
        if f < best[0]:
            best = (f, z, m)
    f, z, m = best
    report(z, n, m, delta, args.target)
    np.save(Path(__file__).resolve().parent /
            f"i_nodes_best_{args.target.replace(',', '')}.npy",
            np.concatenate([[m], z]))
    print(f"\n  BEST over m: {f:.10f}")


if __name__ == "__main__":
    main()
