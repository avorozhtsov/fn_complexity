"""Session brief I -- the C_4 question, attacked in the free parametrisation.

An element of the cone C is, up to the projective scaling that the Hilbert
metric quotients out, a pair of unconstrained sequences

    centres      sigma_1 > sigma_2 > ... > sigma_k
    breakpoints  theta_1 < theta_2 < ... < theta_{k-1}

and  y(theta) = log x_j + logaddexp(sigma_j, theta)  on the j-th piece, with
log x_1 = 0 and  log x_{j+1} = log x_j + logaddexp(sigma_j, theta_j)
                                        - logaddexp(sigma_{j+1}, theta_j).

d(a,b) = osc of y_b - y_a over the FINITE set {breakpoints} u {-inf, +inf},
because y_b - y_a is monotone between consecutive breakpoints.

    python research/realizability/i_c4.py [--k K] [--restarts N]
"""
from __future__ import annotations

import argparse
import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from optimizers import differential_evolution, pattern_search  # noqa: E402

LOG2 = math.log(2.0)


# ---------------------------------------------------------------------------
# a cone element in the free parametrisation
# ---------------------------------------------------------------------------

class Elt:
    __slots__ = ("sig", "th", "lx")

    def __init__(self, sig, th):
        self.sig = np.asarray(sig, float)
        self.th = np.asarray(th, float)
        k = len(self.sig)
        lx = np.zeros(k)
        for j in range(k - 1):
            lx[j + 1] = (lx[j]
                         + np.logaddexp(self.sig[j], self.th[j])
                         - np.logaddexp(self.sig[j + 1], self.th[j]))
        self.lx = lx

    def y(self, theta):
        """log Phi at the given theta (vectorised)."""
        theta = np.atleast_1d(np.asarray(theta, float))
        j = np.searchsorted(self.th, theta, side="right")
        return self.lx[j] + np.logaddexp(self.sig[j], theta)

    @property
    def y_minus(self):          # log R
        return float(self.sig[0])

    @property
    def y_plus(self):           # log Lambda  (value of y(theta) - theta at +inf)
        return float(self.lx[-1])

    @property
    def sigma(self):
        return self.y_minus - self.y_plus


def make(z, k):
    """z (length 2k-1) -> Elt.  Unconstrained: gaps are exp()'d."""
    sig = np.empty(k)
    sig[0] = z[0]
    for j in range(1, k):
        sig[j] = sig[j - 1] - math.exp(np.clip(z[j], -25.0, 6.0))
    th = np.empty(k - 1)
    if k > 1:
        th[0] = z[k]
        for j in range(1, k - 1):
            th[j] = th[j - 1] + math.exp(np.clip(z[k + j], -25.0, 6.0))
    return Elt(sig, th)


def dist(a: Elt, b: Elt):
    pts = np.unique(np.concatenate([a.th, b.th]))
    if len(pts):
        v = b.y(pts) - a.y(pts)
        hi, lo = float(v.max()), float(v.min())
    else:
        hi, lo = -math.inf, math.inf
    e0 = b.y_minus - a.y_minus
    e1 = b.y_plus - a.y_plus
    hi = max(hi, e0, e1)
    lo = min(lo, e0, e1)
    return hi - lo


def parts(a: Elt, b: Elt):
    pts = np.unique(np.concatenate([a.th, b.th]))
    v = b.y(pts) - a.y(pts) if len(pts) else np.array([])
    e0 = b.y_minus - a.y_minus
    e1 = b.y_plus - a.y_plus
    hi = max(float(v.max()) if len(v) else -math.inf, e0, e1)
    lo = min(float(v.min()) if len(v) else math.inf, e0, e1)
    return dict(d=hi - lo, A=0.5 * (hi + lo), P=hi - max(e0, e1),
                Q=min(e0, e1) - lo, dsigma=b.sigma - a.sigma)


def dmat(elts):
    n = len(elts)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = dist(elts[i], elts[j])
    return D


# ---------------------------------------------------------------------------
# targets and the distortion objective
# ---------------------------------------------------------------------------

def graph_metric(n, edges):
    INF = 10 ** 6
    D = np.full((n, n), INF, dtype=float)
    np.fill_diagonal(D, 0)
    for i, j in edges:
        D[i, j] = D[j, i] = 1
    for kk in range(n):
        D = np.minimum(D, D[:, kk][:, None] + D[kk][None, :])
    return D


C4 = graph_metric(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
K4 = 1.0 - np.eye(4)
C5 = graph_metric(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
K23 = graph_metric(5, [(i, 2 + j) for i in range(2) for j in range(3)])


def distortion(z, n, k, delta, iu):
    elts = [make(z[i * (2 * k - 1):(i + 1) * (2 * k - 1)], k) for i in range(n)]
    D = dmat(elts)
    ratios = D[iu] / delta[iu]
    if ratios.min() <= 1e-14:
        return 1e6
    return float(ratios.max() / ratios.min())


def run(delta, k, seed, maxiter=600, restarts=6, pop=16):
    n = delta.shape[0]
    iu = np.triu_indices(n, 1)
    dim = n * (2 * k - 1)
    bounds = []
    for _ in range(n):
        bounds.append((-3.0, 3.0))               # sigma_1
        bounds += [(-6.0, 2.0)] * (k - 1)        # log gaps of sigma
        if k > 1:
            bounds.append((-6.0, 6.0))           # theta_1
            bounds += [(-6.0, 2.0)] * (k - 2)    # log gaps of theta
    assert len(bounds) == dim
    best_z, best = None, math.inf
    for t in range(restarts):
        z, f = differential_evolution(distortion, bounds, args=(n, k, delta, iu),
                                      seed=seed + 5171 * t, maxiter=maxiter,
                                      popsize=pop, F=(0.3, 1.2), CR=0.9)
        for step in (0.5, 0.15, 0.04, 0.01, 2e-3, 5e-4, 1e-4, 2e-5, 5e-6, 1e-6,
                     2e-7, 5e-8):
            z, f = pattern_search(distortion, z, args=(n, k, delta, iu),
                                  step=step, min_step=1e-14, maxiter=200000,
                                  bounds=bounds)
        if f < best:
            best_z, best = z, f
    return best_z, best


def report(z, n, k, delta, name):
    iu = np.triu_indices(n, 1)
    elts = [make(z[i * (2 * k - 1):(i + 1) * (2 * k - 1)], k) for i in range(n)]
    D = dmat(elts)
    ratios = D[iu] / delta[iu]
    sc = ratios.min()
    print(f"\n  --- {name}, k = {k}, distortion = "
          f"{ratios.max()/ratios.min():.9f} ---")
    print("  d matrix / target ratio:")
    for i in range(n):
        print("   " + "  ".join(f"{D[i, j]:9.6f}" for j in range(n)))
    print("  ratios d/delta: " + "  ".join(f"{r:.6f}" for r in ratios))
    print(f"  scale (min ratio) = {sc:.6f}")
    for i, e in enumerate(elts):
        print(f"   elt {i}: sigma = "
              + ", ".join(f"{v:+.5f}" for v in e.sig)
              + " | theta = " + ", ".join(f"{v:+.5f}" for v in e.th)
              + f" | sigma_a = {e.sigma:+.5f}")
    print("  pairwise P, Q, dsigma:")
    for i, j in zip(*iu):
        p = parts(elts[i], elts[j])
        print(f"   ({i},{j}) d={p['d']:.6f}  P={p['P']:.6f}  Q={p['Q']:.6f}  "
              f"dsigma={p['dsigma']:+.6f}  d-|dsigma|={p['d']-abs(p['dsigma']):.6f}")
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", type=int, nargs="*", default=[2, 3, 4, 5, 6, 8])
    ap.add_argument("--restarts", type=int, default=6)
    ap.add_argument("--target", default="C_4")
    args = ap.parse_args()
    delta = {"C_4": C4, "K_4": K4, "C_5": C5, "K_2,3": K23}[args.target]
    n = delta.shape[0]
    best_overall = (math.inf, None, None)
    for k in args.ks:
        z, f = run(delta, k, seed=9001 + 37 * k, restarts=args.restarts)
        print(f"  {args.target}  k={k}  distortion = {f:.10f}", flush=True)
        if f < best_overall[0]:
            best_overall = (f, z, k)
    f, z, k = best_overall
    report(z, n, k, delta, args.target)
    np.save(Path(__file__).resolve().parent /
            f"i_c4_best_{args.target.replace(',', '')}.npy",
            np.concatenate([[k], z]))
    print(f"\n  BEST over k: {f:.10f}")


if __name__ == "__main__":
    main()
