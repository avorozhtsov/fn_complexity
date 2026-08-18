"""Session brief I -- the two sharp-constant questions, asked inside the cone.

Brief G left open:

  (i)  is  eps = P + Q  bounded by  log(1+e^{-Delta})  rather than
       2 log(1+e^{-Delta})?  Every signature search returned 0.9889*log 2 at
       Delta = 0, always with one of P, Q exactly zero.

  (ii) is the maximum triangle curl  (log 2)/2  rather than  3(log 2)/2 ?

Both are suprema over the achievable set, hence -- since the tropical cone C is
the projective closure of that set and d is continuous along the tropical limit
-- suprema over C.  Inside C the objective is an exact finite maximum, so the
hill-climb is over an honest smooth-ish landscape rather than over atoms.

    python research/realizability/i_constants.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import i_cone as T  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

LOG2 = math.log(2.0)


def parts(a: T.Trop, b: T.Trop):
    bs = T.candidates(a, b)
    v = T.phi(a, b, bs)
    e0 = math.log(b.R) - math.log(a.R)
    e1 = math.log(b.Lam) - math.log(a.Lam)
    hi = max(float(v.max()), e0, e1)
    lo = min(float(v.min()), e0, e1)
    P = hi - max(e0, e1)
    Q = min(e0, e1) - lo
    return P, Q, hi - lo, 0.5 * (hi + lo), abs(b.sigma - a.sigma)


def _mk(z, k):
    z = np.asarray(z, float).reshape(-1, k, 2)
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


def neg_eps(z, k, dtarget=None):
    try:
        a, b = _mk(z, k)
    except ValueError:
        return 1e3
    P, Q, _, _, dsig = parts(a, b)
    if dtarget is not None and abs(dsig - dtarget) > 1e-9:
        return 1e3 + abs(dsig - dtarget)
    return -(P + Q)


def neg_minPQ(z, k):
    try:
        a, b = _mk(z, k)
    except ValueError:
        return 1e3
    P, Q, _, _, _ = parts(a, b)
    return -min(P, Q)


def neg_D(z, k):
    try:
        a, b = _mk(z, k)
    except ValueError:
        return 1e3
    P, Q, _, _, _ = parts(a, b)
    return -abs(P - Q) / 2.0


def neg_curl(z, k):
    try:
        fs = _mk(z, k)
    except ValueError:
        return 1e3
    A = T.amatrix(fs)
    return -abs(A[0, 1] + A[1, 2] + A[2, 0])


def climb(fun, nfun, k, seed, args=(), maxiter=500, restarts=6, xmax=8.0):
    bounds = [(0.0, xmax)] * (nfun * k * 2)
    best_z, best = None, math.inf
    for t in range(restarts):
        z, f = differential_evolution(fun, bounds, args=args, seed=seed + 911 * t,
                                      maxiter=maxiter, popsize=14,
                                      F=(0.3, 1.2), CR=0.9)
        for step in (0.5, 0.1, 0.02, 4e-3, 8e-4, 1.6e-4, 3e-5, 6e-6, 1e-6, 2e-7):
            z, f = pattern_search(fun, z, args=args, step=step, min_step=1e-13,
                                  maxiter=80000, bounds=bounds)
        if f < best:
            best_z, best = z, f
    return best_z, best


def main():
    print("=== (i) sup of eps = P + Q over the cone, at Delta = 0 by "
          "construction ===")
    print(f"  {'k lines':>8} {'sup eps':>14} {'/log2':>9} {'P':>12} {'Q':>12} "
          f"{'|dsigma|':>10}")
    for k in (2, 3, 4, 5, 6):
        z, f = climb(neg_eps, 2, k, seed=31 + k, args=(k,))
        a, b = _mk(z, k)
        P, Q, d, A, dsig = parts(a, b)
        print(f"  {k:>8} {-f:14.9f} {-f/LOG2:9.5f} {P:12.8f} {Q:12.8f} "
              f"{dsig:10.6f}", flush=True)

    print("\n=== is min(P,Q) always 0?  maximise min(P,Q) over the cone ===")
    print(f"  {'k lines':>8} {'max min(P,Q)':>16}")
    for k in (2, 3, 4, 5, 6):
        z, f = climb(neg_minPQ, 2, k, seed=77 + k, args=(k,))
        print(f"  {k:>8} {-f:16.10f}", flush=True)

    print("\n=== (ii) sup |D| and sup |triangle curl| over the cone ===")
    print(f"  {'k lines':>8} {'sup |D|':>14} {'/(log2/2)':>11}")
    for k in (2, 3, 4, 6):
        z, f = climb(neg_D, 2, k, seed=131 + k, args=(k,))
        print(f"  {k:>8} {-f:14.9f} {-f/(LOG2/2):11.5f}", flush=True)
    print(f"  {'k lines':>8} {'sup |curl|':>14} {'/(log2/2)':>11} "
          f"{'/(3log2/2)':>11}")
    for k in (2, 3, 4, 6):
        z, f = climb(neg_curl, 3, k, seed=181 + k, args=(k,), restarts=8)
        print(f"  {k:>8} {-f:14.9f} {-f/(LOG2/2):11.5f} {-f/(1.5*LOG2):11.5f}",
              flush=True)


if __name__ == "__main__":
    main()
