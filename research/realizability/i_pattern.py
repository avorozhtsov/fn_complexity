"""Session brief I -- the forced 4-node pattern for C_4, and its feasibility.

WHY FOUR NODES.  Suppose  d = s.C_4  is realised, with the 4-cycle 1-2-3-4.
Let beta+ / beta- be a maximiser / minimiser of phi_13 and gamma+ / gamma- of
phi_24.  Because  phi_13 = phi_12 + phi_23 = phi_14 + phi_43  and each summand
has oscillation exactly s while phi_13 has 2s, beta+ maximises and beta-
minimises EVERY one of phi_12, phi_23, phi_14, phi_43 simultaneously; likewise
gamma+/gamma- for phi_21, phi_14, phi_23, phi_34.  Reading the six values off
gives, modulo an additive constant per point (the projective freedom) and a
common profile kappa (a function of the node),

     y_a(P_k) = kappa_k + c_a + s * T[a][k],
     T = [[0,0,0,0],[1,0,0,1],[2,1,0,1],[1,1,0,0]],
     columns (P_1,...,P_4) = (beta+, gamma+, beta-, gamma-),

and the four points are pairwise distinct.  This module tests, for each of the
24 possible orderings of the four points along theta, whether that pattern is
realisable by a cone element with breakpoints only at those four nodes -- i.e.
with S_a constant on each cell.  On cell k the increment of y_a is

     u_{a,k} = delta_{L_k}( S_{a,k} - theta_k ),
     delta_L(v) = softplus(L - v) - softplus(-v),  a decreasing bijection
                  R -> (0, L),

so u_{a,k} must lie in (0, L_k) and S_{a,k} = theta_k + delta_{L_k}^{-1}(u_{a,k})
must be nonincreasing in k.  The two outer cells cost nothing: taking
S_{a,0} = +inf (a fibre of size 1, slope 0) and S_{a,m} = -inf (multiplicity 1)
makes phi(-inf) = phi(theta_1) and phi(+inf) = phi(theta_m).

    python research/realizability/i_pattern.py
"""
from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from optimizers import differential_evolution, pattern_search  # noqa: E402

T0 = np.array([[0, 0, 0, 0],
               [1, 0, 0, 1],
               [2, 1, 0, 1],
               [1, 1, 0, 0]], dtype=float)


def sp(t):
    return np.logaddexp(0.0, t)


def delta_L(v, L):
    return sp(L - v) - sp(-v)


def delta_inv(u, L):
    """v with delta_L(v) = u, u in (0, L).  Bisection; delta_L is decreasing."""
    lo, hi = -1.0, 1.0
    while delta_L(lo, L) < u:
        lo -= 2.0 + abs(lo)
        if lo < -1e4:
            return -1e4
    while delta_L(hi, L) > u:
        hi += 2.0 + abs(hi)
        if hi > 1e4:
            return 1e4
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if delta_L(mid, L) > u:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-14 * (1 + abs(mid)):
            break
    return 0.5 * (lo + hi)


def feasibility(z, T):
    """z = (log L1, log L2, log L3, rho1, rho2, rho3, log s).

    Returns (-s, violation).  violation = 0 iff the pattern is realised.
    """
    L = np.exp(np.clip(z[0:3], -12.0, 6.0))
    rho = z[3:6]
    s = math.exp(np.clip(z[6], -12.0, 3.0))
    theta = np.array([0.0, L[0], L[0] + L[1], L[0] + L[1] + L[2]])
    viol = 0.0
    S = np.zeros((4, 3))
    for k in range(3):
        u = rho[k] + s * (T[:, k + 1] - T[:, k])
        for a in range(4):
            if u[a] <= 1e-12:
                viol += (1e-12 - u[a]) + 1e-9
            elif u[a] >= L[k] - 1e-12:
                viol += (u[a] - L[k] + 1e-12) + 1e-9
            else:
                S[a, k] = theta[k] + delta_inv(u[a], L[k])
    if viol > 0:
        return 1e3 + viol
    for a in range(4):
        for k in range(2):
            if S[a, k] < S[a, k + 1]:
                viol += S[a, k + 1] - S[a, k]
    if viol > 0:
        return 1e3 + viol
    return -s


def best_for_order(perm, seed=0, restarts=2, maxiter=160):
    T = T0[:, list(perm)]
    b = [(-4.0, 4.0)] * 3 + [(-3.0, 3.0)] * 3 + [(-9.0, 1.0)]
    best_z, best = None, math.inf
    for t in range(restarts):
        z, f = differential_evolution(feasibility, b, args=(T,),
                                      seed=seed + 613 * t, maxiter=maxiter,
                                      popsize=16, F=(0.3, 1.2), CR=0.9)
        for step in (0.4, 0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5,
                     3e-6, 1e-6):
            z, f = pattern_search(feasibility, z, args=(T,), step=step,
                                  min_step=1e-12, maxiter=100000, bounds=b)
        if f < best:
            best_z, best = z, f
    return best_z, best


def main():
    names = "ABCD"          # A = beta+, B = gamma+, C = beta-, D = gamma-
    print("  ordering of (beta+, gamma+, beta-, gamma-) along theta,")
    print("  and the largest scale s for which the forced 4-node pattern is")
    print("  realisable by a cone element  (feasible <=> s > 0):\n")
    print(f"  {'theta-order':>14} {'max s':>14}")
    any_ok = False
    for perm in itertools.permutations(range(4)):
        z, f = best_for_order(perm, seed=11 + 97 * sum(p * i for i, p in enumerate(perm)))
        lbl = "".join(names[p] for p in perm)
        if f < 0:
            any_ok = True
            print(f"  {lbl:>14} {-f:14.8f}   FEASIBLE", flush=True)
        else:
            print(f"  {lbl:>14} {'infeasible':>14}", flush=True)
    print("\n  any feasible ordering:", any_ok)


if __name__ == "__main__":
    main()
