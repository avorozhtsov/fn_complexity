#!/usr/bin/env python3
"""The gauge decomposition of the exchange matrix, and what it says about cycles.

Write L(a,b) = -log C(a->b), which by Theorem 1 of
`paper_finite_fields_maps/docs/exchange_positivity_and_weil.md` equals
sup_beta g(beta) for g = u_b - u_a, u_a = log log Z_a.  Split it into its
symmetric and antisymmetric parts,

    L = S + A,      S(a,b) = (L_ab + L_ba)/2,   A(a,b) = (L_ab - L_ba)/2.

Then, because sup(-g) = -inf g,

    S = (max g - min g)/2 = d/2        the half-RANGE of g   (the metric)
    A = (max g + min g)/2              the MIDRANGE of g     (a 1-form)

so the two halves of the exchange matrix are the range and the midrange of one
and the same function.  Four consequences, each checked below.

1.  a < b  (i.e. C(a->b) < C(b->a))  <=>  A(a,b) > 0.
    The comparison is nothing but the sign pattern of the antisymmetric part.

2.  A triangle is a strict 3-cycle  <=>  |curl A| = sum |A| over its edges and
    min |A| > 0, where curl A = A(a,b) + A(b,c) + A(c,a).  The ratio
    |curl A| / sum|A| lies in [0,1] and equals 1 exactly on cycles, which makes
    it a smooth search objective where sign-hunting gives none.

3.  If both infima of a pair are attained at an endpoint then A is EXACT:

        A(a,b) = psi(b) - psi(a),    psi = (1/2) log[ log(#fibers) * log(max) ]

    i.e. psi = (1/2) log phi with phi the potential of the endpoint-regime
    theorem.  An exact 1-form has zero curl, so no cycle can live inside the
    endpoint regime -- which is that theorem, re-derived as "zero curvature".
    A cycle is therefore a nonzero holonomy of A, and |A - dpsi| measures how
    far a pair is from the regime that forbids it.

4.  No arbitrage (cycle products of M at most 1, in both orientations) says
    exactly that sum S >= |sum A| around every loop: the metric part dominates
    the flux part.

Finally the Gibbs family e^{-t d} is examined at t = 1/2, where it equals the
Szegedy discriminant sqrt(M o M^T) of a quantum walk.  The minimal five-point
negative-type certificate is already not positive semidefinite there.

    python research/m_and_e_and_a_c/gauge_decomposition.py
"""
from __future__ import annotations

import itertools
import math
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from fn_complexity import exchange_rate, exchange_rate_result  # noqa: E402
import t1_2_common as T  # noqa: E402

CERT5 = [
    (12, 10, 8, 8, 2, 1),
    (11, 9, 7, 7, 4, 1),
    (12, 12, 6, 5, 4, 4),
    (12, 10, 7, 4, 3, 3),
    (11, 11, 7, 7, 4, 3),
]
CYCLE = [(6, 3, 3), (7, 2, 1), (6, 5, 1)]


def L(a, b):
    return -math.log(exchange_rate(a, b))


def parts(a, b):
    lab, lba = L(a, b), L(b, a)
    return 0.5 * (lab + lba), 0.5 * (lab - lba)


def psi(a):
    return 0.5 * math.log(math.log(len(a)) * math.log(max(a)))


def endpoint_pair(a, b):
    ends = (0.0, math.inf)
    return (exchange_rate_result(a, b).beta in ends
            and exchange_rate_result(b, a).beta in ends)


def curl_ratio(triangle):
    a, b, c = triangle
    edges = [parts(a, b)[1], parts(b, c)[1], parts(c, a)[1]]
    return abs(sum(edges)) / sum(abs(e) for e in edges), min(abs(e) for e in edges)


def check_range_midrange(families, beta_max=600.0, points=400_001):
    betas = np.concatenate([np.linspace(0.0, beta_max, points)[1:], [1e12]])
    worst = 0.0
    for family in families:
        for a, b in itertools.combinations(family, 2):
            g = T.u_values(b, betas) - T.u_values(a, betas)
            s_grid, a_grid = 0.5 * (g.max() - g.min()), 0.5 * (g.max() + g.min())
            s, ant = parts(a, b)
            worst = max(worst, abs(s - s_grid), abs(ant - a_grid))
    return worst


def check_exactness(pool):
    endpoint, interior = [], []
    for a, b in itertools.combinations(sorted(set(pool)), 2):
        err = abs(parts(a, b)[1] - (psi(b) - psi(a)))
        (endpoint if endpoint_pair(a, b) else interior).append(err)
    return endpoint, interior


def random_pool(n=60, seed=11):
    rng = random.Random(seed)
    pool = set()
    while len(pool) < n:
        r = rng.randint(3, 7)
        s = tuple(sorted((rng.randint(1, 14) for _ in range(r)), reverse=True))
        if len(set(s)) > 1:
            pool.add(s)
    return sorted(pool)


def main():
    print("1. S is the half-range and A the midrange of g = u_b - u_a")
    worst = check_range_midrange([CYCLE, CERT5])
    print(f"   max |solver - beta-grid| over S and A: {worst:.3e}"
          "   (grid resolution, not a discrepancy)\n")

    print("2. |curl A| / sum|A| = 1 exactly on a strict 3-cycle")
    triangles = {
        "known 3-cycle": CYCLE,
        "control, endpoint regime": [(10, 6), (8, 8, 1, 1, 1, 1), (9, 9, 2)],
        "control, mixed": [(6, 3, 3), (7, 2, 1), (9, 9, 2)],
        "cert5 (1,2,3)": CERT5[:3],
        "cert5 (1,3,5)": [CERT5[0], CERT5[2], CERT5[4]],
    }
    for name, tri in triangles.items():
        ratio, margin = curl_ratio(tri)
        print(f"   {name:<24} ratio = {ratio:.10f}   min|A| = {margin:.3e}")
    print()

    print("3. in the endpoint regime A is exact, with potential psi = (1/2) log phi")
    endpoint, interior = check_exactness(random_pool())
    print(f"   both infima at an endpoint : {len(endpoint):4d} pairs, "
          f"max |A - dpsi| = {max(endpoint):.3e}")
    print(f"   at least one interior      : {len(interior):4d} pairs, "
          f"max = {max(interior):.3e}, median = {sorted(interior)[len(interior) // 2]:.3e}\n")

    print("4. no arbitrage: sum S >= |sum A| around every loop")
    tot_s = tot_a = 0.0
    for a, b in zip(CYCLE, CYCLE[1:] + CYCLE[:1]):
        s, ant = parts(a, b)
        tot_s += s
        tot_a += ant
    print(f"   the 3-cycle: sum S = {tot_s:.6f} >= |sum A| = {abs(tot_a):.6f}  "
          f"-> {tot_s >= abs(tot_a)}\n")

    print("5. the Gibbs family e^{-t d} at the Szegedy point t = 1/2")
    families = T.build_families()
    families["cert5"] = CERT5
    for name in ("cert5", "cert13", "greedy25", "greedy30", "greedy40"):
        if name not in families:
            continue
        d = T.distance_matrix(families[name])
        lam = float(np.linalg.eigvalsh(np.exp(-0.5 * d))[0])
        star = T.psd_threshold(d)
        star = "-" if star is None else f"{star:.6f}"
        print(f"   {name:>9} n={len(families[name]):>3}  "
              f"lambda_min(sqrt(M o M^T)) = {lam: .6e}   t* = {star}")
    print("\n   NOTE: the minimal certificate has t* = 1.0918, not the 0.124 quoted in\n"
          "   FINDINGS T1.3 -- that value belongs to the superseded 13-point family.")


if __name__ == "__main__":
    main()
