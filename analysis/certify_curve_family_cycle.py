#!/usr/bin/env python3
"""A strict exchange 3-cycle among genus-2 curve families over F_101.

This answers the central question of `research/session_briefs/B_cycles_among_curve_families.md`
affirmatively.  The three resources are fiber signatures of hyperelliptic
pencils ``f_i : A^2 -> A^1``, ``f_i(x,y) = y^2 - P_i(x)``, whose fibers over
``c in F_101`` are the genus-2 curves ``y^2 = P_i(x) + c``:

    P_1 = x^5 + 70x^4 + 28x^3 + 15x^2 + 11x + 31
    P_2 = x^5 + 42x^4 + 32x^3 + 74x^2 + 96x + 60
    P_3 = x^5 + 72x^4 + 21x^3 +  2x^2 +  6x + 57

Writing ``a < b`` for ``C(a->b) < C(b->a)``, the comparison runs

    f_1 > f_2 > f_3 > f_1.

Nothing is hard-coded except the polynomials: the signatures are recomputed here
by counting points, so the script is reproducible without the search RNG.

WHY IT IS NOT FORBIDDEN BY THE ENDPOINT THEOREM.  If both rates of a pair are
attained at an endpoint then ``a < b <=> phi(a) < phi(b)`` with
``phi = log(#fibers) log(max fiber)``, and a total preorder has no cycles.  Here
four of the six rates are attained at INTERIOR beta (21.0, 93.2, 38.1, 24.8) and
two at ``beta = inf``, so the hypothesis fails on exactly the edges it must.
Addendum 1 to brief B pushes the scalar one order further, to
``phi~ = M - ((3-2sqrt2)/2) m2``, and predicts that a cycle needs exact
``(M, m2)`` degeneracy because ``m2 = nu(P)/q - 1`` is an integer condition.
That is exactly what happened, without being designed for:

    f_1 and f_2 have max fiber 123 and m2 = 0.851485 BOTH, i.e. nu(P) = 187 for
    both, hence identical phi and identical phi~;
    f_3 has max 122 and m2 = 0.990099, so phi~(f_3) < phi~(f_1) = phi~(f_2).

So one edge (f_3 > f_1) contradicts phi~, one edge (f_2 > f_3) agrees with it,
and the third is a phi~ TIE broken by the interior alone -- the "one or two
phi-violating edges, never zero, never three" pattern of the parent brief.

CERTIFICATION.  Each ``A(a,b) = (log C(b->a) - log C(a->b))/2`` is computed three
ways: the package solver, the independent grid solver of
`research/m_and_e_and_a_c/t2_2_common.py`, and a 40-digit mpmath scan.  The
smallest margin is 1.16e-4, six orders above the 1e-10 tie floor of the briefs.

CAVEAT, and it belongs in any write-up.  This is a cycle among SIGNATURES.  The
signature merges families the geometry separates, so a cycle among signatures is
not yet a cycle among the underlying pencils in any stronger sense; settling that
is brief C's question, not this script's.

    python analysis/certify_curve_family_cycle.py
"""
from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "research" / "m_and_e_and_a_c"))

from fn_complexity import exchange_rate_result  # noqa: E402
import t2_2_common as T2  # noqa: E402

Q = 101
POLYS = {
    "f1": [31, 11, 15, 28, 70, 1],
    "f2": [60, 96, 74, 32, 42, 1],
    "f3": [57, 6, 2, 21, 72, 1],
}
CYCLE = ("f1", "f2", "f3")
mp.mp.dps = 40


def signature(coeffs: list[int], q: int = Q) -> tuple[int, ...]:
    """Fiber sizes of f(x,y) = y^2 - P(x) over F_q, by direct point count."""
    x = np.arange(q)
    px = np.zeros(q, dtype=np.int64)
    for c in reversed(coeffs):
        px = (px * x + c) % q
    squares = np.zeros(q, dtype=np.int64)
    np.add.at(squares, (np.arange(q) ** 2) % q, 1)   # #{y : y^2 = v}
    counts = np.zeros(q, dtype=np.int64)
    for value in px:                                  # fiber over c gets P(x)+c
        counts += np.roll(squares, int(value))
    return tuple(int(v) for v in np.sort(counts)[::-1])


def a_solver(a, b) -> float:
    return 0.5 * (-math.log(exchange_rate_result(a, b).rate)
                  + math.log(exchange_rate_result(b, a).rate))


def a_grid(a, b) -> float:
    sa, sb = T2.sig_from_counts(np.array(a)), T2.sig_from_counts(np.array(b))
    return 0.5 * (-math.log(T2.rate(sa, sb)[0]) + math.log(T2.rate(sb, sa)[0]))


def a_mpmath(a, b, points: int = 4000, beta_max: float = 2000.0) -> float:
    """40-digit scan; log Z is smooth so a dense grid brackets the infimum."""
    def log_z(sig, beta):
        return mp.log(mp.fsum(mp.power(v, beta) for v in sig))
    betas = [mp.mpf(beta_max) * mp.mpf(i) / points for i in range(1, points + 1)]
    betas += [mp.mpf(j) / 200 for j in range(1, 200)]   # resolve small beta too
    ratios = [log_z(a, t) / log_z(b, t) for t in betas]
    lo = min(ratios + [mp.log(len(a)) / mp.log(len(b)),
                       mp.log(max(a)) / mp.log(max(b))])
    hi = max(ratios + [mp.log(len(a)) / mp.log(len(b)),
                       mp.log(max(a)) / mp.log(max(b))])
    # C(a->b) = lo and C(b->a) = 1/hi, so A = (-log lo - log hi)/2.
    return float(-(mp.log(lo) + mp.log(hi)) / 2)


def main() -> None:
    sigs = {name: signature(p) for name, p in POLYS.items()}
    print(f"three genus-2 pencils y^2 = P(x) + c over F_{Q}\n")
    for name, s in sigs.items():
        a_c = [Q - n for n in s]
        m2 = sum(v * v for v in a_c) / Q ** 2
        phi = math.log(len(s)) * math.log(max(s))
        phi_t = max(-v for v in a_c) - ((3 - 2 * math.sqrt(2)) / 2) * m2
        assert sum(s) == Q * Q, "signature must total q^2"
        print(f"  {name}: fibers={len(s)}  max={max(s)}  min={min(s)}  "
              f"sum={sum(s)}=q^2   m2={m2:.6f}   phi={phi:.9f}   phi~={phi_t:.6f}")

    print("\n  phi ties f1 and f2 (both max 123): the endpoint invariant cannot")
    print("  separate them, which is why the interior gets to decide.\n")

    print(f"  {'edge':<10} {'A (solver)':>14} {'A (grid)':>14} {'A (mpmath)':>14} "
          f"{'argmin betas':>26}  verdict")
    ok = True
    for u, v in zip(CYCLE, CYCLE[1:] + CYCLE[:1]):
        a, b = sigs[u], sigs[v]
        s1, s2, s3 = a_solver(a, b), a_grid(a, b), a_mpmath(a, b)
        r1, r2 = exchange_rate_result(a, b), exchange_rate_result(b, a)
        betas = f"{r1.beta:.4g} / {r2.beta:.4g}"
        rel = max(abs(s1 - s2), abs(s1 - s3))
        ok &= (s1 < 0) == (s2 < 0) == (s3 < 0)
        print(f"  {u}->{v:<7} {s1:>14.6e} {s2:>14.6e} {s3:>14.6e} {betas:>26}"
              f"  {u} {'<' if s1 > 0 else '>'} {v}   (spread {rel:.1e})")

    edges = [a_solver(sigs[u], sigs[v])
             for u, v in zip(CYCLE, CYCLE[1:] + CYCLE[:1])]
    ratio = abs(sum(edges)) / sum(abs(e) for e in edges)
    print(f"\n  all three edges agree in sign across all three methods: {ok}")
    print(f"  |curl A| / sum|A| = {ratio:.12f}   (1 exactly iff a strict 3-cycle)")
    print(f"  smallest margin  = {min(abs(e) for e in edges):.3e}"
          f"   against a tie floor of 1e-10")
    print(f"  cycle: {' > '.join(CYCLE)} > {CYCLE[0]}")

    print("\n  cycle products, both orientations (no arbitrage requires <= 1):")
    for direction, order in (("forward", CYCLE), ("reverse", CYCLE[::-1])):
        product = 1.0
        for u, v in zip(order, order[1:] + order[:1]):
            product *= exchange_rate_result(sigs[u], sigs[v]).rate
        print(f"    {direction}: {product:.12f}")


if __name__ == "__main__":
    main()
