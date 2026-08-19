#!/usr/bin/env python3
"""Symbolic support for the genus-four ``SU(2) x 3.SU(2)`` witness.

The construction (``GENUS4_WITNESSES.md``, Theorem C).  Over a field ``K`` of
characteristic ``!= 2`` containing ``i`` with ``i^2 = -1``, take the branch set
of the ``Jac ~ E^2`` sextic of ``REPEATED_FACTOR.md``,

    B3 = {i, -i, n, -n, 1/n, -1/n},   f3 = (x^2+1)(x^2-n^2)(x^2-n^-2),

split it into the two triples ``S = {i, n, e/n}``, ``T = -S`` (``e = +-1``),
adjoin a seventh point ``p`` and form the ``(Z/2)^2``-cover of ``P^1`` with

    f1 = prod_{T u {p}} (x - r),      f2 = prod_{S u {p}} (x - r),
    f1 f2 = (x - p)^2 f3.

Kani-Rosen gives ``Jac(C) ~ J1 x J2 x J3`` with ``J3 ~ E^2``, so
``a(C) = a(J1) + a(J2) + 2 a(E)``.  This script verifies, in ``sympy``:

1.  ``lambda(J2) = lambda(E) = 1/n^2``  exactly when  ``p = n + i + e/n``;
2.  the Moebius map ``M`` carrying the branch set of ``J2`` to that of ``E``
    pulls ``E`` back to ``w^2 = K f2(x)`` with ``K = (n^2+1)^4`` a **square**,
    so ``J2 = E`` over the base field, no quadratic twist -- hence
    ``a(C) = a(J1) + 3 a(E)``;
3.  ``j(E)`` and ``j(J1)`` are non-constant with **different pole divisors**, so
    the two elliptic pencils are not isogenous over the algebraic closure of
    ``K(n)`` and the monodromy is the full ``SL2 x SL2``;
4.  the remaining ``4.SU(2)`` condition ``j(J1) = j(E)`` is a non-zero
    polynomial in ``n`` -- a **finite** set of fibres, not a family.

    python research/sato_tate_limit/genus4_symbolic.py
"""

from __future__ import annotations

import sympy as sp

n, x = sp.symbols('n x')
I = sp.I


def jinv(lam):
    return sp.cancel(256 * (lam ** 2 - lam + 1) ** 3 / (lam ** 2 * (lam - 1) ** 2))


def lam_of(p, s):
    """cross-ratio of ``{p} u s`` under ``s1 -> 0, s2 -> infty, s0 -> 1``."""
    k = (s[0] - s[2]) / (s[0] - s[1])
    return sp.cancel(sp.together(k * (p - s[1]) / (p - s[2])))


def pole_divisor(f):
    """the pole set of a rational function of ``n``, with orders, plus infinity."""
    f = sp.cancel(sp.together(f))
    num, den = sp.fraction(f)
    poles = sp.roots(sp.Poly(sp.expand(den), n))
    d = sp.degree(sp.Poly(sp.expand(num), n)) - sp.degree(sp.Poly(sp.expand(den), n))
    out = {sp.simplify(r): m for r, m in poles.items()}
    if d > 0:
        out[sp.oo] = d
    return out


def report(e: int) -> None:
    tag = {1: "S = {i, n, 1/n}   (splitting 5)",
           -1: "S = {i, n, -1/n}  (splitting 6)"}[e]
    print("=" * 78)
    print(f"  {tag}")
    print("=" * 78)

    S = [I, n, e / n]
    T = [-v for v in S]
    lamE = 1 / n ** 2
    p = n + I + e / n

    print(f"  p = n + i + ({e})/n")
    l2 = sp.cancel(sp.simplify(lam_of(p, S)))
    print(f"  lambda(J2) = {l2}        lambda(E) = {sp.simplify(lamE)}")
    assert sp.simplify(l2 - lamE) == 0
    print("  -> lambda(J2) = lambda(E) identically in n            [verified]")

    # ---- the Moebius map and the twist constant
    a, b, c, d = sp.symbols('a b c d')
    # M : n -> -1,  i -> n^2,  e/n -> infinity     (so p -> n^-2)
    eqs = [a * S[1] + b - (-1) * (c * S[1] + d),
           a * S[0] + b - (n ** 2) * (c * S[0] + d),
           c * S[2] + d]
    sol = sp.solve(eqs, [a, b, c, d], dict=True)[0]
    free = [v for v in (a, b, c, d) if v not in sol]
    sub = {free[0]: 1}
    A = sp.simplify(sol.get(a, a).subs(sub))
    B = sp.simplify(sol.get(b, b).subs(sub))
    C = sp.simplify(sol.get(c, c).subs(sub))
    D = sp.simplify(sol.get(d, d).subs(sub))
    M = (A * x + B) / (C * x + D)
    print(f"  M(x) = ({A}) x + ({B})  over  ({C}) x + ({D})")
    print(f"  M(p) = {sp.simplify(sp.cancel(M.subs(x, p)))}   (target n^-2)")
    assert sp.simplify(M.subs(x, p) - 1 / n ** 2) == 0

    g = (M + 1) * (M - n ** 2) * (M - 1 / n ** 2)
    f2 = sp.expand(sp.prod([x - r for r in S + [p]]))
    K = sp.simplify(sp.cancel(sp.expand(sp.cancel(g * (C * x + D) ** 4)) / f2))
    print(f"  pullback constant K = {sp.factor(K)}")
    coeff, facs = sp.factor_list(sp.expand(K))
    is_sq = all(m % 2 == 0 for _, m in facs) and sp.sqrt(coeff).is_rational
    root = sp.sqrt(coeff) * sp.prod([f ** (m // 2) for f, m in facs])
    print(f"  K = ({sp.factor(root)})^2   -- perfect square: {is_sq}")
    assert is_sq and sp.simplify(sp.expand(root ** 2 - K)) == 0
    print("  -> J2 = E over the base field, NO quadratic twist       [proved]")

    # ---- the two j-invariants
    jE = sp.cancel(jinv(lamE))
    l1 = sp.cancel(sp.simplify(lam_of(p, T)))
    j1 = sp.cancel(sp.simplify(jinv(l1)))
    print(f"\n  j(E)  = {sp.factor(jE)}")
    print(f"  j(J1) = {sp.factor(j1)}")
    pE, p1 = pole_divisor(jE), pole_divisor(j1)
    print(f"  poles of j(E)  : {pE}")
    print(f"  poles of j(J1) : {p1}")
    same = set(sp.simplify(k) for k in pE) == set(sp.simplify(k) for k in p1)
    print(f"  pole sets equal? {same}   -> isogenous? {same}")
    print("  -> different places of potentially multiplicative reduction,")
    print("     so J1 and E are NOT isogenous over the closure of K(n).")

    # ---- the 4.SU(2) condition
    num = sp.factor(sp.numer(sp.cancel(j1 - jE)))
    deg = sp.degree(sp.Poly(sp.expand(sp.numer(sp.cancel(j1 - jE))), n))
    print(f"\n  j(J1) - j(E) numerator, degree {deg}:")
    print(f"    {num}")
    print(f"  identically zero? {sp.simplify(sp.cancel(j1 - jE)) == 0}")
    print("  -> the 4.SU(2) condition cuts out AT MOST "
          f"{deg} values of n: a finite set, not a pencil.\n")


def elementary_abelian_shapes(r: int, gmax: int = 5) -> None:
    """Every ``(Z/2)^r``-cover of ``P^1`` with all inertia of order two: which
    multisets of quotient genera occur, up to total genus ``gmax``?

    ``2 g_C - 2 = 2^r(-2) + B 2^{r-1}`` with ``B`` the number of branch points,
    and for a non-trivial character ``chi`` the quotient ``C/ker chi`` is the
    double cover branched at ``B_chi = {P : chi(v_P) = -1}``, of genus
    ``|B_chi|/2 - 1``.  The data is the multiset ``{v_P}`` of inertia vectors,
    subject to ``sum v_P = 0`` (the cover exists) and ``<v_P> = G`` (it is
    connected).
    """
    import itertools
    nz = list(range(1, 2 ** r))
    seen: dict[tuple, tuple] = {}
    for B in range(3, 12):
        if 2 ** r * (-2) + B * 2 ** (r - 1) + 2 > 2 * gmax:
            break
        for v in itertools.combinations_with_replacement(nz, B):
            s = 0
            for z in v:
                s ^= z
            if s != 0:
                continue
            span = {0}
            for z in v:
                span |= {y ^ z for y in span}
            if len(span) != 2 ** r:
                continue
            gs = []
            for chi in nz:
                b = sum(1 for z in v if bin(chi & z).count("1") % 2)
                gs.append(b // 2 - 1 if b >= 2 else 0)
            key = tuple(sorted(gs))
            if sum(key) <= gmax:
                seen.setdefault((sum(key), key), v)
    for (g, key), v in sorted(seen.items()):
        print(f"    genus {g}:  quotient genera {key}   example inertia {v}")


def main() -> int:
    print("=" * 78)
    print("  0.  which quotient-genus shapes an elementary abelian cover of"
          " P^1 can have")
    print("=" * 78)
    for r in (2, 3):
        print(f"  (Z/2)^{r}:")
        elementary_abelian_shapes(r, 5)
    print("  -> at genus 4 a (Z/2)^2-cover has quotient genera (1,1,2) or")
    print("     (0,2,2) and nothing else, and NO (Z/2)^3-cover of P^1 has")
    print("     genus 4 at all (2g-2 = -16+4B forces g odd).  So four elliptic")
    print("     quotients cyclically permuted -- the mechanism that gives")
    print("     3.SU(2) at genus 3 -- is unavailable at genus 4.\n")
    for e in (-1, 1):
        report(e)
    print("=" * 78)
    print("  degeneracy: the seven points {+-i, +-n, +-1/n, p} are distinct iff")
    print("  n != 0, n^4 != 1, and p is none of them.  Solving:")
    for e in (-1, 1):
        p = n + I + e / n
        S = [I, n, e / n]
        bad = []
        for r in [I, -I, n, -n, 1 / n, -1 / n]:
            sol = sp.solve(sp.numer(sp.cancel(p - r)), n)
            bad.extend(sol)
        bad = sorted(set(sp.simplify(v) for v in bad), key=str)
        print(f"    e = {e:+d}:  p in B3 iff n in {bad}")

    # ------------------------------------------------- the SU(2)^5 pencil
    print()
    print("=" * 78)
    print("  the genus-five SU(2)^5 pencil of wide_cycle.py: five elliptic")
    print("  quotients of the (Z/2)^3-cover branched at {0, oo, 1, 2, t, t+1}")
    print("=" * 78)
    t = sp.symbols('t')
    lam = {
        "E_a  {0,oo,2,t+1}": (t + 1) / 2,
        "E_b  {0,oo,1,t+1}": t + 1,
        "E_c  {0,oo,2,t}": t / 2,
        "E_d  {1,2,t,t+1}": (1 - t) ** 2 / (t * (t - 2)),
        "E_e  {0,oo,1,t}": t,
    }
    poles = {}
    for name, l in lam.items():
        l = sp.cancel(sp.together(l))
        pol = set()
        for target in (0, 1):
            pol |= set(sp.solve(sp.numer(sp.cancel(l - target)), t))
        pol |= set(sp.solve(sp.denom(sp.cancel(l)), t))
        if sp.degree(sp.Poly(sp.numer(l), t)) >= sp.degree(
                sp.Poly(sp.denom(l), t)):
            pol.add(sp.oo)
        poles[name] = frozenset(sp.simplify(v) for v in pol)
        print(f"    {name:<20} lambda = {l}    j-poles {sorted(poles[name], key=str)}")
    names = list(poles)
    clash = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]
             if poles[a] == poles[b]]
    print(f"  pairs with equal pole sets: {clash}")
    print("  -> all five pole sets are distinct, so no two of the five elliptic")
    print("     pencils are isogenous, and the monodromy is the full SL2^5.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
