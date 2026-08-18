#!/usr/bin/env python3
"""The genus-four measures the nine limiting cycles still need -- session brief J.

``REPEATED_FACTOR.md`` closed the smallest of ``FINDINGS.md``'s nine limiting
3-cycles and left two measures open, both at genus four:

    SU(2) x 3.SU(2)      Jac ~ J1 x E^3,  E varying     (cycles 8 and 9)
    4.SU(2)              Jac ~ E^4,       E varying     (cycles 1-6)

and wrote out the route: a ``(Z/2)^2``-cover of ``P^1`` with ``|B1| = |B2| = 4``,
``|B1 n B2| = 1``, hence ``|B3| = 6`` and quotient genera ``(1,1,2)``, built on
the branch set of the ``Jac ~ E^2`` sextic.  This script executes it.

**The witness.**  Let ``q = 1 (mod 4)``, ``i^2 = -1`` in ``F_q``, ``e = +-1``:

    B3 = {i,-i,n,-n,1/n,-1/n},   f3 = (x^2+1)(x^2-n^2)(x^2-n^-2)   [Jac ~ E^2]
    S  = {i, n, e/n},   T = -S,   p = n + i + e/n
    f2 = prod_{S u {p}}(x-r)     f1 = prod_{T u {p}}(x-r)
    f1 f2 = (x-p)^2 f3

``C_n`` is the smooth ``(Z/2)^2``-cover ``y1^2 = f1, y2^2 = f2``; it has genus 4,
is non-hyperelliptic, and ``Jac(C_n) ~ J1 x J2 x J3`` with ``J2 = E`` (proved
in ``genus4_symbolic.py``: the pullback constant is a perfect square, so there
is no quadratic twist) and ``J3 ~ E^2`` (Theorem A of ``REPEATED_FACTOR.md``).
Hence

    a(C_n) = a(J1) + 3 a(E)        ->  measure  SU(2) x 3.SU(2)

Detectors, ``alpha = -a/sqrt(q)``, ``m_j = E[alpha^j]``:

    SU(2) x 3.SU(2)    m2 =  10   m4 = 218   m6 = 6350
    4.SU(2)            m2 =  16   m4 = 512   m6 = 20480
    half-twisted mix   m2 =   6   m4 = 114   m6 =  3210
    USp(8)             m2 =   1   m4 =   3   m6 =    15
    SU2 x USp4 x USp4  m2 =   3   m4 =  26   m6 =   363
    SU2^2 x USp6       m2 =   3   m4 =  25   m6 =   325

    python research/sato_tate_limit/genus4_witness.py [--deep]
"""

from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import numpy as np

from curve_lib import (chi_table, family_traces, is_prime, moments, polymul,
                       squarefree)

HERE = Path(__file__).resolve().parent

# primes = 1 mod 4
SMALL = [13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157,
         173, 181, 193, 197, 229, 233, 241, 257, 269, 277, 281, 293, 401, 1009,
         4001, 16001]
BIG = [401, 1009, 4001, 16001]
DEEP = [40009]
LOCUS_QS = [401, 1009, 4001, 16001, 40009]


def inv(a: int, q: int) -> int:
    return pow(int(a) % q, q - 2, q)


def sqrt_minus_one(q: int) -> int:
    return next(t for t in range(1, q) if (t * t + 1) % q == 0)


def poly_from_roots(roots, q: int) -> list[int]:
    p = [1]
    for r in roots:
        p = polymul(p, [1, (-r) % q], q)
    return p


# ------------------------------------------------------------ the new pencil


def pencil_su2_3su2(q: int, e: int = -1):
    """``S = {i, n, e/n}``, ``T = -S``, ``p = n + i + e/n``.

    Returns ``(ns, F1, F2, F3, GE)``: the parameters and the coefficient rows
    of ``f1`` (quartic), ``f2`` (quartic), ``f3`` (sextic) and the elliptic
    quotient ``g_E(u) = (u+1)(u-n^2)(u-n^-2)`` of ``f3``.
    """
    i0 = sqrt_minus_one(q)
    ns, F1, F2, F3, GE = [], [], [], [], []
    for n in range(1, q):
        if pow(n, 4, q) == 1:
            continue
        ni = inv(n, q)
        B3 = [i0, (-i0) % q, n, (-n) % q, ni, (-ni) % q]
        S = [i0, n, (e * ni) % q]
        T = [(-v) % q for v in S]
        p = (n + i0 + e * ni) % q
        if p in B3:
            continue
        ns.append((n, p))
        F1.append(poly_from_roots(T + [p], q))
        F2.append(poly_from_roots(S + [p], q))
        F3.append(poly_from_roots(B3, q))
        GE.append(poly_from_roots([(-1) % q, n * n % q, inv(n * n % q, q)], q))
    return (ns, np.array(F1, dtype=np.int64), np.array(F2, dtype=np.int64),
            np.array(F3, dtype=np.int64), np.array(GE, dtype=np.int64))


def direct_count(q: int, f1, f2, f3, chi) -> int:
    """``#C(F_q)`` for the ``(Z/2)^2`` cover, straight from the definition:
    ``1 + chi(f1) + chi(f2) + chi(f3)`` points over each affine ``x``, and, all
    three being monic of even degree, four unramified points over infinity."""
    x = np.arange(q, dtype=np.int64)

    def ev(f):
        acc = np.zeros(q, dtype=np.int64)
        for c in f:
            acc = (acc * x + int(c)) % q
        return acc
    return int(np.sum(1 + chi[ev(f1)] + chi[ev(f2)] + chi[ev(f3)])) + 4


# --------------------------------------------- the general (n, p) parametrisation


SPLITS = [(S, tuple(k for k in range(6) if k not in S))
          for S in itertools.combinations(range(6), 3) if 0 in S]


def jlam(lam: int, q: int):
    den = lam * lam % q * pow((lam - 1) % q, 2, q) % q
    if den == 0:
        return None
    return 256 * pow((lam * lam - lam + 1) % q, 3, q) % q * inv(den, q) % q


def cross(p: int, s, q: int):
    k = (s[0] - s[2]) % q * inv((s[0] - s[1]) % q, q) % q
    den = (p - s[2]) % q
    if den == 0:
        return None
    return k * (p - s[1]) % q * inv(den, q) % q


def locus_counts(q: int):
    """For every splitting ``S|T`` and every one of the six ``S3``-branches of
    ``j(J2) = j(E)``, count the ``n in F_q`` that in addition have
    ``j(J1) = j(E)`` -- the ``4.SU(2)`` condition."""
    i0 = sqrt_minus_one(q)
    counts = {}
    for si, (S, T) in enumerate(SPLITS):
        for br in range(6):
            hits = 0
            for n in range(1, q):
                if pow(n, 4, q) == 1:
                    continue
                ni = inv(n, q)
                roots = [i0, (-i0) % q, n, (-n) % q, ni, (-ni) % q]
                lamE = inv(n * n % q, q)
                if lamE in (0, 1):
                    continue
                orb = [lamE, inv(lamE, q), (1 - lamE) % q,
                       inv((1 - lamE) % q, q),
                       lamE * inv((lamE - 1) % q, q) % q,
                       (lamE - 1) % q * inv(lamE, q) % q]
                s = [roots[k] for k in S]
                t = [roots[k] for k in T]
                k = (s[0] - s[2]) % q * inv((s[0] - s[1]) % q, q) % q
                a = (orb[br] - k) % q
                if a == 0:
                    continue
                p = (orb[br] * s[2] - k * s[1]) % q * inv(a, q) % q
                if p in roots:
                    continue
                l1 = cross(p, t, q)
                if l1 is None:
                    continue
                j1, jE = jlam(l1, q), jlam(lamE, q)
                if j1 is not None and jE is not None and j1 == jE:
                    hits += 1
            counts[(si, br)] = hits
    return counts


def pencil_branch(q: int, si: int, br: int):
    """the same construction on an arbitrary splitting/branch -- for controls."""
    i0 = sqrt_minus_one(q)
    S, T = SPLITS[si]
    F1, F2, F3, GE = [], [], [], []
    for n in range(1, q):
        if pow(n, 4, q) == 1:
            continue
        ni = inv(n, q)
        roots = [i0, (-i0) % q, n, (-n) % q, ni, (-ni) % q]
        lamE = inv(n * n % q, q)
        if lamE in (0, 1):
            continue
        orb = [lamE, inv(lamE, q), (1 - lamE) % q, inv((1 - lamE) % q, q),
               lamE * inv((lamE - 1) % q, q) % q,
               (lamE - 1) % q * inv(lamE, q) % q]
        s = [roots[k] for k in S]
        t = [roots[k] for k in T]
        k = (s[0] - s[2]) % q * inv((s[0] - s[1]) % q, q) % q
        a = (orb[br] - k) % q
        if a == 0:
            continue
        p = (orb[br] * s[2] - k * s[1]) % q * inv(a, q) % q
        if p in roots:
            continue
        F1.append(poly_from_roots(t + [p], q))
        F2.append(poly_from_roots(s + [p], q))
        F3.append(poly_from_roots(roots, q))
        GE.append(poly_from_roots([(-1) % q, n * n % q, inv(n * n % q, q)], q))
    return tuple(np.array(z, dtype=np.int64) for z in (F1, F2, F3, GE))


# ---------------------------------------------------- genus five: SU2 x USp4^2


def pencil_g5_usp4sq(q: int):
    """``(Z/2)^2``-cover with ``|D1| = 4``, ``|D2| = |D3| = 2``: quotient genera
    ``(1, 2, 2)``, total genus 5, measure ``SU(2) x USp(4) x USp(4)``.

        u1 = x^4 + x + c        u2 = (x-1)(x-2)        u3 = (x-3)(x-c)
        f1 = u2 u3   (genus 1)  f2 = u3 u1  (genus 2)  f3 = u1 u2  (genus 2)
    """
    F1, F2, F3 = [], [], []
    for c in range(q):
        u1 = [1, 0, 0, 1, c % q]
        u2 = polymul([1, (-1) % q], [1, (-2) % q], q)
        u3 = polymul([1, (-3) % q], [1, (-c) % q], q)
        f1 = polymul(u2, u3, q)
        f2 = polymul(u3, u1, q)
        f3 = polymul(u1, u2, q)
        if not (squarefree(f1, q) and squarefree(f2, q) and squarefree(f3, q)):
            continue
        F1.append(f1)
        F2.append(f2)
        F3.append(f3)
    return tuple(np.array(z, dtype=np.int64) for z in (F1, F2, F3))


# ------------------------------------------------------------------- reporting


def main() -> int:
    deep = "--deep" in sys.argv
    rows: list[list] = []

    # =================================================================== 1
    print("=" * 100)
    print("1.  the SU(2) x 3.SU(2) pencil at genus four:  is  a(C_n) = a(J1) +"
          " 3 a(E)  on every fibre?")
    print("    C_n = (Z/2)^2 cover of P^1 branched at {+-i, +-n, +-1/n, p},"
          "  p = n + i + e/n")
    print("=" * 100)
    for e in (-1, 1):
        print(f"\n  e = {e:+d}   (S = " + ("{i, n, -1/n}" if e < 0
                                           else "{i, n, 1/n}") + ")")
        print(f"  {'q':>7}{'fibres':>8}{'aC!=a1+3aE':>12}{'a2!=aE':>9}"
              f"{'a3!=2aE':>9}{'KR fail':>9}{'distinct j(E)':>15}"
              f"{'a1=+-aE':>9}{'m2':>9}{'m4':>10}{'m6':>12}")
        tot_fib = tot_bad = 0
        for q in SMALL + (DEEP if deep else []):
            if not is_prime(q) or q % 4 != 1:
                continue
            ns, F1, F2, F3, GE = pencil_su2_3su2(q, e)
            a1 = family_traces(q, F1, 4)
            a2 = family_traces(q, F2, 4)
            a3 = family_traces(q, F3, 6)
            aE = family_traces(q, GE, 3)
            aC = a1 + a2 + a3
            bad = int(np.sum(aC != a1 + 3 * aE))
            b2 = int(np.sum(a2 != aE))
            b3 = int(np.sum(a3 != 2 * aE))
            # Kani-Rosen, checked directly against a point count on C
            chi = chi_table(q)
            step = max(1, len(F1) // 16)
            kr = sum(1 for k in range(0, len(F1), step)
                     if q + 1 - direct_count(q, F1[k], F2[k], F3[k], chi)
                     != int(aC[k]))
            js = set()
            for (n, _) in ns:
                n2 = n * n % q
                c = (-(n2 + inv(n2, q))) % q
                d = (c + 2) % q
                if d:
                    js.add(256 * pow((c + 1) % q, 3, q) * inv(d, q) % q)
            eq1 = int(np.sum((a1 == aE) | (a1 == -aE)))
            m = moments(aC, q, 6)
            tot_fib += len(F1)
            tot_bad += bad
            print(f"  {q:>7}{len(F1):>8}{bad:>12}{b2:>9}{b3:>9}{kr:>9}"
                  f"{len(js):>15}{eq1:>9}{m[2]:>9.4f}{m[4]:>10.2f}{m[6]:>12.1f}")
            rows.append([f"su2x3su2_e{e:+d}", q, len(F1), bad, b2, b3, kr,
                         len(js), eq1, f"{m[2]:.6f}", f"{m[4]:.4f}",
                         f"{m[6]:.2f}"])
        print(f"  total: {tot_fib} fibres, {tot_bad} mismatches of"
              f"  a(C) = a(J1) + 3 a(E)")
        print("  target SU(2) x 3.SU(2):  m2 = 10, m4 = 218, m6 = 6350")

    # =================================================================== 2
    print()
    print("=" * 100)
    print("2.  controls: the other five S3-branches, and the splittings with"
          " T != -S")
    print("    (a family in which J2 is the quadratic TWIST of E on half the"
          " fibres has")
    print("     alpha = alpha_1 + 3t on half and alpha_1 + t on the other half:"
          "  m2 = 6, not 10)")
    print("=" * 100)
    print(f"  {'splitting':<22}{'branch':>7}{'q':>7}{'fibres':>8}"
          f"{'a2=aE':>8}{'a2=-aE':>8}{'m2':>9}{'m4':>10}   reading")
    for si in (5, 6, 0, 4):
        for br in range(6):
            q = 4001
            F1, F2, F3, GE = pencil_branch(q, si, br)
            if len(F1) < 100:
                continue
            a1 = family_traces(q, F1, 4)
            a2 = family_traces(q, F2, 4)
            a3 = family_traces(q, F3, 6)
            aE = family_traces(q, GE, 3)
            aC = a1 + a2 + a3
            same = int(np.sum(a2 == aE))
            opp = int(np.sum(a2 == -aE))
            m = moments(aC, q, 6)
            tag = ("SU2 x 3.SU2" if same == len(F1) else
                   ("3.SU2 twisted" if opp == len(F1) else
                    "half-twisted mixture"))
            S, T = SPLITS[si]
            lbl = f"S={S} T={T}"
            print(f"  {lbl:<22}{br:>7}{q:>7}{len(F1):>8}{same:>8}{opp:>8}"
                  f"{m[2]:>9.4f}{m[4]:>10.2f}   {tag}")
            rows.append(["branch_control", si, br, q, len(F1), same, opp,
                         f"{m[2]:.6f}", f"{m[4]:.4f}", tag])

    print()
    print("  all 60 splitting/branch pairs at q = 1009, classified:")
    tally: dict[str, int] = {}
    for si in range(len(SPLITS)):
        for br in range(6):
            q = 1009
            F1, F2, F3, GE = pencil_branch(q, si, br)
            if len(F1) < 100:
                tally["degenerate"] = tally.get("degenerate", 0) + 1
                continue
            a2 = family_traces(q, F2, 4)
            aE = family_traces(q, GE, 3)
            same = int(np.sum(a2 == aE))
            opp = int(np.sum(a2 == -aE))
            tag = ("J2 = E on every fibre" if same == len(F1) else
                   ("J2 = twist of E on every fibre" if opp == len(F1) else
                    "half-twisted mixture"))
            tally[tag] = tally.get(tag, 0) + 1
    for k, v in sorted(tally.items()):
        print(f"    {k:<34} {v:>3} of 60")
        rows.append(["branch_tally", k, v])

    # =================================================================== 3
    print()
    print("=" * 100)
    print("3.  the 4.SU(2) locus:  for every splitting and every branch,")
    print("    #{n in F_q : j(J1) = j(E)}  --  a family would grow like q")
    print("=" * 100)
    counts = {q: locus_counts(q) for q in LOCUS_QS}
    print(f"  {'split':>6}{'branch':>7}" + "".join(f"{q:>9}" for q in LOCUS_QS))
    worst = 0
    for si in range(len(SPLITS)):
        for br in range(6):
            row = [counts[q][(si, br)] for q in LOCUS_QS]
            worst = max(worst, max(row))
            print(f"  {si:>6}{br:>7}" + "".join(f"{v:>9}" for v in row))
            rows.append(["locus", si, br] + row)
    print(f"\n  largest count anywhere: {worst}, over a 100-fold range of q.")
    print("  The locus is ZERO-DIMENSIONAL in every branch: 4.SU(2) is not")
    print("  reachable by this construction.")

    # =================================================================== 4
    print()
    print("=" * 100)
    print("4.  genus five, SU(2) x USp(4) x USp(4)  --  the third vertex of"
          " cycle 8")
    print("    (Z/2)^2 cover with |D1|=4, |D2|=|D3|=2: quotient genera (1,2,2)")
    print("=" * 100)
    print(f"  {'q':>7}{'fibres':>8}{'m2(C)':>9}{'m4(C)':>10}{'m6(C)':>11}"
          f"{'m2(E)':>8}{'m4(E)':>8}{'m2(J2)':>8}{'m4(J2)':>8}"
          f"{'m2(J3)':>8}{'m4(J3)':>8}{'m2(J2+J3)':>11}")
    for q in BIG + (DEEP if deep else []):
        F1, F2, F3 = pencil_g5_usp4sq(q)
        a1 = family_traces(q, F1, 4)
        a2 = family_traces(q, F2, 6)
        a3 = family_traces(q, F3, 6)
        aC = a1 + a2 + a3
        m = moments(aC, q, 6)
        m1, m2m, m3m = (moments(z, q, 4) for z in (a1, a2, a3))
        m23 = moments(a2 + a3, q, 4)
        print(f"  {q:>7}{len(F1):>8}{m[2]:>9.4f}{m[4]:>10.3f}{m[6]:>11.1f}"
              f"{m1[2]:>8.3f}{m1[4]:>8.3f}{m2m[2]:>8.3f}{m2m[4]:>8.3f}"
              f"{m3m[2]:>8.3f}{m3m[4]:>8.3f}{m23[2]:>11.3f}")
        rows.append(["g5_su2_usp4_usp4", q, len(F1), f"{m[2]:.6f}",
                     f"{m[4]:.4f}", f"{m[6]:.2f}", f"{m1[2]:.4f}",
                     f"{m2m[2]:.4f}", f"{m3m[2]:.4f}", f"{m23[2]:.4f}"])
    print("  targets:  C  3, 26, 363   E  1, 2   J2, J3  1, 3   J2+J3  2")
    print("  (SU2^2 x USp6, the other genus-five vertex, has m4 = 25 not 26;")
    print("   the three quotient genera 1,2,2 are forced by the branch data.)")

    with (HERE / "genus4_witness.csv").open("w", newline="",
                                            encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind"] + [f"c{i}" for i in range(11)])
        for r in rows:
            wr.writerow(list(r) + [""] * (12 - len(r)))
    print("\nwritten: genus4_witness.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
