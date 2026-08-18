#!/usr/bin/env python3
"""The wide limiting 3-cycle at ``alpha_max <= 14``, and its four witnesses.

``FINDINGS.md`` searched the symplectic cone only to ``alpha_max <= 12`` (genus
six) and concluded that **every** limiting 3-cycle needs a Jacobian with a
repeated isogeny factor -- "the multiplicity-free sub-cone is transitive, 7308
oriented triangles, 0 cycles".  That is an artefact of the cutoff.  Pushing the
same search to ``alpha_max <= 14`` (``genus4_cycles.py --wide``) produces

    SU(2)^5           <  USp(14)  <  USp(6) x USp(6)  <  SU(2)^5
      genus 5            genus 7        genus 6            margin 4.0116e-2

with **no repeated factor anywhere**, and

    SU(2) x 2.SU(2) x USp(4)  <  USp(14)  <  USp(6) x USp(6)
      genus 5                     genus 7       genus 6      margin 6.1583e-2

which needs only the multiplicity two of ``REPEATED_FACTOR.md``'s Theorem A.
Both margins are far wider than the ``1.206e-2`` of the cycle that session
witnessed, so the implied ``q_0`` is much smaller.

This script exhibits and verifies all four vertices.

**SU(2)^5, genus five** -- the ``(Z/2)^3``-cover of ``P^1`` branched at the six
points ``0, infty, 1, 2, t, t+1`` with inertia vectors

    0 |-> e1,  infty |-> e1,  1 |-> e2,  2 |-> e1+e2,  t |-> e3,  t+1 |-> e1+e3

(their sum is 0 and they span, so the cover is connected).  Of the seven
non-trivial characters, five have a four-point branch locus and two have a
two-point one, so the quotient genera are ``(1,1,1,1,1,0,0)`` and

    Jac(C_t) ~ E_a x E_b x E_c x E_d x E_e ,
    E_a : y^2 = x(x-2)(x-t-1)          E_b : y^2 = x(x-1)(x-t-1)
    E_c : y^2 = x(x-2)(x-t)            E_d : y^2 = (x-1)(x-2)(x-t)(x-t-1)
    E_e : y^2 = x(x-1)(x-t)            (the Legendre pencil)

Their ``lambda``-invariants are ``(t+1)/2, t+1, t/2, (1-t)^2/(t(t-2)), t``, whose
``j``-maps have the five **distinct** pole sets ``{-1,1,inf}``, ``{-1,0,inf}``,
``{0,2,inf}``, ``{0,1,2,inf}``, ``{0,1,inf}``; different places of potentially
multiplicative reduction means no two of the five are isogenous, so by Goursat
the monodromy is the full ``SL2^5`` and the measure is ``SU(2)^5``, ``m2 = 5``.

**SU(2) x 2.SU(2) x USp(4), genus five** -- ``(Z/2)^2``-cover with
``|D1| = 4``, ``|D2| = |D3| = 2``:

    D1 = {i,-i,1/m,-1/m}   D2 = {1,2}   D3 = {m,-m}
    f1 = (x-1)(x-2)(x^2-m^2)                    genus 1
    f2 = (x^2+1)(x^2-m^2)(x^2-m^-2)             genus 2, Jac ~ E^2  (Theorem A)
    f3 = (x^2+1)(x^2-m^-2)(x-1)(x-2)            genus 2

All three polynomials are rational for ``m in F_q``: no congruence on ``q``.

**USp(14), genus seven** -- ``y^2 = x^15 + x + c``.
**USp(6) x USp(6), genus six** -- ``y^2 = (x^2)^7 + x^2 + c``, the even split.

    python research/sato_tate_limit/wide_cycle.py [--deep]
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

from curve_lib import (chi_table, compose_square, family_traces, is_prime,
                       moments, polymul, squarefree)

HERE = Path(__file__).resolve().parent

BIG = [401, 1009, 4001, 16001]
DEEP = [40009]


def inv(a: int, q: int) -> int:
    return pow(int(a) % q, q - 2, q)


def pr(roots, q: int) -> list[int]:
    p = [1]
    for r in roots:
        p = polymul(p, [1, (-r) % q], q)
    return p


# ------------------------------------------------------------- SU(2)^5, g = 5

PARAMS = {
    "a=2, b=t, c=t+1": lambda t, q: (2 % q, t % q, (t + 1) % q),
    "a=t, b=t^2+1, c=t+3": lambda t, q: (t % q, (t * t + 1) % q, (t + 3) % q),
    "a=t, b=t+1, c=t+2  [TRAP]": lambda t, q: (t % q, (t + 1) % q, (t + 2) % q),
}


def pencil_su2_5(q: int, f):
    """the five elliptic quotients of the ``(Z/2)^3``-cover."""
    out = [[], [], [], [], []]
    ts = []
    for t in range(q):
        a, b, c = f(t, q)
        if len({0, 1, a, b, c}) != 5:
            continue
        polys = [pr([0, a, c], q), pr([0, 1, c], q), pr([0, a, b], q),
                 pr([1, a, b, c], q), pr([0, 1, b], q)]
        if not all(squarefree(v, q) for v in polys):
            continue
        ts.append(t)
        for k in range(5):
            out[k].append(polys[k])
    return ts, [np.array(z, dtype=np.int64) for z in out]


BRANCH_V = {"0": 1, "oo": 1, "1": 2, "a": 3, "b": 4, "c": 5}


def par(u: int, v: int) -> int:
    return bin(u & v).count("1") % 2


def branch_polys(q: int, a: int, b: int, c: int):
    """for each of the eight characters, the monic polynomial whose roots are
    the FINITE branch points of that quotient, and whether infinity is one."""
    pts = [(0, 1), (1, 2), (a, 3), (b, 4), (c, 5)]        # (point, inertia)
    out = {}
    for m in range(8):
        roots = [x for x, v in pts if par(m, v)]
        out[m] = (pr(roots, q), par(m, 1))                 # infinity has v = 1
    return out


def direct_su2_5(q: int, a: int, b: int, c: int, chi) -> int:
    """``#C(F_q)`` of the smooth projective (Z/2)^3-cover, counted from the
    local structure of the normalisation and nothing else:

    * away from the branch locus the fibre is the product of the three double
      covers of a basis, ``prod_i (1 + chi(f_i(x)))``;
    * at a branch point with inertia ``<v>`` there are ``8/2 = 4`` geometric
      points, and they are rational exactly when ``f_psi(x)`` is a square for
      every ``psi`` in the annihilator ``v^perp`` (a group of order four), the
      count being 4 or 0;
    * at infinity, whose inertia is also of order two, the same rule applies
      and all the relevant polynomials are monic of even degree, giving 4.
    """
    F = branch_polys(q, a, b, c)
    basis = [F[1][0], F[2][0], F[4][0]]
    x = np.arange(q, dtype=np.int64)

    def ev(f):
        acc = np.zeros(q, dtype=np.int64)
        for cc in f:
            acc = (acc * x + int(cc)) % q
        return acc
    vals = {m: ev(F[m][0]) for m in range(1, 8)}
    prod = (1 + chi[vals[1]]) * (1 + chi[vals[2]]) * (1 + chi[vals[4]])
    branch = {0: 1, 1: 2, a: 3, b: 4, c: 5}
    total = 0
    for xx in range(q):
        if xx in branch:
            v = branch[xx]
            perp = [m for m in range(1, 8) if par(m, v) == 0]
            total += 4 if all(chi[int(vals[m][xx])] == 1 for m in perp) else 0
        else:
            total += int(prod[xx])
    return total + 4


def su2_5_direct_check(q: int, f, chi, nsample: int = 8) -> int:
    """compare ``a(C) = sum_chi a(J_chi)`` against that direct count."""
    bad = 0
    step = max(1, q // nsample)
    for t in range(0, q, step):
        a, b, c = f(t, q)
        if len({0, 1, a, b, c}) != 5:
            continue
        quots = [pr([0, a, c], q), pr([0, 1, c], q), pr([0, a, b], q),
                 pr([1, a, b, c], q), pr([0, 1, b], q)]
        if not all(squarefree(v, q) for v in quots):
            continue
        degs = [3, 3, 3, 4, 3]
        s = sum(int(family_traces(q, np.array([p], dtype=np.int64), d)[0])
                for p, d in zip(quots, degs))
        if q + 1 - direct_su2_5(q, a, b, c, chi) != s:
            bad += 1
    return bad


# --------------------------------------------- SU(2) x 2.SU(2) x USp(4), g = 5


def pencil_su2_2su2_usp4(q: int):
    F1, F2, F3, GE = [], [], [], []
    ms = []
    for m in range(1, q):
        if pow(m, 4, q) == 1:
            continue
        m2 = m * m % q
        mi2 = inv(m2, q)
        u1 = polymul([1, 0, 1], [1, 0, (-mi2) % q], q)
        u2 = polymul([1, (-1) % q], [1, (-2) % q], q)
        u3 = [1, 0, (-m2) % q]
        f1 = polymul(u2, u3, q)
        f2 = polymul(u3, u1, q)
        f3 = polymul(u1, u2, q)
        if not (squarefree(f1, q) and squarefree(f2, q) and squarefree(f3, q)):
            continue
        ms.append(m)
        F1.append(f1)
        F2.append(f2)
        F3.append(f3)
        GE.append(pr([(-1) % q, m2, mi2], q))
    return ms, [np.array(z, dtype=np.int64) for z in (F1, F2, F3, GE)]


# ------------------------------------------------------- USp(14), USp(6)^2


def pencil_usp14(q: int):
    F = []
    for c in range(q):
        p = [1] + [0] * 13 + [1, c % q]
        if squarefree(p, q):
            F.append(p)
    return np.array(F, dtype=np.int64)


def pencil_usp6_usp6(q: int):
    G, H, O = [], [], []
    for c in range(q):
        g = [1] + [0] * 5 + [1, c % q]          # u^7 + u + c
        if not squarefree(g, q):
            continue
        G.append(compose_square(g))
        H.append(g)
        O.append(g + [0])
    return tuple(np.array(z, dtype=np.int64) for z in (G, H, O))


def main() -> int:
    deep = "--deep" in sys.argv
    qs = BIG + (DEEP if deep else [])
    rows: list[list] = []

    print("=" * 100)
    print("1.  SU(2)^5 at genus five:  the (Z/2)^3-cover branched at"
          "  {0, oo, 1, a, b, c}")
    print("    quotient genera (1,1,1,1,1,0,0);  a(C) = sum of five elliptic"
          " traces")
    print("=" * 100)
    print(f"  {'parametrisation':<28}{'q':>7}{'fibres':>8}{'m2':>9}{'m4':>10}"
          f"{'m6':>11}{'linked pairs':>14}{'KR fail':>9}")
    for name, f in PARAMS.items():
        for q in qs:
            if not is_prime(q):
                continue
            ts, FS = pencil_su2_5(q, f)
            AS = [family_traces(q, F, 4 if k == 3 else 3)
                  for k, F in enumerate(FS)]
            aC = sum(AS)
            m = moments(aC, q, 6)
            linked = [(i, j) for i in range(5) for j in range(i + 1, 5)
                      if np.all(AS[i] == AS[j]) or np.all(AS[i] == -AS[j])]
            kr = su2_5_direct_check(q, f, chi_table(q), 12)
            print(f"  {name:<28}{q:>7}{len(ts):>8}{m[2]:>9.4f}{m[4]:>10.2f}"
                  f"{m[6]:>11.1f}{str(linked):>14}{kr:>9}")
            rows.append(["su2_5", name, q, len(ts), f"{m[2]:.6f}",
                         f"{m[4]:.4f}", f"{m[6]:.2f}", str(linked), kr])
        print()
    print("  target SU(2)^5:  m2 = 5, m4 = 70, m6 = 1525")
    print("  (the third row is a trap: with a = t, b = t+1, c = t+2 the")
    print("   lambda-invariants t and (t+1)/t are in the same S3-orbit, so two")
    print("   of the five quotients are isomorphic and the measure is")
    print("   SU(2)^3 x 2.SU(2) with m2 = 7, not SU(2)^5 with m2 = 5.)")

    print()
    print("=" * 100)
    print("2.  SU(2) x 2.SU(2) x USp(4) at genus five")
    print("=" * 100)
    print(f"  {'q':>7}{'fibres':>8}{'m2':>9}{'m4':>10}{'m6':>11}"
          f"{'a(J2)!=2a(E)':>14}{'m2(E1)':>9}{'m2(J2)':>9}{'m2(J3)':>9}")
    for q in qs:
        if not is_prime(q):
            continue
        ms, (F1, F2, F3, GE) = pencil_su2_2su2_usp4(q)
        a1 = family_traces(q, F1, 4)
        a2 = family_traces(q, F2, 6)
        a3 = family_traces(q, F3, 6)
        aE = family_traces(q, GE, 3)
        aC = a1 + a2 + a3
        m = moments(aC, q, 6)
        bad = int(np.sum(a2 != 2 * aE))
        print(f"  {q:>7}{len(ms):>8}{m[2]:>9.4f}{m[4]:>10.2f}{m[6]:>11.1f}"
              f"{bad:>14}{moments(a1, q, 2)[2]:>9.4f}"
              f"{moments(a2, q, 2)[2]:>9.4f}{moments(a3, q, 2)[2]:>9.4f}")
        rows.append(["su2_2su2_usp4", q, len(ms), f"{m[2]:.6f}",
                     f"{m[4]:.4f}", f"{m[6]:.2f}", bad])
    print("  target:  m2 = 6, m4 = 91, m6 = 2034;  blocks 1, 4, 1")

    print()
    print("=" * 100)
    print("3.  USp(14) at genus seven and USp(6) x USp(6) at genus six")
    print("=" * 100)
    print(f"  {'q':>7}{'USp14 m2':>10}{'m4':>9}{'USp6^2 m2':>11}{'m4':>9}"
          f"{'half m2':>9}{'half m2':>9}{'split fail':>11}")
    for q in qs:
        if not is_prime(q):
            continue
        a = family_traces(q, pencil_usp14(q), 15)
        m = moments(a, q, 4)
        G, H, O = pencil_usp6_usp6(q)
        aG = family_traces(q, G, 14)
        aH = family_traces(q, H, 7)
        aO = family_traces(q, O, 8)
        mg = moments(aG, q, 4)
        print(f"  {q:>7}{m[2]:>10.4f}{m[4]:>9.3f}{mg[2]:>11.4f}{mg[4]:>9.3f}"
              f"{moments(aH, q, 2)[2]:>9.4f}{moments(aO, q, 2)[2]:>9.4f}"
              f"{int(np.sum(aG != aH + aO)):>11}")
        rows.append(["usp14_usp6sq", q, f"{m[2]:.6f}", f"{m[4]:.4f}",
                     f"{mg[2]:.6f}", f"{mg[4]:.4f}",
                     int(np.sum(aG != aH + aO))])
    print("  targets:  USp(14) 1, 3    USp(6)xUSp(6) 2, 12    halves 1, 1")

    with (HERE / "wide_cycle.csv").open("w", newline="",
                                        encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind"] + [f"c{i}" for i in range(8)])
        for r in rows:
            wr.writerow(list(r) + [""] * (9 - len(r)))
    print("\nwritten: wide_cycle.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
