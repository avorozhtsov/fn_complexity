#!/usr/bin/env python3
"""The one arithmetic input of ``FINDINGS.md``: families with ``Jac ~ A^k``.

Brief F reduced the whole ``q -> infinity`` non-transitivity question to one
existence statement (``FINDINGS.md``, "Open"):

    an explicit one-parameter family of curves over ``F_q`` whose Jacobian is
    isogenous to ``A^k`` for some ``k >= 2``, with ``A`` varying.

The construction that settles it is

    C :  y^2 = f(x),   f **even** and **palindromic** of degree 2m, m **odd**,

for which ``Jac(C) ~ J x J`` with ``J`` the Jacobian of the even half
``v^2 = g(u)``, ``f(x) = g(x^2)``.  Proof in ``REPEATED_FACTOR.md``; the
mechanism is that ``iota : (x,y) -> (1/x, y/x^m)`` conjugates the even
involution ``sigma : (x,y) -> (-x,y)`` to ``sigma . (hyperelliptic)`` exactly
when ``m`` is odd, so the two elliptic-type quotients are isomorphic over the
base field.

    m = 3   f = (x^2+1)(x^4 + c x^2 + 1)          Jac ~ E^2       2.SU(2)
    m = 5   degree 10                             Jac ~ A_2^2     2.USp(4)
    m = 7   degree 14                             Jac ~ A_3^2     2.USp(6)

and, choosing the half ``J`` to be a rationally split abelian surface,

    m = 5, J ~ E_1 x E_2   Jac ~ E_1^2 x E_2^2    2.SU(2) x 2.SU(2)

which is the genus-four vertex of the smallest cycle of ``FINDINGS.md``.

Detectors, with ``alpha_c = -a_c/sqrt(q)`` and ``m_j = E[alpha^j]``:

    USp(2g), g >= 2       m2 = 1  m4 =  3
    SU(2)                 m2 = 1  m4 =  2   m6 =    5
    SU(2) x SU(2)         m2 = 2  m4 = 10   m6 =   70
    Jac ~ E^2             m2 = 4  m4 = 32   m6 =  320
    2.USp(4)              m2 = 4  m4 = 48   m6 =  896
    2.USp(6)              m2 = 4  m4 = 48   m6 =  960
    E_1^2 x E_2^2         m2 = 8  m4 = 160  m6 = 4480
    swap of SU(2)^2       m2 = 1  m4 =  5   m6 =   35

plus the exact fibre-by-fibre identity ``a_c(C) = 2 a_c(J)``, which is the
decisive check: a moment can be faked by a coincidence, an identity on every
fibre of every prime cannot.

    python research/sato_tate_limit/repeated_factor.py
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

from curve_lib import (compose_square, even_palindromic, family_traces,
                       is_prime, moments, polymul, squarefree)

HERE = Path(__file__).resolve().parent

SMALL = [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
         79, 83, 89, 97, 101, 211, 401, 1009, 2003, 4001, 8009, 16001]
BIG = [401, 1009, 4001, 16001]

# The genus-four vertex is the slowest of the three to settle -- its `m2` is
# still 7.52 at q = 4001 -- so it gets an optional deeper sweep.  ``--deep``
# adds these; the cost is O(q^2) per family, a few minutes at q = 64007.
DEEP = [32003, 64007]


# --------------------------------------------------------------- the pencils


def pencil_E2(q: int):
    """``y^2 = (x^2+1)(x^4 + c x^2 + 1)`` and its elliptic quotient
    ``E_c : v^2 = (u+1)(u^2 + c u + 1)``."""
    cs, F, G = [], [], []
    for c in range(q):
        g = [1, (c + 1) % q, (c + 1) % q, 1]
        if not squarefree(g, q):
            continue
        cs.append(c)
        G.append(g)
        F.append(compose_square(g))
    return (np.array(cs), np.array(F, dtype=np.int64),
            np.array(G, dtype=np.int64))


def pencil_palindromic(q: int, m: int, params):
    """``y^2 = f(x)``, ``f`` even+palindromic of degree ``2m``, and its even
    half ``v^2 = g(u)`` of degree ``m``."""
    F, G = [], []
    for c in range(q):
        f = even_palindromic(m, params(c))
        g = f[::2]
        if not squarefree(f, q):
            continue
        F.append(f)
        G.append(g)
    return np.array(F, dtype=np.int64), np.array(G, dtype=np.int64)


def pencil_E1sq_E2sq(q: int):
    """Genus four, ``Jac ~ E_1^2 x E_2^2``  --  measure ``2.SU2 x 2.SU2``.

    ``g(u) = (u+1)(u^2 - R u + 1)(u^2 - S u + 1)`` is the palindromic quintic
    with roots ``-1, r, 1/r, s, 1/s``; ``R = r + 1/r``, ``S = s + 1/s``.  The
    genus-two half ``v^2 = g(u)`` is bielliptic with the involution
    ``theta : u -> (r u + beta)/(u - r)``, ``beta = r - s(1+r)``, which pairs
    ``(infty, r) (-1, s) (1/r, 1/s)``; that forces ``r s^2 - r s + 1 - r = 0``,
    i.e. ``r(5r-4)`` a square, and the split is defined over the base field --
    rather than over a quadratic extension -- exactly when the fixed points of
    ``theta`` are rational, i.e. when ``(r+1)(r-s)`` is a square.  Both conics
    are rational; solving them gives the single parameter ``w``:

        r = -1 / (w^4 - 3w^2 + 1),      s = 2 - w^2.
    """
    F, G, ws = [], [], []
    for w in range(q):
        w2 = (w * w) % q
        den = (w2 * w2 - 3 * w2 + 1) % q
        s = (2 - w2) % q
        if den == 0 or s == 0:
            continue
        r = (-pow(den, q - 2, q)) % q
        R = (r + pow(r, q - 2, q)) % q
        S = (s + pow(s, q - 2, q)) % q
        g = polymul(polymul([1, 1], [1, (-R) % q, 1], q), [1, (-S) % q, 1], q)
        if not squarefree(g, q):
            continue
        ws.append(w)
        G.append(g)
        F.append(compose_square(g))
    return (np.array(ws), np.array(F, dtype=np.int64),
            np.array(G, dtype=np.int64))


def pencil_E1sq_E2sq_swap(q: int):
    """The **trap**: the same bielliptic pairing without the rationality of the
    fixed points of ``theta``.  Pairing ``(infty, r)(1, -r)(s, -s)`` in the
    coordinate where ``rho -> -rho`` -- the split is then defined only over a
    quadratic extension, the geometric monodromy contains the swap, and the
    measure is ``2 x [1/2 (SU2 * SU2) + 1/2 delta_0]``, not
    ``2.SU2 x 2.SU2``: ``m2 = 4`` instead of ``8``."""
    F, G = [], []
    for k in range(q):
        k2 = (k * k) % q
        d1 = ((k2 - 1) * (k2 + 3)) % q
        d2 = pow((k2 - 1) % q, 2, q)
        if d1 == 0 or d2 == 0:
            continue
        R = (2 * (k2 * k2 + 2 * k2 + 5)) % q * pow(d1, q - 2, q) % q
        S = (2 * (k2 * k2 + 6 * k2 + 1)) % q * pow(d2, q - 2, q) % q
        g = polymul(polymul([1, 1], [1, (-R) % q, 1], q), [1, (-S) % q, 1], q)
        if not squarefree(g, q):
            continue
        G.append(g)
        F.append(compose_square(g))
    return np.array(F, dtype=np.int64), np.array(G, dtype=np.int64)


def pencil_g5(q: int):
    """Genus five, ``Jac ~ E_1 x E_2 x A_3``  --  measure ``SU2^2 x USp(6)``.

    ``y^2 = f(x^2)`` with ``f(u) = ((u-1)^2 - 1)((u-1)^2 - 2)((u-1)^2 - c)``.
    The genus-two half ``v^2 = f(u)`` is even in ``u-1``, so it splits
    rationally into two elliptic curves; the genus-three half ``w^2 = u f(u)``
    does **not** inherit the symmetry, because ``u -> 2-u`` does not fix
    ``{0, infty}``, and keeps full ``USp(6)`` monodromy.  Compare
    ``pencil_g5_trap``.
    """
    F, H, O = [], [], []
    for c in range(q):
        z2 = [1, (-2) % q, 1]
        f = [1]
        for t in (1, 2, c):
            f = polymul(f, [z2[0], z2[1], (z2[2] - t) % q], q)
        if not squarefree(f, q):
            continue
        H.append(f)
        O.append(f + [0])
        F.append(compose_square(f))
    return (np.array(F, dtype=np.int64), np.array(H, dtype=np.int64),
            np.array(O, dtype=np.int64))


def pencil_g5_trap(q: int):
    """The other **trap**: ``y^2 = F(x^2)`` with ``F`` a *palindromic* sextic.
    The genus-two half does split rationally, but the branch set of the
    genus-three half is ``{0, infty} u roots(F)``, which ``u -> 1/u`` also
    preserves, so it splits too and the measure is ``SU2^3 x USp(4)``,
    ``m2 = 4``, not ``SU2^2 x USp(6)``, ``m2 = 3``."""
    F = []
    for c in range(q):
        f = [1, 1, 3, c, 3, 1, 1]
        if not squarefree(f, q):
            continue
        F.append(compose_square(f))
    return np.array(F, dtype=np.int64)


def pencil_generic(q: int, deg: int):
    F = []
    for c in range(q):
        f = [1] + [0] * (deg - 2) + [1, c]
        if squarefree(f, q):
            F.append(f)
    return np.array(F, dtype=np.int64)


def pencil_even_split(q: int, n: int):
    F = []
    for c in range(q):
        g = [1] + [0] * (n - 2) + [1, c]
        if squarefree(g, q):
            F.append(compose_square(g))
    return np.array(F, dtype=np.int64)


def pencil_E3(q: int):
    """Genus three, ``Jac ~ E^3``  --  measure ``3.SU(2)``, needs ``q = 1 (3)``.

    The ``(Z/2)^2``-cover of ``P^1`` branched at ``{a, 1, za, z, z^2 a, z^2}``,
    ``z`` a primitive cube root of unity.  Its three intermediate quotients are
    the elliptic curves ``y^2 = f_i`` with

        B_1 = {a, 1, za, z},  B_2 = {za, z, z^2a, z^2},  B_3 = {a, 1, z^2a, z^2}

    and ``u -> z u`` carries ``B_1 -> B_3 -> B_2 -> B_1``, scaling the quartic
    by ``z``.  ``z = g^{(q-1)/3}`` is always a square when ``3 | q-1``, so the
    three quotients are isomorphic **over ``F_q``** and ``a_c(C) = 3 a_c(E)``.
    Returns the three quartic families; ``Jac(C) ~ J_1 x J_2 x J_3`` by
    Kani-Rosen, so ``a_C = a_1 + a_2 + a_3``.
    """
    z = None
    for g in range(2, q):
        t = pow(g, (q - 1) // 3, q)
        if t != 1:
            z = t
            break
    z2 = z * z % q
    R = ([], [], [])
    for a in range(q):
        sets = ((a, 1, z * a % q, z), (z * a % q, z, z2 * a % q, z2),
                (a, 1, z2 * a % q, z2))
        polys = []
        ok = True
        for S in sets:
            p = [1]
            for r in S:
                p = polymul(p, [1, (-r) % q, ], q)
            if not squarefree(p, q):
                ok = False
                break
            polys.append(p)
        if not ok:
            continue
        for i in range(3):
            R[i].append(polys[i])
    return tuple(np.array(x, dtype=np.int64) for x in R)


# ------------------------------------------------------------------ reports


def line(name: str, q: int, a: np.ndarray, target: str, extra: str = "") -> list:
    m = moments(a, q, 6)
    print(f"  {name:<40}{q:>7}{len(a):>8}{m[2]:>9.4f}{m[4]:>10.3f}"
          f"{m[6]:>12.2f}   {target}{extra}")
    return [name, q, len(a), f"{m[2]:.6f}", f"{m[4]:.4f}", f"{m[6]:.3f}",
            target]


HEAD = (f"  {'family':<40}{'q':>7}{'fibres':>8}{'m2':>9}{'m4':>10}"
        f"{'m6':>12}   target")


def main() -> int:
    rows: list[list] = []

    # ------------------------------------------------------------- part 1
    print("=" * 100)
    print("1.  y^2 = (x^2+1)(x^4 + c x^2 + 1):  is  a_c(C) = 2 a_c(E)  on"
          " every fibre?")
    print("=" * 100)
    print(f"  {'q':>7}{'fibres':>9}{'mismatches':>12}{'distinct j':>12}"
          f"{'m2':>9}{'m4':>10}{'m6':>12}")
    total_fib = 0
    total_bad = 0
    for q in SMALL:
        if not is_prime(q):
            continue
        cs, F, G = pencil_E2(q)
        aC = family_traces(q, F, 6)
        aE = family_traces(q, G, 3)
        bad = int(np.sum(aC != 2 * aE))
        total_fib += len(cs)
        total_bad += bad
        js = set()
        for c in cs.tolist():
            den = (c + 2) % q
            if den:
                js.add(256 * pow((c + 1) % q, 3, q) * pow(den, q - 2, q) % q)
        m = moments(aC, q, 6)
        print(f"  {q:>7}{len(cs):>9}{bad:>12}{len(js):>12}"
              f"{m[2]:>9.4f}{m[4]:>10.3f}{m[6]:>12.2f}")
        rows.append(["E2_identity", q, len(cs), bad, len(js),
                     f"{m[2]:.6f}", f"{m[4]:.4f}", f"{m[6]:.3f}"])
    print(f"\n  {total_fib} fibres over {len([q for q in SMALL if is_prime(q)])}"
          f" primes, {total_bad} mismatches.")
    print("  target 2.SU(2):  m2 = 4, m4 = 32, m6 = 320")
    print("  j(E_c) = 256 (c+1)^3 / (c+2) is non-constant, so the family is"
          " NOT isotrivial.")

    # ------------------------------------------------------------- part 2
    print()
    print("=" * 100)
    print("2.  even + palindromic of degree 2m, m odd:  a_c(C) = 2 a_c(J)"
          " and the measure 2.USp(2h)")
    print("=" * 100)
    print(f"  {'family':<40}{'q':>7}{'fibres':>8}{'mismatch':>10}{'m2':>9}"
          f"{'m4':>10}{'m6':>12}   target")
    specs = [
        ("m=3  (genus 2)  -> 2.SU(2)", 3, lambda c: (c,), 6, 3,
         "4, 32, 320"),
        ("m=5  (genus 4)  -> 2.USp(4)", 5, lambda c: (c, 3), 10, 5,
         "4, 48, 896"),
        ("m=7  (genus 6)  -> 2.USp(6)", 7, lambda c: (c, 3, 5), 14, 7,
         "4, 48, 960"),
    ]
    for label, m, par, degf, degg, target in specs:
        for q in BIG:
            F, G = pencil_palindromic(q, m, par)
            aC = family_traces(q, F, degf)
            aJ = family_traces(q, G, degg)
            bad = int(np.sum(aC != 2 * aJ))
            mm = moments(aC, q, 6)
            print(f"  {label:<40}{q:>7}{len(F):>8}{bad:>10}{mm[2]:>9.4f}"
                  f"{mm[4]:>10.3f}{mm[6]:>12.2f}   {target}")
            rows.append([f"palindromic_m{m}", q, len(F), bad,
                         f"{mm[2]:.6f}", f"{mm[4]:.4f}", f"{mm[6]:.3f}",
                         target])

    # ------------------------------------------------------------- part 3
    print()
    print("=" * 100)
    print("3.  the three vertices of the smallest cycle of FINDINGS.md")
    print("=" * 100)
    print(HEAD)
    for q in BIG + (DEEP if "--deep" in sys.argv else []):
        ws, F, G = pencil_E1sq_E2sq(q)
        aC = family_traces(q, F, 10)
        aJ = family_traces(q, G, 5)
        bad = int(np.sum(aC != 2 * aJ))
        rows.append(line("genus 4  2.SU2 x 2.SU2", q, aC,
                         "8, 160, 4480", f"   [a_C=2a_J off on {bad}]"))
        rows.append(line("   its half J ~ E_1 x E_2", q, aJ, "2, 10, 70"))
    print()
    for q in BIG:
        F, H, O = pencil_g5(q)
        rows.append(line("genus 5  SU2^2 x USp(6)", q,
                         family_traces(q, F, 12), "3, 25, 325"))
        rows.append(line("   its genus-2 half", q, family_traces(q, H, 6),
                         "2, 10, 70"))
        rows.append(line("   its genus-3 half", q, family_traces(q, O, 7),
                         "1, 3, 15"))
    print()
    for q in BIG:
        rows.append(line("genus 6  USp(12)  y^2=x^13+x+c", q,
                         family_traces(q, pencil_generic(q, 13), 13),
                         "1, 3, 15"))

    # ------------------------------------------------------- part 3b: k = 3
    print()
    print("=" * 100)
    print("3b.  multiplicity three:  genus 3 with Jac ~ E^3  (needs q = 1 mod 3)")
    print("=" * 100)
    print(f"  {'family':<40}{'q':>7}{'fibres':>8}{'a1=a2=a3':>10}{'m2':>9}"
          f"{'m4':>10}{'m6':>12}   target")
    for q in [211, 601, 1009, 4021, 16033]:
        if not is_prime(q) or q % 3 != 1:
            continue
        A1, A2, A3 = pencil_E3(q)
        a1 = family_traces(q, A1, 4)
        a2 = family_traces(q, A2, 4)
        a3 = family_traces(q, A3, 4)
        same = int(np.sum((a1 == a2) & (a2 == a3)))
        aC = a1 + a2 + a3
        mm = moments(aC, q, 6)
        print(f"  {'(Z/2)^2 cover, z-symmetric  [g3]':<40}{q:>7}{len(a1):>8}"
              f"{same:>10}{mm[2]:>9.4f}{mm[4]:>10.3f}{mm[6]:>12.2f}   "
              f"3.SU(2): 9, 162, 3645")
        rows.append(["E3_genus3", q, len(a1), same, f"{mm[2]:.6f}",
                     f"{mm[4]:.4f}", f"{mm[6]:.3f}", "9, 162, 3645"])

    # ------------------------------------------------------------- part 4
    print()
    print("=" * 100)
    print("4.  the two traps -- families that look right and are not")
    print("=" * 100)
    print(HEAD)
    for q in BIG:
        F, G = pencil_E1sq_E2sq_swap(q)
        rows.append(line("swap split (irrational theta fixed pts)", q,
                         family_traces(q, F, 10), "4, 80, 2240  (NOT 8)"))
        rows.append(line("   its half (swap measure)", q,
                         family_traces(q, G, 5), "1, 5, 35  (NOT 2)"))
    print()
    for q in BIG:
        rows.append(line("y^2=F(x^2), F palindromic sextic", q,
                         family_traces(q, pencil_g5_trap(q), 12),
                         "4, 45  (SU2^3 x USp4, NOT 3)"))

    # ------------------------------------------------------------- part 5
    print()
    print("=" * 100)
    print("5.  controls")
    print("=" * 100)
    print(HEAD)
    for q in BIG:
        rows.append(line("generic genus 2  y^2=x^5+x+c", q,
                         family_traces(q, pencil_generic(q, 5), 5),
                         "1, 3, 14"))
    for q in BIG:
        rows.append(line("generic genus 4  y^2=x^9+x+c", q,
                         family_traces(q, pencil_generic(q, 9), 9),
                         "1, 3, 15"))
    for q in BIG:
        rows.append(line("split  y^2=x^6+x^2+c", q,
                         family_traces(q, pencil_even_split(q, 3), 6),
                         "2, 10, 70"))

    # ------------------------------------------------------------- part 6
    print()
    print("=" * 100)
    print("6.  the pencil as an honest map  Phi : A^2 -> A^1  (a map of finite"
          " sets, which is all")
    print("    the signature framework needs):  Phi = (y^2 - (x^2+1)(x^4+1)) /"
          " (x^2(x^2+1))")
    print("    off the base locus, and the base locus -- x in {0} u {+-i}, "
          "exactly q or 3q points --")
    print("    distributed one or three to a fibre.  Then sum_c N_c = q^2 on"
          " the nose.")
    print("=" * 100)
    from curve_lib import chi_table
    print(f"  {'q':>7}{'base pts/fibre':>16}{'sum N_c - q^2':>16}"
          f"{'sum a_c':>10}{'a_map - a_proj':>18}{'max |d alpha|':>15}")
    for q in [101, 401, 1009, 4001, 16001]:
        ch = chi_table(q)
        x = np.arange(q, dtype=np.int64)
        base = (x == 0)
        if pow(q - 1, (q - 1) // 2, q) == 1:            # -1 is a square
            i0 = next(t for t in range(q) if (t * t + 1) % q == 0)
            base |= (x == i0) | (x == (q - i0) % q)
        good = ~base
        pad = int(base.sum())                           # points per fibre
        x2 = (x * x) % q
        x4 = (x2 * x2) % q
        x6 = (x4 * x2) % q
        b0 = (x6 + x4 + x2 + 1) % q                     # (x^2+1)(x^4+1)
        b1 = (x4 + x2) % q                              # x^2 (x^2+1)
        tot = 0
        N = np.empty(q, dtype=np.int64)
        for c in range(q):
            fv = (b0 + c * b1) % q
            N[c] = int(np.sum(1 + ch[fv[good]])) + pad
        tot = int(N.sum())
        amap = q - N
        cs, F, _ = pencil_E2(q)
        aproj = family_traces(q, F, 6)
        shifts = amap[cs] - aproj
        st = sorted(set(shifts.tolist()))
        print(f"  {q:>7}{pad:>16}{tot - q * q:>16}{int(amap.sum()):>10}"
              f"{str(st):>18}{max(abs(shifts)) / math.sqrt(q):>15.5f}")
        rows.append(["map_presentation", q, pad, tot - q * q,
                     int(amap.sum()), str(st),
                     f"{max(abs(shifts)) / math.sqrt(q):.6f}", ""])
    print("  a_map - a_proj is a single constant at every q, so alpha shifts"
          " uniformly by O(q^-1/2)")
    print("  and the limiting measure is unchanged.")

    with (HERE / "repeated_factor.csv").open("w", newline="",
                                             encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["family", "q", "fibres", "a", "b", "c", "d", "e"])
        for r in rows:
            wr.writerow(list(r) + [""] * (8 - len(r)))
    print("\nwritten: repeated_factor.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
