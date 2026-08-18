#!/usr/bin/env python3
"""The cycle search restricted to what an explicit pencil now witnesses.

``symplectic_search.py`` found nine 3-cycles among the 71 symplectic measures
with ``alpha_max <= 12`` and reported that **all nine** need a block ``k.G``
with ``k >= 2`` -- a Jacobian isogenous to ``A^k``.  ``repeated_factor.py``
exhibits such families.  This script asks what that buys:

* how many of the nine survive when multiplicities are capped at ``k <= 2``,
  which is what the even+palindromic construction supplies;
* and it re-verifies the surviving cycle with a **second, independent**
  implementation of the Sato-Tate cumulant generating function -- the Bessel
  Toeplitz-minus-Hankel determinant

      E[e^{tau tr}]_{USp(2N)} = det( I_{i-j}(2 tau) - I_{i+j}(2 tau) )_{i,j=1..N}

  carried at 40+ working digits in ``mpmath`` all the way through the midrange,
  with the extrema located by golden section on ``log tau`` in ``mpf``
  arithmetic, so no float64 ever touches the headline numbers.

    python research/sato_tate_limit/witness_search.py
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import numpy as np
from mpmath import besseli, mp, mpf

import st_lib as S

HERE = Path(__file__).resolve().parent
TAU = S.tau_grid(1e-4, 1e5, 1201)
CAP = 12.0


# ----------------------------------------------- independent mpmath K_USp(2N)


def K_usp(rank: int, tau) -> "mpf":
    """``log E[e^{tau tr}]`` over ``USp(2 rank)``, Bessel determinant form."""
    t = mpf(tau)
    if t == 0:
        return mpf(0)
    saved = mp.dps
    mp.dps = 50 + int(rank * (rank - 1) * math.log10(max(2.0 * float(t), 10.0)))
    try:
        z = 2 * t
        cache: dict[int, mpf] = {}

        def I(n: int) -> mpf:
            n = abs(n)
            if n not in cache:
                cache[n] = besseli(n, z)
            return cache[n]

        m = mp.matrix(rank, rank)
        for i in range(1, rank + 1):
            for j in range(1, rank + 1):
                m[i - 1, j - 1] = I(i - j) - I(i + j)
        return mp.log(mp.det(m))
    finally:
        mp.dps = saved


def K_measure(blocks, tau) -> "mpf":
    """``blocks`` is a tuple ``(rank, multiplicity)``; the trace is
    ``sum_i k_i tr(g_i)`` with the ``g_i`` independent."""
    return sum((K_usp(r, mpf(k) * mpf(tau)) for r, k in blocks), mpf(0))


def alpha_max(blocks) -> int:
    return sum(2 * r * k for r, k in blocks)


def psi(blocks, tau) -> "mpf":
    return K_measure(blocks, tau) / mpf(tau)


def midrange_mp(A, B, lo=mpf("1e-6"), hi=mpf("1e6"), n=241):
    """``mid = (sup + inf)/2`` of ``Psi_A - Psi_B`` over ``[0, infty]``, with a
    coarse ``mpf`` scan followed by golden-section polish of both extrema."""
    saved = mp.dps
    mp.dps = 60
    try:
        lg0, lg1 = mp.log(lo), mp.log(hi)
        xs = [lg0 + (lg1 - lg0) * mpf(i) / (n - 1) for i in range(n)]

        def f(s):
            t = mp.e ** s
            return psi(A, t) - psi(B, t)

        vals = [f(s) for s in xs]
        imax = max(range(n), key=lambda i: vals[i])
        imin = min(range(n), key=lambda i: vals[i])

        def polish(i, sign):
            a = xs[max(i - 1, 0)]
            b = xs[min(i + 1, n - 1)]
            phi = (mp.sqrt(5) - 1) / 2
            x1, x2 = b - phi * (b - a), a + phi * (b - a)
            f1, f2 = sign * f(x1), sign * f(x2)
            for _ in range(120):
                if f1 > f2:
                    b, x2, f2 = x2, x1, f1
                    x1 = b - phi * (b - a)
                    f1 = sign * f(x1)
                else:
                    a, x1, f1 = x1, x2, f2
                    x2 = a + phi * (b - a)
                    f2 = sign * f(x2)
                if b - a < mpf("1e-30"):
                    break
            return sign * max(f1, f2), mp.e ** (x1 if f1 > f2 else x2)

        hv, ht = polish(imax, mpf(1))
        lv, lt = polish(imin, mpf(-1))
        end = mpf(alpha_max(A) - alpha_max(B))
        sup, sup_at = (hv, ht) if hv >= max(mpf(0), end) else \
            ((end, mp.inf) if end >= 0 else (mpf(0), mpf(0)))
        inf, inf_at = (lv, lt) if lv <= min(mpf(0), end) else \
            ((end, mp.inf) if end <= 0 else (mpf(0), mpf(0)))
        return (sup + inf) / 2, sup, inf, sup_at, inf_at
    finally:
        mp.dps = saved


# --------------------------------------------------- the float64 cone search


def factors(cap: float, kmax: int) -> list[S.Factor]:
    out = []
    for g in (1, 2, 3, 4, 5, 6):
        for k in range(1, min(kmax, int(cap // (2 * g))) + 1):
            out.append(S.Factor("SU2" if g == 1 else f"USp{2 * g}", k, 0.0))
    return [f for f in out if f.alpha_max <= cap + 1e-9]


def products(cap: float, facs: list[S.Factor]) -> list[S.Measure]:
    facs = sorted(facs, key=lambda f: (f.alpha_max, f.group, f.multiplicity))
    out: list[S.Measure] = []

    def rec(start: int, used: float, chosen: list[S.Factor]) -> None:
        if chosen:
            out.append(S.Measure(tuple(chosen)))
        for i in range(start, len(facs)):
            f = facs[i]
            if used + f.alpha_max <= cap + 1e-9:
                rec(i, used + f.alpha_max, chosen + [f])

    rec(0, 0.0, [])
    seen: dict[str, S.Measure] = {}
    for m in out:
        seen.setdefault(m.label, m)
    return list(seen.values())


def mid_matrix(psis: np.ndarray, amax: np.ndarray) -> np.ndarray:
    n = psis.shape[0]
    mid = np.zeros((n, n))
    for i in range(n):
        d = psis[i][None, :] - psis
        end = amax[i] - amax
        hi = np.maximum(np.maximum(d.max(axis=1), 0.0), end)
        lo = np.minimum(np.minimum(d.min(axis=1), 0.0), end)
        mid[i] = 0.5 * (hi + lo)
    np.fill_diagonal(mid, 0.0)
    return mid


def cycles_of(mid: np.ndarray, tol: float = 1e-6):
    n = mid.shape[0]
    strict = mid < -tol
    out = []
    total = 0
    for i, j, k in itertools.combinations(range(n), 3):
        for a, b, c in ((i, j, k), (i, k, j)):
            total += 1
            if strict[a, b] and strict[b, c] and strict[c, a]:
                out.append((a, b, c, min(-mid[a, b], -mid[b, c], -mid[c, a])))
    return total, out


def main() -> int:
    rows: list[list] = []
    print("=" * 92)
    print("the cone search with the multiplicity cap the construction supplies")
    print("=" * 92)
    for kmax in (1, 2, 3, 4, 6):
        lib = products(CAP, factors(CAP, kmax))
        lib.sort(key=lambda m: (m.alpha_max, m.variance, m.tail, m.label))
        amax = np.array([m.alpha_max for m in lib])
        psis = np.array([m.Psi(TAU) for m in lib])
        mid = mid_matrix(psis, amax)
        total, cyc = cycles_of(mid)
        seen = set()
        uniq = []
        for a, b, c, marg in cyc:
            key = frozenset((a, b, c))
            if key not in seen:
                seen.add(key)
                uniq.append((a, b, c, marg))
        best = max((m for *_, m in uniq), default=float("nan"))
        print(f"\n  multiplicity <= {kmax}:  {len(lib):>3} measures, "
              f"{total:>6} oriented triangles, {len(uniq):>2} distinct "
              f"3-cycles"
              + (f", widest margin {best:.4e}" if uniq else ""))
        for a, b, c, marg in sorted(uniq, key=lambda r: -r[3]):
            print(f"      {lib[a].label:<26} < {lib[b].label:<12} < "
                  f"{lib[c].label:<20} margin {marg:.4e}   genus "
                  f"{amax[a] / 2:g}/{amax[b] / 2:g}/{amax[c] / 2:g}")
            rows.append(["cap", kmax, lib[a].label, lib[b].label,
                         lib[c].label, f"{marg:.6e}"])
        rows.append(["summary", kmax, len(lib), total, len(uniq),
                     f"{best:.6e}" if uniq else ""])

    # ------------------------------------------------- the witnessed cycle
    print()
    print("=" * 92)
    print("the surviving cycle, recomputed at 40+ digits with an independent")
    print("Bessel determinant  E[e^{tau tr}] = det(I_{i-j}(2tau) - I_{i+j}(2tau))")
    print("=" * 92)

    # cross-check the independent implementation against st_lib first
    print("\n  cross-check against st_lib (float64), tau = 0.7, 3, 40:")
    for name, rank in (("SU2", 1), ("USp4", 2), ("USp6", 3), ("USp12", 6)):
        for t in (0.7, 3.0, 40.0):
            a = float(S.group_K(name, np.array([t]))[0])
            b = float(K_usp(rank, t))
            print(f"    {name:<6} tau={t:<5g}  st_lib {a:>18.12f}   "
                  f"Bessel {b:>18.12f}   diff {abs(a - b):.3e}")
            rows.append(["crosscheck", name, t, f"{a:.12f}", f"{b:.12f}",
                         f"{abs(a - b):.3e}"])

    TRI = [
        ("2.SU2 x 2.SU2", ((1, 2), (1, 2))),
        ("USp12",         ((6, 1),)),
        ("SU2 x SU2 x USp6", ((1, 1), (1, 1), (3, 1))),
    ]
    print(f"\n  {'edge':<44}{'midrange (40 digits)':>44}")
    margins = []
    for i in range(3):
        (na, A), (nb, B) = TRI[i], TRI[(i + 1) % 3]
        mid, sup, inf, sat, iat = midrange_mp(A, B)
        mp.dps = 42
        print(f"  {na + '  ->  ' + nb:<44}{mp.nstr(mid, 40):>44}")
        print(f"      sup {mp.nstr(sup, 25)} at tau = {mp.nstr(sat, 12)}")
        print(f"      inf {mp.nstr(inf, 25)} at tau = {mp.nstr(iat, 12)}")
        margins.append(-mid)
        rows.append(["cycle40", na, nb, mp.nstr(mid, 40), mp.nstr(sup, 25),
                     mp.nstr(inf, 25), mp.nstr(sat, 12), mp.nstr(iat, 12)])
    strict = all(m > 0 for m in margins)
    print(f"\n  all three midranges negative -- strict 3-cycle: {strict}")
    print(f"  smallest margin  {mp.nstr(min(margins), 30)}")
    rows.append(["cycle40-margin", mp.nstr(min(margins), 30), strict,
                 "", "", "", "", ""])

    with (HERE / "witness_search.csv").open("w", newline="",
                                            encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind", "a", "b", "c", "d", "e", "f", "g"])
        for r in rows:
            wr.writerow(list(r) + [""] * (8 - len(r)))
    print("\nwritten: witness_search.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
