#!/usr/bin/env python3
"""What the new genus-four witness buys, at 40 digits -- session brief J.

Part 1 recomputes, with the **independent** Bessel Toeplitz-minus-Hankel
determinant of ``witness_search.py``,

    E[e^{tau tr}]_{USp(2N)} = det( I_{i-j}(2 tau) - I_{i+j}(2 tau) )_{i,j=1..N}

the two limiting 3-cycles of ``FINDINGS.md`` whose genus-four vertex is
``SU(2) x 3.SU(2)`` -- cycles 8 and 9 of the nine, and the only two that
``genus4_witness.py`` unlocks:

    SU2 x 3.SU2  <  USp12  <  SU2 x USp4 x USp4  <  SU2 x 3.SU2
    SU2 x 3.SU2  <  USp12  <  SU2^2 x USp6      <  SU2 x 3.SU2

Part 2 is the optional second front of the brief: the cone search with the
isogeny multiplicity capped at ``k``, pushed from ``alpha_max <= 12`` (genus 6)
to ``alpha_max <= 16`` (genus 8), to see whether a *wider* cycle -- and hence a
smaller ``q_0`` -- lives at higher genus.

    python research/sato_tate_limit/genus4_cycles.py [--wide]
"""

from __future__ import annotations

import csv
import itertools
import math
import sys
from pathlib import Path

import numpy as np
from mpmath import besseli, mp, mpf

import st_lib as S
from witness_search import K_usp, midrange_mp

HERE = Path(__file__).resolve().parent

TRIANGLES = {
    "cycle 8": [("SU2 x 3.SU2", ((1, 1), (1, 3))),
                ("USp12", ((6, 1),)),
                ("SU2 x USp4 x USp4", ((1, 1), (2, 1), (2, 1)))],
    "cycle 9": [("SU2 x 3.SU2", ((1, 1), (1, 3))),
                ("USp12", ((6, 1),)),
                ("SU2 x SU2 x USp6", ((1, 1), (1, 1), (3, 1)))],
    "cycle 7 (control)": [("2.SU2 x 2.SU2", ((1, 2), (1, 2))),
                          ("USp12", ((6, 1),)),
                          ("SU2 x SU2 x USp6", ((1, 1), (1, 1), (3, 1)))],
    "wide, multiplicity <= 2": [("SU2 x 2.SU2 x USp4", ((1, 1), (1, 2), (2, 1))),
                                ("USp14", ((7, 1),)),
                                ("USp6 x USp6", ((3, 1), (3, 1)))],
    "wide, MULTIPLICITY-FREE": [("SU2^5", ((1, 1),) * 5),
                                ("USp14", ((7, 1),)),
                                ("USp6 x USp6", ((3, 1), (3, 1)))],
}


# --------------------------------------------------- ranks 7 and 8 for st_lib


def extend_library() -> None:
    for g in (7, 8):
        name = f"USp{2 * g}"
        if name in S.GROUPS:
            continue
        S.GROUPS[name] = S.Group(
            name=name, kind="sp", rank=g, epsilon=0, alpha_max=2.0 * g,
            variance=1.0, tail=S.edge_exponent(g, 0.5), realisable=True,
            note="big monodromy, generic hyperelliptic pencil")


def K_ranks(tau: float, maxrank: int) -> list[float]:
    """``log E[e^{tau tr}]`` over ``USp(2r)``, ``r = 1..maxrank``, from ONE
    table of Bessel functions ``I_0(2 tau) .. I_{2 maxrank}(2 tau)``."""
    saved = mp.dps
    mp.dps = 40 + int(maxrank * (maxrank - 1)
                      * math.log10(max(2.0 * tau, 10.0)))
    try:
        z = 2 * mpf(tau)
        I = [besseli(k, z) for k in range(2 * maxrank + 1)]
        out = []
        for r in range(1, maxrank + 1):
            m = mp.matrix(r, r)
            for i in range(1, r + 1):
                for j in range(1, r + 1):
                    m[i - 1, j - 1] = I[abs(i - j)] - I[i + j]
            out.append(float(mp.log(mp.det(m))))
        return out
    finally:
        mp.dps = saved


def wide_search(cap_alpha: float, kmax: int, npts: int, rows: list) -> None:
    """products of ``USp(2g)`` blocks with multiplicity <= kmax and
    alpha_max <= cap_alpha; cycles of the midrange comparison."""
    maxrank = int(cap_alpha // 2)
    tau = np.exp(np.linspace(math.log(1e-4), math.log(1e5), npts))
    grid = sorted(set([float(t) for t in tau]
                      + [float(k * t) for t in tau for k in range(2, kmax + 1)]))
    table: dict[float, list[float]] = {}
    for t in grid:
        table[t] = K_ranks(t, maxrank)

    facs = [(g, k) for g in range(1, maxrank + 1)
            for k in range(1, kmax + 1) if 2 * g * k <= cap_alpha]
    facs.sort(key=lambda f: (2 * f[0] * f[1], f))

    meas: list[tuple[tuple, ...]] = []

    def rec(start: int, used: float, chosen: list) -> None:
        if chosen:
            meas.append(tuple(chosen))
        for idx in range(start, len(facs)):
            g, k = facs[idx]
            if used + 2 * g * k <= cap_alpha + 1e-9:
                rec(idx, used + 2 * g * k, chosen + [(g, k)])

    rec(0, 0.0, [])
    lab = {}
    uniq = {}
    for m in meas:
        s = " x ".join(("SU2" if g == 1 else f"USp{2*g}") if k == 1
                       else f"{k}.{'SU2' if g == 1 else f'USp{2*g}'}"
                       for g, k in m)
        uniq.setdefault(s, m)
    keys = sorted(uniq)
    lib = [uniq[s] for s in keys]
    amax = np.array([sum(2 * g * k for g, k in m) for m in lib], dtype=float)
    psis = np.empty((len(lib), npts))
    for i, m in enumerate(lib):
        acc = np.zeros(npts)
        for g, k in m:
            acc += np.array([table[float(k * t)][g - 1] for t in tau])
        psis[i] = acc / tau

    n = len(lib)
    mid = np.zeros((n, n))
    for i in range(n):
        d = psis[i][None, :] - psis
        end = amax[i] - amax
        hi = np.maximum(np.maximum(d.max(axis=1), 0.0), end)
        lo = np.minimum(np.minimum(d.min(axis=1), 0.0), end)
        mid[i] = 0.5 * (hi + lo)
    np.fill_diagonal(mid, 0.0)

    strict = mid < -1e-6
    seen, cyc, total = set(), [], 0
    for i, j, k in itertools.combinations(range(n), 3):
        for a, b, c in ((i, j, k), (i, k, j)):
            total += 1
            if strict[a, b] and strict[b, c] and strict[c, a]:
                key = frozenset((a, b, c))
                if key not in seen:
                    seen.add(key)
                    cyc.append((a, b, c,
                                min(-mid[a, b], -mid[b, c], -mid[c, a])))
    best = max((m for *_, m in cyc), default=float("nan"))
    print(f"\n  alpha_max <= {cap_alpha:g}, multiplicity <= {kmax}: "
          f"{n} measures, {total} oriented triangles, {len(cyc)} distinct "
          f"3-cycles" + (f", widest margin {best:.4e}" if cyc else ""))
    rows.append(["wide-summary", cap_alpha, kmax, n, total, len(cyc),
                 f"{best:.6e}" if cyc else ""])
    for a, b, c, marg in sorted(cyc, key=lambda r: -r[3])[:12]:
        print(f"      {keys[a]:<26} < {keys[b]:<20} < {keys[c]:<26} "
              f"margin {marg:.4e}   genus "
              f"{amax[a]/2:g}/{amax[b]/2:g}/{amax[c]/2:g}")
        rows.append(["wide-cycle", cap_alpha, kmax, keys[a], keys[b], keys[c],
                     f"{marg:.6e}"])
    if len(cyc) > 12:
        print(f"      ... and {len(cyc) - 12} more")


def main() -> int:
    rows: list[list] = []

    print("=" * 96)
    print("1.  cross-check of the two independent implementations of K")
    print("=" * 96)
    for name, rank in (("SU2", 1), ("USp4", 2), ("USp6", 3), ("USp12", 6)):
        for t in (0.7, 3.0, 40.0):
            a = float(S.group_K(name, np.array([t]))[0])
            b = float(K_usp(rank, t))
            print(f"    {name:<6} tau={t:<5g}  st_lib {a:>18.12f}   "
                  f"Bessel {b:>18.12f}   diff {abs(a - b):.3e}")
            rows.append(["crosscheck", name, t, f"{a:.12f}", f"{b:.12f}",
                         f"{abs(a - b):.3e}"])

    print()
    print("=" * 96)
    print("2.  the cycles the new genus-four witness unlocks, at 40 digits")
    print("=" * 96)
    for label, tri in TRIANGLES.items():
        print(f"\n  {label}")
        margins = []
        for i in range(3):
            (na, A), (nb, B) = tri[i], tri[(i + 1) % 3]
            midv, sup, inf, sat, iat = midrange_mp(A, B)
            mp.dps = 42
            print(f"    {na + '  ->  ' + nb:<46}{mp.nstr(midv, 40):>46}")
            print(f"        sup {mp.nstr(sup, 25)} at tau = {mp.nstr(sat, 12)}")
            print(f"        inf {mp.nstr(inf, 25)} at tau = {mp.nstr(iat, 12)}")
            margins.append(-midv)
            rows.append(["cycle40", label, na, nb, mp.nstr(midv, 40),
                         mp.nstr(sup, 25), mp.nstr(inf, 25),
                         mp.nstr(sat, 12), mp.nstr(iat, 12)])
        strict = all(m > 0 for m in margins)
        print(f"    strict 3-cycle: {strict};  smallest margin "
              f"{mp.nstr(min(margins), 30)}")
        rows.append(["cycle40-margin", label, mp.nstr(min(margins), 30),
                     strict])

    if "--wide" in sys.argv:
        print()
        print("=" * 96)
        print("3.  the second front: the capped cone search pushed to genus 8")
        print("=" * 96)
        extend_library()
        print("  K_ranks (one Bessel table, all ranks) against K_usp:")
        worst = 0.0
        for t in (0.7, 3.0, 40.0, 1000.0):
            v = K_ranks(t, 8)
            for r in (1, 2, 3, 6):
                worst = max(worst, abs(v[r - 1] - float(K_usp(r, t))))
        print(f"    worst absolute difference over 16 probes: {worst:.3e}")
        rows.append(["K_ranks-crosscheck", f"{worst:.3e}"])
        for cap in (12.0, 14.0, 16.0):
            for kmax in (1, 2, 3):
                wide_search(cap, kmax, 601, rows)

    with (HERE / "genus4_cycles.csv").open("w", newline="",
                                           encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind"] + [f"c{i}" for i in range(8)])
        for r in rows:
            wr.writerow(list(r) + [""] * (9 - len(r)))
    print("\nwritten: genus4_cycles.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
