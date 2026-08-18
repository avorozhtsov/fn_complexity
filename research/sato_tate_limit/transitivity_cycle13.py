#!/usr/bin/env python3
"""The genus-13 same-genus 3-cycle, verified.

``transitivity_dominance.py`` finds the first 3-cycle among **multiplicity-free**
symplectic products at genus 13:

    USp14 x USp4^2 x SU2^2   <   USp10 x USp8 x USp4^2   <   USp12 x USp8 x SU2^3   <  ...
       (7,2,2,1,1)                   (5,4,2,2)                  (6,4,1,1,1)

All three have ``alpha_max = 26``, none has a repeated isogeny factor, so the
obstruction that brief F's cross-genus cycles all needed -- a Jacobian
isogenous to ``A^k`` with ``k >= 2`` -- is absent here.

This script re-verifies it independently of the genus-15 sweep:

* only the ranks actually used (``g in {1,2,4,5,6,7}``) enter, so the working
  precision is far higher relative to the matrix size;
* the grid is wider and finer (``10^-6 .. 10^8``, 3001 points against
  ``10^-4 .. 10^5``, 1201);
* the interior extrema are polished by golden section on ``log tau``;
* the tail beyond the grid is bounded analytically, ``|D| <= |dt| log(tau)/tau``
  eventually, and checked numerically to be decaying;
* the ``USp(2g)`` blocks are re-validated against their known cumulants
  (``kappa_2 = 1`` for all ``g``; ``kappa_4 = -1`` for ``g = 1`` and ``0`` for
  ``g >= 2``).

    python research/sato_tate_limit/transitivity_cycle13.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from mpmath import mp, mpf

import kappa_lib as KL

HERE = Path(__file__).resolve().parent
GRID = np.geomspace(1e-6, 1e8, 3001)
EXTRA = 40

CYCLES = [
    ("genus 13, margin 8.0e-3",
     [(7, 2, 2, 1, 1), (5, 4, 2, 2), (6, 4, 1, 1, 1)]),
    ("genus 13, margin 4.3e-4",
     [(8, 1, 1, 1, 1, 1), (4, 4, 4, 1), (6, 3, 2, 2)]),
]


def label(lam):
    return " x ".join(f"USp{2 * p}" if p > 1 else "SU2" for p in lam)


class Prod:
    def __init__(self, lam):
        self.lam = tuple(lam)
        self.alpha_max = mpf(2 * sum(lam))
        self._z = {}
        for g in set(lam):
            dps = KL.working_dps(g + 1, 1.0) + EXTRA
            h = KL.hankels(0.0, g, dps)
            with mp.workdps(dps):
                self._z[g] = mp.log(h[g])

    def K(self, tau):
        t = mpf(tau)
        out = mpf(0)
        for g in self.lam:
            dps = KL.working_dps(g + 1, float(t)) + EXTRA
            h = KL.hankels(t, g, dps)
            with mp.workdps(dps):
                out += mp.log(h[g]) - self._z[g]
        return out

    def Psi(self, tau):
        return self.K(tau) / mpf(tau)


def golden(f, lo, hi, maximise=True, iters=70):
    phi = (mp.sqrt(5) - 1) / 2
    a, b = mp.log(mpf(lo)), mp.log(mpf(hi))
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = f(mp.e ** c), f(mp.e ** d)
    for _ in range(iters):
        if (fc > fd) if maximise else (fc < fd):
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = f(mp.e ** c)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = f(mp.e ** d)
    x = mp.e ** ((a + b) / 2)
    return x, f(x)


def main() -> int:
    mp.dps = 50
    print("=" * 78)
    print("0.  cumulant check on the blocks used")
    print("=" * 78)
    # K(tau) = k2 tau^2/2 + k4 tau^4/24 + O(tau^6): two samples give both
    h = mpf("1e-3")
    print(f"  {'g':>3}{'kappa_2 (must be 1)':>28}"
          f"{'kappa_4 (-1 at g=1, else 0)':>32}")
    for g in (1, 2, 3, 4, 5, 6, 7, 8):
        q = Prod((g,))
        f1, f2 = q.K(h), q.K(2 * h)
        k2 = (16 * f1 - f2) / (6 * h ** 2)
        k4 = 2 * (f2 - 4 * f1) / h ** 4
        print(f"  {g:>3}{mp.nstr(k2, 20):>28}{mp.nstr(k4, 10):>32}")

    rows = []
    kappa_needed = sorted({p for _, cyc in CYCLES for lam in cyc for p in lam})
    print(f"\n  ranks in play: {kappa_needed}")
    kap, _ = KL.kappa_and_b(GRID, max(kappa_needed))

    for name, cyc in CYCLES:
        print("\n" + "=" * 78)
        print(name)
        print("=" * 78)
        prods = [Prod(l) for l in cyc]
        Ks = [sum(kap[p - 1] for p in l) for l in cyc]
        worst_grid = -np.inf
        worst_pol = None
        for i in range(3):
            u, v = prods[i], prods[(i + 1) % 3]
            Ku, Kv = Ks[i], Ks[(i + 1) % 3]
            d = (Ku - Kv) / GRID
            dd = np.concatenate([[0.0], d, [0.0]])
            hi_g, lo_g = float(dd.max()), float(dd.min())
            mid_g = 0.5 * (hi_g + lo_g)
            worst_grid = max(worst_grid, mid_g)
            ihi, ilo = int(dd.argmax()), int(dd.argmin())

            def D(t, u=u, v=v):
                return u.Psi(t) - v.Psi(t)

            hi_p, lo_p = mpf(0), mpf(0)
            thi = tlo = "endpoint"
            if 0 < ihi <= GRID.size:
                j = ihi - 1
                a = GRID[max(j - 1, 0)]
                b = GRID[min(j + 1, GRID.size - 1)]
                t, val = golden(D, a, b, maximise=True)
                if val > hi_p:
                    hi_p, thi = val, t
            if 0 < ilo <= GRID.size:
                j = ilo - 1
                a = GRID[max(j - 1, 0)]
                b = GRID[min(j + 1, GRID.size - 1)]
                t, val = golden(D, a, b, maximise=False)
                if val < lo_p:
                    lo_p, tlo = val, t
            mid_p = (hi_p + lo_p) / 2
            worst_pol = mid_p if worst_pol is None else max(worst_pol, mid_p)
            print(f"  edge {label(cyc[i])}")
            print(f"    ->  {label(cyc[(i + 1) % 3])}")
            print(f"    mid  grid  {mid_g:+.10f}")
            print(f"    mid  polished  {mp.nstr(mid_p, 30)}")
            print(f"      sup {mp.nstr(hi_p, 16)} at tau = "
                  f"{thi if isinstance(thi, str) else mp.nstr(thi, 10)}")
            print(f"      inf {mp.nstr(lo_p, 16)} at tau = "
                  f"{tlo if isinstance(tlo, str) else mp.nstr(tlo, 10)}")
            # tail check
            tails = [mpf("1e5"), mpf("1e6"), mpf("1e7"), mpf("1e8")]
            tv = [D(t) for t in tails]
            print("      tail |D| at tau = 1e5,1e6,1e7,1e8: "
                  + ", ".join(mp.nstr(abs(x), 5) for x in tv))
            assert all(abs(tv[k + 1]) < abs(tv[k]) for k in range(3)), \
                "tail not decaying"
            rows.append([name, label(cyc[i]), label(cyc[(i + 1) % 3]),
                         f"{mid_g:.10f}", mp.nstr(mid_p, 30),
                         mp.nstr(hi_p, 16), mp.nstr(lo_p, 16),
                         thi if isinstance(thi, str) else mp.nstr(thi, 10),
                         tlo if isinstance(tlo, str) else mp.nstr(tlo, 10)])
        print(f"\n  margin (grid)     {-worst_grid:.10f}")
        print(f"  margin (polished) {mp.nstr(-worst_pol, 20)}")
        print(f"  verdict: {'3-CYCLE' if worst_pol < 0 else 'transitive'}")
        rows.append([name, "MARGIN", "", f"{-worst_grid:.10f}",
                     mp.nstr(-worst_pol, 20), "", "", "", ""])

    with (HERE / "transitivity_cycle13.csv").open("w", newline="",
                                                  encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["cycle", "mu", "nu", "mid_grid", "mid_polished",
                     "sup", "inf", "tau_sup", "tau_inf"])
        wr.writerows(rows)
    print("\nwritten: transitivity_cycle13.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
