#!/usr/bin/env python3
"""Independent verification of the limiting 3-cycles.

Each of the three midranges is recomputed three ways:

* on the search grid, ``tau`` geometric from ``1e-4`` to ``1e5``, 3001 points;
* on an independent grid, ``1e-5`` to ``1e6``, 1501 points -- different range,
  different spacing, no shared node except by accident;
* with the two extrema located by golden section on ``log tau`` rather than read
  off a grid, so the reported ``sup`` and ``inf`` are exact to ``1e-12`` and the
  grid discretisation drops out entirely.

Also reported: the contact temperatures, which say where on the ``tau`` axis the
comparison is decided, and (for the reader's orientation) the ``(alpha_max, m2,
t)`` coordinates of each measure.

    python research/sato_tate_limit/verify_cycle.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

import st_lib as S

HERE = Path(__file__).resolve().parent

E = S.Factor("SU2")


def M(*factors: S.Factor) -> S.Measure:
    return S.Measure(tuple(factors))


CYCLES = {
    "widest symplectic": [
        M(S.Factor("SU2", 4)),                      # genus 4, Jac ~ E^4
        M(S.Factor("USp6"), S.Factor("USp6")),      # genus 6, Jac ~ A_3 x A_3'
        M(E, S.Factor("USp4", 2)),                  # genus 5, Jac ~ E' x A^2
    ],
    "multiplicity two only": [
        M(S.Factor("SU2", 2), S.Factor("SU2", 2)),  # genus 4, Jac ~ E_1^2 x E_2^2
        M(S.Factor("USp12"),),                      # genus 6, generic pencil
        M(E, E, S.Factor("USp6")),                  # genus 5, E_3 x E_4 x A_3
    ],
    "two vertices multiplicity-free": [
        M(S.Factor("SU2", 4)),                      # genus 4, Jac ~ E^4
        M(S.Factor("USp4"), S.Factor("USp8")),      # genus 6, Jac ~ A_2 x A_4
        M(E, E, E, E, E),                           # genus 5, five elliptic
    ],
    "torus (not realisable, for contrast)": [
        M(S.Factor("U1", 3)),
        M(E, E, S.Factor("SU2", 2)),
        M(S.Factor("U1"), S.Factor("U1", 1, 0.5), S.Factor("U1", 1, 0.5),
          S.Factor("U1", 1, 0.5)),
    ],
}


def extrema(a: S.Measure, b: S.Measure, tau: np.ndarray):
    d = a.K(tau) / tau - b.K(tau) / tau
    logtau = np.log(tau)

    def value(s: float) -> float:
        t = np.array([math.exp(s)])
        return float((a.K(t) - b.K(t))[0] / t[0])

    def polish(i: int, sign: float) -> tuple[float, float]:
        lo, hi = logtau[max(i - 1, 0)], logtau[min(i + 1, logtau.size - 1)]
        phi = (math.sqrt(5.0) - 1.0) / 2.0
        x1, x2 = hi - phi * (hi - lo), lo + phi * (hi - lo)
        f1, f2 = sign * value(x1), sign * value(x2)
        for _ in range(80):
            if f1 > f2:
                hi, x2, f2 = x2, x1, f1
                x1 = hi - phi * (hi - lo)
                f1 = sign * value(x1)
            else:
                lo, x1, f1 = x1, x2, f2
                x2 = lo + phi * (hi - lo)
                f2 = sign * value(x2)
            if hi - lo < 1e-13:
                break
        s = x1 if f1 > f2 else x2
        return sign * max(f1, f2), math.exp(s)

    hi_v, hi_t = polish(int(np.argmax(d)), +1.0)
    lo_v, lo_t = polish(int(np.argmin(d)), -1.0)
    end = a.alpha_max - b.alpha_max
    sup, sup_at = (hi_v, hi_t) if hi_v >= max(0.0, end) else \
        ((end, math.inf) if end >= 0.0 else (0.0, 0.0))
    inf, inf_at = (lo_v, lo_t) if lo_v <= min(0.0, end) else \
        ((end, math.inf) if end <= 0.0 else (0.0, 0.0))
    return 0.5 * (sup + inf), sup, inf, sup_at, inf_at


def main() -> int:
    grid_a = S.tau_grid(1e-4, 1e5, 3001)
    grid_b = S.tau_grid(1e-5, 1e6, 1501)
    rows = []
    for name, tri in CYCLES.items():
        print("=" * 78)
        print(name)
        print("=" * 78)
        for m in tri:
            print(f"  {m.label:<34} genus {m.alpha_max / 2:>3g}   "
                  f"alpha_max {m.alpha_max:>4g}   m2 {m.variance:>5g}   "
                  f"t {m.tail:>5g}")
        print()
        print(f"  {'edge':<52}{'grid A':>13}{'grid B':>13}{'polished':>13}")
        margins = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            a, b = tri[i], tri[j]

            def grid_mid(g):
                d = a.K(g) / g - b.K(g) / g
                dd = np.concatenate([[0.0], d, [a.alpha_max - b.alpha_max]])
                return 0.5 * (float(dd.max()) + float(dd.min()))

            ma, mb = grid_mid(grid_a), grid_mid(grid_b)
            mp_, sup, inf, sup_at, inf_at = extrema(a, b, grid_a)
            print(f"  {a.label + ' -> ' + b.label:<52}"
                  f"{ma:>13.8f}{mb:>13.8f}{mp_:>13.8f}")
            print(f"      sup {sup:+.8f} at tau = {sup_at:.4f};   "
                  f"inf {inf:+.8f} at tau = {inf_at:.4f}")
            margins.append(-mp_)
            rows.append([name, a.label, b.label, f"{ma:.10f}", f"{mb:.10f}",
                         f"{mp_:.10f}", f"{sup:.10f}", f"{inf:.10f}",
                         f"{sup_at:.6g}", f"{inf_at:.6g}"])
        strict = all(m > 0 for m in margins)
        print(f"\n  all three midranges negative -- strict 3-cycle: {strict}")
        print(f"  smallest margin {min(margins):.6e}\n")

    with (HERE / "verify_cycle.csv").open("w", newline="",
                                          encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["cycle", "from", "to", "mid_grid_A", "mid_grid_B",
                     "mid_polished", "sup", "inf", "argsup_tau", "arginf_tau"])
        wr.writerows(rows)
    print("written: verify_cycle.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
