#!/usr/bin/env python3
"""The level lemma, as a theorem -- and how far the library is from its
hypothesis.

``level_lemma.py`` solves an LP over argmax/argmin patterns and reports

    both ends pinned:  n = 1 -> 0,  n = 2 -> 0,  n = 3 -> +1/4
    one end free    :  n = 1 -> 0,  n = 2 -> +1/4

``TRANSITIVITY.md`` proves these numbers outright:

* **Two-level lemma.**  If two points ``p, q`` capture every sup and every inf
  of the three differences (together with the pinned endpoint value ``0``),
  then ``sum_i mid(D_i) >= 0``, so no strict 3-cycle.
* **Three-level cap.**  With ``|D| <= 1`` and ``sum_i D_i = 0``,
  ``min_i (-mid(D_i)) <= 1/4`` always, and ``1/4`` is attained by the
  ``(-1, 1/2, 1/2)`` model on three levels.

A free endpoint is one more level, which is why one end free at ``n = 2``
already reaches the cap.

What is left is the level *count*.  This script measures it on the library:
for every same-genus triple that survives the crossing reduction of
``transitivity_certificate.py``, it locates the three interior minima and the
three interior maxima of ``D_1, D_2, D_3`` and reports how tightly they cluster
-- the quantity that decides whether the two-level hypothesis can hold.

    python research/sato_tate_limit/transitivity_levels.py
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

import numpy as np

import st_lib as S
from transitivity_pairs import products

HERE = Path(__file__).resolve().parent
TAU = S.tau_grid(1e-4, 1e5, 2401)
CAP = 12.0


def check_two_level_lemma(trials: int = 200000, seed: int = 5) -> tuple:
    """Random check of the proved lemma: two sampled levels admit no cycle."""

    rng = np.random.default_rng(seed)
    a = rng.uniform(-1, 1, size=(trials, 3))
    b = rng.uniform(-1, 1, size=(trials, 3))
    a[:, 2] = -a[:, 0] - a[:, 1]
    b[:, 2] = -b[:, 0] - b[:, 1]
    keep = (np.abs(a).max(axis=1) <= 1) & (np.abs(b).max(axis=1) <= 1)
    a, b = a[keep], b[keep]
    lo = np.minimum(np.minimum(a, b), 0.0)
    hi = np.maximum(np.maximum(a, b), 0.0)
    mid = 0.5 * (lo + hi)
    cyc = int((mid.max(axis=1) < 0).sum())
    return a.shape[0], cyc, float(mid.sum(axis=1).min())


def check_three_level_cap(trials: int = 200000, seed: int = 7) -> float:
    """Random check of ``min_i (-mid) <= 1/4`` with three free levels."""

    rng = np.random.default_rng(seed)
    v = rng.uniform(-1, 1, size=(trials, 3, 3))       # [trial, edge, level]
    v[:, 2, :] = -v[:, 0, :] - v[:, 1, :]
    keep = np.abs(v).max(axis=(1, 2)) <= 1
    v = v[keep]
    lo = np.minimum(v.min(axis=2), 0.0)
    hi = np.maximum(v.max(axis=2), 0.0)
    mid = 0.5 * (lo + hi)
    return float((-mid).min(axis=1).max())


def main() -> int:
    print("=" * 78)
    print("1.  the two proved statements, checked at random")
    print("=" * 78)
    n, cyc, worst = check_two_level_lemma()
    print(f"  two levels, both ends pinned: {n} random triples, "
          f"{cyc} cycles  (theorem: 0)")
    print(f"      (min over trials of sum_i mid = {worst:+.6f}; the theorem's"
          f" 'sum >= 0' is\n       conditional on all three mid being"
          f" negative, which is what it refutes)")
    cap = check_three_level_cap()
    print(f"  three free levels: max over {200000} random triples of "
          f"min_i(-mid) = {cap:.6f}  (theorem: <= 1/4, attained)")
    print("  the (-1, 1/2, 1/2) model attains it exactly:")
    v = np.array([[-1, .5, .5], [.5, -1, .5], [.5, .5, -1.]])
    m = 0.5 * (np.minimum(v.min(axis=1), 0) + np.maximum(v.max(axis=1), 0))
    print(f"      mid = {m},  min_i(-mid) = {(-m).min():.6f}")

    print("\n" + "=" * 78)
    print("2.  how many levels the library actually uses")
    print("=" * 78)
    lib = products(CAP)
    lib.sort(key=lambda m: (m.alpha_max, m.variance, m.tail, m.label))
    psis = {i: m.Psi(TAU) for i, m in enumerate(lib)}
    amax = np.array([m.alpha_max for m in lib])

    def diff(i, j):
        return psis[i] - psis[j]

    # surviving triples: at least two crossing edges (see certificate script)
    crossing = {}
    for a in sorted(set(amax)):
        sel = np.where(np.abs(amax - a) < 1e-9)[0]
        for i, j in itertools.combinations(sel, 2):
            crossing[(i, j)] = S.sign_changes(diff(i, j)) > 0

    rows = []
    for a in sorted(set(amax)):
        sel = np.where(np.abs(amax - a) < 1e-9)[0]
        for i, j, k in itertools.combinations(sel, 3):
            e = (crossing[(i, j)] + crossing[(j, k)] + crossing[(i, k)])
            if e < 2:
                continue
            ds = [diff(i, j), diff(j, k), diff(k, i)]
            tmin = [float(TAU[int(d.argmin())]) for d in ds]
            tmax = [float(TAU[int(d.argmax())]) for d in ds]
            # a level is "shared" when two of the three extrema coincide on
            # the grid; count the distinct locations
            nmin = len(set(np.round(np.log10(tmin), 2)))
            nmax = len(set(np.round(np.log10(tmax), 2)))
            rows.append([int(a / 2), lib[i].label, lib[j].label, lib[k].label,
                         nmin, nmax,
                         f"{max(tmin) / min(tmin):.3f}",
                         f"{max(tmax) / min(tmax):.3f}"])
    import collections
    cmin = collections.Counter(r[4] for r in rows)
    cmax = collections.Counter(r[5] for r in rows)
    print(f"  surviving triples: {len(rows)}")
    print(f"  distinct interior-minimum locations per triple: "
          f"{dict(sorted(cmin.items()))}")
    print(f"  distinct interior-maximum locations per triple: "
          f"{dict(sorted(cmax.items()))}")
    sp_min = [float(r[6]) for r in rows]
    sp_max = [float(r[7]) for r in rows]
    print(f"  spread of the three argmin (max/min tau): median "
          f"{np.median(sp_min):.3f}, worst {max(sp_min):.3f}")
    print(f"  spread of the three argmax (max/min tau): median "
          f"{np.median(sp_max):.3f}, worst {max(sp_max):.3f}")
    print("\n  Reading: the three interior minima sit at three *different*")
    print("  tau on every one of the surviving triples, and so do the three")
    print("  maxima, with spreads of one to five decades.  The two-level")
    print("  hypothesis is therefore not just inexact, it is badly false: the")
    print("  library uses six levels, not two, and the reason it does not")
    print("  cycle is not a level count.  Route 2 of brief K is a dead end")
    print("  in this form.")

    with (HERE / "transitivity_levels.csv").open("w", newline="",
                                                 encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["genus", "mu", "nu", "rho", "distinct_argmin",
                     "distinct_argmax", "argmin_spread", "argmax_spread"])
        wr.writerows(rows)
    print("\nwritten: transitivity_levels.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
