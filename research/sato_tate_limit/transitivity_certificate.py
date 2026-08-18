#!/usr/bin/env python3
"""Same-genus transitivity of the full symplectic library, reduced to a
handful of triples.

Two facts do almost all the work.

**(P1) Pointwise domination is transitive.**  If ``Psi_mu <= Psi_nu``
everywhere then ``mid(Psi_mu - Psi_nu) <= 0``; and ``<=`` between functions is
a partial order.  Hence *at most one edge of a 3-cycle can be a non-crossing
edge*: two non-crossing edges of a 3-cycle are consecutive, and chaining their
dominations contradicts the third edge.  (``FINDINGS.md`` says a 3-cycle needs
**all three** pairs to cross; that is one edge too strong -- see the
corrections section of ``TRANSITIVITY.md``.)

**(P2) Same-genus differences cross at most once.**  Sign changes of
``Psi_mu - Psi_nu`` are sign changes of ``M_mu - M_nu = int e^{tau a} d(mu-nu)``;
the kernel ``e^{tau a}`` is totally positive, so by Schoenberg's
variation-diminishing property that number is at most the number of sign
changes of ``mu - nu`` in ``a``, and both measures being symmetric halves it.
On the library it is 0 or 1 on every one of the 765 same-genus pairs.

Together: 718 of the 765 pairs are pointwise dominations, only 47 cross, and a
3-cycle must use at least two of those 47.  This script enumerates the triples
that survive and evaluates them.

    python research/sato_tate_limit/transitivity_certificate.py
"""

from __future__ import annotations

import csv
import itertools
from collections import defaultdict
from pathlib import Path

import numpy as np

import st_lib as S
from transitivity_pairs import products

HERE = Path(__file__).resolve().parent
TAU = S.tau_grid(1e-4, 1e5, 1201)
CAP = 12.0


def main() -> int:
    lib = products(CAP)
    lib.sort(key=lambda m: (m.alpha_max, m.variance, m.tail, m.label))
    amax = np.array([m.alpha_max for m in lib])
    var = np.array([m.variance for m in lib])
    tail = np.array([m.tail for m in lib])
    psis = np.array([m.Psi(TAU) for m in lib])

    print("=" * 78)
    print("1.  crossings inside a genus")
    print("=" * 78)
    mid: dict[tuple[int, int], float] = {}
    cross: dict[tuple[int, int], int] = {}
    by_genus: dict[float, list[int]] = defaultdict(list)
    for i, m in enumerate(lib):
        by_genus[m.alpha_max].append(i)
    npairs = ncross = 0
    tie_rows = []
    for a, sel in sorted(by_genus.items()):
        for i, j in itertools.combinations(sel, 2):
            d = psis[i] - psis[j]
            dd = np.concatenate([[0.0], d, [0.0]])
            m = 0.5 * (float(dd.max()) + float(dd.min()))
            sc = S.sign_changes(d)
            mid[(i, j)], mid[(j, i)] = m, -m
            cross[(i, j)] = cross[(j, i)] = sc
            npairs += 1
            ncross += sc > 0
            dm, dt = var[i] - var[j], tail[i] - tail[j]
            como = dm * dt > 0
            if como != (sc > 0):
                tie_rows.append([lib[i].label, lib[j].label, dm, dt, sc,
                                 f"{m:.8f}"])
    print(f"  same-genus pairs {npairs};  crossing {ncross};  "
          f"pointwise dominations {npairs - ncross}")
    print(f"  maximum number of sign changes on any pair: "
          f"{max(cross.values())}")
    print(f"\n  pairs where 'crossing <=> (m2, t) comonotone' is not decided "
          f"by the criterion: {len(tie_rows)}")
    for r in tie_rows:
        print(f"      {r[0]:<28}{r[1]:<28} dm2={r[2]:>3g} dt={r[3]:>6g} "
              f"crossings={r[4]}  mid={r[5]}")

    print("\n" + "=" * 78)
    print("2.  triples that survive (P1): at least two crossing edges")
    print("=" * 78)
    rows = []
    total_triples = surviving = 0
    worst = np.inf
    for a, sel in sorted(by_genus.items()):
        cnt = 0
        for i, j, k in itertools.combinations(sel, 3):
            total_triples += 1
            e = cross[(i, j)] + cross[(j, k)] + cross[(i, k)]
            if e >= 2:
                cnt += 1
                surviving += 1
                v = min(max(mid[(i, j)], mid[(j, k)], mid[(k, i)]),
                        max(mid[(i, k)], mid[(k, j)], mid[(j, i)]))
                worst = min(worst, v)
                rows.append([a / 2, lib[i].label, lib[j].label, lib[k].label,
                             e, f"{v:.8f}"])
        ntri = len(sel) * (len(sel) - 1) * (len(sel) - 2) // 6
        print(f"  genus {int(a / 2)}:  {len(sel):>3} measures, {ntri:>6} "
              f"unordered triples, {cnt:>4} with >= 2 crossing edges")
    print(f"\n  total unordered same-genus triples : {total_triples}")
    print(f"  surviving (P1)                     : {surviving}")
    print(f"  of these, cyclic                   : "
          f"{sum(1 for r in rows if float(r[5]) < 0)}")
    print(f"  worst (least positive) certificate : {worst:.8f}")

    print("\n  the ten tightest surviving triples:")
    rows.sort(key=lambda r: float(r[5]))
    print(f"  {'g':>2}  {'mu':<26}{'nu':<26}{'rho':<26}{'edges':>6}"
          f"{'max mid':>11}")
    for r in rows[:10]:
        print(f"  {int(r[0]):>2}  {r[1]:<26}{r[2]:<26}{r[3]:<26}{r[4]:>6}"
              f"{float(r[5]):>11.6f}")

    with (HERE / "transitivity_certificate.csv").open("w", newline="",
                                                      encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["genus", "mu", "nu", "rho", "crossing_edges",
                     "best_max_mid"])
        wr.writerows(rows)
    print("\nwritten: transitivity_certificate.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
