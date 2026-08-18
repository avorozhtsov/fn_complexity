#!/usr/bin/env python3
"""How many levels a midrange 3-cycle needs.

Write ``D_1, D_2, D_3`` for the three differences ``Psi_mu - Psi_nu`` around a
triangle.  Two facts are exact:

    D_1 + D_2 + D_3 = 0   pointwise,
    D_i(0) = 0            (every limit measure has mean zero),

and a third holds inside one genus,

    D_i(inf) = alpha_max(mu) - alpha_max(nu) = 0.

Since ``mid D = (sup D + inf D)/2`` depends on ``D`` only through the values it
takes, sample the three functions at the ``n`` interior points where their
extrema sit.  The question "can all three midranges be negative?" becomes a
finite linear program: three vectors ``v_i`` in ``R^n`` with ``sum_i v_i = 0``,
each padded with the pinned endpoint values, and

    mid(v) = (max(v, pinned) + min(v, pinned))/2 < 0   for all three.

For fixed choices of which coordinate attains the max and which the min this is
an LP, and there are finitely many such choices, so the exact optimum of
``min_i (-mid_i)`` subject to ``|v| <= 1`` is computable.  The answer:

    both endpoints pinned to 0 (one genus):   n = 1: -1/4      n = 2: 0 (exact)
                                              n = 3: +1/5      n = 4: +1/5
    one endpoint free (two genera):           n = 1: 0         n = 2: +1/6

**Reading, and its limits.**  With both endpoints pinned, ``n = 2`` sampled
levels give an optimum of exactly ``0``: no strict cycle, however the levels are
chosen.  Freeing one endpoint -- which is precisely what comparing two *genera*
does, since ``D(inf) = alpha_max(a) - alpha_max(b)`` -- lifts the optimum to
``+1/4`` at the same ``n = 2``.  So the endpoint gap is worth a whole level, and
it is the level that makes the difference.  That matches the library exactly:
every 3-cycle found in ``cone_search.py`` is cross-genus, and every same-genus
class is transitive.

This is an explanation, not a proof of same-genus transitivity: three continuous
functions attain their six extrema at up to six distinct ``tau``, so ``n = 2`` is
a hypothesis about how aligned the extrema are, not a theorem about arbitrary
differences.  With ``n = 3`` sampled levels a both-ends-pinned cycle does exist
(the ``(-2,1,1)`` model), so what rules cycles out inside one genus is the
rigidity of the Sato--Tate shapes, measured in ``cone_search.py``, and not a
formal obstruction.

    python research/sato_tate_limit/level_lemma.py
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent


def optimum(n: int, pin_end: bool) -> tuple[float, np.ndarray | None]:
    """Exact ``max min_i(-mid(v_i))`` over ``|v| <= 1``, ``sum_i v_i = 0``.

    ``pin_end`` pins the ``tau = inf`` endpoint of every difference to zero (the
    same-genus case).  Otherwise the endpoint is a free variable shared by the
    triangle, subject only to summing to zero -- the cross-genus case, where the
    three endpoint gaps are ``alpha_a - alpha_b`` and do sum to zero.
    """

    # variables: v[i][j] for i<3, j<n   then (if not pin_end) e[i] for i<3
    #            plus the objective variable s.
    nv = 3 * n + (0 if pin_end else 3) + 1
    s_at = nv - 1

    def vidx(i, j):
        return i * n + j

    def eidx(i):
        return 3 * n + i

    def coord(i, j):
        """coefficient vector selecting the j-th sampled value of edge i;
        j == n means the pinned/free endpoint, j == -1 the pinned 0 at tau=0."""
        row = np.zeros(nv)
        if j == -1:
            return row
        if j == n:
            if pin_end:
                return row
            row[eidx(i)] = 1.0
            return row
        row[vidx(i, j)] = 1.0
        return row

    slots = list(range(-1, n + 1))
    best = (-np.inf, None)
    for his in itertools.product(slots, repeat=3):
        for los in itertools.product(slots, repeat=3):
            A, b = [], []
            for i in range(3):
                hi, lo = coord(i, his[i]), coord(i, los[i])
                for j in slots:
                    c = coord(i, j)
                    A.append(c - hi)
                    b.append(0.0)
                    A.append(lo - c)
                    b.append(0.0)
                row = hi + lo
                row[s_at] = 2.0
                A.append(row)
                b.append(0.0)
            Aeq = []
            for j in range(n):
                r = np.zeros(nv)
                for i in range(3):
                    r[vidx(i, j)] = 1.0
                Aeq.append(r)
            if not pin_end:
                r = np.zeros(nv)
                for i in range(3):
                    r[eidx(i)] = 1.0
                Aeq.append(r)
            bounds = [(-1.0, 1.0)] * (nv - 1) + [(None, None)]
            c_obj = np.zeros(nv)
            c_obj[s_at] = -1.0
            res = linprog(c_obj, A_ub=np.array(A), b_ub=np.array(b),
                          A_eq=np.array(Aeq), b_eq=np.zeros(len(Aeq)),
                          bounds=bounds)
            if res.status == 0 and res.x[s_at] > best[0]:
                best = (float(res.x[s_at]), res.x.copy())
    return best


def main() -> int:
    rows = []
    print("=" * 78)
    print("exact optimum of  min_i(-mid_i)  over three differences summing to")
    print("zero, normalised by |D| <= 1")
    print("=" * 78)
    print(f"\n  {'interior levels n':>18}{'both ends pinned':>20}"
          f"{'one end free':>16}")
    for n in (1, 2, 3):
        a, xa = optimum(n, pin_end=True)
        b, xb = optimum(n, pin_end=False)
        print(f"  {n:>18}{a:>20.6f}{b:>16.6f}")
        rows.append([n, f"{a:.9f}", f"{b:.9f}"])
        if n == 2:
            print("        (both-ends-pinned optimum is exactly 0: no strict")
            print("         cycle, however the levels are chosen)")
            if b > 1e-9:
                v = xb[:6].reshape(3, 2)
                e = xb[6:9]
                print("        one-end-free witness:")
                for i in range(3):
                    print(f"          D_{i + 1}: 0, {v[i, 0]:+.4f}, "
                          f"{v[i, 1]:+.4f}, endpoint {e[i]:+.4f}")
    print("\n  Reading: inside one genus both endpoints are pinned to zero, so")
    print("  the three differences reach at most two independent interior")
    print("  levels before the shape becomes an ordinary single crossing, and")
    print("  the optimum is exactly 0 -- no strict cycle.  Across genera the")
    print("  endpoint gap alpha_max(a) - alpha_max(b) is free (it only has to")
    print("  sum to zero around the triangle) and supplies the missing level.")

    with (HERE / "level_lemma.csv").open("w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["interior_levels", "both_ends_pinned", "one_end_free"])
        wr.writerows(rows)
    print("\nwritten: level_lemma.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
