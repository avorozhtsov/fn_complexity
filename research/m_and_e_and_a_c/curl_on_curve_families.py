#!/usr/bin/env python3
"""The curl of A on curve families: a potential-free search statistic for cycles.

Brief B searches for cycles by looking for pairs where the computed comparison
disagrees with a scalar `phi`; addendum 1 sharpens `phi` to
`phi~ = M - ((3-2sqrt2)/2) m2` and observes that this is STILL a total order, so
the search must aim one order deeper.  That is a ladder: every time the
expansion is pushed further, a new scalar appears and has to be subtracted.

The gauge decomposition (see `gauge_decomposition.py`) ends the ladder.  Write

    L(a,b) = -log C(a->b) = S(a,b) + A(a,b),
    S = (L + L^T)/2 = d/2,     A = (L - L^T)/2,

so that `a < b  <=>  A(a,b) > 0`.  For an antisymmetric edge function on a
COMPLETE graph, the triangle sums

    curl A (a,b,c) = A(a,b) + A(b,c) + A(c,a)

all vanish exactly iff A = d psi for some potential psi.  Hence:

  * every scalar invariant, at every order of any expansion, contributes exactly
    zero to curl A.  phi, phi~, and whatever comes next are all annihilated
    without being computed;
  * a strict 3-cycle forces |curl A| = sum|A| over the triangle, so the ratio
    r = |curl A| / sum|A| lies in [0,1] and equals 1 exactly on a cycle;
  * max over triangles of |curl A| is a direct, expansion-free measure of how far
    the comparison is from being given by ANY scalar.

This script computes r and max|curl A| on pools of curve families over F_q and
reports them against the solver's precision floor.

    python research/m_and_e_and_a_c/curl_on_curve_families.py
"""
from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))

import ffmaps  # noqa: E402
import t2_2_common as C  # noqa: E402

FLOOR = 1e-10  # differences below this are ties (parent brief)


def pool(q: int, count: int, degree: int, seed: int = 5):
    """Distinct hyperelliptic fiber signatures y^2 = P(x) + c over F_q."""
    rng = np.random.default_rng(seed)
    seen, out = set(), []
    for _coeffs, counts in ffmaps.random_hyperelliptic(q, count, degree, rng):
        key = tuple(int(v) for v in np.sort(counts)[::-1])
        if key in seen:
            continue
        seen.add(key)
        out.append(C.sig_from_counts(counts, name=f"h{len(out)}"))
    return out


def antisymmetric(sigs):
    n = len(sigs)
    a = np.zeros((n, n))
    beta_arg = {}
    for i, j in itertools.combinations(range(n), 2):
        cij, bij = C.rate(sigs[i], sigs[j])
        cji, bji = C.rate(sigs[j], sigs[i])
        a[i, j] = 0.5 * (-math.log(cij) + math.log(cji))
        a[j, i] = -a[i, j]
        beta_arg[(i, j)] = (bij, bji)
    return a, beta_arg


def scan(a: np.ndarray):
    """Best (largest r) triangle and the largest |curl A| anywhere."""
    n = a.shape[0]
    best_r, best_curl = (0.0, None), (0.0, None)
    for i, j, k in itertools.combinations(range(n), 3):
        edges = (a[i, j], a[j, k], a[k, i])
        total = abs(sum(edges))
        denom = sum(abs(e) for e in edges)
        if denom == 0.0:
            continue
        if total / denom > best_r[0]:
            best_r = (total / denom, (i, j, k, min(abs(e) for e in edges)))
        if total > best_curl[0]:
            best_curl = (total, (i, j, k))
    return best_r, best_curl


def main() -> None:
    print("curl A on hyperelliptic families y^2 = P(x) + c over F_q")
    print("r = |curl A| / sum|A|;  r = 1 exactly on a strict 3-cycle\n")
    print(f"{'q':>5} {'g':>3} {'n':>4} {'triangles':>10} "
          f"{'max r':>12} {'max |curl A|':>14} {'median |A|':>12} {'margin':>10}")
    for q, degree in ((101, 5), (211, 5), (211, 7), (503, 5), (1009, 5)):
        sigs = pool(q, 90, degree)
        if len(sigs) < 3:
            print(f"{q:>5} {(degree - 1) // 2:>3} {len(sigs):>4}   too few distinct signatures")
            continue
        a, _ = antisymmetric(sigs)
        (r, r_info), (curl, _) = scan(a)
        off = np.abs(a[np.triu_indices_from(a, 1)])
        print(f"{q:>5} {(degree - 1) // 2:>3} {len(sigs):>4} "
              f"{math.comb(len(sigs), 3):>10} {r:>12.6f} {curl:>14.3e} "
              f"{np.median(off):>12.3e} {r_info[3]:>10.3e}")

    print(f"\nprecision floor for a rate difference: {FLOOR:.0e}")
    print("Read the table as follows.  max r is how close any triangle comes to")
    print("closing; r = 1 would be a cycle.  max |curl A| is the size of the")
    print("expansion-free residual -- everything a scalar invariant CANNOT explain,")
    print("at every order at once.  Compare it to the floor: a residual below")
    print("1e-10 is not a small violation, it is no measurement at all.")


if __name__ == "__main__":
    main()
