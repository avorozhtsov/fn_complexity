#!/usr/bin/env python3
"""How low in the genus can a same-``alpha_max`` 3-cycle be pushed?

``transitivity_mixtures.py`` finds same-genus 3-cycles among *uniform*
2-mixtures of multiplicity-free symplectic products, the first at genus 8.  The
mixing weights are a free parameter, so this script takes the closest triples of
each genus and optimises the three weights, to find the lowest genus at which
the class of mean-zero measures with a given ``alpha_max`` stops being totally
ordered by the midrange.

Each vertex is ``p mu_a + (1-p) mu_b`` for two partitions ``a, b`` of the genus;
the objective is the cycle margin ``-max_i mid_i`` over the three weights.

    python research/sato_tate_limit/transitivity_freeweight.py
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

import kappa_lib as KL
from transitivity_dominance import partitions
from transitivity_mixtures import mid_matrix, mid_of, mixture_K

HERE = Path(__file__).resolve().parent
TAU = np.geomspace(1e-4, 1e5, 1201)
TAU_B = np.geomspace(1e-5, 1e6, 1501)
GLO, GHI = 4, 9
NCAND = 40


def margin(p, ing, Ks, tau) -> float:
    """cycle margin of the triple with weights ``p``; positive means a cycle."""

    K = [mixture_K(Ks[a], Ks[b], float(w)) for (a, b), w in zip(ing, p)]
    psi = [k / tau for k in K]
    m = [mid_of(psi[0] - psi[1]), mid_of(psi[1] - psi[2]),
         mid_of(psi[2] - psi[0])]
    return -max(m)


def main() -> int:
    kappa, _ = KL.kappa_and_b(TAU, 12)
    kappa_b, _ = KL.kappa_and_b(TAU_B, 12)
    rows = []
    print("=" * 78)
    print("free-weight 2-mixtures: the lowest genus that cycles")
    print("=" * 78)
    print(f"  {'genus':>5}{'candidates':>12}{'best margin':>16}"
          f"{'cycle?':>9}   winner")
    for genus in range(GLO, GHI + 1):
        lams = partitions(genus)
        labels = ["".join(map(str, lam)) for lam in lams]
        Ks = [sum(kappa[p - 1] for p in lam) for lam in lams]
        Ks_b = [sum(kappa_b[p - 1] for p in lam) for lam in lams]
        n = len(lams)

        # ingredient pairs: every unordered pair of partitions, plus the
        # degenerate "pair" (a, a) which is the pure product mu_a
        ing_all = [(i, i) for i in range(n)]
        ing_all += list(itertools.combinations(range(n), 2))
        mKs = [mixture_K(Ks[a], Ks[b], 0.5) for a, b in ing_all]
        mid = mid_matrix(mKs, TAU)

        cands = []
        N = len(ing_all)
        for i, j, k in itertools.combinations(range(N), 3):
            for a, b, c in ((i, j, k), (i, k, j)):
                cands.append((max(mid[a, b], mid[b, c], mid[c, a]),
                              (a, b, c)))
        cands.sort()
        cands = cands[:NCAND]

        best = (-np.inf, None, None)
        for _, trip in cands:
            ing = [ing_all[t] for t in trip]
            f = lambda p: -margin(p, ing, Ks, TAU)
            res = differential_evolution(f, [(0.0, 1.0)] * 3, seed=17,
                                         maxiter=60, popsize=12, tol=1e-10,
                                         polish=True)
            pol = minimize(f, res.x, method="Nelder-Mead",
                           options={"xatol": 1e-12, "fatol": 1e-14,
                                    "maxiter": 4000})
            p = pol.x if pol.fun < res.fun else res.x
            val = -min(pol.fun, res.fun)
            if val > best[0]:
                best = (val, ing, np.clip(p, 0.0, 1.0))
        val, ing, p = best
        # independent grid
        valb = margin(p, ing, Ks_b, TAU_B)
        names = [labels[a] if a == b else f"{p_i:.6f}*{labels[a]}"
                 f" + {1 - p_i:.6f}*{labels[b]}"
                 for (a, b), p_i in zip(ing, p)]
        print(f"  {genus:>5}{len(cands):>12}{val:>16.8f}"
              f"{('YES' if val > 0 else 'no'):>9}")
        for nm in names:
            print(f"            {nm}")
        print(f"            grid B margin {valb:+.8f}")
        rows.append([genus, f"{val:.10g}", f"{valb:.10g}",
                     val > 0] + names)

    with (HERE / "transitivity_freeweight.csv").open("w", newline="",
                                                     encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["genus", "margin_gridA", "margin_gridB", "is_cycle",
                     "vertex1", "vertex2", "vertex3"])
        wr.writerows(rows)
    print("\nwritten: transitivity_freeweight.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
