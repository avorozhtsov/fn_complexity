#!/usr/bin/env python3
"""E3.4 -- brief C's claim, with the numbers in it.

Brief C proposed: *no real-valued invariant of such families is compatible with
asymptotic conversion*.  With the flux in hand the statement can be made exact
in three parts, and this script measures all three.

1. **No scalar reproduces the comparison, and that is a theorem, not a fit.**
   ``a < b <=> A(a,b) > 0`` would follow from a scalar ``psi`` only if
   ``sign A = sign d psi``; a strict 3-cycle makes that impossible.  A *certified
   lower bound* on how many ordered pairs any scalar must get wrong is the size
   of a set of edge-disjoint 3-cycles: distinct edge-disjoint cycles need
   distinct reversed edges.  Greedy packing gives such a set.

2. **How much a scalar does explain.**  The least-squares potential
   ``psi_opt`` gets a measured fraction of ordered pairs right.  That is a lower
   bound on the best achievable (minimum feedback arc set is NP-hard) and it
   brackets the answer together with (1).

3. **The energy split**, which is the ``L^2`` version of the same statement.

Writes ``e3_brief_c.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import flux
import pools
from e2_potential import COARSE, FINE, statistics_of


def edge_disjoint_cycles(A: np.ndarray, tie: float = flux.TIE, limit: int = 4_000_000):
    """A greedy set of edge-disjoint strict 3-cycles."""

    n = A.shape[0]
    P = A > tie
    used = np.zeros((n, n), dtype=bool)
    chosen = []
    for i in range(n):
        for j in range(i + 1, n):
            if used[i, j] or abs(A[i, j]) <= tie:
                continue
            for k in range(j + 1, n):
                if used[i, k] or used[j, k] or abs(A[i, k]) <= tie or abs(A[j, k]) <= tie:
                    continue
                forward = P[i, j] and P[j, k] and P[k, i]
                backward = P[j, i] and P[k, j] and P[i, k]
                if forward or backward:
                    used[i, j] = used[j, i] = True
                    used[i, k] = used[k, i] = True
                    used[j, k] = used[k, j] = True
                    chosen.append((i, j, k))
                    break
    return chosen


def main() -> None:
    rows = []
    for q in (11, 13):
        S, _ = pools.arithmetic_pool(q)
        A = flux.flux_matrix(S, flux.beta_grid(q, COARSE, FINE))
        h = flux.hodge(A)
        n = len(S)
        pairs = n * (n - 1) // 2
        cycles, triangles = flux.cycle_count(A)
        packing = edge_disjoint_cycles(A)
        st = statistics_of(S, q)
        psi_end = st["half_log_phi"]
        G = psi_end[None, :] - psi_end[:, None]
        mask = np.abs(A) > flux.TIE
        phi_right = float((np.sign(A[mask]) == np.sign(G[mask])).mean())
        rec = {
            "q": q,
            "n": n,
            "unordered_pairs": pairs,
            "ties": h["ties"],
            "grad_energy_pct": 100 * h["grad_energy"],
            "curl_energy_pct": 100 * h["curl_energy"],
            "psi_opt_order_pct": 100 * h["order_agreement"],
            "phi_order_pct": 100 * phi_right,
            "cycles": cycles,
            "triangles": triangles,
            "edge_disjoint_cycles": len(packing),
            "min_wrong_pairs_lower_bound": len(packing),
            "min_wrong_pct_lower_bound": 100 * len(packing) / pairs,
        }
        rows.append(rec)
        print(f"q = {q}, n = {n}: {pairs} unordered pairs, {h['ties']} ties")
        print(f"  energy       gradient {rec['grad_energy_pct']:.3f}%   "
              f"curl {rec['curl_energy_pct']:.3f}%")
        print(f"  order        psi_opt  {rec['psi_opt_order_pct']:.3f}%   "
              f"phi = (1/2)log phi  {rec['phi_order_pct']:.3f}%")
        print(f"  cycles       {cycles} strict 3-cycles of {triangles} triangles; "
              f"{len(packing)} of them edge-disjoint")
        print(f"  so ANY scalar misorders at least {len(packing)} of the {pairs} "
              f"pairs ({rec['min_wrong_pct_lower_bound']:.4f}%), and the best one "
              f"found misorders {100 - rec['psi_opt_order_pct']:.3f}%")
    with open("e3_brief_c.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
