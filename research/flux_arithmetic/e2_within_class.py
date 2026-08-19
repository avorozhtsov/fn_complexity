#!/usr/bin/env python3
"""E2, within a `phi`-class -- where all the content of the flux actually is.

A global regression of ``psi_opt`` on ``(1/2) log phi`` scores ``R^2 = 0.99`` and
says almost nothing, because at a fixed ``q`` every pool member has ``q`` fibers,
so ``phi`` is a monotone function of the single integer ``max_c N_c`` and takes
only about ``2 g sqrt(q)`` values.  Between classes ``phi`` decides; *inside* a
class it is exactly tied, and the curve-family session's census shows that inside
a class not one rate is attained at an endpoint.  So the within-class flux is the
whole of what a scalar has to explain and cannot.

This script restricts the flux to each class of equal largest fiber, redoes the
Hodge split there, and asks what the within-class potential is.  It also settles
the sign of the ``m2`` correction, which brief B's addendum got backwards.

Writes ``e2_within_class.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import flux
import pools
from e2_potential import COARSE, FINE, order_agreement, r2, statistics_of

MODELS = [
    ("m2", ["m2"]),
    ("m3", ["m3"]),
    ("log mult", ["log_mult"]),
    ("log minN", ["log_minN"]),
    ("m2, m3", ["m2", "m3"]),
    ("m2, m3, m4", ["m2", "m3", "m4"]),
    ("m2, m3, m4, log mult", ["m2", "m3", "m4", "log_mult"]),
    ("m2, m3, m4, log mult, log minN", ["m2", "m3", "m4", "log_mult", "log_minN"]),
]


def main() -> None:
    rows = []
    for q in (11, 13, 17, 19):
        if q <= 13:
            S, _ = pools.arithmetic_pool(q)
            tag = "exhaustive"
        else:
            S = pools.sampled_pool(q, draws=12000)
            if len(S) > 1200:  # the flux matrix costs O(n^2 * grid)
                rng = np.random.default_rng(q)
                S = S[np.sort(rng.choice(len(S), 1200, replace=False))]
            tag = "sampled"
        A = flux.flux_matrix(S, flux.beta_grid(q, COARSE, FINE))
        st = statistics_of(S, q)
        maxN = st["maxN"]
        print(f"\n=== q={q} ({tag}), n={len(S)}, "
              f"{len(np.unique(maxN))} classes of equal largest fiber ===")

        # aggregate over classes: stack the within-class problems
        agg = {name: [0.0, 0.0] for name, _ in MODELS}
        tot_pairs = tot_agree_m2 = 0
        curl_num = curl_den = 0.0
        psi_all, key_all = [], []
        used = 0
        for value in np.unique(maxN):
            idx = np.flatnonzero(maxN == value)
            if len(idx) < 3:
                continue
            used += 1
            B = A[np.ix_(idx, idx)]
            h = flux.hodge(B)
            curl_num += np.linalg.norm(B - (h["psi"][None, :] - h["psi"][:, None])) ** 2
            curl_den += np.linalg.norm(B) ** 2
            sub = {k: v[idx] for k, v in st.items()}
            psi = h["psi"]
            psi_all.append(psi)
            key_all.append(np.full(len(idx), value))
            for name, cols in MODELS:
                X = np.column_stack([np.ones(len(idx))] + [sub[c] for c in cols])
                beta, *_ = np.linalg.lstsq(X, psi, rcond=None)
                res = psi - X @ beta
                agg[name][0] += float((res**2).sum())
                agg[name][1] += float(((psi - psi.mean()) ** 2).sum())
            # the addendum's rule, restricted to where it acts
            mask = np.abs(B) > flux.TIE
            G = sub["m2"][None, :] - sub["m2"][:, None]
            live = mask & (np.abs(G) > 0)
            tot_pairs += int(live.sum())
            # "the larger m2 precedes" means psi decreasing in m2
            tot_agree_m2 += int((np.sign(B[live]) == -np.sign(G[live])).sum())

        curl = math.sqrt(curl_num / curl_den)
        print(f"  {used} classes with 3 or more members")
        print(f"  within-class curl/|A| = {curl:.6f}   "
              f"(whole pool: {flux.hodge(A)['curl_frac']:.6f})")
        print(f"  the addendum's rule 'the larger m2 precedes' is right on "
              f"{tot_agree_m2/tot_pairs*100:.1f}% of the {tot_pairs} ordered "
              f"within-class pairs it decides")
        print("  regressions of the WITHIN-CLASS potential")
        for name, _ in MODELS:
            value = 1 - agg[name][0] / agg[name][1]
            rows.append({"q": q, "pool": tag, "scope": "within-class",
                         "model": name, "R2": value})
            print(f"    R2 = {value:.6f}   {name}")
        rows.append({"q": q, "pool": tag, "scope": "within-class",
                     "model": "[curl/|A| within class]", "R2": curl})
        rows.append({"q": q, "pool": tag, "scope": "within-class",
                     "model": "[larger m2 precedes: fraction right]",
                     "R2": tot_agree_m2 / tot_pairs})

        # the sign of the m2 coefficient, pooled over classes
        num = den = 0.0
        for psi, key in zip(psi_all, key_all):
            idx = np.flatnonzero(maxN == key[0])
            x = st["m2"][idx] - st["m2"][idx].mean()
            num += float(x @ (psi - psi.mean()))
            den += float(x @ x)
        print(f"  pooled within-class slope of psi on m2: {num/den:+.6f} "
              f"(the addendum requires it to be negative)")
        rows.append({"q": q, "pool": tag, "scope": "within-class",
                     "model": "[pooled slope of psi on m2]", "R2": num / den})

    with open("e2_within_class.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
