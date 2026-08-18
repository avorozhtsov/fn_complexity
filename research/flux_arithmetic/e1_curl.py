#!/usr/bin/env python3
"""The triangle statistic ``r = |curl A| / sum|A|`` on the arithmetic pools and
on the matched controls, for reconciliation with
``research/m_and_e_and_a_c/curl_on_curve_families.py``.

That script introduced ``r`` as an expansion-free search statistic and is right
to: for an antisymmetric ``A`` on a *complete* graph, all triangle sums vanish
iff ``A = d psi``, so every scalar invariant contributes exactly zero to
``curl A`` at every order.  (Proof: if every triangle sum vanishes, fix a base
point ``o`` and set ``psi(x) = A(o, x)``; the triangle ``(o, x, y)`` gives
``A(x, y) = psi(y) - psi(x)``.  The converse is immediate.)

It sampled ~90 signatures per field at ``q = 101 ... 1009``, where the pool is a
sample.  This script runs the same statistic on the *exhaustive* pools at
``q = 11, 13`` -- where 132 and 1475 strict 3-cycles are certified, so
``max r = 1`` is known in advance -- and on the matched controls, so the two
measurements can be compared.

Writes ``e1_curl.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import e1_hodge
import flux
import pools


def curl_stats(A: np.ndarray) -> dict:
    """max/median ``|curl A|``, max ``r``, and the tail of the ``r`` histogram."""

    n = A.shape[0]
    absA = np.abs(A)
    max_curl = 0.0
    max_r = 0.0
    hi = 0
    curl_sum = 0.0
    curl_sq = 0.0
    count = 0
    bins = np.zeros(11, dtype=np.int64)
    for i in range(n - 2):
        j = np.arange(i + 1, n)
        # curl(i, j, k) = A[i,j] + A[j,k] + A[k,i] = A[i,j] - A[i,k] + A[j,k]
        block = A[i, j][:, None] - A[i, j][None, :] + A[np.ix_(j, j)]
        denom = absA[i, j][:, None] + absA[i, j][None, :] + absA[np.ix_(j, j)]
        mask = np.triu(np.ones_like(block, dtype=bool), 1)
        c = np.abs(block[mask])
        d = denom[mask]
        ratio = np.where(d > 0, c / np.maximum(d, 1e-300), 0.0)
        max_curl = max(max_curl, float(c.max(initial=0.0)))
        max_r = max(max_r, float(ratio.max(initial=0.0)))
        hi += int((ratio > 0.9).sum())
        curl_sum += float(c.sum())
        curl_sq += float((c**2).sum())
        count += int(c.size)
        bins += np.bincount(np.minimum((ratio * 10).astype(int), 10), minlength=11)
    return {
        "triangles": count,
        "max_curl": max_curl,
        "rms_curl": math.sqrt(curl_sq / count) if count else float("nan"),
        "mean_curl": curl_sum / count if count else float("nan"),
        "max_r": max_r,
        "r_above_0.9": hi,
        "r_hist": bins.tolist(),
    }


def main() -> None:
    rows = []
    for q in (11, 13):
        betas = flux.beta_grid(q, e1_hodge.COARSE, e1_hodge.FINE)
        S, _ = pools.arithmetic_pool(q)
        entries = [("arithmetic", S)]
        for kind in e1_hodge.CONTROLS:
            entries.append((f"control:{kind}",
                            pools.control_pool(q, len(S), kind, reference=S, seed=q)))
        for name, pool in entries:
            A = flux.flux_matrix(pool, betas)
            cs = curl_stats(A)
            off = ~np.eye(len(A), dtype=bool)
            cs.update(q=q, pool=name, n=len(pool),
                      median_absA=float(np.median(np.abs(A[off]))),
                      cycles=flux.cycle_count(A)[0])
            rows.append(cs)
            print(f"q={q:3d} {name:20s} n={cs['n']:4d}  max r={cs['max_r']:.6f}  "
                  f"max|curl|={cs['max_curl']:.3e}  rms|curl|={cs['rms_curl']:.3e}  "
                  f"median|A|={cs['median_absA']:.3e}  cycles={cs['cycles']}")
    with open("e1_curl.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
