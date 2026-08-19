#!/usr/bin/env python3
"""E1 supplement.

Two jobs.

1. The arithmetic-versus-matched-control comparison at ``q = 17, 19, 23, 29``,
   where the enumeration is sampled rather than exhaustive, at a fixed pool size
   so the (already measured, and negligible) ``n``-dependence cannot confound
   the reading.

2. **What actually sets the curl fraction.**  Pool every arithmetic pool and
   every control at every ``q``, record the spread statistics of the trace
   vectors alongside the curl fraction, and regress.  If the arithmetic pools
   are ordinary residuals of that regression, the curl is a function of the
   spread and not of the arithmetic.

Writes ``e1_supplement.csv``.
"""

from __future__ import annotations

import csv
import math
import statistics

import numpy as np

import e1_hodge
import flux
import pools

N = 280
REPS = 3
QS = (11, 13, 17, 19, 23, 29)


def spread(S: np.ndarray, q: int) -> dict:
    a = q - S.astype(float)
    m2 = (a**2).sum(axis=1) / q**2
    maxN = S.max(axis=1).astype(float)
    return {
        "mean_m2": float(m2.mean()),
        "sd_m2": float(m2.std()),
        "sd_logloger": float(np.std(np.log(np.log(maxN)))),
        "mean_alpha_max": float(((maxN - q) / math.sqrt(q)).mean()),
        "sd_alpha_max": float(((maxN - q) / math.sqrt(q)).std()),
    }


def main() -> None:
    rows = []
    for q in QS:
        betas = flux.beta_grid(q, e1_hodge.COARSE, e1_hodge.FINE)
        if q <= 13:
            full, _ = pools.arithmetic_pool(q)
        else:
            full = pools.sampled_pool(q, draws=12000)
        rng = np.random.default_rng(q)
        print(f"\nq={q}: pool has {len(full)} signatures; using n={N}")
        for kind in ("arithmetic",) + e1_hodge.CONTROLS:
            vals, cycs = [], []
            for rep in range(REPS):
                if kind == "arithmetic":
                    S = full[np.sort(rng.choice(len(full), N, replace=False))]
                else:
                    S = pools.control_pool(q, N, kind, reference=full, seed=31 * rep + q)
                row = e1_hodge.stats(S, betas)
                row.update(spread(S, q))
                row.update(pool=kind, q=q, replicate=rep)
                rows.append(row)
                vals.append(row["curl_frac"])
                cycs.append(row["cycles"])
            print(f"  {kind:18s} curl={statistics.mean(vals):.6f} "
                  f"+-{statistics.pstdev(vals):.6f}   "
                  f"cycles={statistics.mean(cycs):7.1f}   "
                  f"m2={rows[-1]['mean_m2']:.3f}+-{rows[-1]['sd_m2']:.3f}")

    with open("e1_supplement.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # --- what sets the curl fraction?
    print("\n" + "=" * 74)
    print("regression of curl/|A| on the spread statistics, over all "
          f"{len(rows)} pools (arithmetic and control alike)")
    y = np.array([r["curl_frac"] for r in rows])
    names = ["mean_m2", "sd_m2", "sd_alpha_max", "sd_logloger"]
    X = np.column_stack([np.ones(len(y))] + [np.array([r[k] for r in rows]) for k in names])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    R2 = 1 - resid.var() / y.var()
    sigma = resid.std()
    print(f"  R2 = {R2:.4f}, residual sd = {sigma:.5f}")
    for name, coefficient in zip(["const"] + names, beta):
        print(f"    {name:14s} {coefficient:+.6f}")
    print("\n  where the arithmetic pools sit, in residual standard deviations:")
    for q in QS:
        z = [resid[i] / sigma for i, r in enumerate(rows)
             if r["pool"] == "arithmetic" and r["q"] == q]
        print(f"    q={q:3d}   " + "  ".join(f"{v:+.2f}" for v in z))
    others = [resid[i] / sigma for i, r in enumerate(rows) if r["pool"] != "arithmetic"]
    print(f"    controls: mean {statistics.mean(others):+.2f}, "
          f"sd {statistics.pstdev(others):.2f}, "
          f"range {min(others):+.2f}..{max(others):+.2f}")


if __name__ == "__main__":
    main()
