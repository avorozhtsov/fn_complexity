#!/usr/bin/env python3
"""E1 -- the Hodge decomposition of the flux on arithmetic pools, against
properly matched random controls.

Three blocks, all written to CSV.

``e1_smoke.csv``     the calibration run of ``research/realizability/
                     tournament_seed.py`` reproduced through this pipeline,
                     plus its seed-to-seed spread.
``e1_pools.csv``     the full pools: arithmetic ``F_11`` (296) and ``F_13``
                     (698), and matched random controls at the same ``q`` and
                     the same ``n``.
``e1_scaling.csv``   the same statistics as a function of ``n``, because the
                     curl fraction of a pool of size ``n`` is not comparable
                     with that of a pool of size ``m``.

Run from ``research/flux_arithmetic``.
"""

from __future__ import annotations

import csv
import math
import random
import statistics

import numpy as np

import flux
import pools

COARSE, FINE = 2400, 13600  # grid convergence: |dA| < 5e-9, see FINDINGS


def stats(signatures, betas, want_cycles=True) -> dict:
    A = flux.flux_matrix(signatures, betas)
    h = flux.hodge(A)
    n = len(A)
    off = ~np.eye(n, dtype=bool)
    out = {
        "n": n,
        "grad_frac": h["grad_frac"],
        "curl_frac": h["curl_frac"],
        "curl_energy": h["curl_energy"],
        "order_agreement": h["order_agreement"],
        "ties": h["ties"],
        "min_abs_A": float(np.abs(A[off]).min()),
        "sd_psi": float(np.std(h["psi"])),
        "rms_A": float(np.sqrt((A[off] ** 2).mean())),
    }
    if want_cycles:
        c, t = flux.cycle_count(A)
        out["cycles"] = c
        out["triangles"] = t
        out["cycle_frac"] = c / t if t else float("nan")
    return out


# --------------------------------------------------------------- smoke test


def smoke() -> list[dict]:
    """The random-integer calibration of brief E, and its variance."""

    rows = []
    betas = np.geomspace(1e-4, 4000.0, 20000)
    for seed in range(11, 21):
        rnd = random.Random(seed)

        def draw():
            k = rnd.randint(2, 7)
            return tuple(sorted((rnd.randint(1, 40) for _ in range(k)), reverse=True))

        pool = sorted({draw() for _ in range(300)})
        for n in (8, 16, 24):
            sub = rnd.sample(pool, n)
            row = stats(sub, betas)
            row["seed"] = seed
            row["pool"] = "random-integer (tournament_seed)"
            rows.append(row)
    return rows


# ------------------------------------------------------------- the controls

CONTROLS = ("loose", "m2matched", "marginal", "sigshuffle", "maxmatched")


def replicates(n: int) -> int:
    return 8 if n <= 100 else 4 if n <= 300 else 2


def main() -> None:
    smoke_rows = smoke()
    with open("e1_smoke.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(smoke_rows[0]))
        w.writeheader()
        w.writerows(smoke_rows)
    print("smoke test -- random integer signatures, 10 seeds")
    for n in (8, 16, 24):
        sel = [r for r in smoke_rows if r["n"] == n]
        c = [r["curl_frac"] for r in sel]
        cyc = [r["cycles"] for r in sel]
        print(
            f"  n={n:3d}  curl/|A| = {min(c):.4f}..{max(c):.4f} "
            f"(median {statistics.median(c):.4f})   3-cycles {min(cyc)}..{max(cyc)}"
        )

    pool_rows: list[dict] = []
    scaling_rows: list[dict] = []
    for q in (11, 13):
        betas = flux.beta_grid(q, COARSE, FINE)
        S, _ = pools.arithmetic_pool(q)
        n_full = len(S)
        row = stats(S, betas)
        row.update(pool="arithmetic", q=q, replicate=0)
        pool_rows.append(row)
        print(
            f"\nq={q} arithmetic  n={row['n']}  grad={row['grad_frac']:.6f} "
            f"curl={row['curl_frac']:.6f}  cycles={row['cycles']}  "
            f"agree={row['order_agreement']:.5f}"
        )

        for kind in CONTROLS:
            vals, cycs = [], []
            for rep in range(5):
                C = pools.control_pool(q, n_full, kind, reference=S, seed=1000 * rep + q)
                if len(C) < n_full:
                    print(f"  (control {kind}: only {len(C)} distinct)")
                row = stats(C, betas)
                row.update(pool=f"control:{kind}", q=q, replicate=rep)
                pool_rows.append(row)
                vals.append(row["curl_frac"])
                cycs.append(row["cycles"])
            print(
                f"  control {kind:11s} curl={statistics.mean(vals):.6f} "
                f"+-{statistics.pstdev(vals):.6f}   cycles={statistics.mean(cycs):.1f}"
            )

        # --- n-scaling, arithmetic subsamples against the same controls
        rng = np.random.default_rng(q)
        sizes = [k for k in (8, 16, 24, 50, 100, 200, 296, 500, 698) if k <= n_full]
        for n in sizes:
            for rep in range(replicates(n)):
                idx = rng.choice(n_full, n, replace=False)
                row = stats(S[idx], betas)
                row.update(pool="arithmetic", q=q, replicate=rep)
                scaling_rows.append(row)
            for kind in CONTROLS:
                for rep in range(replicates(n)):
                    C = pools.control_pool(q, n, kind, reference=S, seed=97 * rep + 7 * n + q)
                    row = stats(C, betas)
                    row.update(pool=f"control:{kind}", q=q, replicate=rep)
                    scaling_rows.append(row)

    with open("e1_pools.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pool_rows[0]))
        w.writeheader()
        w.writerows(pool_rows)
    with open("e1_scaling.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(scaling_rows[0]))
        w.writeheader()
        w.writerows(scaling_rows)

    print("\nn-scaling of curl/|A|  (mean over 8 replicates)")
    for q in (11, 13):
        print(f"  q={q}")
        sizes = sorted({r["n"] for r in scaling_rows if r["q"] == q})
        header = "    n     " + "".join(f"{k:>14s}" for k in ["arithmetic"] + list(CONTROLS))
        print(header)
        for n in sizes:
            cells = []
            for key in ["arithmetic"] + [f"control:{k}" for k in CONTROLS]:
                sel = [r["curl_frac"] for r in scaling_rows
                       if r["q"] == q and r["n"] == n and r["pool"] == key]
                cells.append(f"{statistics.mean(sel):14.4f}" if sel else " " * 14)
            print(f"  {n:5d}" + "".join(cells))


if __name__ == "__main__":
    main()
