"""G: numerical verification of the structure theorem and its corollaries.

The theorem proved in FINDINGS.md says, with R = log r, Lam = log max,
sigma = log(R/Lam), s = log beta:

  (S1)  max(R, beta*Lam) <= F(beta) <= R + beta*Lam
  (S2)  U(s) = log Lam + max(sigma, s) + w(s),  0 <= w <= log(1+e^{-|s-sigma|})
  (S3)  U' in (0,1); hence w is 1-Lipschitz in s
  (S4)  w increases on s < sigma and decreases on s > sigma: unimodal, peak
        exactly at s = sigma; w(peak) = log(log Z(tau)/log r)
  (S5)  d(a,b) = |sigma_b - sigma_a| + P + Q,   P,Q >= 0
  (S6)  A(a,b) = psi(b) - psi(a) + D,   D = (P-Q)/2,  |D| < (log 2)/2

This script checks every one of them on a random integer pool and on an
adversarial pool of two-scale signatures, and writes structure.csv.

    python research/realizability/structure.py
"""
from __future__ import annotations

import csv
import itertools
import math
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

HERE = Path(__file__).resolve().parent


def random_pool(n, seed, rmax=7, vmax=40):
    rng = random.Random(seed)
    pool = set()
    while len(pool) < n:
        k = rng.randint(2, rmax)
        t = tuple(sorted((rng.randint(1, vmax) for _ in range(k)), reverse=True))
        if t[0] > 1:
            pool.add(t)
    return [C.Sig.of(t) for t in sorted(pool)]


def adversarial_pool():
    """Two-scale and flat signatures: the extremal shapes of the theorem."""
    out = []
    for r in (2, 3, 5, 10, 30, 100, 1000, 10 ** 4, 10 ** 8):
        for M in (2, 3, 7, 40, 1000):
            if r > 1:
                out.append(C.Sig.compressed((M, 1), (1, r - 1)))   # thin
            out.append(C.Sig.compressed((M,), (r,)))               # flat
    for r in (4, 8, 16, 64, 10 ** 6):
        for M in (5, 50):
            out.append(C.Sig.compressed((M, 1), (r // 2, r // 2)))
    return out


def check_pointwise(sigs, grid):
    """S1-S4 pointwise, returns the worst violations."""
    v = {"S1lo": 0.0, "S1hi": 0.0, "wlo": 0.0, "env": 0.0,
         "Ulo": 0.0, "Uhi": 0.0, "lip": 0.0, "uni": 0.0, "peak": 0.0}
    beta = np.exp(grid)
    h = grid[1] - grid[0]
    for a in sigs:
        F = a.F(beta)
        v["S1lo"] = max(v["S1lo"], float(np.max(np.maximum(a.R, beta * a.Lam) - F)))
        v["S1hi"] = max(v["S1hi"], float(np.max(F - (a.R + beta * a.Lam))))
        w = a.w(grid)
        v["wlo"] = max(v["wlo"], float(np.max(-w)))
        v["env"] = max(v["env"], float(np.max(w - C.envelope(grid, a.sigma))))
        U = a.U(grid)
        dU = np.diff(U) / h
        v["Ulo"] = max(v["Ulo"], float(np.max(-dU)))
        v["Uhi"] = max(v["Uhi"], float(np.max(dU - 1.0)))
        dw = np.diff(w) / h
        v["lip"] = max(v["lip"], float(np.max(np.abs(dw)) - 1.0))
        # unimodality: dw >= 0 left of sigma, <= 0 right of sigma.  The one
        # interval straddling sigma is skipped: w has a slope jump of -1 there
        # (U is smooth, max(sigma,s) is not), so its difference quotient mixes
        # the two branches and is not a difference quotient of either.
        left = grid[1:] <= a.sigma
        right = grid[:-1] >= a.sigma
        if left.any():
            v["uni"] = max(v["uni"], float(np.max(-dw[left])))
        if right.any():
            v["uni"] = max(v["uni"], float(np.max(dw[right])))
        # peak value
        wpeak = float(a.w(np.array([a.sigma]))[0])
        v["peak"] = max(v["peak"], abs(wpeak - float(w.max())))
    return v


def check_pairwise(sigs, grid, limit=4000, seed=0):
    rng = random.Random(seed)
    pairs = list(itertools.combinations(range(len(sigs)), 2))
    if len(pairs) > limit:
        pairs = rng.sample(pairs, limit)
    worst = {"S5": 0.0, "Pneg": 0.0, "Qneg": 0.0, "maxD": 0.0, "maxeps": 0.0,
             "maxP": 0.0, "maxQ": 0.0}
    argD = argE = None
    for i, j in pairs:
        p = C.parts(sigs[i], sigs[j], grid)
        worst["S5"] = max(worst["S5"], abs(p["d"] - (abs(p["dsigma"]) + p["P"] + p["Q"])))
        worst["Pneg"] = max(worst["Pneg"], -p["P"])
        worst["Qneg"] = max(worst["Qneg"], -p["Q"])
        if abs(p["D"]) > worst["maxD"]:
            worst["maxD"], argD = abs(p["D"]), (sigs[i], sigs[j])
        if p["eps"] > worst["maxeps"]:
            worst["maxeps"], argE = p["eps"], (sigs[i], sigs[j])
        worst["maxP"] = max(worst["maxP"], p["P"])
        worst["maxQ"] = max(worst["maxQ"], p["Q"])
    return worst, argD, argE


def sharpness_ladder():
    """a = (r,1,...,1) with r atoms against b = flat (s,...,s): D -> log2/2."""
    rows = []
    b = C.Sig.compressed((5,), (5,))
    for r in (10, 100, 10 ** 3, 10 ** 4, 10 ** 6, 10 ** 9, 10 ** 12,
              10 ** 20, 10 ** 40, 10 ** 80):
        a = C.Sig.compressed((float(r), 1.0), (1.0, float(r) - 1.0))
        p = C.parts(a, b, C.make_grid(-30.0, 220.0, 0.002))
        rows.append((r, p["P"], p["Q"], p["D"], p["d"]))
    return rows


def main():
    grid = C.make_grid(-25.0, 25.0, 0.002)
    pools = {
        "random(2-7 atoms, 1-40)": random_pool(120, 11),
        "random(2-12 atoms, 1-400)": random_pool(120, 12, rmax=12, vmax=400),
        "adversarial two-scale/flat": adversarial_pool(),
    }
    print("=== S1-S4, pointwise (positive number = violation) ===")
    for name, pool in pools.items():
        v = check_pointwise(pool, grid)
        print(f"  {name}  (n={len(pool)})")
        print(f"    F >= max(R,bL): {v['S1lo']:+.2e}   F <= R+bL: {v['S1hi']:+.2e}")
        print(f"    w >= 0        : {v['wlo']:+.2e}   w <= envelope: {v['env']:+.2e}")
        print(f"    U' >= 0       : {v['Ulo']:+.2e}   U' <= 1     : {v['Uhi']:+.2e}")
        print(f"    |w'| <= 1     : {v['lip']:+.2e}   unimodal    : {v['uni']:+.2e}")
        print(f"    peak at sigma : {v['peak']:.2e} (|w(sigma)-max w|)")

    print("\n=== S5, S6, pairwise ===")
    rows = []
    for name, pool in pools.items():
        worst, argD, argE = check_pairwise(pool, grid)
        print(f"  {name}")
        print(f"    |d - (|dsigma| + P + Q)| <= {worst['S5']:.2e}"
              f"    P>=0: {-worst['Pneg']:+.2e}   Q>=0: {-worst['Qneg']:+.2e}")
        print(f"    max P = {worst['maxP']:.6f}   max Q = {worst['maxQ']:.6f}"
              f"   (bound log 2 = {C.LOG2:.6f})")
        print(f"    max |D| = {worst['maxD']:.6f}   (bound log2/2 = {C.LOG2/2:.6f})"
              f"   at {argD[0].atoms if argD else None} / {argD[1].atoms if argD else None}")
        print(f"    max eps = {worst['maxeps']:.6f}   (bound 2log2 = {2*C.LOG2:.6f})"
              f"   at {argE[0].atoms if argE else None} / {argE[1].atoms if argE else None}")
        rows.append({"pool": name, "n": len(pool), "S5_resid": worst["S5"],
                     "max_P": worst["maxP"], "max_Q": worst["maxQ"],
                     "max_absD": worst["maxD"], "max_eps": worst["maxeps"]})

    print("\n=== sharpness of |D| < log2/2 ===")
    print("  a = (r,1,...,1) with r atoms;  b = (5,5,5,5,5) flat")
    print(f"  {'r':>8} {'P':>12} {'Q':>12} {'D':>12} {'d':>12}")
    sharp = []
    for r, P, Q, D, d in sharpness_ladder():
        print(f"  {r:>8} {P:12.8f} {Q:12.8f} {D:12.8f} {d:12.8f}")
        sharp.append({"r": r, "P": P, "Q": Q, "D": D, "d": d})
    print(f"  limit as r -> infinity: P -> log 2 = {C.LOG2:.8f}, Q -> 0, "
          f"D -> log2/2 = {C.LOG2/2:.8f}")

    with (HERE / "structure.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["table", "key", "n", "S5_resid", "max_P", "max_Q", "max_absD", "max_eps"])
        for row in rows:
            wr.writerow(["pool", row["pool"], row["n"], f"{row['S5_resid']:.3e}",
                         f"{row['max_P']:.9f}", f"{row['max_Q']:.9f}",
                         f"{row['max_absD']:.9f}", f"{row['max_eps']:.9f}"])
        wr.writerow([])
        wr.writerow(["ladder", "r", "P", "Q", "D", "d"])
        for row in sharp:
            wr.writerow(["ladder", row["r"], f"{row['P']:.9f}", f"{row['Q']:.9f}",
                         f"{row['D']:.9f}", f"{row['d']:.9f}"])
    print(f"\nwrote {HERE / 'structure.csv'}")


if __name__ == "__main__":
    main()
