"""G3 -- the rigidity inequality for the curl fraction.

Brief E measured that ||curl||/||A|| is a function of the spread of the pool
(R^2 = 0.93) and is flat in n from 8 to 698.  The structure theorem explains
that and turns it into an inequality.  Since A = d(psi) + D and the Hodge
gradient projection annihilates no gradient, the residual of A is the residual
of D, so ||curl|| <= ||D||_F, and D is bounded three ways:

  (B1) computable, exact:   |D(a,b)| <= (d(a,b) - |sigma_a - sigma_b|)/2
  (B2) a priori from sigma: |D(a,b)| <= (1/2) log(1 + e^{-|sigma_a-sigma_b|})
  (B3) universal:           |D(a,b)| <= (log 2)/2

Each gives  ||curl||/||A|| <= ||bound||_F / ||A||_F,  with no n anywhere: the
bounds are per-edge, which is exactly why the measured statistic is flat in n.
B3 in the form  ||curl||/||A|| <= (log2)/2 / rms(A)  is the promised
"function of the spread": rms(A) is dominated by the spread of psi.

    python research/realizability/g3_bound.py
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "flux_arithmetic"))
import common as C  # noqa: E402
import pools as P  # noqa: E402

HERE = Path(__file__).resolve().parent
GRID = C.make_grid(-20.0, 20.0, 0.004)


def matrices(sigs, grid=GRID):
    tab = np.vstack([s.U(grid) for s in sigs])
    eR = np.array([math.log(s.R) for s in sigs])
    eL = np.array([math.log(s.Lam) for s in sigs])
    n = len(sigs)
    A = np.zeros((n, n))
    d = np.zeros((n, n))
    for i in range(n):
        diff = tab - tab[i]
        hi = np.maximum(diff.max(axis=1), np.maximum(eR - eR[i], eL - eL[i]))
        lo = np.minimum(diff.min(axis=1), np.minimum(eR - eR[i], eL - eL[i]))
        A[i] = 0.5 * (hi + lo)
        d[i] = hi - lo
    np.fill_diagonal(A, 0.0)
    np.fill_diagonal(d, 0.0)
    return A, d


def report(name, sigs, rows):
    n = len(sigs)
    A, d = matrices(sigs)
    sg = np.array([s.sigma for s in sigs])
    ps = np.array([s.psi for s in sigs])
    ell = np.abs(sg[None, :] - sg[:, None])
    psi = -A.mean(axis=1)
    G = psi[None, :] - psi[:, None]
    nA = np.linalg.norm(A)
    ratio = np.linalg.norm(A - G) / nA
    off = ~np.eye(n, dtype=bool)
    rms = nA / math.sqrt(n * (n - 1))

    D_true = A - (ps[None, :] - ps[:, None])
    b1 = 0.5 * (d - ell)
    b2 = 0.5 * np.log1p(np.exp(-ell))
    np.fill_diagonal(b2, 0.0)
    b3 = np.where(off, C.LOG2 / 2, 0.0)
    viol1 = float(np.max(np.abs(D_true)[off] - b1[off]))
    viol2 = float(np.max(np.abs(D_true)[off] - b2[off]))

    B1 = np.linalg.norm(b1) / nA
    B2 = np.linalg.norm(b2) / nA
    B3 = np.linalg.norm(b3) / nA
    print(f"  {name:<26} n={n:<4} ratio={ratio:.4f}  rms|A|={rms:.4f}  "
          f"sd(psi)={ps.std():.4f}  sd(sigma)={sg.std():.4f}")
    print(f"      bounds  B1={B1:8.4f}  B2={B2:8.4f}  B3={B3:8.4f}   "
          f"(B3 = (log2/2)/rms|A| = {C.LOG2/2/rms:.4f})")
    print(f"      |D| <= B1 violation {viol1:+.2e}   |D| <= B2 violation {viol2:+.2e}")
    rows.append([name, n, ratio, rms, ps.std(), sg.std(), B1, B2, B3, viol1, viol2])
    return ratio, B1


def main():
    rows = []
    rng = np.random.default_rng(11)

    print("=== short random integer signatures (2-7 entries, values 1-40) ===")
    sys.path.insert(0, str(HERE))
    from gpools import integer_pool
    big = integer_pool(400, seed=11)
    for n in (8, 16, 24, 48, 96):
        sel = [big[i] for i in rng.choice(len(big), n, replace=False)]
        report(f"random short, n={n}", sel, rows)

    print("\n=== q-entry near-flat signatures (brief E shape) ===")
    for q in (11, 13, 17):
        pool = P.arithmetic_pool(q)
        arr = pool[0] if isinstance(pool, tuple) else pool
        arr = np.asarray(arr)
        idx = rng.choice(len(arr), min(60, len(arr)), replace=False)
        sigs = [C.Sig.of(tuple(int(v) for v in arr[i] if v > 0)) for i in idx]
        report(f"arithmetic q={q}", sigs, rows)
        ctrl = P.control_pool(q, 60, "marginal", reference=arr, seed=7)
        carr = np.asarray(ctrl[0] if isinstance(ctrl, tuple) else ctrl)
        sigs = [C.Sig.of(tuple(int(v) for v in row if v > 0)) for row in carr]
        report(f"control q={q}", sigs, rows)

    print("\n=== focused pools (r and M fixed: sigma and psi exactly tied) ===")
    from g1_mine import bucket
    for (r, M) in ((6, 12), (5, 12), (7, 12)):
        sigs = bucket(r, M, limit=60, seed=1)
        report(f"bucket r={r},M={M}", sigs, rows)

    with (HERE / "g3_bound.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["pool", "n", "curl_fraction", "rms_absA", "sd_psi",
                     "sd_sigma", "bound_B1", "bound_B2", "bound_B3",
                     "B1_violation", "B2_violation"])
        for row in rows:
            wr.writerow([row[0], row[1]] + [f"{v:.6g}" for v in row[2:]])
    print(f"\nwrote {HERE/'g3_bound.csv'}")


if __name__ == "__main__":
    main()
