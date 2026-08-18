"""G3 -- the global Hodge residual ||curl||/||A||, hill-climbed honestly.

The residual fraction is scale-free in A, so a family with A ~ 0 can report a
residual fraction of 1 while carrying no comparison at all.  Every maximum
below is therefore reported together with ||A|| and with min |A|, and only
configurations whose smallest |A| clears the 1e-10 tie threshold by four
orders of magnitude are accepted.

Search space: n signatures, r fibers each, log-atoms in [0,5]; that is n*r
real parameters against the n(n-1)/2 numbers of A, and one parameter is pure
scale.

    python research/realizability/g3_hodge.py [nmax]
"""
from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import realize as R  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

HERE = Path(__file__).resolve().parent
FLOOR = 1e-6          # required smallest |A| in an accepted configuration


def build(x, n, r):
    return R._sigs_from(x, n, r)


def stats(x, n, r):
    A = R._A_matrix(build(x, n, r))
    na = float(np.linalg.norm(A))
    psi = -A.mean(axis=1)
    G = psi[None, :] - psi[:, None]
    iu = np.triu_indices(n, 1)
    return A, na, float(np.linalg.norm(A - G) / na if na > 0 else 0.0), \
        float(np.abs(A[iu]).min())


def objective(x, n, r):
    """-residual fraction, with a hard floor on the smallest |A|."""
    try:
        A, na, ratio, mn = stats(x, n, r)
    except ValueError:
        return 1e3
    if mn < FLOOR:
        return 1e3 * (1.0 + (FLOOR - mn) / FLOOR)
    return -ratio


def soft(x, n, r):
    """A feasibility-first surrogate: get every |A| above the floor."""
    try:
        A, na, ratio, mn = stats(x, n, r)
    except ValueError:
        return 1e3
    iu = np.triu_indices(n, 1)
    feas = float(np.mean(np.tanh(np.abs(A[iu]) / FLOOR)))
    return -(feas + ratio)


def climb(n, r, seed, maxiter=250, restarts=3):
    bounds = [(0.0, 5.0)] * (n * r)
    best = (-math.inf, None)
    for k in range(restarts):
        x, _ = differential_evolution(soft, bounds, args=(n, r), seed=seed + 131 * k,
                                      maxiter=maxiter, popsize=12, F=(0.3, 1.2), CR=0.9)
        x, _ = pattern_search(soft, x, args=(n, r), step=0.2, min_step=1e-6,
                              maxiter=15000, bounds=bounds)
        if objective(x, n, r) < 1e2:
            for step in (0.15, 0.02, 2e-3, 2e-4):
                x, _ = pattern_search(objective, x, args=(n, r), step=step,
                                      min_step=1e-9, maxiter=15000, bounds=bounds)
        f = objective(x, n, r)
        if f < 1e2 and -f > best[0]:
            best = (-f, x)
    return best


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    rows = []
    print("=== hill-climbed max ||curl||/||A||  (accepted only if min|A| > 1e-6) ===")
    print(f"  {'n':>3} {'r':>3} {'ratio':>10} {'||A||':>11} {'min|A|':>11} "
          f"{'3-cycles':>12} {'psi spread':>11}")
    for n in (3, 4, 5, 6, 8, 10, 12, 16, 24):
        if n > nmax:
            break
        r = max(3, min(6, n))
        t0 = time.time()
        v, x = climb(n, r, seed=500 + n, maxiter=200 if n <= 8 else 120,
                     restarts=3 if n <= 8 else 2)
        if x is None:
            print(f"  {n:>3} {r:>3}  no admissible configuration found")
            rows.append([n, r, "", "", "", "", ""])
            continue
        A, na, ratio, mn = stats(x, n, r)
        sigs = build(x, n, r)
        cyc = len(C.three_cycles(A))
        spread = max(s.psi for s in sigs) - min(s.psi for s in sigs)
        print(f"  {n:>3} {r:>3} {ratio:10.6f} {na:11.3e} {mn:11.3e} "
              f"{cyc:>5}/{math.comb(n,3):<6} {spread:11.3e}   {time.time()-t0:.0f}s")
        rows.append([n, r, ratio, na, mn, cyc, spread])

    with (HERE / "g3_hodge.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["n", "r", "max_residual_fraction", "normA", "min_absA",
                     "three_cycles", "psi_spread"])
        for row in rows:
            wr.writerow(row)
    print(f"\nwrote {HERE/'g3_hodge.csv'}")


if __name__ == "__main__":
    main()
