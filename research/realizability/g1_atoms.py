"""G1 -- how many fibers does a signature need before cycles appear?

A signature with r fibers is r real parameters, one of which is removed by the
exact reparametrisation symmetry  a -> a^p  (beta -> p beta), which fixes both
d and A.  So an r-fiber signature carries r-1 essential parameters and a
triangle carries 3r-3 of them, against 3 sign conditions.  The naive count
therefore predicts cycles from r = 2 on.  This script tests that.

r = 1 is settled by hand: u_a(s) = s + log(log a), so u_b - u_a is the
constant log(log b / log a), the flow is exact and the tournament is a total
order.  r = 2 is the first open case; it is swept densely below.

    python research/realizability/g1_atoms.py
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
from optimizers import differential_evolution, pattern_search  # noqa: E402

HERE = Path(__file__).resolve().parent


def two_fiber_grid(nx=70, nd=70, xlo=0.25, xhi=4.0):
    """A dense grid of two-fiber signatures (e^x, e^{x-delta})."""
    sigs, params = [], []
    for x in np.linspace(xlo, xhi, nx):
        for f in np.linspace(0.0, 1.0, nd):
            delta = f * x
            sigs.append(C.Sig.from_logs([x, x - delta]))
            params.append((x, delta))
    return sigs, params


def A_matrix(sigs, grid):
    tab = np.vstack([s.U(grid) for s in sigs])
    eR = np.array([math.log(s.R) for s in sigs])
    eL = np.array([math.log(s.Lam) for s in sigs])
    n = len(sigs)
    A = np.empty((n, n), dtype=np.float64)
    for i in range(n):
        diff = tab - tab[i]
        hi = np.maximum(diff.max(axis=1), np.maximum(eR - eR[i], eL - eL[i]))
        lo = np.minimum(diff.min(axis=1), np.minimum(eR - eR[i], eL - eL[i]))
        A[i] = 0.5 * (hi + lo)
    np.fill_diagonal(A, 0.0)
    return A


def count_cycles(A, tol=0.0):
    S = A > tol
    n = A.shape[0]
    total = 0
    example = None
    for i in range(n):
        out = np.flatnonzero(S[i])
        inn = np.flatnonzero(S[:, i])
        if len(out) == 0 or len(inn) == 0:
            continue
        sub = S[np.ix_(out, inn)]
        c = int(sub.sum())
        total += c
        if c and example is None:
            a, b = np.argwhere(sub)[0]
            example = (i, int(out[a]), int(inn[b]))
    return total // 3, example


def neg_cycle_margin_r(v, r):
    """-(min of the three A's around the triangle) for r-fiber signatures."""
    try:
        ss = []
        for i in range(3):
            xs = np.sort(np.clip(v[r * i:r * (i + 1)], 0.0, 6.0))[::-1]
            if xs[0] <= 1e-9:
                return 1e3
            ss.append(C.Sig.from_logs(xs))
    except ValueError:
        return 1e3
    g = C.make_grid(-14.0, 14.0, 0.005)
    As = [C.parts(a, b, g)["A"] for a, b in zip(ss, ss[1:] + ss[:1])]
    return -min(As)


def soft_cycle_r(v, r, kappa=2e-3):
    try:
        ss = []
        for i in range(3):
            xs = np.sort(np.clip(v[r * i:r * (i + 1)], 0.0, 6.0))[::-1]
            if xs[0] <= 1e-9:
                return 1e3
            ss.append(C.Sig.from_logs(xs))
    except ValueError:
        return 1e3
    g = C.make_grid(-14.0, 14.0, 0.005)
    As = [C.parts(a, b, g)["A"] for a, b in zip(ss, ss[1:] + ss[:1])]
    return -float(sum(math.tanh(t / kappa) for t in As))


def hunt_cycle(r, seeds=8, maxiter=250):
    bounds = [(0.0, 6.0)] * (3 * r)
    best = (-math.inf, None)
    for k in range(seeds):
        x, _ = differential_evolution(soft_cycle_r, bounds, args=(r,), seed=11 * k + r,
                                      maxiter=maxiter, popsize=16, F=(0.3, 1.2), CR=0.9)
        x, _ = pattern_search(soft_cycle_r, x, args=(r,), step=0.2, min_step=1e-7,
                              maxiter=8000, bounds=bounds)
        m = -neg_cycle_margin_r(x, r)
        if m > 0:
            for step in (0.1, 0.01, 1e-3):
                x, _ = pattern_search(neg_cycle_margin_r, x, args=(r,), step=step,
                                      min_step=1e-9, maxiter=8000, bounds=bounds)
            m = -neg_cycle_margin_r(x, r)
        if m > best[0]:
            best = (m, x)
    return best


def main():
    rows = []
    grid = C.make_grid(-11.0, 11.0, 0.01)

    print("=== r = 2: dense sweep of the whole two-fiber family ===")
    print("  (the family is 2-dimensional: (x, delta) with atoms e^x >= e^{x-delta};")
    print("   one dimension is pure scale, so the essential family is 1-dimensional)")
    for nx, nd, xlo, xhi in ((60, 60, 0.25, 4.0), (50, 50, 0.05, 1.2),
                             (50, 50, 1.0, 12.0)):
        sigs, params = two_fiber_grid(nx, nd, xlo, xhi)
        t0 = time.time()
        A = A_matrix(sigs, grid)
        cyc, ex = count_cycles(A)
        iu = np.triu_indices(len(sigs), 1)
        ties = int((np.abs(A[iu]) < 1e-10).sum())
        print(f"  x in [{xlo},{xhi}], {len(sigs)} signatures, "
              f"{math.comb(len(sigs),3):,} triples: {cyc} three-cycles, "
              f"{ties} ties  ({time.time()-t0:.0f}s)")
        rows.append(["r=2 sweep", f"x in [{xlo},{xhi}]", len(sigs),
                     math.comb(len(sigs), 3), cyc, ties])

    print("\n=== r = 2: is A exact on the two-fiber family? ===")
    sigs, params = two_fiber_grid(40, 40, 0.25, 4.0)
    A = A_matrix(sigs, grid)
    psi = -A.mean(axis=1)
    G = psi[None, :] - psi[:, None]
    print(f"  Hodge residual ||A - grad||/||A|| = "
          f"{np.linalg.norm(A - G)/np.linalg.norm(A):.3e}   "
          f"max |A - grad| = {np.abs(A - G).max():.3e}")
    rows.append(["r=2 hodge", "", len(sigs), "",
                 float(np.linalg.norm(A - G) / np.linalg.norm(A)),
                 float(np.abs(A - G).max())])

    print("\n=== directed-search for a 3-cycle at each fiber count r ===")
    print(f"  {'r':>3} {'best min A around triangle':>28}")
    for r in (2, 3, 4):
        t0 = time.time()
        m, x = hunt_cycle(r)
        ok = m > 0
        print(f"  {r:>3} {m:28.3e}   cycle found = {ok}   {time.time()-t0:.0f}s")
        if ok:
            xs = [np.sort(np.clip(x[r * i:r * (i + 1)], 0, 6))[::-1] for i in range(3)]
            for t in xs:
                print("        " + str([f"{math.exp(u):.5f}" for u in t]))
        rows.append(["cycle hunt", f"r={r}", m, ok, "", ""])

    with (HERE / "g1_atoms.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["what", "range", "signatures", "triples", "cycles/value", "extra"])
        for row in rows:
            wr.writerow(row)
    print(f"\nwrote {HERE/'g1_atoms.csv'}")


if __name__ == "__main__":
    main()
