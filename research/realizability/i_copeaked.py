"""Session brief I -- the co-peaked-bump relaxation does NOT obstruct C_4.

Brief I proposed reducing the small-scale case to: "which 4-point metrics are
oscillations of differences of co-peaked unimodal bumps?", and expected a
linear relation forced by co-peakedness.  There is none.  The class

    W = { w : R -> [0, inf) : w(+-inf) = 0, w unimodal with peak at 0,
                              w 1-Lipschitz, w(t) <= E(t) = log(1+e^{-|t|}) }

is exactly what FINDINGS Theorem A gives for signatures with equal sigma, and
it realises s.C_4 EXACTLY for every s <= log(2)/3, by the explicit four bumps
built below.  So the relaxation is strictly weaker than the truth and cannot
prove anything about C_4.

(The same relaxation also over-shoots the sharp defect constant: within W the
supremum of osc(w_b - w_a) is 2h* = 0.96242 with h* the root of
h = log(1+e^{-h}), whereas for genuine signatures the supremum is log 2 --
Theorem 3 of OBSTRUCTION.md.  Both facts are checked here.)

    python research/realizability/i_copeaked.py
"""
from __future__ import annotations

import itertools
import math

import numpy as np

LOG2 = math.log(2.0)
C4 = np.array([[0., 1., 2., 1.],
               [1., 0., 1., 2.],
               [2., 1., 0., 1.],
               [1., 2., 1., 0.]])


def E(t):
    return np.log1p(np.exp(-np.abs(t)))


# node abscissae on the left half, and the value pattern (in units of varsigma)
NODES = np.array([-1.5, -0.9, -0.6, 0.0])
VALS = np.array([[0., 2., 2., 3.],
                 [0., 1., 2., 2.],
                 [1., 1., 2., 3.],
                 [0., 1., 1., 3.]])


def bumps(v, grid):
    """The four bumps of the construction, evaluated on ``grid``."""
    out = []
    for a in range(4):
        y = VALS[a] * v
        h = y[-1]
        w = np.empty_like(grid)
        left = grid <= NODES[0]
        w[left] = y[0] * np.log1p(np.exp(grid[left])) / math.log1p(
            math.exp(NODES[0]))
        mid = (grid > NODES[0]) & (grid <= 0.0)
        w[mid] = np.interp(grid[mid], NODES, y)
        right = grid > 0.0
        w[right] = h * np.log1p(np.exp(-grid[right])) / LOG2
        out.append(np.minimum(w, E(grid)))     # keep it under the envelope
    return out


def check(v, grid):
    ws = bumps(v, grid)
    rep = {}
    rep["max w"] = max(float(w.max()) for w in ws)
    rep["envelope slack"] = min(float((E(grid) - w).min()) for w in ws)
    lip = 0.0
    uni = 0.0
    for w in ws:
        lip = max(lip, float(np.abs(np.diff(w) / np.diff(grid)).max()))
        k = int(np.argmax(w))
        uni = max(uni, float(max(np.diff(w[:k + 1]).min(initial=0.0) * -1,
                                 np.diff(w[k:]).max(initial=0.0))))
        rep["peak at"] = float(grid[k])
    rep["Lipschitz const"] = lip
    rep["unimodality violation"] = uni
    D = np.zeros((4, 4))
    for i, j in itertools.combinations(range(4), 2):
        phi = ws[j] - ws[i]
        D[i, j] = D[j, i] = float(phi.max() - phi.min())
    iu = np.triu_indices(4, 1)
    r = D[iu] / C4[iu]
    rep["d"] = D
    rep["distortion"] = float(r.max() / r.min())
    rep["scale"] = float(r.min())
    return rep


def relaxed_eps_sup(grid):
    """sup osc(w_b - w_a) inside W, attained by an explicit pair."""
    lo, hi = 0.0, LOG2
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if m < math.log1p(math.exp(-m)):
            lo = m
        else:
            hi = m
    h = 0.5 * (lo + hi)
    tent = np.maximum(h - np.abs(grid), 0.0)
    plat = np.minimum(np.full_like(grid, h), E(grid))
    wa = np.where(grid <= 0.0, plat, tent)      # plateau left, tent right
    wb = np.where(grid <= 0.0, tent, plat)      # tent left, plateau right
    assert (wa <= E(grid) + 1e-15).all() and (wb <= E(grid) + 1e-15).all()
    phi = wb - wa
    return h, float(phi.max() - phi.min())


def main():
    grid = np.linspace(-40.0, 40.0, 800001)
    print("=== s.C_4 inside the co-peaked-bump relaxation ===")
    print(f"  {'varsigma':>10} {'scale':>12} {'distortion':>16} "
          f"{'max w':>10} {'Lip':>8} {'env slack':>11}")
    for v in (0.05, 0.10, 0.15, 0.20, LOG2 / 3):
        rep = check(v, grid)
        print(f"  {v:10.5f} {rep['scale']:12.8f} {rep['distortion']:16.10f} "
              f"{rep['max w']:10.6f} {rep['Lipschitz const']:8.4f} "
              f"{rep['envelope slack']:11.3e}", flush=True)
    rep = check(0.15, grid)
    print("\n  d matrix at varsigma = 0.15:")
    for i in range(4):
        print("   " + "  ".join(f"{rep['d'][i, j]:10.7f}" for j in range(4)))
    print(f"  unimodality violation {rep['unimodality violation']:.2e}, "
          f"peak at {rep['peak at']:.3f}")

    h, e = relaxed_eps_sup(grid)
    print("\n=== the relaxation also over-shoots the defect constant ===")
    print(f"  h* (root of h = log(1+e^-h)) = {h:.10f}")
    print(f"  osc(w_b - w_a) for the extremal relaxed pair = {e:.10f}")
    print(f"  = 2h* = {2*h:.10f}   vs the true supremum log 2 = {LOG2:.10f}")


if __name__ == "__main__":
    main()
