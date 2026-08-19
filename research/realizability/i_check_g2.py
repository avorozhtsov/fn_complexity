"""Session brief I -- audit of FINDINGS Sec. 4.3's C_4 distortion 1.255692.

g2_metrics.distortion evaluates d on a fixed grid of step 0.01 in s = log beta
and takes the RAW grid max/min, with no parabolic refinement (unlike
common.matrices).  That always UNDER-estimates d, by a different amount for
each pair, so the reported distortion is not an upper bound for the
configuration it returns.  Here the same search is re-run and the returned
configuration is re-scored with the certified Lipschitz-bracket extrema of
common.certified_extrema.

    python research/realizability/i_check_g2.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import g2_metrics as G  # noqa: E402
import realize as R  # noqa: E402


def certified_distortion(sigs, delta):
    n = len(sigs)
    iu = np.triu_indices(n, 1)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mx, mn, _, _ = C.certified_extrema(sigs[i], sigs[j], step=2e-4)
            D[i, j] = D[j, i] = mx - mn
    r = D[iu] / delta[iu]
    return float(r.max() / r.min()), D


def main():
    delta = G.targets()["C_4"]
    n = 4
    print(f"  {'r':>3} {'grid distortion':>18} {'certified distortion':>22} "
          f"{'bias':>10}")
    best = (math.inf, None, None)
    for r in (3, 4, 6):
        x, f = G.realise_metric(delta, r, seed=17 + n, maxiter=250)
        sigs = R._sigs_from(x, n, r)
        cf, D = certified_distortion(sigs, delta)
        print(f"  {r:>3} {f:18.6f} {cf:22.6f} {cf-f:10.2e}", flush=True)
        if cf < best[0]:
            best = (cf, sigs, D)
    cf, sigs, D = best
    print(f"\n  best certified distortion over r = 3,4,6 : {cf:.9f}")
    print("  d matrix:")
    for i in range(n):
        print("   " + "  ".join(f"{D[i, j]:10.7f}" for j in range(n)))
    for i, s in enumerate(sigs):
        print(f"   sig {i}: log-atoms = "
              + ", ".join(f"{v:.5f}" for v in s.xs)
              + f"  mults = {s.mults}  sigma = {s.sigma:+.5f}")


if __name__ == "__main__":
    main()
