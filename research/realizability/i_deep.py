"""Session brief I -- deep search for C_4 in the exact finite (node) model.

Basin hopping on top of i_nodes: DE, long compass polish, then repeated
perturb-and-repolish from the incumbent.  Seeded, among others, from the
scale functions S(theta) = log((F - beta F')/F') of brief G's own best C_4
signature family, so the search provably starts no worse than FINDINGS Sec 4.3.

    python research/realizability/i_deep.py --target C_4 --ms 4 5 6 8
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import i_nodes as N  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

STEPS = (0.6, 0.25, 0.1, 0.04, 0.015, 6e-3, 2.5e-3, 1e-3, 4e-4, 1.5e-4, 6e-5,
         2.5e-5, 1e-5, 4e-6, 1.5e-6, 6e-7, 2.5e-7, 1e-7, 4e-8, 1.5e-8, 6e-9)

# brief G's best C_4 signature family (r = 6), re-derived in i_check_g2.py
G_BEST = ([(5.00000, 3.46833, 1.83864, 1.46401, 0.00000),
           (5.00000, 3.95612, 3.24713, 2.23097, 0.00000),
           (4.35282, 4.19924, 3.42504, 2.86000, 2.54824, 1.94151),
           (4.31133, 3.08399, 1.64404, 1.38879, 1.01713, 0.38547)],
          [(1, 1, 1, 1, 2), (2, 1, 1, 1, 1), (1,) * 6, (1,) * 6])


def scale_function(s: C.Sig, th):
    b = np.exp(np.atleast_1d(np.asarray(th, float)))
    x, m = s.x, s.m
    top = b[:, None] * x[None, :]
    w = np.exp(top - top.max(axis=1, keepdims=True)) * m[None, :]
    Z = w.sum(axis=1)
    F = np.log(Z) + top.max(axis=1)
    Fp = (w * x[None, :]).sum(axis=1) / Z
    return np.log((F - b * Fp) / Fp)


def seed_from_signatures(m, lo=-3.5, hi=3.0):
    sigs = [C.Sig.compressed([math.exp(v) for v in xx], mm)
            for xx, mm in zip(*G_BEST)]
    th = np.linspace(lo, hi, m)
    mids = np.concatenate([[th[0] - 1.0], (th[:-1] + th[1:]) / 2, [th[-1] + 1.0]])
    S = np.vstack([scale_function(s, mids) for s in sigs])
    S = np.minimum.accumulate(S, axis=1)
    return pack(th, S)


def pack(th, S):
    m = len(th)
    n = S.shape[0]
    z = [th[0]] + [math.log(max(th[k] - th[k - 1], 1e-12)) for k in range(1, m)]
    for a in range(n):
        z.append(S[a, 0])
        for k in range(1, m + 1):
            z.append(math.log(max(S[a, k - 1] - S[a, k], 1e-12)))
    return np.array(z)


def polish(z, n, m, delta, iu, b):
    f = N.distortion(z, n, m, delta, iu)
    for step in STEPS:
        z, f = pattern_search(N.distortion, z, args=(n, m, delta, iu), step=step,
                              min_step=1e-16, maxiter=500000, bounds=b)
    return z, f


def deep(delta, m, seed, de_restarts=6, hops=40, maxiter=700, use_seed=True):
    n = delta.shape[0]
    iu = np.triu_indices(n, 1)
    b = [(-6.0, 6.0)] + [(-8.0, 3.0)] * (m - 1)
    for _ in range(n):
        b += [(-6.0, 8.0)] + [(-10.0, 3.0)] * m
    rng = np.random.default_rng(seed)
    best_z, best = None, math.inf
    if use_seed and n == 4:
        z0 = np.clip(seed_from_signatures(m), [q[0] for q in b], [q[1] for q in b])
        best_z, best = polish(z0, n, m, delta, iu, b)
    for t in range(de_restarts):
        z, _ = differential_evolution(N.distortion, b, args=(n, m, delta, iu),
                                      seed=seed + 7717 * t, maxiter=maxiter,
                                      popsize=12, F=(0.3, 1.2), CR=0.9)
        z, f = polish(z, n, m, delta, iu, b)
        if f < best:
            best_z, best = z, f
    for h in range(hops):
        amp = 0.6 * (0.85 ** (h % 12))
        z = best_z + rng.normal(scale=amp, size=len(best_z))
        z = np.clip(z, [q[0] for q in b], [q[1] for q in b])
        z, f = polish(z, n, m, delta, iu, b)
        if f < best - 1e-12:
            best_z, best = z, f
    return best_z, best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ms", type=int, nargs="*", default=[4, 5, 6, 8])
    ap.add_argument("--target", default="C_4")
    ap.add_argument("--hops", type=int, default=40)
    ap.add_argument("--de", type=int, default=6)
    args = ap.parse_args()
    delta = N.TARGETS[args.target]
    n = delta.shape[0]
    best = (math.inf, None, None)
    for m in args.ms:
        z, f = deep(delta, m, seed=6101 + 71 * m, de_restarts=args.de,
                    hops=args.hops)
        print(f"  {args.target}  m={m}  distortion = {f:.12f}", flush=True)
        if f < best[0]:
            best = (f, z, m)
    f, z, m = best
    N.report(z, n, m, delta, args.target)
    np.save(Path(__file__).resolve().parent /
            f"i_deep_best_{args.target.replace(',', '')}.npy",
            np.concatenate([[m], z]))
    print(f"\n  BEST over m: {f:.12f}")


if __name__ == "__main__":
    main()
