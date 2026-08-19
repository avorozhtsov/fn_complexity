"""Session brief I -- numerical verification of the four theorems.

T1  U is strictly convex in s = log beta:   U'' = nu(1-nu) + beta^2 F''/F > 0.
T2  the scale function  T(beta) = (F - beta F')/F'  is positive and strictly
    decreasing, and  U'(s) = beta/(beta + T);  hence
    sign (U_b - U_a)'  =  sign (T_a - T_b).
T3  the sharp metric bound  d(a,b) <= |dsigma| + log(1 + e^{-|dsigma|})
    (half of FINDINGS Corollary A2).
T4  the sharp scale bound for C_4:  any realisation needs  s <= log(golden).

    python research/realizability/i_verify.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import gpools  # noqa: E402

PHI = (1 + math.sqrt(5)) / 2


def pieces(s: C.Sig, th):
    """F, F', F'', and the CANCELLATION-FREE  G = F - beta F' > 0.

    With Lam the top log-atom and y_i = Lam - x_i >= 0,
        F     = beta Lam + log S,      S = sum m_i e^{-beta y_i}
        F'    = Lam - (sum m_i y_i e^{-beta y_i}) / S
        G     = F - beta F' = log S + beta (sum m_i y_i e^{-beta y_i}) / S
    and both summands of G are >= 0, so no cancellation occurs.
    """
    b = np.exp(np.asarray(th, float))
    x, m = s.x, s.m
    Lam = float(x[0])
    y = Lam - x
    w = np.exp(-np.multiply.outer(b, y)) * m[None, :]
    S = w.sum(axis=1)
    F = b * Lam + np.log(S)
    p = w / S[:, None]
    my = (p * y[None, :]).sum(axis=1)
    Fp = (p * x[None, :]).sum(axis=1)            # cancellation-free (x >= 0)
    Fpp = (p * (x[None, :] - Fp[:, None]) ** 2).sum(axis=1)
    G = np.log(S) + b * my                       # = F - beta F',  both terms >=0
    return b, F, Fp, Fpp, G


def main():
    pools = (gpools.integer_pool(70, seed=11)
             + gpools.integer_pool(70, seed=12, rmax=12, vmax=400)
             + gpools.integer_pool(70, seed=13, rmax=20, vmax=10 ** 5)
             + [C.Sig.compressed([10 ** 6, 1], [1, 10 ** 8]),
                C.Sig.compressed([2, 1], [1, 10 ** 12]),
                C.Sig.of((7,) * 40),
                C.Sig.of((10 ** 5, 3, 3, 2, 2, 2, 1, 1))])
    th = np.linspace(-25.0, 25.0, 4001)

    print("=== T1  U'' > 0 in s = log beta  (proved) ===")
    worst = math.inf
    for s in pools:
        b, F, Fp, Fpp, G = pieces(s, th)
        nu = b * Fp / F
        d2 = nu * (G / F) + b ** 2 * Fpp / F
        worst = min(worst, float(d2.min()))
    print(f"  min over {len(pools)} signatures x {len(th)} points of "
          f"U''  =  {worst:.6e}   (>= 0; zero only by underflow at |s| = 25)")
    tm = np.linspace(-6.0, 6.0, 2001)
    w2 = math.inf
    for sg in pools:
        b, F, Fp, Fpp, G = pieces(sg, tm)
        w2 = min(w2, float((b * Fp / F * (G / F) + b ** 2 * Fpp / F).min()))
    print(f"  min over the same signatures on |s| <= 6 : {w2:.6e}")

    print("\n=== T2  scale function T = (F - beta F')/F' decreasing; "
          "U' = beta/(beta+T) ===")
    worst_dec, worst_id = 0.0, 0.0
    for s in pools:
        b, F, Fp, Fpp, G = pieces(s, th)
        T = G / Fp
        den = np.maximum(T[:-1], 1e-300)
        worst_dec = max(worst_dec, float(np.max(np.diff(T) / den)))
        worst_id = max(worst_id, float(np.max(np.abs(b * Fp / F - b / (b + T)))))
    print(f"  max RELATIVE increase of T along the grid : {worst_dec:.3e} "
          f" (<= 0)")
    print(f"  max |U' - beta/(beta+T)|          : {worst_id:.3e}")

    print("\n=== T3  d <= |dsigma| + log(1 + e^{-|dsigma|})  (proved, sharp) ===")
    worst = -math.inf
    worst_pair = None
    n = 0
    for i in range(0, len(pools), 1):
        for j in range(i + 1, len(pools)):
            a, bsig = pools[i], pools[j]
            dd, _ = C.d_and_A(a, bsig)
            D = abs(a.sigma - bsig.sigma)
            slack = dd - (D + math.log1p(math.exp(-D)))
            n += 1
            if slack > worst:
                worst, worst_pair = slack, (i, j)
    print(f"  {n} pairs; worst violation = {worst:.3e}  "
          f"(negative means the bound holds)")

    print("\n  sharpness ladder: a_r = (r, 1, ..., 1) against a flat b")
    print(f"  {'r':>8} {'dsigma':>10} {'d':>12} {'bound':>12} {'ratio':>9}")
    for k in (2, 4, 8, 20, 60, 200):
        a = C.Sig.compressed([10 ** k, 1], [1, 10 ** k])      # sigma = 0
        for lam in (0.0, 0.5, 1.5):
            bsig = C.Sig.compressed([int(round(math.exp(k * math.log(10)
                                                        * math.exp(-lam))))],
                                    [10 ** k])
            dd, _ = C.d_and_A(a, bsig)
            D = abs(a.sigma - bsig.sigma)
            bd = D + math.log1p(math.exp(-D))
            print(f"  {k:>8} {D:10.5f} {dd:12.8f} {bd:12.8f} {dd/bd:9.5f}")

    print("\n=== T4  the sharp large-scale bound for C_4 ===")
    print(f"  s <= log(golden ratio) = {math.log(PHI):.10f}   "
          f"(from the sharp T3)")
    lo, hi = 0.0, 3.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid <= 2 * math.log1p(math.exp(-mid)):
            lo = mid
        else:
            hi = mid
    print(f"  s <= {lo:.10f}                     "
          f"(from FINDINGS' unsharpened 2 log(1+e^-D))")


if __name__ == "__main__":
    main()
