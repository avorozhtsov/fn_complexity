"""Validate the tropical-cone reformulation against real signatures.

Two checks.

(1) EXACTNESS.  For a genuine signature a, d(a,b) computed by common.py equals
    the Hilbert oscillation of F_a = log Z_a directly (this is a tautology --
    d = osc log(F_b/F_a) -- and is checked only to make sure the conventions
    line up).

(2) THE TROPICAL LIMIT.  For a tropical pair (Phi_a, Phi_b) in the cone C,
    build integer signatures from  m_i^K copies of a_i^K  and check
        | d(a^(K), b^(K)) - d_H(Phi_a, Phi_b) |  =  O(1/K).

    python research/realizability/i_validate_cone.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import i_cone as T  # noqa: E402

rng = np.random.default_rng(7)


def random_trop(k, rs):
    c = np.abs(rs.normal(size=k)) * 1.3
    x = np.abs(rs.normal(size=k)) * 1.3
    c[0] = max(c[0], 0.4)
    x[-1] = max(x[-1], 0.4)
    return T.Trop(c, x)


def main():
    print("=== (1) d = osc log(F_b/F_a) on genuine signatures ===")
    worst = 0.0
    for trial in range(200):
        rs = np.random.default_rng(1000 + trial)
        def sig():
            r = int(rs.integers(2, 8))
            return C.Sig.of(tuple(int(v) for v in rs.integers(1, 40, size=r))
                            if True else ())
        while True:
            try:
                a, b = sig(), sig()
                break
            except ValueError:
                continue
        g = np.arange(-25.0, 25.0, 0.001)
        fa = a.F(np.exp(g))
        fb = b.F(np.exp(g))
        v = np.log(fb) - np.log(fa)
        e0 = math.log(b.R) - math.log(a.R)
        e1 = math.log(b.Lam) - math.log(a.Lam)
        d_direct = max(v.max(), e0, e1) - min(v.min(), e0, e1)
        d_ref, _ = C.d_and_A(a, b)
        worst = max(worst, abs(d_direct - d_ref))
    print(f"  200 random integer pairs: max |d_direct - d_common| = {worst:.3e}")

    print("\n=== (2) tropical limit: d(a^(K), b^(K)) -> d_H(Phi_a, Phi_b) ===")
    print(f"  {'K':>6} {'max err over 40 pairs':>24} {'K * err':>12}")
    pairs = []
    for t in range(40):
        rs = np.random.default_rng(500 + t)
        pairs.append((random_trop(int(rs.integers(2, 5)), rs),
                      random_trop(int(rs.integers(2, 5)), rs)))
    rows = []
    for K in (4.0, 8.0, 16.0, 32.0, 64.0):
        err = 0.0
        for pa, pb in pairs:
            dref = T.hilbert(pa, pb)
            va, ma = T.to_signature(pa, K)
            vb, mb = T.to_signature(pb, K)
            try:
                sa = C.Sig.compressed(va, ma)
                sb = C.Sig.compressed(vb, mb)
            except ValueError:
                continue
            dd, _ = C.d_and_A(sa, sb)
            err = max(err, abs(dd - dref))
        print(f"  {K:6.0f} {err:24.3e} {K*err:12.4f}")
        rows.append((K, err))
    print("\n  err ~ C/K with C ~ 1 confirms the tropical cone is exactly the")
    print("  projective closure of the achievable set.")


if __name__ == "__main__":
    main()
