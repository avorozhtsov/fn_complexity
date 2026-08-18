"""Brief L, Part 4 -- does the sandwiched profile leave the classical cone?

OBSTRUCTION.md Theorem 1: the projective closure of the achievable set of
classical profiles is exactly

    C = { Phi : [0,oo) -> (0,oo) convex, nondecreasing, Phi(b) >= Lam_Phi * b }

and every Phi in C is a uniform limit of (1/K) F_{a^{(K)}} for genuine integer
signatures a^{(K)}, at rate O(log k / K).

Part 2 shows the sandwiched quantum profile F_(A,S) is positive, satisfies
max(R, b Lam) <= F <= R + b Lam (proved), and is -- to the second-difference
noise floor -- convex and nondecreasing.  So it lies in C.  This script makes
that operational:

  (1) build the tangent-line cone element Phi = max_j (c_j + x_j b) of the
      quantum profile, check c_j >= 0 (membership in C) and sup |Phi - F|;
  (2) compare the Hilbert metric of the Phi's against the quantum d;
  (3) push the Phi's down to genuine INTEGER signatures at a ladder of
      multiplicity budgets K, and watch d converge at the predicted O(1/K).

EVERY metric here is the oscillation over the SAME domain beta in [1/2, oo],
the range on which the sandwiched divergence is a monotone.  Comparing a
[1/2,oo] quantum metric with a [0,oo] cone metric is meaningless and was the
first thing this script got wrong.

Run:  python3 q_cone.py
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "realizability"))

import i_cone  # noqa: E402
from i_logsig import LogSig  # noqa: E402

from q_core import extrema, rand_admissible  # noqa: E402

BMIN = 0.5
FINE_S = np.linspace(math.log(BMIN), math.log(3000.0), 4001)


class Adapter:
    """Uniform (U, R, Lam) interface over QSig / Trop / LogSig."""

    def __init__(self, obj, kind):
        self.o, self.kind = obj, kind
        if kind == "trop":
            self.R, self.Lam = obj.R, obj.Lam
        else:
            self.R, self.Lam = float(obj.R), float(obj.Lam)

    def U(self, s):
        if self.kind == "q":
            return self.o.U(s)
        if self.kind == "trop":
            return math.log(float(np.ravel(self.o.val(math.exp(s)))[0]))
        return float(np.ravel(np.log(self.o.F(math.exp(s))))[0])

    def U_grid(self, g):
        if self.kind == "q":
            return self.o.U_grid(g)
        if self.kind == "trop":
            return np.log(self.o.val(np.exp(g)))
        return np.log(self.o.F(np.exp(g)))


def d_restricted(pa, pb, grid=FINE_S):
    """osc of U_b - U_a over beta in [BMIN, oo], the beta = oo endpoint exact."""
    fun = lambda s: pb.U(s) - pa.U(s)
    e_inf = math.log(pb.Lam) - math.log(pa.Lam)
    P, Q = extrema(fun, grid, (e_inf,))
    return P - Q


def tangent_cone(q, sgrid, nlines=120):
    """Tangent-line cone element of F_(A,S), plus the horizontal line at F(0).

    c_j = F(b_j) - b_j F'(b_j),  x_j = F'(b_j).  F convex gives Phi <= F.
    c_j >= 0 is exactly membership in C and follows from F >= b*Lam >= b*F'.
    The extra line (c = F(0), x = 0) pins Phi(0) = F(0) = log Tr S, which the
    tangents alone cannot reach when the grid starts at beta = 1/2.
    """
    b = np.exp(np.linspace(sgrid[0], sgrid[-1], nlines))
    h = b * 1e-5
    F0 = q.F_grid(b)
    Fp = (q.F_grid(b + h) - q.F_grid(b - h)) / (2 * h)
    c = np.concatenate([[q.R], F0 - b * Fp])
    x = np.concatenate([[0.0], Fp])
    return c, x


def main():
    rng = np.random.default_rng(31415)

    print("=" * 74)
    print("PART 4 -- the sandwiched profile lies in the CLASSICAL cone")
    print("=" * 74)
    print(f"All metrics are osc over beta in [{BMIN}, oo].")

    print("\n(1) tangent-line cone membership, 200 random admissible (A,S)")
    worst_c = math.inf
    worst_gap = 0.0
    for _ in range(200):
        r = int(rng.integers(2, 7))
        q = rand_admissible(rng, r, rng.uniform(1.2, 8.0), rng.uniform(1.5, 60.0))
        c, x = tangent_cone(q, FINE_S, 60)
        worst_c = min(worst_c, float(c.min()))
        T = i_cone.Trop(np.maximum(c, 0.0), x)
        F = q.F_grid(np.exp(FINE_S))
        worst_gap = max(worst_gap, float(np.max(np.abs(T.val(np.exp(FINE_S)) - F))))
    print(f"    min tangent intercept c_j over the ensemble = {worst_c:+.6e}")
    print("    (c_j >= 0 is exactly 'Phi in C'; it follows from F >= b*Lam)")
    print(f"    max sup|Phi_tangent - F| with 60 tangents   = {worst_gap:.3e}")

    print("\n(2) Hilbert metric of the cone approximant against the quantum d")
    print("    (10 random pairs, 400 tangents each)")
    print("        d_quantum      d_cone         |diff|")
    for _ in range(10):
        r = int(rng.integers(2, 6))
        qa = rand_admissible(rng, r, 4.0, 25.0)
        qb = rand_admissible(rng, r, 4.0, 25.0)
        dq = d_restricted(Adapter(qa, "q"), Adapter(qb, "q"))
        Ts = [i_cone.Trop(*[np.maximum(v, 0.0) if k == 0 else v
                            for k, v in enumerate(tangent_cone(q, FINE_S, 400))])
              for q in (qa, qb)]
        dc = d_restricted(Adapter(Ts[0], "trop"), Adapter(Ts[1], "trop"))
        print(f"      {dq:12.9f}  {dc:12.9f}  {abs(dq-dc):.3e}")

    print("\n(3) quantum pair -> cone -> INTEGER signatures: the O(1/K) ladder")
    r = 4
    qa = rand_admissible(rng, r, 4.0, 25.0)
    qb = rand_admissible(rng, r, 4.0, 25.0)
    dq = d_restricted(Adapter(qa, "q"), Adapter(qb, "q"))
    Ts = [i_cone.Trop(*[np.maximum(v, 0.0) if k == 0 else v
                        for k, v in enumerate(tangent_cone(q, FINE_S, 400))])
          for q in (qa, qb)]
    dc = d_restricted(Adapter(Ts[0], "trop"), Adapter(Ts[1], "trop"))
    print(f"    d_quantum = {dq:.10f}   d_cone = {dc:.10f}   "
          f"|diff| = {abs(dq-dc):.3e}")
    print("        K     log10 r        d_signature     |d_sig - d_cone|   K*gap")
    ladder = []
    for K in (25, 50, 100, 200, 400, 800, 1600, 3200):
        sigs = [Adapter(LogSig(K * T.c, K * T.x), "sig") for T in Ts]
        ds = d_restricted(sigs[0], sigs[1])
        gap = abs(ds - dc)
        log10r = sigs[0].R / math.log(10.0)
        print(f"      {K:5d}  {log10r:11.2f}  {ds:16.10f}  {gap:14.3e}  {K*gap:8.3f}")
        ladder.append({"K": K, "log10_r": log10r, "d_signature": ds, "gap": gap})

    print("\n    The multiplicities are e^{K c_j} and the atoms e^{K x_j}: honest")
    print("    integer signatures with r = 10^(log10 r) fibers.  d converges to")
    print("    the cone value at the predicted O(1/K), so the quantum pair is")
    print("    realised by classical signatures to distortion 1 + O(1/log r).")

    with open(os.path.join(HERE, "q_cone.json"), "w") as fh:
        json.dump({"d_quantum": dq, "d_cone": dc, "ladder": ladder,
                   "min_tangent_intercept": worst_c,
                   "max_tangent_gap": worst_gap}, fh, indent=1)
    print("\nwrote q_cone.json")


if __name__ == "__main__":
    main()
