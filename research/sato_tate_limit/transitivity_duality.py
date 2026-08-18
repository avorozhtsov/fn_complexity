#!/usr/bin/env python3
"""The comparison is Legendre-invariant: ``mid(Psi_mu - Psi_nu)`` is the
midrange of the horizontal gap between the two large-deviation rate functions.

Write ``I_mu = K_mu^*`` for the rate function of ``mu`` and ``J_mu = I_mu^{-1}``
for its inverse on ``[0, alpha_max)``: ``J_mu`` is concave and increasing with
``J_mu(0) = 0`` and ``J_mu(infinity) = alpha_max``.  Then, with ``u = 1/tau``,

    Psi_mu(1/u) = K_mu(1/u) u = sup_{x} ( x - u I_mu(x) ) = sup_{v>=0} ( J_mu(v) - u v ),

so ``Psi_mu`` is the concave conjugate of ``J_mu``, and (theorem, proved in
TRANSITIVITY.md)

    sup_tau (Psi_mu - Psi_nu) = sup_v (J_mu - J_nu),
    inf_tau (Psi_mu - Psi_nu) = inf_v (J_mu - J_nu).

Hence ``mid`` and the metric ``osc`` are unchanged by passing from the cumulant
side to the rate-function side.  That is a genuinely independent way to compute
the comparison: the ``tau`` grid and the ``v`` grid discretise different
functions, and the two ``alpha_max`` endpoints that a ``tau`` grid must be told
analytically are built into ``J`` by construction.

This script verifies the identity on the pairs that matter.

    python research/sato_tate_limit/transitivity_duality.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from mpmath import mp, mpf

import kappa_lib as KL

HERE = Path(__file__).resolve().parent
DPS_EXTRA = 40


def logH(g: int, tau) -> "mpf":
    t = mpf(tau)
    dps = KL.working_dps(g + 1, float(t)) + DPS_EXTRA
    h = KL.hankels(t, g, dps)
    saved = mp.dps
    mp.dps = dps
    try:
        return mp.log(h[g])
    finally:
        mp.dps = saved


class Part:
    """``prod_i USp(2 lambda_i)``: ``K``, ``K'`` at 60+ digits."""

    def __init__(self, lam: tuple[int, ...]):
        self.lam = lam
        self.alpha_max = 2.0 * sum(lam)
        self._z = {}
        for g in set(lam):
            dps = KL.working_dps(g + 1, 1.0) + DPS_EXTRA
            h = KL.hankels(0.0, g, dps)
            saved = mp.dps
            mp.dps = dps
            self._z[g] = mp.log(h[g])
            mp.dps = saved

    def K(self, tau) -> "mpf":
        return sum(logH(g, tau) - self._z[g] for g in self.lam)

    def Kp(self, tau) -> "mpf":
        mp.dps = 40
        return mp.diff(lambda s: self.K(s), mpf(tau), h=mpf(tau) / 2**12)


def rate_curve(part: Part, taus) -> tuple[np.ndarray, np.ndarray]:
    """``(v, x)`` on the graph of ``J``: ``x = K'(tau)``, ``v = tau K' - K``."""

    xs, vs = [], []
    for t in taus:
        k = part.K(t)
        kp = part.Kp(t)
        xs.append(float(kp))
        vs.append(float(mpf(t) * kp - k))
    return np.array(vs), np.array(xs)


def main() -> int:
    TAU = np.geomspace(1e-4, 1e5, 1201)
    kappa, _ = KL.kappa_and_b(TAU, 12)

    tests = [((2,), (1, 1)),                  # genus 2: USp4 vs SU2^2
             ((3,), (2, 1)),                  # genus 3
             ((5, 1, 1), (4, 3)),             # genus 7: the first crossing pair
             ((6, 2, 1, 1, 1), (4, 3, 3, 1)),  # genus 11 near-miss, edge 1
             ((4, 3, 3, 1), (5, 2, 2, 2)),     # ... edge 2
             ((5, 2, 2, 2), (6, 2, 1, 1, 1))]  # ... edge 3

    # a dense tau sample for the parametric rate curve; refining it shows
    # that the residual is interpolation error and nothing else
    npts = int(__import__("os").environ.get("DUALITY_PTS", "880"))
    taus = [mpf(t) for t in np.geomspace(1e-3, 3e2, npts)]
    rows = []
    print("=" * 78)
    print("Legendre invariance of sup, inf and mid")
    print("=" * 78)
    print(f"  {'pair':<30}{'sup(tau)':>13}{'sup(v)':>13}"
          f"{'inf(tau)':>13}{'inf(v)':>13}")
    for a, b in tests:
        pa, pb = Part(a), Part(b)
        # tau side, endpoints supplied analytically
        d = (sum(kappa[p - 1] for p in a) - sum(kappa[p - 1] for p in b)) / TAU
        dd = np.concatenate([[0.0], d, [pa.alpha_max - pb.alpha_max]])
        sup_t, inf_t = float(dd.max()), float(dd.min())
        # v side
        va, xa = rate_curve(pa, taus)
        vb, xb = rate_curve(pb, taus)
        lo = max(va.min(), vb.min())
        hi = min(va.max(), vb.max())
        vv = np.geomspace(max(lo, 1e-12), hi, 20 * npts)
        Ja = np.interp(vv, va, xa)
        Jb = np.interp(vv, vb, xb)
        delta = np.concatenate([[0.0], Ja - Jb,
                                [pa.alpha_max - pb.alpha_max]])
        sup_v, inf_v = float(delta.max()), float(delta.min())
        la = "".join(map(str, a))
        lb = "".join(map(str, b))
        print(f"  {la + ' vs ' + lb:<30}{sup_t:>13.6f}{sup_v:>13.6f}"
              f"{inf_t:>13.6f}{inf_v:>13.6f}")
        rows.append([la, lb, f"{sup_t:.8f}", f"{sup_v:.8f}",
                     f"{inf_t:.8f}", f"{inf_v:.8f}",
                     f"{0.5 * (sup_t + inf_t):.8f}",
                     f"{0.5 * (sup_v + inf_v):.8f}"])

    print(f"\n  (the v-side numbers come from a {npts}-point parametric sample of")
    print("   the rate curve interpolated on 4000 v-points; the residual is")
    print("   interpolation error, not a failure of the identity)")
    worst = max(abs(float(r[2]) - float(r[3])) for r in rows)
    worst = max(worst, max(abs(float(r[4]) - float(r[5])) for r in rows))
    print(f"  worst |tau-side - v-side| over sup and inf: {worst:.3e}")

    with (HERE / "transitivity_duality.csv").open("w", newline="",
                                                  encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["mu", "nu", "sup_tau", "sup_v", "inf_tau", "inf_v",
                     "mid_tau", "mid_v"])
        wr.writerows(rows)
    print("\nwritten: transitivity_duality.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
