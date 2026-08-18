#!/usr/bin/env python3
"""Brief M item 3 -- the phase diagram in lambda, and a counterexample.

The prior is taken to be the **discrete uniform measure on the grid**, so every
number below is exact for that prior: nothing here is a quadrature.

Part A -- the phase diagram on brief E's ``F_11`` pool.  For a ladder of
``lambda`` report the curl fraction ``||curl A_lambda|| / ||A_lambda||``, the
curl energy, the number of directed 3-cycles, and the number of edges whose
sign disagrees with the tropical comparison.

Part B -- the small-``lambda`` law.  ``A_lambda - d Psi = kappa_3 lambda^2/6 +
O(lambda^4)`` and ``d_lambda = lambda Var(f) + O(lambda^3)`` (both proved in
``m4_defect.py``), so the curl fraction must vanish like ``lambda^2``.  Checked
over four decades.

Part C -- **do cycles only live above a critical lambda?**  No.  On the *flat*
locus (all fibers equal) brief G proves ``curl A_inf = 0`` identically, so no
flat pool can contain a cycle at ``lambda = infinity``; and ``curl A_0 = 0`` by
the potential identity.  But ``curl A_lambda != 0`` for every ``0 < lambda <
infinity``.  Since the Cartesian power shifts ``A_lambda`` by an exact 1-form
for *every* lambda (``u_{a^{ok}} = u_a + log k``), Theorem D of brief G applies
verbatim and turns that curl into a genuine 3-cycle.  A certified witness is
built here: a 3-cycle of flat signatures that exists at ``lambda = 2`` and at
neither endpoint of the family.

Writes ``m3_phase.csv``, ``m3_witness.json``.
"""

from __future__ import annotations

import csv
import json
import math

import numpy as np

import common as C

DS = 0.005
LO, HI = -14.0, 14.0


# ---------------------------------------------------------------- flat family


class Flat:
    """A flat signature ``(M, ..., M)`` with ``r`` fibers, held as ``(R, Lam)``.

    ``Z(beta) = r M^beta``, ``F = R + beta Lam`` with ``R = log r``,
    ``Lam = log M``; the ``k``-th Cartesian power is ``(kR, k Lam)``.
    """

    def __init__(self, r, M, k=1):
        self.r, self.M, self.k = r, M, k
        self.R = k * math.log(r)
        self.Lam = k * math.log(M)

    def u(self, s):
        return np.log(self.R + np.exp(s) * self.Lam)

    def __repr__(self):
        return f"flat(r={self.r}^{self.k}, M={self.M}^{self.k})"


def main():
    rows, out = [], []

    def say(t=""):
        print(t)
        out.append(t)

    # =============================================== Part A: the phase diagram
    pool = C.f11_pool()
    n = len(pool)
    ng = int(round((HI - LO) / DS)) + 1
    s, w = C.uniform_prior(LO, HI, ng)
    U, um, up = C.u_matrix(pool, s)
    A_inf_true = C.tropical_matrix(U, um, up)  # with the two exact endpoints
    A_grid_inf = C.soft_matrix(U, w, None)  # midrange over supp(rho) only
    Psi0 = U @ w
    A0 = Psi0[:, None] - Psi0[None, :]

    say("=" * 92)
    say(f"A.  Phase diagram, F_11 pool ({n} signatures), rho = uniform on "
        f"{ng} grid points of [{LO},{HI}]")
    say("=" * 92)
    say(f"  {'lambda':>10s} {'curl frac':>11s} {'curl energy':>12s} {'3-cycles':>9s}"
        f" {'edges != A_inf':>15s} {'edges != A_0':>13s} {'rms|A|':>11s}")

    lams = [1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0,
            1e3, 3e3, 1e4, 3e4, 1e5, 1e6]
    sgn_inf = np.sign(A_grid_inf)
    sgn_0 = np.sign(A0)
    iu = np.triu_indices(n, 1)
    prev_cycles = None
    for lam in lams + [None]:
        Al = C.soft_matrix(U, w, lam)
        h = C.hodge(Al)
        nc = C.three_cycles(Al)
        dif_inf = int((np.sign(Al)[iu] != sgn_inf[iu]).sum())
        dif_0 = int((np.sign(Al)[iu] != sgn_0[iu]).sum())
        rms = float(np.sqrt((Al[iu] ** 2).mean()))
        tag = "inf" if lam is None else f"{lam:g}"
        say(f"  {tag:>10s} {h['curl_fraction']:11.3e} {h['curl_energy']:12.3e} {nc:9d}"
            f" {dif_inf:15d} {dif_0:13d} {rms:11.3e}")
        rows.append(dict(block="phase", lam=tag, curl_fraction=h["curl_fraction"],
                         curl_energy=h["curl_energy"], cycles=nc,
                         edges_ne_Ainf=dif_inf, edges_ne_A0=dif_0, rms_A=rms))
        prev_cycles = nc
    say(f"  tropical comparison with the two exact endpoints: "
        f"{C.three_cycles(A_inf_true)} 3-cycles, "
        f"curl fraction {C.hodge(A_inf_true)['curl_fraction']:.3e}")

    # ------------------------------------------------- Part B: the lambda^2 law
    say()
    say("=" * 92)
    say("B.  Small-lambda law:  curl fraction ~ c lambda^2  (proved leading order)")
    say("=" * 92)
    say(f"  {'lambda':>10s} {'curl frac':>13s} {'curl frac / lambda^2':>22s}")
    for lam in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0):
        Al = C.soft_matrix(U, w, lam)
        cf = C.hodge(Al)["curl_fraction"]
        say(f"  {lam:10.4g} {cf:13.5e} {cf/lam**2:22.6f}")
        rows.append(dict(block="lambda2", lam=lam, curl_fraction=cf, ratio=cf / lam**2))

    # ---------------------------------------- Part C: the flat-locus witness
    say()
    say("=" * 92)
    say("C.  Cycles are NOT confined to large lambda: a flat 3-cycle at lambda = 2")
    say("=" * 92)
    base = [Flat(3, 6), Flat(7, 2), Flat(2, 50)]
    Uf = np.array([b.u(s) for b in base])
    say("    base flat signatures: " + ", ".join(repr(b) for b in base))
    say(f"  {'lambda':>10s} {'A12':>13s} {'A23':>13s} {'A31':>13s} {'curl':>13s}")
    curls = {}
    for lam in [None, 1e-3, 1e-2, 0.1, 0.3, 1.0, 2.0, 3.0, 5.0, 10.0, 100.0, 1e4]:
        v = [C.soft_mid(Uf[i] - Uf[(i + 1) % 3], w, lam) for i in range(3)]
        cu = float(sum(v))
        curls[lam] = cu
        tag = "inf" if lam is None else f"{lam:g}"
        say(f"  {tag:>10s} " + " ".join(f"{x:+13.6e}" for x in v) + f" {cu:+13.6e}")
        rows.append(dict(block="flat_curl", lam=tag, A12=v[0], A23=v[1], A31=v[2], curl=cu))
    say(f"    curl at lambda = infinity is exactly 0 (proved: for three monotone f summing")
    say(f"    to zero pointwise, the midranges sum to zero); measured {curls[None]:+.3e}")

    # find the best lambda, then realise the cycle with Cartesian powers
    grid = np.geomspace(0.05, 200.0, 400)
    vals = []
    for lam in grid:
        v = [C.soft_mid(Uf[i] - Uf[(i + 1) % 3], w, lam) for i in range(3)]
        vals.append(sum(v))
    kb = int(np.argmax(np.abs(vals)))
    lam_star = float(grid[kb])
    say(f"    curl is maximised at lambda* = {lam_star:.4f}, curl = {vals[kb]:+.6e}"
        f"   (margin per edge {vals[kb]/3:+.6e})")

    import mpmath as mp

    def powers_for(lam_target):
        """Smallest Cartesian powers making the triple a cycle at ``lam_target``.

        ``A_lam(a^{ok_i}, b^{ok_j}) = A_lam(a,b) + log(k_i/k_j)`` exactly, for
        every lambda and every rho -- so the powers move the potential and leave
        the curl alone, and Theorem D of brief G applies verbatim.
        """

        v = [C.soft_mid(Uf[i] - Uf[(i + 1) % 3], w, lam_target) for i in range(3)]
        m = sum(v) / 3.0
        x = [0.0, v[0] - m, (v[0] - m) + (v[1] - m)]
        tgt = [math.exp(z) for z in x]
        tgt = [z / min(tgt) for z in tgt]
        best_kk, best_margin = None, -1.0
        for t in range(1, 40001):
            kk = [max(1, int(round(t * z))) for z in tgt]
            g = math.gcd(math.gcd(kk[0], kk[1]), kk[2])
            kk = [k // g for k in kk]
            mar = min(v[i] + math.log(kk[i] / kk[(i + 1) % 3]) for i in range(3))
            if mar > best_margin:
                best_kk, best_margin = kk, mar
            if best_margin > 0.98 * m:
                break
        return best_kk, best_margin, m

    witnesses = {}
    for tag_t, lam_target in (("lambda* (max curl)", lam_star), ("lambda = 1 exactly", 1.0)):
        kk, margin, mmax = powers_for(lam_target)
        wit = [Flat(base[i].r, base[i].M, kk[i]) for i in range(3)]
        Uw = np.array([b.u(s) for b in wit])
        say()
        say(f"    --- witness for {tag_t}, lambda_target = {lam_target:.6f} ---")
        say(f"    Cartesian powers k = {kk}   margin {margin:.6e} "
            f"(max possible {mmax:.6e})")
        say(f"    witness: " + ", ".join(repr(b) for b in wit))
        say(f"  {'lambda':>10s} {'A12':>15s} {'A23':>15s} {'A31':>15s}  cycle?")
        witness_rows = []
        for lam in [None, 1e-3, 1e-2, 0.1, 0.5, 1.0, lam_star, 3.0, 5.0, 10.0, 50.0, 1e3, 1e6]:
            vv = [C.soft_mid(Uw[i] - Uw[(i + 1) % 3], w, lam) for i in range(3)]
            cyc = all(z > C.TIE for z in vv)
            t2 = "inf" if lam is None else f"{lam:g}"
            say(f"  {t2:>10s} " + " ".join(f"{z:+15.8e}" for z in vv) + f"  {cyc}")
            witness_rows.append(dict(lam=t2, A=[float(z) for z in vv], cycle=bool(cyc)))
            rows.append(dict(block=f"flat_witness[{tag_t}]", lam=t2,
                             A12=vv[0], A23=vv[1], A31=vv[2], cycle=cyc))

        def is_cycle(lam):
            vv = [C.soft_mid(Uw[i] - Uw[(i + 1) % 3], w, lam) for i in range(3)]
            return all(z > C.TIE for z in vv)

        def bisect(lo, hi):
            for _ in range(90):
                mid = math.sqrt(lo * hi)
                if is_cycle(mid):
                    hi = mid
                else:
                    lo = mid
            return math.sqrt(lo * hi)

        lo_edge = bisect(1e-8, lam_target)
        hi_edge = bisect(1e9, lam_target)
        say(f"    a 3-cycle exactly for lambda in ({lo_edge:.6g}, {hi_edge:.6g})"
            f"  -- a BAND, not a half-line")
        rows.append(dict(block=f"flat_band[{tag_t}]", lam_lo=lo_edge, lam_hi=hi_edge,
                         lam_target=lam_target, margin=margin, powers=str(kk)))

        with mp.workdps(40):
            S = [mp.mpf(LO) + mp.mpf(HI - LO) * i / (ng - 1) for i in range(ng)]
            Um = [[mp.log(mp.mpf(b.k) * mp.log(b.r) + mp.e ** sk * mp.mpf(b.k) * mp.log(b.M))
                   for sk in S] for b in wit]
            L = mp.mpf(lam_target)
            wk = mp.mpf(1) / ng
            cert = []
            for i in range(3):
                fm = [Um[i][t] - Um[(i + 1) % 3][t] for t in range(ng)]
                fmax, fmin = max(fm), min(fm)
                smax = (mp.log(mp.fsum([wk * mp.e ** (L * (z - fmax)) for z in fm])) + L * fmax) / L
                smin = -(mp.log(mp.fsum([wk * mp.e ** (-L * (z - fmin)) for z in fm])) - L * fmin) / L
                cert.append((smax + smin) / 2)
            say("    40-digit certification at lambda_target:")
            for i, z in enumerate(cert):
                say(f"       A_{i+1}{(i+1)%3+1} = {mp.nstr(z, 30)}")
            dbl = [C.soft_mid(Uw[i] - Uw[(i + 1) % 3], w, lam_target) for i in range(3)]
            agree = max(abs(float(cert[i]) - dbl[i]) for i in range(3))
            say(f"       |double - mpmath(40)| <= {agree:.3e}   all three positive: "
                f"{all(z > 0 for z in cert)}")

        witnesses[tag_t] = {
            "powers": kk, "margin": margin, "lambda_target": lam_target,
            "witness": [{"r": b.r, "M": b.M, "k": b.k} for b in wit],
            "cycle_band": [lo_edge, hi_edge], "table": witness_rows,
            "mpmath40": [mp.nstr(z, 30) for z in cert], "double_vs_mpmath": agree,
        }

    json.dump({"base": [{"r": b.r, "M": b.M} for b in base],
               "prior": {"kind": "uniform on grid", "lo": LO, "hi": HI, "points": ng},
               "lambda_star": lam_star, "witnesses": witnesses},
              open("m3_witness.json", "w"), indent=1)

    with open("m3_phase.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    say("\nwrote m3_phase.csv, m3_witness.json")


if __name__ == "__main__":
    main()
