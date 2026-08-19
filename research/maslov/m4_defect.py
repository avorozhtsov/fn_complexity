#!/usr/bin/env python3
"""Brief M item 4 -- how the structure theorem deforms in lambda.

Write ``K(t) = log INT e^{t (f - INT f rho)} rho`` for the centred cumulant
generating function of ``f = u_a - u_b`` under the prior.  Then, *identically*,

    A_lambda - (Psi(a) - Psi(b)) = (K(lambda) - K(-lambda)) / (2 lambda)
                                 = sum_{j odd >= 3} kappa_j lambda^{j-1} / j!
                                 = kappa_3 lambda^2 / 6 + O(lambda^4)
    S_lambda = (softmax - softmin)/2
             = (K(lambda) + K(-lambda)) / (2 lambda)
             = lambda Var_rho(f) / 2 + O(lambda^3)

and since ``K(t) >= 0`` for all ``t`` (Jensen),

    |A_lambda - d Psi|  <=  S_lambda          for every rho and every lambda.

That is the exact deformation of brief G's defect bound: at ``lambda = infinity``
it reads ``|D| <= d/2``, the weak form of (B1).  Both sides vanish at
``lambda = 0``, the left one twice as fast.

What does **not** deform is the sharp universal constant ``(log 2)/2``.  It comes
from the structure theorem, which controls the *midrange* against the two
asymptotic values of ``f``; a soft midrange at finite ``lambda`` is a
``rho``-weighted interior average, and it can sit anywhere between the two
asymptotes, so ``|A_lambda - d psi|`` grows like ``|sigma_a - sigma_b|/2``,
without bound.  Demonstrated on flat pairs below.

What *does* deform is the potential-free statement.  Around any directed cycle
of ``A_lambda`` the potential cancels, so

    |curl A_lambda| = sum_e |A_lambda(e)|  <=  sum_e S_lambda(e),

the soft-metric perimeter -- the exact ``lambda``-deformation of brief D(d).
The sharp constant is then measured by search.

Writes ``m4_defect.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import common as C

RNG = np.random.default_rng(11)


def cumulants(f, w, order=5):
    """Centred cumulants of ``f`` under the discrete prior ``w``."""

    m = float((w * f).sum())
    c = f - m
    mu = [float((w * c**k).sum()) for k in range(order + 1)]
    k2 = mu[2]
    k3 = mu[3]
    k4 = mu[4] - 3 * mu[2] ** 2
    k5 = mu[5] - 10 * mu[3] * mu[2]
    return m, k2, k3, k4, k5


def u_from_logs(x, s):
    """``u`` of a signature whose log-entries are ``x`` (reals >= 0)."""

    z = np.outer(np.exp(s), np.asarray(x, dtype=float))
    mx = z.max(axis=1)
    return np.log(mx + np.log(np.exp(z - mx[:, None]).sum(axis=1)))


def main():
    rows, out = [], []

    def say(t=""):
        print(t)
        out.append(t)

    s, w = C.uniform_prior(-14.0, 14.0, 5601)
    pool = C.f11_pool()
    n = len(pool)
    U, um, up = C.u_matrix(pool, s)
    Psi0 = U @ w

    # ---------------------------------------------------- 4a. the two identities
    say("=" * 88)
    say("4a.  The cumulant identities (identities; the numbers only check the code)")
    say("=" * 88)
    say(f"  {'lambda':>9s} {'max |A_lam - dPsi - odd(K)/2lam|':>34s} "
        f"{'max |S_lam - even(K)/2lam|':>28s} {'max |A-dPsi|/S_lam':>20s}")
    idx = list(range(0, n, 7))
    for lam in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0):
        e1 = e2 = ratio = 0.0
        for i in idx:
            for j in idx:
                if i >= j:
                    continue
                f = U[i] - U[j]
                mean = float((w * f).sum())
                c = f - mean
                Kp = math.log(float((w * np.exp(np.clip(lam * c, -700, 700))).sum()))
                Km = math.log(float((w * np.exp(np.clip(-lam * c, -700, 700))).sum()))
                Al = C.soft_mid(f, w, lam)
                Sl = 0.5 * (C.softmax(f, w, lam) - C.softmin(f, w, lam))
                e1 = max(e1, abs((Al - (Psi0[i] - Psi0[j])) - (Kp - Km) / (2 * lam)))
                e2 = max(e2, abs(Sl - (Kp + Km) / (2 * lam)))
                if Sl > 0:
                    ratio = max(ratio, abs(Al - (Psi0[i] - Psi0[j])) / Sl)
        say(f"  {lam:9.4g} {e1:34.3e} {e2:28.3e} {ratio:20.12f}")
        rows.append(dict(block="identity", lam=lam, err_odd=e1, err_even=e2, max_ratio=ratio))
    say("  the last column is the proved bound |A_lam - dPsi| <= S_lam, never exceeded")

    # ---------------------------------------------------- 4b. leading order
    say()
    say("=" * 88)
    say("4b.  Leading order:  (A_lam - dPsi)/lam^2 -> kappa_3/6 ,  S_lam/lam -> Var/2")
    say("=" * 88)
    say(f"  {'pair':>14s} {'lambda':>9s} {'(A-dPsi)/lam^2':>18s} {'kappa_3/6':>14s}"
        f" {'S_lam/lam':>14s} {'Var/2':>14s}")
    for (i, j) in ((0, 40), (17, 200), (48, 64)):
        f = U[i] - U[j]
        mean, k2, k3, k4, k5 = cumulants(f, w)
        for lam in (1e-3, 1e-2, 1e-1):
            Al = C.soft_mid(f, w, lam)
            Sl = 0.5 * (C.softmax(f, w, lam) - C.softmin(f, w, lam))
            say(f"  {str((i,j)):>14s} {lam:9.4g} {(Al-(Psi0[i]-Psi0[j]))/lam**2:18.10e}"
                f" {k3/6:14.6e} {Sl/lam:14.10f} {k2/2:14.10f}")
            rows.append(dict(block="leading", pair=str((i, j)), lam=lam,
                             defect_over_lam2=(Al - (Psi0[i] - Psi0[j])) / lam**2,
                             kappa3_over_6=k3 / 6, S_over_lam=Sl / lam, var_over_2=k2 / 2))

    # ------------------------------------- 4c. (log 2)/2 does not deform
    say()
    say("=" * 88)
    say("4c.  The sharp constant (log 2)/2 does NOT deform: |A_lam - d psi| is unbounded")
    say("=" * 88)
    say("     flat pairs a = (M_a,...) with r_a fibers; psi = 1/2 log(log r . log M);")
    say("     at lambda = infinity brief G gives |A - d psi| = 0 exactly on the flat locus.")
    say(f"  {'|sigma_a-sigma_b|':>18s} {'lam=0':>12s} {'lam=1':>12s} {'lam=10':>12s}"
        f" {'lam=1e4':>12s} {'lam=inf':>12s}  {'(log2)/2':>10s}")
    sw, ww = C.uniform_prior(-60.0, 60.0, 24001)   # wide enough that both asymptotes
    for (Ra, La, Rb, Lb) in [(1.0, 1.0, 1.0, math.e), (1.0, 1.0, 1.0, math.e**3),
                             (1.0, 1.0, 1.0, math.e**6), (1.0, 1.0, 1.0, math.e**10),
                             (1.0, 1.0, 1.0, math.e**16)]:                # are inside supp(rho)
        ua = np.log(Ra + np.exp(sw) * La)
        ub = np.log(Rb + np.exp(sw) * Lb)
        f = ua - ub
        psia = 0.5 * math.log(Ra * La)
        psib = 0.5 * math.log(Rb * Lb)
        dsig = abs(math.log(Ra / La) - math.log(Rb / Lb))
        vals = []
        for lam in (None, 1.0, 10.0, 1e4):
            vals.append(C.soft_mid(f, ww, lam) - (psia - psib))
        v0 = float((ww * f).sum()) - (psia - psib)
        say(f"  {dsig:18.4f} {v0:12.5f} {vals[1]:12.5f} {vals[2]:12.5f} {vals[3]:12.5f}"
            f" {vals[0]:12.5f}  {math.log(2)/2:10.5f}")
        rows.append(dict(block="no_deform", dsigma=dsig, lam0=v0, lam1=vals[1],
                         lam10=vals[2], lam1e4=vals[3], laminf=vals[0]))
    say("     the lam = inf column is 0 to machine precision (brief G, proved);")
    say()
    say("     same flat pairs, but a prior concentrated near beta = 1 (uniform on [-2,2]),")
    say("     where the growth saturates the proved envelope |A_lam - d psi| <= |dsigma|/2 + log 2:")
    say(f"  {'|sigma_a-sigma_b|':>18s} {'lam=0':>12s} {'lam=1':>12s} {'lam=10':>12s}"
        f" {'lam=inf':>12s}  {'envelope':>10s}")
    sn, wn = C.uniform_prior(-2.0, 2.0, 20001)
    for D in (1.0, 3.0, 6.0, 10.0, 16.0, 24.0):
        fn = np.log(1.0 + np.exp(sn)) - np.log(1.0 + np.exp(sn + D))
        fw = np.log(1.0 + np.exp(sw)) - np.log(1.0 + np.exp(sw + D))
        dpsi = 0.0 - 0.5 * D   # psi_a = 0, psi_b = D/2
        v = [float((wn * fn).sum()) - dpsi] + [C.soft_mid(fn, wn, L) - dpsi for L in (1.0, 10.0)]
        v.append(C.soft_mid(fw, ww, None) - dpsi)
        say(f"  {D:18.4f} {v[0]:12.5f} {v[1]:12.5f} {v[2]:12.5f} {v[3]:12.5f}"
            f"  {0.5*D + math.log(2):10.4f}")
        rows.append(dict(block="no_deform_narrow", dsigma=D, lam0=v[0], lam1=v[1],
                         lam10=v[2], laminf=v[3], envelope=0.5 * D + math.log(2)))
    say("     (log 2)/2 = 0.34657 is already exceeded at |dsigma| = 3, at lambda = 1.")

    # ------------------------------------- 4d. the cycle bound, deformed
    say()
    say("=" * 88)
    say("4d.  |curl A_lam| <= sum_e S_lam(e)  (proved), and how tight")
    say("=" * 88)
    Uc = [C.u_of(a, s) for a in C.CYCLE]
    say(f"  {'lambda':>9s} {'|curl A_lam|':>15s} {'sum S_lam':>13s} {'slack':>9s}"
        f"  {'3(log2)/2':>10s}")
    for lam in (0.01, 0.1, 1.0, 10.0, 100.0, 300.0, 1e3, 1e4, None):
        cu, sm = 0.0, 0.0
        for i in range(3):
            f = Uc[i] - Uc[(i + 1) % 3]
            cu += C.soft_mid(f, w, lam)
            sm += (0.5 * (f.max() - f.min()) if lam is None
                   else 0.5 * (C.softmax(f, w, lam) - C.softmin(f, w, lam)))
        tag = "inf" if lam is None else f"{lam:g}"
        say(f"  {tag:>9s} {abs(cu):15.6e} {sm:13.6e} {sm/abs(cu):9.2f}  {3*math.log(2)/2:10.5f}")
        rows.append(dict(block="cycle_bound", lam=tag, curl=abs(cu), sumS=sm, slack=sm / abs(cu)))

    # ------------------------------------- 4e. the sharp constant, searched
    say()
    say("=" * 88)
    say("4e.  mean_e |A_lam| around a directed 3-cycle ( = |curl A_lam|/3 ): a WEAK search.")
    say("     14 restarts of a compass search over 3-fiber triples.  At lambda = infinity it")
    say("     reaches only 0.159 of brief G's own hill-climb optimum (0.34274), so every")
    say("     entry is a lower bound.  The decisive measurement is m4b_cycle_strength.py.")
    say("=" * 88)

    ss, ws = C.uniform_prior(-14.0, 14.0, 1401)

    def curl_of(x, lam):
        us = [u_from_logs(v, ss) for v in x]
        return sum(C.soft_mid(us[i] - us[(i + 1) % 3], ws, lam) for i in range(3))

    say(f"  {'lambda':>9s} {'max |curl|/3':>14s} {'vs (log2)/2':>12s}  best triple (log-entries)")
    for lam in (1.0, 3.0, 10.0, 100.0, 1e4, None):
        best, bx = 0.0, None
        for _ in range(14):
            x = [np.sort(RNG.uniform(0.0, 4.0, 3))[::-1].copy() for _ in range(3)]
            for v in x:
                v[0] = max(v[0], 0.7)
            cur = abs(curl_of(x, lam))
            step, it = 0.8, 0
            while step > 3e-3 and it < 120:
                it += 1
                moved = False
                for t in range(3):
                    for u_ in range(3):
                        for d in (+step, -step):
                            y = [v.copy() for v in x]
                            y[t][u_] = max(0.0, y[t][u_] + d)
                            c = abs(curl_of(y, lam))
                            if c > cur + 1e-14:
                                x, cur, moved = y, c, True
                if not moved:
                    step *= 0.5
            if cur > best:
                best, bx = cur, [v.copy() for v in x]
        tag = "inf" if lam is None else f"{lam:g}"
        say(f"  {tag:>9s} {best/3:14.8f} {best/3/(math.log(2)/2):12.4f}  "
            + " | ".join("(" + ",".join(f"{math.exp(z):.3f}" for z in v) + ")" for v in bx))
        rows.append(dict(block="sharp", lam=tag, max_mean_abs_A=best / 3,
                         over_log2_half=best / 3 / (math.log(2) / 2),
                         triple=str([[float(math.exp(z)) for z in v] for v in bx])))

    with open("m4_defect.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    say("\nwrote m4_defect.csv")


if __name__ == "__main__":
    main()
