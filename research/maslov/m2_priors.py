#!/usr/bin/env python3
"""Brief M item 2 -- is there a canonical prior?

Three strands.

**(a) The endpoint prior is brief E's potential.**  ``Psi_rho(a) = INT u_a rho``
with ``rho = (delta_{-inf} + delta_{+inf})/2`` is
``(log log r + log log M)/2 = (1/2) log(log r . log M) = (1/2) log phi``,
the brief-D(c)/E endpoint potential -- *exactly*, not approximately.  So
``1/2 log phi`` is already a member of the family and its 82.494 % of ordered
pairs is the number to beat.

**(b) The energy gauge.**  ``a -> a^(t) = (a_i^t)`` gives
``u_{a^(t)}(s) = u_a(s + log t)``: the gauge acts on the ``s``-line by
translation.  Hence ``A_inf`` (and the whole exchange geometry) is gauge
*invariant*, while ``Psi_rho`` is only gauge *equivariant*:
``Psi_rho(a^(t)) = Psi_{rho shifted}(a)``.  A canonical prior would have to be a
translation-invariant probability measure on ``R``.  There is none.  This is a
proof that no ``rho`` is canonical; the numbers below say how much it costs.

**(c) How good can Psi_rho be?**  Order agreement with the tropical comparison
on brief E's ``F_11`` pool, for: the endpoint prior, single temperatures
``delta_s`` (the fixed-``beta`` orders of item 5), uniform/Gaussian priors, and
an *optimised* prior found by mirror descent plus a Frank-Wolfe polish on the
exact 0-1 count.  Benchmarks: ``1/2 log phi`` 82.494 %, HodgeRank ``psi_opt``
98.314 %, and the hard ceiling 99.904 % (>= 42 pairs must be misordered).

Writes ``m2_priors.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import common as C

RNG = np.random.default_rng(20260818)


def agreement_from_psi(psi, A, tie=C.TIE):
    return C.order_agreement(psi, A, tie)


def main():
    rows, out = [], []

    def say(t=""):
        print(t)
        out.append(t)

    pool = C.f11_pool()
    n = len(pool)
    s, w_unif = C.uniform_prior(-14.0, 14.0, 14001)
    U, um, up = C.u_matrix(pool, s)
    A = C.tropical_matrix(U, um, up)
    hod = C.hodge(A)
    iu = np.triu_indices(n, 1)
    npairs = iu[0].size

    say("=" * 78)
    say(f"F_11 pool: {n} signatures, {npairs} unordered pairs")
    say(f"  grad energy {hod['grad_energy']:.6f}   curl energy {hod['curl_energy']:.6f}"
        f"   3-cycles {C.three_cycles(A)}")
    say("=" * 78)

    # ---- (a) the endpoint prior IS 1/2 log phi -------------------------------
    psi_end = np.array([C.psi_endpoint(a) for a in pool])
    psi_end_family = 0.5 * (um + up)  # Psi of rho = (delta_-inf + delta_+inf)/2
    say()
    say("(a) endpoint prior  rho = (delta_{-inf} + delta_{+inf})/2")
    say(f"    max |Psi_rho - (1/2) log phi| = {np.abs(psi_end_family - psi_end).max():.3e}"
        "   (an identity, shown only to confirm the code)")
    a_end = agreement_from_psi(psi_end, A)
    say(f"    order agreement: {a_end:.5f}")
    rows.append(dict(block="endpoint", prior="0.5*(delta_-inf+delta_+inf) = 1/2 log phi",
                     agreement=a_end))

    # ---- (b) the energy gauge ------------------------------------------------
    say()
    say("(b) energy gauge  a -> (a_i^t):  u_{a^(t)}(s) = u_a(s + log t)")
    for t in (2, 3, 5):
        gp = [tuple(int(x) ** t for x in a) for a in pool]
        Ug, umg, upg = C.u_matrix(gp, s)
        Ag = C.tropical_matrix(Ug, umg, upg)
        dev = float(np.abs(Ag - A).max())
        psi_g = np.array([C.psi_endpoint(a) for a in gp])
        say(f"    t={t}:  max|A_inf(a^t) - A_inf(a)| = {dev:.3e}"
            f"    agreement of 1/2 log phi on the gauged pool: "
            f"{agreement_from_psi(psi_g, Ag):.5f}")
        rows.append(dict(block="gauge", t=t, max_dev_A=dev,
                         agreement=agreement_from_psi(psi_g, Ag)))

    # ---- (c) fixed-temperature orders (delta_s) ------------------------------
    say()
    say("(c) single temperatures  rho = delta_s  (the fixed-beta orders)")
    scan = np.linspace(-6.0, 10.0, 641)
    agr = np.array([agreement_from_psi(C.u_matrix(pool, np.array([x]))[0][:, 0], A) for x in scan])
    kbest = int(agr.argmax())
    say(f"    best single temperature: s = {scan[kbest]:+.4f} (beta = {math.exp(scan[kbest]):.4g})"
        f"   agreement {agr[kbest]:.5f}")
    say(f"    s = -inf (order by r):   {agreement_from_psi(um, A):.5f}")
    say(f"    s = +inf (order by M):   {agreement_from_psi(up, A):.5f}")
    for x in (-4, -2, -1, 0, 1, 2, 3, 4, 6, 8):
        k = int(np.abs(scan - x).argmin())
        rows.append(dict(block="delta_s", s=float(scan[k]), agreement=float(agr[k])))
    say("    agreement(s):  " + "  ".join(
        f"s={x:+g}:{agr[int(np.abs(scan-x).argmin())]:.4f}" for x in (-4, -2, 0, 2, 4, 6, 8)))
    rows.append(dict(block="delta_s_best", s=float(scan[kbest]), agreement=float(agr[kbest])))

    # ---- (c2) named priors ---------------------------------------------------
    say()
    say("(c2) named priors")
    named = []
    for lo, hi in [(-4, 4), (-8, 8), (-14, 14), (-24, 24), (0, 6), (-2, 4)]:
        named.append((f"uniform s in [{lo},{hi}]", *C.uniform_prior(lo, hi, 14001)))
    for mu, sd in [(0.0, 1.0), (0.0, 2.0), (0.0, 4.0), (1.0, 2.0), (2.0, 2.0)]:
        named.append((f"gaussian mu={mu} sd={sd}", *C.gaussian_prior(mu, sd, 14001)))
    # uniform in beta on (0, B]:  density prop. to e^s
    for B in (10.0, 100.0, 1000.0):
        ss = np.linspace(-14.0, math.log(B), 14001)
        ww = np.exp(ss)
        ww /= ww.sum()
        named.append((f"uniform in beta on (0,{B:g}]", ss, ww))
    # the pool's own sigma distribution
    sig = np.array([C.sigma_of(a) for a in pool])
    ss = np.linspace(sig.min() - 1e-6, sig.max() + 1e-6, 4001)
    ww = np.histogram(sig, bins=np.r_[ss - (ss[1] - ss[0]) / 2, ss[-1] + (ss[1] - ss[0]) / 2])[0].astype(float)
    ww /= ww.sum()
    named.append(("empirical sigma of the pool", ss, ww))
    for name, sp, wp in named:
        Up, _, _ = C.u_matrix(pool, sp)
        psi = Up @ wp
        a_ = agreement_from_psi(psi, A)
        say(f"    {name:32s} {a_:.5f}")
        rows.append(dict(block="named_prior", prior=name, agreement=a_))

    # ---- (c3) optimised prior -----------------------------------------------
    say()
    say("(c3) optimised prior on a 561-atom grid  s in [-6, 12]")
    sg = np.linspace(-6.0, 12.0, 561)
    Ug, _, _ = C.u_matrix(pool, sg)
    G = sg.size
    sgn = np.sign(A)
    mask = (np.abs(A) > C.TIE).astype(float)

    def count(rho):
        psi = Ug @ rho
        d = psi[:, None] - psi[None, :]
        ok = mask * (np.sign(d) == sgn) * (np.abs(d) > 0)
        return float(ok[iu].sum()) / npairs

    def mirror(rho, tau, steps, lr):
        for _ in range(steps):
            psi = Ug @ rho
            d = psi[:, None] - psi[None, :]
            m = sgn * d / tau
            g = -sgn * mask / (1.0 + np.exp(m))  # d softplus(-m)/d d
            gpsi = g.sum(axis=1) - g.sum(axis=0)
            grad = gpsi @ Ug
            grad = grad - grad.mean()
            sc = np.abs(grad).max()
            if sc == 0:
                break
            rho = rho * np.exp(-lr * grad / sc)
            rho /= rho.sum()
        return rho

    best_rho, best = np.ones(G) / G, None
    best = count(best_rho)
    for tau in (3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        cand = mirror(best_rho.copy(), tau, 400, 0.5)
        if count(cand) > best:
            best_rho, best = cand, count(cand)
    say(f"    after mirror descent: {best:.5f}")
    # Frank-Wolfe polish directly on the 0-1 count
    improved = True
    sweeps = 0
    while improved and sweeps < 8:
        improved = False
        sweeps += 1
        for k in range(G):
            for al in (0.5, 0.25, 0.1, 0.05, 0.02, 0.01):
                cand = (1 - al) * best_rho
                cand[k] += al
                c = count(cand)
                if c > best + 1e-12:
                    best_rho, best, improved = cand, c, True
        say(f"    sweep {sweeps}: {best:.5f}")
    supp = np.argsort(best_rho)[::-1][:8]
    say(f"    optimised prior agreement: {best:.5f}")
    say("    top atoms (s, weight): "
        + ", ".join(f"({sg[k]:+.2f}, {best_rho[k]:.3f})" for k in supp if best_rho[k] > 1e-3))
    rows.append(dict(block="optimised", agreement=best,
                     support=str([(round(float(sg[k]), 3), round(float(best_rho[k]), 4))
                                  for k in supp if best_rho[k] > 1e-3])))
    np.save("m2_best_prior.npy", np.vstack([sg, best_rho]))

    # ---- benchmarks ----------------------------------------------------------
    say()
    say("=" * 78)
    say("benchmarks on the same pool")
    say("=" * 78)
    a_opt = agreement_from_psi(hod["psi"], A)
    say(f"    1/2 log phi   (= endpoint prior)      {a_end:.5f}")
    say(f"    best single temperature delta_s       {agr[kbest]:.5f}")
    say(f"    best optimised Psi_rho                {best:.5f}")
    say(f"    HodgeRank psi_opt (not of Psi_rho form) {a_opt:.5f}")
    say(f"    ceiling from 42 edge-disjoint cycles   {1 - 42 / npairs:.5f}")
    rows.append(dict(block="summary", half_log_phi=a_end, best_delta_s=float(agr[kbest]),
                     best_optimised=best, psi_opt=a_opt, ceiling=1 - 42 / npairs))

    with open("m2_priors.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    say("\nwrote m2_priors.csv, m2_best_prior.npy")


if __name__ == "__main__":
    main()
