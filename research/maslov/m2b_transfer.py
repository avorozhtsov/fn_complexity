#!/usr/bin/env python3
"""Brief M item 2, continued -- does the optimised prior transfer?

``m2_priors.py`` fits a prior to the ``F_11`` pool and evaluates it there.  That
is in-sample and proves nothing about "which ``rho`` is best".  Here the prior
is fitted on one pool and scored on another (``F_11`` <-> ``F_13``), alongside
the fixed, unfitted priors.

Benchmarks (brief E, reproduced independently in ``m2_priors.py``):
``psi_opt`` 98.314 % / 98.651 %, ``1/2 log phi`` 82.494 % / 83.316 % at
``q = 11 / 13``.

Writes ``m2b_transfer.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import common as C

GRID_LO, GRID_HI, GRID_N = -6.0, 12.0, 561
DS = 0.005
LO, HI = -14.0, 14.0


def build(q):
    pool = C._pencil_pool(q)
    ng = int(round((HI - LO) / DS)) + 1
    s, w = C.uniform_prior(LO, HI, ng)
    U, um, up = C.u_matrix(pool, s)
    A = C.tropical_matrix(U, um, up)
    sg = np.linspace(GRID_LO, GRID_HI, GRID_N)
    Ug, _, _ = C.u_matrix(pool, sg)
    return dict(q=q, pool=pool, n=len(pool), A=A, Ug=Ug, sg=sg,
                psi_opt=C.hodge(A)["psi"],
                psi_phi=np.array([C.psi_endpoint(a) for a in pool]),
                hod=C.hodge(A))


def counter(P):
    n = P["n"]
    iu = np.triu_indices(n, 1)
    sgn = np.sign(P["A"])
    mask = (np.abs(P["A"]) > C.TIE).astype(float)
    npairs = iu[0].size

    def count(rho):
        psi = P["Ug"] @ rho
        d = psi[:, None] - psi[None, :]
        ok = mask * (np.sign(d) == sgn) * (np.abs(d) > 0)
        return float(ok[iu].sum()) / npairs

    return count, sgn, mask, iu, npairs


def fit_prior(P, seed_rho=None):
    count, sgn, mask, iu, npairs = counter(P)
    Ug = P["Ug"]
    G = Ug.shape[1]

    def mirror(rho, tau, steps, lr):
        for _ in range(steps):
            psi = Ug @ rho
            d = psi[:, None] - psi[None, :]
            m = np.clip(sgn * d / tau, -60.0, 60.0)
            g = -sgn * mask / (1.0 + np.exp(m))
            gpsi = g.sum(axis=1) - g.sum(axis=0)
            grad = gpsi @ Ug
            grad -= grad.mean()
            sc = np.abs(grad).max()
            if sc == 0:
                break
            rho = rho * np.exp(-lr * grad / sc)
            rho /= rho.sum()
        return rho

    best = np.ones(G) / G if seed_rho is None else seed_rho.copy()
    bv = count(best)
    for tau in (3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        cand = mirror(best.copy(), tau, 400, 0.5)
        if count(cand) > bv:
            best, bv = cand, count(cand)
    for _ in range(6):
        improved = False
        for k in range(G):
            for al in (0.5, 0.25, 0.1, 0.05, 0.02, 0.01):
                cand = (1 - al) * best
                cand[k] += al
                c = count(cand)
                if c > bv + 1e-12:
                    best, bv, improved = cand, c, True
        if not improved:
            break
    return best, bv


def main():
    rows, out = [], []

    def say(t=""):
        print(t)
        out.append(t)

    P11, P13 = build(11), build(13)
    for P in (P11, P13):
        say(f"F_{P['q']}: {P['n']} signatures, grad energy {P['hod']['grad_energy']:.6f}, "
            f"curl energy {P['hod']['curl_energy']:.6f}, 3-cycles {C.three_cycles(P['A'])}")

    named = {}
    for lo, hi in [(-8, 8), (0, 6), (-2, 4)]:
        named[f"uniform[{lo},{hi}]"] = C.uniform_prior(lo, hi, 14001)
    named["gaussian mu=2 sd=2"] = C.gaussian_prior(2.0, 2.0, 14001)
    named["gaussian mu=0 sd=4"] = C.gaussian_prior(0.0, 4.0, 14001)

    say()
    say("=" * 84)
    say("fixed (unfitted) priors, scored on both pools")
    say("=" * 84)
    say(f"  {'prior':28s} {'F_11':>10s} {'F_13':>10s}")
    for name, (sp, wp) in named.items():
        vals = []
        for P in (P11, P13):
            Up, _, _ = C.u_matrix(P["pool"], sp)
            vals.append(C.order_agreement(Up @ wp, P["A"]))
        say(f"  {name:28s} {vals[0]:10.5f} {vals[1]:10.5f}")
        rows.append(dict(block="named", prior=name, f11=vals[0], f13=vals[1]))
    for name, key in (("1/2 log phi (canonical)", "psi_phi"), ("HodgeRank psi_opt", "psi_opt")):
        vals = [C.order_agreement(P[key], P["A"]) for P in (P11, P13)]
        say(f"  {name:28s} {vals[0]:10.5f} {vals[1]:10.5f}")
        rows.append(dict(block="benchmark", prior=name, f11=vals[0], f13=vals[1]))

    say()
    say("=" * 84)
    say("fitted priors -- in sample and transferred")
    say("=" * 84)
    r11, v11 = fit_prior(P11)
    r13, v13 = fit_prior(P13)
    c11, c13 = counter(P11)[0], counter(P13)[0]
    say(f"  fitted on F_11:  in sample {v11:.5f}   transferred to F_13 {c13(r11):.5f}")
    say(f"  fitted on F_13:  in sample {v13:.5f}   transferred to F_11 {c11(r13):.5f}")
    rows.append(dict(block="fitted", prior="fit on F_11", f11=v11, f13=c13(r11)))
    rows.append(dict(block="fitted", prior="fit on F_13", f11=c11(r13), f13=v13))
    for tag, rho in (("F_11", r11), ("F_13", r13)):
        sg = P11["sg"]
        top = np.argsort(rho)[::-1][:6]
        say(f"  support of the fit on {tag}: "
            + ", ".join(f"({sg[k]:+.2f}, {rho[k]:.3f})" for k in top if rho[k] > 1e-3))
        rows.append(dict(block="support", prior=f"fit on {tag}",
                         support=str([(round(float(sg[k]), 3), round(float(rho[k]), 4))
                                      for k in top if rho[k] > 1e-3])))
    np.save("m2b_priors.npy", np.vstack([P11["sg"], r11, r13]))

    say()
    say("  ceilings (>= k edge-disjoint 3-cycles must be misordered):")
    say(f"    F_11: {1 - 42/np.triu_indices(P11['n'],1)[0].size:.5f}   "
        f"F_13: {1 - 214/np.triu_indices(P13['n'],1)[0].size:.5f}")

    with open("m2b_transfer.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    say("\nwrote m2b_transfer.csv, m2b_priors.npy")


if __name__ == "__main__":
    main()
