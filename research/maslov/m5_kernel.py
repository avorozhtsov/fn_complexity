#!/usr/bin/env python3
"""Brief M item 5 -- transition kernels, literally.

**The fixed-beta object.**  At inverse temperature ``beta`` a signature ``a`` is
a Gibbs system on its fibers with energies ``eps_i = -log a_i``:

    pi^a_beta(i) = a_i^beta / Z_a(beta),
    E_a(beta) = sum_i pi^a_beta(i) log a_i,
    S_a(beta) = -sum_i pi^a_beta(i) log pi^a_beta(i),
    F_a(beta) = log Z_a(beta) = S_a(beta) + beta E_a(beta).

The transition kernel from ``a``'s fibers to ``b``'s fibers that is compatible
with both Gibbs states and carries no extra information is the source-free
(maximum-entropy) transport plan ``K^{a->b}_beta(i, j) = pi^b_beta(j)``, and the
*rate* attached to it -- the number of copies of ``b`` one copy of ``a`` buys at
that single temperature -- is the ratio of free energies

    c_beta(a -> b) = F_a(beta) / F_b(beta) = exp(u_a(s) - u_b(s)),   s = log beta,

because ``f^{ok}`` fits into ``g^{on}`` asymptotically iff
``k F_f(beta) <= n F_g(beta)`` for every ``beta`` (the paper's Theorem 1).  The
framework's exchange rate is ``C(a->b) = inf_beta c_beta(a->b)``.

**Two facts, both identities.**

  (i)  ``c_beta(a->b) c_beta(b->a) = 1``:  at one temperature the trade is
       reversible, so the fixed-``beta`` comparison is the total order of the
       scalar ``F_.(beta)``.  No cycles at any single temperature.

  (ii) Replace the infimum by the ``lambda``-power mean over a temperature
       prior ``rho``:

           C_lambda(a->b) = ( INT c_beta(a->b)^{-lambda} rho )^{-1/lambda},
           A_lambda(a,b)  = (1/2) log( C_lambda(a->b) / C_lambda(b->a) ).

       This **is** the brief's family, exactly.  ``lambda = infinity`` is the
       infimum (the framework); ``lambda = 1`` is the plain arithmetic average
       of the fixed-``beta`` rates -- the user's literal question; ``lambda = 0``
       is the geometric mean.  Because the geometric mean is the unique power
       mean with ``M_0(1/X) = 1/M_0(X)``, and only for it is
       ``C_0(a->b) C_0(b->a) = 1`` -- the trade becomes reversible again, and a
       reversible exchange rate is a potential.  **That is why cycles die at
       lambda = 0, and it is also why the metric dies with them.**

  Also ``A_lambda`` is even in ``lambda``: the harmonic mean (lambda = -1) and
  the arithmetic mean (lambda = +1) induce the same comparison.

**A genuinely third construction.**  Averaging the fixed-``beta`` *verdicts*
rather than the rates gives the Condorcet rule
``V(a,b) = INT sign(u_a - u_b) rho``, which is not ``A_lambda`` for any
``lambda`` -- it sees only the sign pattern of ``f``, never its size.  It has
cycles of its own.  So the temperature-indexed total orders form a *profile* and
the framework is one aggregation rule on it; ``lambda = 0`` is the scoring
(Borda-like) rule, always transitive; majority is a third rule.

Writes ``m5_kernel.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import common as C


def gibbs(sig, s):
    """``(pi, E, S, F)`` of the Gibbs state at ``beta = e^s``."""

    a = np.array(sig, dtype=float)
    beta = math.exp(s)
    z = beta * np.log(a)
    m = z.max()
    p = np.exp(z - m)
    Z = p.sum()
    p = p / Z
    F = m + math.log(Z)  # = log Z_a(beta)
    E = float((p * np.log(a)).sum())
    S = float(-(p * np.log(p)).sum())
    return p, E, S, F


def main():
    rows, out = [], []

    def say(t=""):
        print(t)
        out.append(t)

    say("=" * 86)
    say("5a.  The Gibbs data and F = S + beta E   (an identity; shown to check the code)")
    say("=" * 86)
    worst = 0.0
    for sig in [(6, 3, 3), (7, 2, 1), (6, 5, 1), (22, 17, 16, 12, 10, 10, 8, 8, 8, 6, 4)]:
        for s in (-3.0, -0.5, 0.0, 1.0, 4.0):
            p, E, S, F = gibbs(sig, s)
            worst = max(worst, abs(F - (S + math.exp(s) * E)))
    say(f"    max |log Z - (S + beta E)| over 20 (signature, beta) = {worst:.3e}")
    rows.append(dict(block="gibbs", max_dev=worst))

    say()
    say("=" * 86)
    say("5b.  c_beta(a->b) c_beta(b->a) = 1: the fixed-beta trade is reversible")
    say("=" * 86)
    pool = C.f11_pool()[:60]
    s, w = C.uniform_prior(-14.0, 14.0, 5601)
    U, um, up = C.u_matrix(pool, s)
    dev = 0.0
    for i in range(0, 60, 7):
        for j in range(0, 60, 11):
            f = U[i] - U[j]
            dev = max(dev, float(np.abs(np.exp(f) * np.exp(-f) - 1.0).max()))
    say(f"    max |c(a->b) c(b->a) - 1| = {dev:.3e}")
    say("    => at a single temperature the comparison is the total order of F_.(beta):")
    say("       no 3-cycles at any beta, for any pool, by construction.")
    rows.append(dict(block="reversible", max_dev=dev))

    say()
    say("=" * 86)
    say("5c.  A_lambda is the lambda-power mean of the fixed-beta rates")
    say("=" * 86)
    say("     C_lambda(a->b) = (INT c_beta(a->b)^-lambda rho)^-1/lambda ;"
        "  A_lambda = 1/2 log(C_lambda(a->b)/C_lambda(b->a))")
    pairs = [(C.CYCLE[0], C.CYCLE[1]), (C.CYCLE[1], C.CYCLE[2]),
             (pool[3], pool[17]), (pool[40], pool[55])]
    say(f"  {'pair':>8s} {'lambda':>9s} {'A_lambda (soft mid)':>21s} {'via power means':>18s}"
        f" {'diff':>10s}")
    for k, (a, b) in enumerate(pairs):
        ua, ub = C.u_of(a, s), C.u_of(b, s)
        f = ua - ub
        for lam in (0.25, 1.0, 4.0, 100.0):
            direct = C.soft_mid(f, w, lam)
            # power means, computed the long way round from the rates themselves
            mx = float(np.max(-lam * f))
            Cab = math.exp(-(math.log(float((w * np.exp(-lam * f - mx)).sum())) + mx) / lam)
            mx2 = float(np.max(lam * f))
            Cba = math.exp(-(math.log(float((w * np.exp(lam * f - mx2)).sum())) + mx2) / lam)
            via = 0.5 * math.log(Cab / Cba)
            say(f"  {k:>8d} {lam:9.4g} {direct:21.15f} {via:18.12f} {abs(direct-via):10.2e}")
            rows.append(dict(block="power_mean", pair=k, lam=lam, direct=direct, via=via,
                             diff=abs(direct - via)))

    say()
    say("     lambda = 1 is literally the arithmetic average of the rates:")
    ok = True
    for a, b in pairs:
        f = C.u_of(a, s) - C.u_of(b, s)
        lhs = C.soft_mid(f, w, 1.0)
        rhs = float((w * np.exp(f)).sum()) - float((w * np.exp(-f)).sum())
        ok &= (lhs > 0) == (rhs > 0)
        say(f"       A_1 = {lhs:+.9f}   E[c(a->b)] - E[c(b->a)] = {rhs:+.9f}   "
            f"signs agree: {(lhs>0)==(rhs>0)}")
        rows.append(dict(block="lambda1", A1=lhs, mean_rate_gap=rhs, signs_agree=(lhs > 0) == (rhs > 0)))

    say()
    say("     evenness in lambda (harmonic vs arithmetic mean give the same comparison):")
    dev = 0.0
    for a, b in pairs:
        f = C.u_of(a, s) - C.u_of(b, s)
        for lam in (0.5, 1.0, 7.0):
            dev = max(dev, abs(C.soft_mid(f, w, lam) - C.soft_mid(f, w, -lam)))
    say(f"       max |A_lambda - A_-lambda| = {dev:.3e}")
    rows.append(dict(block="even", max_dev=dev))

    say()
    say("=" * 86)
    say("5d.  At lambda = 0 the trade is lossless: the metric collapses")
    say("=" * 86)
    say(f"  {'lambda':>10s} {'S_lambda(known cycle edge 1)':>30s} {'S_lambda/lambda':>18s}"
        f"  {'Var_rho(f)/2':>14s}")
    f = C.u_of(C.CYCLE[0], s) - C.u_of(C.CYCLE[1], s)
    var = float((w * f**2).sum() - ((w * f).sum()) ** 2)
    for lam in (1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e4):
        Sl = 0.5 * (C.softmax(f, w, lam) - C.softmin(f, w, lam))
        say(f"  {lam:10.4g} {Sl:30.12e} {Sl/lam:18.10f}  {var/2:14.10f}")
        rows.append(dict(block="lossless", lam=lam, S_lambda=Sl, ratio=Sl / lam, half_var=var / 2))
    say(f"    S_lambda -> (lambda/2) Var_rho(f) as lambda -> 0   [proved: "
        f"S_lambda = (K(lam)+K(-lam))/(2 lam)]")
    say(f"    S_infinity = osc(f)/2 = {0.5*(f.max()-f.min()):.9f}")

    say()
    say("=" * 86)
    say("5e.  Averaging the VERDICTS instead: the Condorcet rule, a third construction")
    say("=" * 86)
    pool = C.f11_pool()
    n = len(pool)
    U, um, up = C.u_matrix(pool, s)
    A_inf = C.tropical_matrix(U, um, up)
    V = np.empty((n, n))
    for i in range(n):
        D = U[i][None, :] - U
        V[i] = (w * np.sign(D)).sum(axis=1)
    V = 0.5 * (V - V.T)
    np.fill_diagonal(V, 0.0)
    say(f"    Condorcet V on the F_11 pool: {C.three_cycles(V)} directed 3-cycles"
        f"  (tropical A_inf has {C.three_cycles(A_inf)})")
    iu = np.triu_indices(n, 1)
    dis = int((np.sign(V)[iu] != np.sign(A_inf)[iu]).sum())
    say(f"    V and A_inf disagree on {dis} of {iu[0].size} pairs "
        f"({dis/iu[0].size:.4%})")
    hv = C.hodge(V)
    say(f"    V curl energy {hv['curl_energy']:.4e}   (A_inf: {C.hodge(A_inf)['curl_energy']:.4e})")
    rows.append(dict(block="condorcet", cycles_V=C.three_cycles(V),
                     cycles_Ainf=C.three_cycles(A_inf), disagreements=dis,
                     curl_energy_V=hv["curl_energy"]))

    # V is not A_lambda for any lambda: find two pairs with equal V, opposite A_lambda
    say()
    say("    V is not A_lambda for any lambda -- V depends on f only through sign(f):")
    Vv, Av = V[iu], None
    found = None
    for lam in (0.5, 2.0, 10.0, 1e3, None):
        Al = C.soft_matrix(U, w, lam)
        av = Al[iu]
        order = np.argsort(Vv)
        for t in range(order.size - 1):
            p, q = order[t], order[t + 1]
            if abs(Vv[p] - Vv[q]) < 1e-12 and np.sign(av[p]) != np.sign(av[q]) \
               and min(abs(av[p]), abs(av[q])) > 1e-6:
                found = (lam, float(Vv[p]), float(av[p]), float(av[q]))
                break
        tag = "inf" if lam is None else f"{lam:g}"
        say(f"       lambda={tag:>6s}: equal-V pairs with opposite A_lambda: "
            f"{'yes ' + str(found[1:]) if found else 'none in this scan'}")
        rows.append(dict(block="V_not_A", lam=tag, witness=str(found)))
        found = None

    with open("m5_kernel.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    say("\nwrote m5_kernel.csv")


if __name__ == "__main__":
    main()
