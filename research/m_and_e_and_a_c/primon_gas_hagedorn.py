#!/usr/bin/env python3
"""What breaks at beta = 1, in resource language: the primon gas is Hagedorn.

Section 3 of `exchange_positivity_and_weil.md` observes that the tensor product
of the prime ladders is the resource ``{1,2,3,...}`` with partition function
``zeta``, and then flags a caveat: the Euler product converges only for
``Re s > 1``, so the critical strip is exactly the region where the monotone
DIVERGES, and no bridge from the operational theory to RH crosses that gap.

The caveat is not a defect of the framework.  This resource is the primon gas
(B. Julia, 1990), whose non-commutative completion is the Bost-Connes
C*-dynamical system, and ``beta = 1`` is its phase transition.  The transition is
of Hagedorn type, and that has an exact meaning in the language of the first
paper, where the exchange rate is the supporting-line construction on the Gibbs
entropy-energy curve.

The thermodynamics is closed-form.  With energy levels ``E_n = log n``:

    Z(beta)  = sum_n n^-beta            = zeta(beta)
    U(beta)  = -zeta'(beta)/zeta(beta)  = sum_n Lambda(n) n^-beta
    S(beta)  = log zeta(beta) + beta U(beta)

so the mean energy of the primon gas is the von Mangoldt Dirichlet series -- the
same object the Weil pairing produces in T1.4, arriving here as an ENERGY rather
than as a spectral remainder.  As ``beta -> 1+``, ``zeta(beta) ~ 1/(beta-1)``
gives

    U ~ 1/(beta-1),   S ~ 1/(beta-1),   S/U -> 1,

which is the Hagedorn signature: the entropy-energy curve is asymptotically the
straight line ``S = beta_H U`` with ``beta_H = 1``.  In the resource language the
supporting line of slope ``1`` touches the curve only at infinity, so:

  * for ``beta > 1`` the Gibbs state exists, the contact point is finite, and the
    rate against any other resource is a well-defined infimum;
  * at ``beta = 1`` there is no contact point and no normalisable Gibbs state --
    infinite entropy at finite energy density.  The infimum defining ``C`` has no
    minimiser, and the resource cannot be quoted;
  * this is exactly why ``C(xi -> P) = C(P -> xi) = 0`` unrestricted, and why the
    companion paper had to impose a temperature window.

So the honest statement is not "the monotone diverges, and that is where RH
lives".  It is: **the zeta resource has a Hagedorn temperature at beta = 1, the
window the companion paper imposes is the physical response to it, and the
theory of what replaces the Gibbs state below it is the KMS classification of
Bost-Connes.**  That is the direction, and it is a programme, not a session.

    python research/m_and_e_and_a_c/primon_gas_hagedorn.py
"""
from __future__ import annotations

import math

import mpmath as mp

mp.mp.dps = 30


def von_mangoldt(n: int) -> float:
    m, d = n, 2
    while d * d <= m:
        if m % d == 0:
            while m % d == 0:
                m //= d
            return math.log(d) if m == 1 else 0.0
        d += 1
    return math.log(n) if n >= 2 else 0.0


def thermodynamics(beta):
    beta = mp.mpf(beta)
    z = mp.zeta(beta)
    u = -mp.diff(mp.zeta, beta) / z
    return z, u, mp.log(z) + beta * u


def main() -> None:
    print("primon gas: states = positive integers, E_n = log n\n")

    print("1. the partition function is zeta and the energy is von Mangoldt")
    terms = 200_000
    for beta in (1.5, 2.0, 3.0):
        z, u, _ = thermodynamics(beta)
        z_sum = sum(mp.mpf(n) ** (-beta) for n in range(1, terms + 1))
        u_sum = sum(mp.mpf(von_mangoldt(n)) * mp.mpf(n) ** (-beta)
                    for n in range(2, terms + 1))
        tail = mp.mpf(terms) ** (1 - beta) / (beta - 1)
        print(f"   beta = {float(beta):<4}  |zeta - sum n^-b| = {float(abs(z - z_sum)):.3e}"
              f"  (tail bound {float(tail):.1e})"
              f"   |(-zeta'/zeta) - sum Lambda(n) n^-b| = {float(abs(u - u_sum)):.3e}")

    print("\n2. the Hagedorn divergence at beta -> 1+:  S/U -> beta_H = 1")
    print(f"   {'beta':>12} {'U':>16} {'S':>16} {'S/U':>12} {'U*(beta-1)':>12}")
    for eps in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8):
        beta = 1 + mp.mpf(eps)
        _, u, s = thermodynamics(beta)
        print(f"   {float(beta):>12.8f} {float(u):>16.6f} {float(s):>16.6f}"
              f" {float(s / u):>12.8f} {float(u * mp.mpf(eps)):>12.8f}")

    print("\n   U*(beta-1) -> 1 and S/U -> 1: the entropy-energy curve approaches the")
    print("   straight line S = U, so the supporting line of slope 1 has no finite")
    print("   contact point.  beta = 1 is a Hagedorn temperature, not a coordinate")
    print("   singularity.")

    print("\n3. what a finite resource does instead, for contrast")
    print("   a truncated Euler factor P(p,K) = {1, p, ..., p^K}:")
    for p, k in ((2, 8), (2, 128), (7, 8)):
        zk = sum(mp.mpf(p) ** (-mp.mpf(1.0) * j) for j in range(k + 1))
        u_k = (sum(j * math.log(p) * mp.mpf(p) ** (-mp.mpf(1.0) * j)
                   for j in range(k + 1)) / zk)
        print(f"     p={p:<3} K={k:<4} at beta = 1:  Z = {float(zk):.6f}"
              f"   U = {float(u_k):.6f}   finite, so the contact point exists")
    print("   Every finite portfolio has a finite Hagedorn-free spectrum; only the")
    print("   infinite tensor product acquires the transition.  This is the same")
    print("   sentence as 'the zeros of zeta are a phenomenon of the infinite tensor")
    print("   product', now with a temperature attached to it.")


if __name__ == "__main__":
    main()
