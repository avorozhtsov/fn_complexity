#!/usr/bin/env python3
"""Stage 5: what a KMS comparison of resources would be, and where it stops.

Checkpoint one (`primon_gas_hagedorn.py`) established that beta = 1 is the
Hagedorn temperature of the primon gas.  This file asks the follow-up: is there a
comparison of resources built from KMS states rather than from partition
functions, and does it continue below beta = 1?

THE REDUCTION THAT MATTERS.  A signature ``a`` is a Hamiltonian: put
``H_a`` with eigenvalues ``-log a_i``, so that

    Tr exp(-beta H_a) = sum_i a_i^beta = Z_a(beta).

Every monotone of the first paper is therefore a quantity **relative to the
trace**, and the exchange rate is a ratio of two of them.  In a factor of type
III there is no trace at all -- and the Bost-Connes KMS state for beta <= 1
generates a type III_1 factor, while the ones for beta > 1 are type I_infinity.
So the beta = 1 transition is a transition in the TYPE of the algebra, and what
dies with it is not the convergence of a series but the whole class of absolute
quantities: trace, density matrix, entropy, partition function, and with them
thermomajorisation and the "second laws" of quantum thermodynamics, which are
statements about eigenvalue vectors of rho relative to a Gibbs state.

What survives in type III is everything RELATIVE: Araki's relative entropy, the
Connes cocycle, the Araki-Masuda L^p spaces and the sandwiched Renyi divergences
``D_alpha(omega || phi)``, which are defined for arbitrary von Neumann algebras
and satisfy data processing for ``alpha >= 1/2``.

    So the comparison does continue below beta = 1, but only in relative form:
    the VALUE of a resource is gone, the RATE between two of them is not.

That is the project's own methodological choice -- exchange rates are primitive,
no scalar value exists -- arriving as a structure theorem rather than a taste.

WHAT THIS SCRIPT CHECKS, all on the beta > 1 side where the algebra is type I and
everything is computable:

  1. the monotones are divergences against the trace,
     D_alpha(rho_a(beta) || tr) = (log Z_a(alpha beta) - alpha log Z_a(beta))/(alpha-1);
  2. between two primon-gas KMS states the Renyi divergence is EXACTLY the
     convexity defect of log zeta along a chord,
     D_alpha(w_b1 || w_b2) = (L(b_alpha) - alpha L(b1) - (1-alpha) L(b2))/(alpha-1),
     L = log zeta,  b_alpha = alpha b1 + (1-alpha) b2;
  3. and therefore the family is TRUNCATED by the Hagedorn wall: it is finite
     exactly while b_alpha > 1, i.e. for alpha < alpha* = (b2-1)/(b2-b1).
     **Only finitely many of the second laws survive above the transition.**

    python research/m_and_e_and_a_c/kms_comparison.py
"""
from __future__ import annotations

import math

import mpmath as mp

mp.mp.dps = 30


def log_z(signature, beta):
    return mp.log(mp.fsum(mp.power(v, beta) for v in signature))


def renyi_against_trace(signature, beta, alpha):
    """D_alpha(rho_a(beta) || tr) for the diagonal state rho_a(beta)."""
    beta, alpha = mp.mpf(beta), mp.mpf(alpha)
    z = mp.fsum(mp.power(v, beta) for v in signature)
    p = [mp.power(v, beta) / z for v in signature]
    return mp.fsum(mp.power(x, alpha) for x in p).__class__(
        mp.log(mp.fsum(mp.power(x, alpha) for x in p)) / (alpha - 1))


def primon_renyi_direct(b1, b2, alpha, terms=400_000):
    """D_alpha between two primon-gas Gibbs states, by direct summation."""
    b1, b2, alpha = mp.mpf(b1), mp.mpf(b2), mp.mpf(alpha)
    z1, z2 = mp.zeta(b1), mp.zeta(b2)
    s = mp.fsum(mp.power(n, -(alpha * b1 + (1 - alpha) * b2))
                for n in range(1, terms + 1))
    return mp.log(s / (mp.power(z1, alpha) * mp.power(z2, 1 - alpha))) / (alpha - 1)


def primon_renyi_closed(b1, b2, alpha):
    """The same, as a convexity defect of log zeta."""
    b1, b2, alpha = mp.mpf(b1), mp.mpf(b2), mp.mpf(alpha)
    chord = alpha * b1 + (1 - alpha) * b2
    if chord <= 1:
        return mp.inf
    lg = mp.log
    return (lg(mp.zeta(chord)) - alpha * lg(mp.zeta(b1))
            - (1 - alpha) * lg(mp.zeta(b2))) / (alpha - 1)


def main() -> None:
    print("1. the first paper's monotones are divergences against the TRACE")
    sig = (12, 10, 8, 8, 2, 1)
    print(f"   signature {sig}")
    for beta, alpha in ((1.0, 2.0), (0.5, 3.0), (2.0, 0.5)):
        direct = renyi_against_trace(sig, beta, alpha)
        closed = (log_z(sig, alpha * beta)
                  - mp.mpf(alpha) * log_z(sig, beta)) / (mp.mpf(alpha) - 1)
        print(f"   beta={beta:<4} alpha={alpha:<4}  D_alpha = {float(direct): .12f}"
              f"   from Z: {float(closed): .12f}   diff {float(abs(direct-closed)):.2e}")
    print("   Every Z_a(beta) is recoverable from this family.  It needs a trace;")
    print("   a type III factor has none, which is the whole obstruction.\n")

    print("2. between primon-gas KMS states, D_alpha is the convexity defect of log zeta")
    print(f"   {'b1':>5} {'b2':>5} {'alpha':>7} {'direct sum':>16} {'closed form':>16} {'diff':>10}")
    for b1, b2, alpha in ((2.0, 3.0, 0.5), (1.5, 4.0, 0.25), (3.0, 2.0, 0.75),
                          (2.5, 2.0, 2.0)):
        d = primon_renyi_direct(b1, b2, alpha)
        c = primon_renyi_closed(b1, b2, alpha)
        print(f"   {b1:>5} {b2:>5} {alpha:>7} {float(d):>16.10f} {float(c):>16.10f}"
              f" {float(abs(d - c)):>10.1e}")
    print("   (direct sums are truncated at 4e5 terms, so the residual is the tail)\n")

    print("3. the Hagedorn wall truncates the alpha-family")
    for b1, b2 in ((1.2, 3.0), (2.0, 5.0), (1.05, 2.0)):
        star = (b2 - 1) / (b2 - b1)
        print(f"   b1={b1:<5} b2={b2:<5}  alpha* = (b2-1)/(b2-b1) = {star:.6f}")
        for alpha in (star - 0.1, star - 0.001, star + 0.001, star + 0.1):
            chord = alpha * b1 + (1 - alpha) * b2
            value = primon_renyi_closed(b1, b2, alpha)
            state = "finite" if mp.isfinite(value) else "DIVERGENT"
            shown = f"{float(value):.6f}" if mp.isfinite(value) else "inf"
            print(f"      alpha={alpha:>9.6f}  chord={chord:>9.6f}  "
                  f"D_alpha = {shown:>12}   {state}")
        print()
    print("   Above the transition the second laws are available only for")
    print("   alpha < alpha*.  The Hagedorn temperature is not merely a boundary in")
    print("   beta: it truncates the family of monotones that governs conversion.")


if __name__ == "__main__":
    main()
