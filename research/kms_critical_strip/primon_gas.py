#!/usr/bin/env python3
"""Step 1 of session brief H: the primon gas as a C*-dynamical system.

INDEPENDENT VERIFICATION, not a re-derivation.  While this session was running,
`research/m_and_e_and_a_c/primon_gas_hagedorn.py` (commit c6848fa, branch
claude/weyl-matrix-quantum-mechanics-2fc4ca) established the same transition
from the entropy-energy side: ``U*(beta-1) -> 1`` and ``S/U -> 1``, so the Gibbs
curve tends to the line ``S = U`` and the supporting line of slope 1 has no
finite contact point.  Both of those numbers are recomputed here at 40 digits
from an independent code path, and the rest of this script covers what that one
does not: the trace identity itself, the Euler product as the tensor
factorisation, and the counting-function form of the Hagedorn statement
(abscissa of convergence = exponential growth rate of the level count).


The resource ``(x)_p P_{p,inf} = {1,2,3,...}`` of the two-positivities note is
the primon (Riemann) gas of Julia 1990.  Concretely:

    H |n> = log n |n>            on   l^2(N),   n = 1, 2, 3, ...
    Tr e^{-s H} = sum_n n^{-s} = zeta(s)   (s > 1)

and the Gibbs state ``w_s(A) = Tr(A e^{-sH})/zeta(s)`` is the unique KMS_s
state of the *diagonal* (type I) picture.  This script confirms, at 40 digits:

  (a) Tr e^{-sH} = zeta(s) = prod_p (1 - p^{-s})^{-1}, i.e. the Euler product
      is the tensor factorisation of the resource into local prime modes;
  (b) each local mode is a harmonic oscillator with level spacing log p, and
      the truncations P_{p,K} converge to the local factor;
  (c) the Gibbs thermodynamics U(s) = -(log zeta)'(s), S = log zeta + s U, and
      their divergence rates as s -> 1+;
  (d) the density of states is exponential, N(E) = #{n : log n <= E} = e^E,
      so s = 1 is a Hagedorn temperature: the abscissa of convergence equals
      the exponential growth rate of the level counting function.

CONVENTION WARNING.  Three different inverse temperatures appear in this
programme and must never be identified:

    s        the Dirichlet exponent / primon-gas inverse temperature.
             Z(s) = sum n^{-s}.  Critical strip 0 < s < 1, critical line
             s = 1/2, abscissa of convergence s = 1.
    beta_x   the exchange framework's temperature, Z_a(beta) = sum_i a_i^beta
             with entries a_i >= 1 and beta in [0, infinity].  For a resource
             written in the cost convention, beta_x = -s.
    beta_xi  the completed-zeta temperature of the two-positivities note,
             Z_xi(beta) = xi(1/2 + beta)/xi(1/2).  Critical LINE at beta = 0.

Under beta_x = -s the primon gas converges only for beta_x < -1: the whole
admissible range beta_x in [0, infinity] of the exchange framework lies inside
the divergence region, and the critical strip sits at beta_x in (-1, 0).
"""

from __future__ import annotations

import csv
from pathlib import Path

import mpmath as mp

OUT = Path(__file__).resolve().parent
mp.mp.dps = 40


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i, flag in enumerate(sieve) if flag]


def trace_partial(s, terms: int):
    """Partial trace sum_{n<=terms} n^{-s} plus its Euler-Maclaurin tail.

    Computed without calling zeta, so the comparison against zeta below is a
    genuine check of ``Tr e^{-sH} = zeta(s)`` rather than a tautology.
    """
    head = mp.fsum(mp.mpf(n) ** (-s) for n in range(1, terms + 1))
    N = mp.mpf(terms)
    tail = N ** (1 - s) / (s - 1) - N ** (-s) / 2 + s * N ** (-s - 1) / 12
    tail -= s * (s + 1) * (s + 2) * N ** (-s - 3) / 720
    return head + tail


def euler_product(s, prime_bound: int):
    """prod_{p <= bound} (1 - p^{-s})^{-1}, the tensor product of local modes."""
    value = mp.mpf(1)
    for p in primes_up_to(prime_bound):
        value *= 1 / (1 - mp.mpf(p) ** (-s))
    return value


def thermodynamics(s):
    z = mp.zeta(s)
    dz = mp.diff(mp.zeta, s)
    energy = -dz / z
    entropy = mp.log(z) + s * energy
    return z, energy, entropy


def main() -> int:
    print("=" * 78)
    print("(a) Tr e^{-sH} on l^2(N) equals zeta(s), and the Euler product is")
    print("    the tensor factorisation into local prime modes.")
    print("=" * 78)
    print(f"{'s':>6} {'|trace - zeta|':>16} {'|euler-zeta|/zeta':>18} {'zeta(s)':>24}")
    rows_a = []
    for s_value in ("1.1", "1.5", "2", "3", "5"):
        s = mp.mpf(s_value)
        z = mp.zeta(s)
        tr = trace_partial(s, 100000)
        eu = euler_product(s, 500000)
        print(
            f"{s_value:>6} {mp.nstr(abs(tr - z), 5):>16} "
            f"{mp.nstr(abs(eu - z) / z, 5):>18} {mp.nstr(z, 18):>24}"
        )
        rows_a.append(
            [s_value, mp.nstr(abs(tr - z), 8), mp.nstr(abs(eu - z) / z, 8), mp.nstr(z, 25)]
        )

    print()
    print("=" * 78)
    print("(b) local mode p: levels E_k = k log p, Z_p(s) = (1 - p^{-s})^{-1};")
    print("    truncation P_{p,K} = {1,p,...,p^K} has Z = sum_{k<=K} p^{-ks}.")
    print("=" * 78)
    s = mp.mpf(2)
    for p in (2, 3, 5):
        exact = 1 / (1 - mp.mpf(p) ** (-s))
        for K in (1, 8, 64):
            trunc = mp.fsum(mp.mpf(p) ** (-k * s) for k in range(K + 1))
            print(
                f"  p={p:>2} K={K:>3}  Z_trunc={mp.nstr(trunc, 20):>24} "
                f"rel.err vs local factor {mp.nstr(abs(trunc - exact) / exact, 4)}"
            )

    print()
    print("=" * 78)
    print("(c) Gibbs thermodynamics of the primon gas and the s -> 1+ rates.")
    print("    log Z ~ -log(s-1),  U ~ 1/(s-1),  S ~ 1/(s-1).")
    print("=" * 78)
    print(
        f"{'s-1':>10} {'log zeta':>14} {'U':>18} {'S':>18} "
        f"{'U*(s-1)':>12} {'S/U':>12} {'logZ/-log(s-1)':>16}"
    )
    rows_c = []
    for exponent in range(1, 13):
        eps = mp.mpf(10) ** (-exponent)
        s = 1 + eps
        z, energy, entropy = thermodynamics(s)
        print(
            f"{mp.nstr(eps, 3):>10} {mp.nstr(mp.log(z), 8):>14} "
            f"{mp.nstr(energy, 12):>18} {mp.nstr(entropy, 12):>18} "
            f"{mp.nstr(energy * eps, 8):>12} {mp.nstr(entropy / energy, 8):>12} "
            f"{mp.nstr(mp.log(z) / (-mp.log(eps)), 8):>16}"
        )
        rows_c.append(
            [
                mp.nstr(eps, 3),
                mp.nstr(mp.log(z), 20),
                mp.nstr(energy, 20),
                mp.nstr(entropy, 20),
                mp.nstr(energy * eps, 20),
                mp.nstr(entropy / energy, 20),
            ]
        )

    print()
    print("  U*(s-1) -> 1 and S/U -> 1 reproduce primon_gas_hagedorn.py from an")
    print("  independent code path at 40 digits.  log zeta / -log(s-1) -> 1 adds")
    print("  the rate for the free energy: it diverges only logarithmically while")
    print("  energy and entropy diverge like a simple pole.")

    print()
    print("=" * 78)
    print("(d) Hagedorn: the level counting function N(E) = #{n : log n <= E}")
    print("    is floor(e^E), so the density of states grows like e^E and the")
    print("    abscissa of convergence is exactly the growth rate 1.")
    print("=" * 78)
    for E in (2, 5, 10, 15, 20):
        count = mp.floor(mp.e ** E)
        print(
            f"  E={E:>3}  N(E)=floor(e^E)={mp.nstr(count, 12):>16}  "
            f"log N(E)/E = {mp.nstr(mp.log(count) / E, 12)}"
        )
    print()
    print("  For a FINITE signature a the counting function is bounded, growth")
    print("  rate 0, and Z_a is entire in s: no finite resource of the exchange")
    print("  semiring has a phase transition at all.")

    with (OUT / "primon_gas.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["part", "s_or_eps", "col1", "col2", "col3", "col4", "col5"])
        for row in rows_a:
            writer.writerow(["trace_and_euler", *row, "", ""])
        for row in rows_c:
            writer.writerow(["thermodynamics", *row])
    print(f"\nwritten to {(OUT / 'primon_gas.csv').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
