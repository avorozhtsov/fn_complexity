#!/usr/bin/env python3
"""Step 4 of session brief H: what the KMS classification actually replaces the
Gibbs state with, and whether the Galois symmetry can do any resource-theoretic
work.

FIRST, A CORRECTION TO THE BRIEF.  Session brief H states the Bost-Connes
phase transition backwards.  It says: "for beta > 1 there is a unique KMS_beta
state ... at beta <= 1 the Gibbs state does not exist; KMS_beta states still do,
they form a simplex, and for beta <= 1 the symmetry group Zhat* acts on them --
spontaneous symmetry breaking".  The theorem (Bost-Connes 1995, Thm 5(a,b);
Connes-Marcolli Thm 3.1) is the other way round, in the direction that
spontaneous symmetry breaking always goes -- symmetry breaks at LOW
temperature:

  * for 0 < beta <= 1 there is a UNIQUE KMS_beta state.  It is a factor state
    of type III_1.  The Zhat* action fixes it, so the symmetry is unbroken.
  * for beta > 1 the extremal KMS_beta states are a free transitive Zhat*-torsor,
    indexed by embeddings of Q^ab in C.  They are type I_infinity and their
    partition function is zeta(beta).

So the Galois-permuted FAMILY of states exists exactly where the exchange
framework already lives (beta > 1), and at beta = 1/2 there is one state and
nothing to compare it with.  The brief's hope -- "if the extremal KMS states
below beta = 1 give a family of comparisons permuted by Gal(Q^ab/Q)" -- is
excluded by the uniqueness half of the theorem, not by any numerical failure.

WHAT IS COMPUTED HERE.

1. The extremal KMS_s states for s > 1, restricted to the group algebra
   C*(Z/q) < C*(Q/Z), are the residue distributions of the zeta measure:

       P_j^{(g)}(s) = zeta(s)^{-1} sum_{n : g n = j mod q} n^{-s}
                    = q^{-s} zeta(s, {g^{-1} j}/q) / zeta(s),

   with g in (Z/q)* the Galois label.  Verified against direct summation.

2. **The Galois action is a relabelling.**  P^{(g)} is the pushforward of
   P^{(1)} by multiplication by g, hence a PERMUTATION of the same probability
   vector.  Every extremal KMS_s state therefore has the same partition
   function, the same energy spectrum, the same entropy -- and every exchange
   monotone Z_a(beta) = sum_i a_i^beta is a function of the multiset of levels
   alone.  A label-blind functional cannot see a relabelling.  Measured.

3. The Zhat*-average of the extremal states is the Ramanujan-sum formula

       F_s(q) = (1/phi(q)) sum_{d | q} d^{1-s} mu(q/d)
              = prod_{p | q} p^{-(k_p - 1)s} (p^{1-s} - 1)/(p - 1),

   in which zeta has cancelled.  This is Bost-Connes' unique KMS_s state for
   0 < s <= 1; at s = infinity it is mu(q)/phi(q).  It is checked to be a
   genuine state at s = 1/2 (positive definite on Z/q for every q tested).

4. **The naive continuation is not a state.**  q^{-s} zeta(s, j/q)/zeta(s)
   continues to s < 1, keeps total mass 1, and acquires NEGATIVE entries.  The
   threshold is located: the first negative entry appears at 1 - s_c ~ 1/q
   (measured q(1 - s_c) -> 1), and at s = 1/2 every modulus q >= 4 already has
   one.  So "regularised Gibbs weights" on the critical strip are a signed
   measure, not a resource.

STANDING OBSTRUCTIONS, unchanged by anything here:
  * atomic measures sum_i delta_{a_i} are not admissible Weil test functions,
    so E is a finite-rank truncation and |Z_a(1/2 + i gamma)| does not decay;
  * the exchange monotone diverges in the critical strip.
The KMS route addresses the second and not the first.
"""

from __future__ import annotations

import csv
from math import gcd
from pathlib import Path

import mpmath as mp

OUT = Path(__file__).resolve().parent
mp.mp.dps = 40


def factorise(n: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    m, d = n, 2
    while d * d <= m:
        while m % d == 0:
            factors[d] = factors.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        factors[m] = factors.get(m, 0) + 1
    return factors


def mobius(n: int) -> int:
    factors = factorise(n)
    if any(e > 1 for e in factors.values()):
        return 0
    return -1 if len(factors) % 2 else 1


def euler_phi(n: int) -> int:
    result = n
    for p in factorise(n):
        result = result // p * (p - 1)
    return result


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def residue_weights(q: int, s, label: int = 1):
    """Extremal KMS_s state, restricted to C*(Z/q), as a vector over j=0..q-1.

    P_j = zeta(s)^{-1} sum_{n >= 1, g n == j (mod q)} n^{-s}, continued in s by
    the Hurwitz zeta function.  ``label`` is the Galois index g in (Z/q)*.
    """
    inverse = pow(label, -1, q)
    weights = []
    for j in range(q):
        residue = (inverse * j) % q
        shift = mp.mpf(q if residue == 0 else residue) / q
        weights.append(q ** (-s) * mp.zeta(s, shift) / mp.zeta(s))
    return weights


def residue_weights_direct(q: int, s, terms: int = 400000):
    """Same vector by direct summation plus an Euler-Maclaurin tail (s > 1)."""
    partial = [mp.mpf(0)] * q
    for n in range(1, terms + 1):
        partial[n % q] += mp.mpf(n) ** (-s)
    N = mp.mpf(terms)
    tail_total = N ** (1 - s) / (s - 1) - N ** (-s) / 2 + s * N ** (-s - 1) / 12
    for j in range(q):
        partial[j] += tail_total / q
    z = mp.zeta(s)
    return [value / z for value in partial]


def symmetric_state(q: int, s):
    """F_s(q) = (1/phi(q)) sum_{d|q} d^{1-s} mu(q/d), the Zhat*-average."""
    if q == 1:
        return mp.mpf(1)
    total = mp.fsum(mp.mpf(d) ** (1 - s) * mobius(q // d) for d in divisors(q))
    return total / euler_phi(q)


def symmetric_state_product(q: int, s):
    """prod_{p|q} p^{-(k-1)s} (p^{1-s} - 1)/(p - 1), the Euler-product form."""
    value = mp.mpf(1)
    for p, k in factorise(q).items():
        value *= mp.mpf(p) ** (-(k - 1) * s) * (mp.mpf(p) ** (1 - s) - 1) / (p - 1)
    return value


def symmetric_fourier(q: int, s):
    """Fourier coefficients over Z/q of m -> F_s(q/gcd(m,q)); a state iff >= 0."""
    values = [symmetric_state(q // gcd(m, q), s) for m in range(q)]
    coefficients = []
    for j in range(q):
        total = mp.mpf(0)
        for m in range(q):
            total += values[m] * mp.cos(-2 * mp.pi * j * m / q)
        coefficients.append(total / q)
    return coefficients


def main() -> int:
    rows = []

    print("=" * 78)
    print("1. Extremal KMS_s states at s > 1 are the residue distributions of")
    print("   the zeta measure.  Hurwitz form vs direct summation.")
    print("=" * 78)
    print(f"{'q':>4} {'s':>5} {'max |Hurwitz - direct|':>26} {'sum P_j - 1':>16}")
    for q in (3, 5, 12):
        for s_value in ("1.5", "2", "3"):
            s = mp.mpf(s_value)
            a = residue_weights(q, s)
            b = residue_weights_direct(q, s)
            worst = max(abs(x - y) for x, y in zip(a, b))
            print(
                f"{q:>4} {s_value:>5} {mp.nstr(worst, 5):>26} "
                f"{mp.nstr(mp.fsum(a) - 1, 5):>16}"
            )
            rows.append(["extremal_vs_direct", q, s_value, mp.nstr(worst, 10),
                         mp.nstr(mp.fsum(a) - 1, 10)])

    print()
    print("=" * 78)
    print("2. The Zhat* action is a relabelling: P^{(g)} is a permutation of")
    print("   P^{(1)}.  All extremal KMS_s states share every label-blind")
    print("   functional, in particular every exchange monotone.")
    print("=" * 78)
    q, s = 12, mp.mpf(2)
    base = residue_weights(q, s)
    print(f"  q = {q}, s = {s}")
    print(f"  g = 1  P = [{', '.join(mp.nstr(x, 8) for x in base)}]")
    for g in sorted(x for x in range(1, q) if gcd(x, q) == 1):
        twisted = residue_weights(q, s, label=g)
        sorted_gap = max(
            abs(x - y) for x, y in zip(sorted(base, key=float), sorted(twisted, key=float))
        )
        pointwise_gap = max(abs(x - y) for x, y in zip(base, twisted))
        print(
            f"  g = {g:>2}  max|sorted P^(g) - sorted P^(1)| = {mp.nstr(sorted_gap, 4):>10}"
            f"   max|P^(g) - P^(1)| = {mp.nstr(pointwise_gap, 6)}"
        )
        rows.append(["galois_permutation", q, g, mp.nstr(sorted_gap, 10),
                     mp.nstr(pointwise_gap, 10)])
    print("  The sorted vectors agree to working precision while the unsorted")
    print("  ones differ by O(1).  Entropy, energy and partition function are")
    print("  symmetric functions of the vector, so they are Galois-invariant;")
    print("  Z_a(beta) = sum_i a_i^beta is one of them.")
    print()
    print("  Proof, one line: the Zhat* action commutes with sigma_t and fixes")
    print("  every isometry mu_n, hence fixes H, hence fixes Tr e^{-sH} = zeta(s)")
    print("  and every functional of the level multiset.  Equivalently the whole")
    print("  Galois orbit lives in the ZERO-ENERGY sector, sigma_t(e(a)) = e(a).")

    print()
    print("=" * 78)
    print("3. The Zhat*-average is the Ramanujan formula; zeta cancels out of it.")
    print("   F_s(q) = (1/phi(q)) sum_{d|q} d^{1-s} mu(q/d).")
    print("=" * 78)
    print(f"{'q':>4} {'s':>5} {'average of extremals':>24} {'F_s(q)':>24} {'diff':>10}")
    for q in (3, 4, 5, 6, 12):
        for s_value in ("1.5", "3"):
            s = mp.mpf(s_value)
            units = [g for g in range(1, q) if gcd(g, q) == 1]
            average = mp.fsum(
                mp.polylog(s, mp.e ** (2j * mp.pi * g / q)) for g in units
            ) / (len(units) * mp.zeta(s))
            closed = symmetric_state(q, s)
            product = symmetric_state_product(q, s)
            print(
                f"{q:>4} {s_value:>5} {mp.nstr(mp.re(average), 18):>24} "
                f"{mp.nstr(closed, 18):>24} {mp.nstr(abs(average - closed), 4):>10}"
            )
            assert abs(closed - product) < mp.mpf(10) ** -30
            rows.append(["symmetric_average", q, s_value, mp.nstr(mp.re(average), 12),
                         mp.nstr(closed, 12)])
    print("   (the Euler-product form agrees with the divisor sum to 1e-30 in")
    print("    every case, asserted in code)")
    print()
    print("   Limits of F_s(q):  s -> infinity gives mu(q)/phi(q); s = 1 gives 0")
    print("   for q > 1 (Haar measure on Zhat, the critical state); s = 0 gives 1.")
    print(f"{'q':>4} {'mu(q)/phi(q)':>16} {'F_20(q)':>18} {'F_1(q)':>12} {'F_0(q)':>10}")
    for q in (2, 3, 4, 5, 6, 10, 12, 30):
        print(
            f"{q:>4} {mp.nstr(mp.mpf(mobius(q)) / euler_phi(q), 10):>16} "
            f"{mp.nstr(symmetric_state(q, mp.mpf(20)), 10):>18} "
            f"{mp.nstr(symmetric_state(q, mp.mpf(1)), 4):>12} "
            f"{mp.nstr(symmetric_state(q, mp.mpf(0)), 4):>10}"
        )
        rows.append(["symmetric_limits", q, mp.nstr(mp.mpf(mobius(q)) / euler_phi(q), 12),
                     mp.nstr(symmetric_state(q, mp.mpf(20)), 12), ""])

    print()
    print("   Is F_{1/2} a state?  Its Fourier coefficients on Z/q must be >= 0.")
    print(f"{'q':>4} {'min Fourier coeff at s=1/2':>30} {'sum':>10}")
    half = mp.mpf(1) / 2
    for q in (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 30, 36, 60):
        coefficients = symmetric_fourier(q, half)
        worst = min(coefficients, key=float)
        print(f"{q:>4} {mp.nstr(worst, 16):>30} {mp.nstr(mp.fsum(coefficients), 6):>10}")
        rows.append(["unique_state_positivity", q, mp.nstr(worst, 16),
                     mp.nstr(mp.fsum(coefficients), 10), ""])
    print("   All nonnegative: the Bost-Connes KMS_{1/2} state is a genuine")
    print("   state, as the theorem says.  It is also the ONLY one.")

    print()
    print("=" * 78)
    print("4. The naive continuation of the Gibbs weights below s = 1 is NOT a")
    print("   state: mass stays 1, entries go negative.")
    print("=" * 78)
    print(f"{'q':>4} {'min_j P_j at s=1/2':>24} {'sum_j P_j':>14} "
          f"{'s_c (first sign change)':>24} {'q(1-s_c)':>10}")
    for q in (2, 3, 4, 5, 8, 12, 20, 50, 100):
        weights = residue_weights(q, half)
        worst = min(weights, key=float)
        total = mp.fsum(weights)

        def least(s_value, q=q):
            return min(residue_weights(q, mp.mpf(s_value)), key=float)

        lo, hi = mp.mpf("0.001"), mp.mpf("0.999")
        critical = mp.mpf("nan")
        if float(least(hi)) > 0 > float(least(lo)):
            for _ in range(140):
                mid = (lo + hi) / 2
                if float(least(mid)) > 0:
                    hi = mid
                else:
                    lo = mid
            critical = (lo + hi) / 2
        print(
            f"{q:>4} {mp.nstr(worst, 16):>24} {mp.nstr(total, 10):>14} "
            f"{mp.nstr(critical, 16):>24} "
            f"{mp.nstr(q * (1 - critical), 6) if critical == critical else '-':>10}"
        )
        rows.append(["continuation_negativity", q, mp.nstr(worst, 16),
                     mp.nstr(total, 12), mp.nstr(critical, 16)])
    print("   q(1 - s_c) -> 1: the first negative weight appears at")
    print("   1 - s_c ~ 1/q, so every modulus q >= 4 already fails at")
    print("   s = 1/2 (q = 3 holds out until s = 0.4313).")
    print("   The 'analytically continued Gibbs state' is a")
    print("   signed measure of total mass 1, and the resource-theoretic reading")
    print("   of a negative weight is a negative multiplicity.")

    with (OUT / "kms_states.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["kind", "q", "arg", "value1", "value2"])
        writer.writerows(rows)
    print(f"\nwritten to {(OUT / 'kms_states.csv').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
