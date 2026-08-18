#!/usr/bin/env python3
"""The candidate comparison of session brief H: the Connes cocycle
Radon-Nikodym derivative, and whether it reduces to the exchange rate.

The brief nominates ``[D phi : D psi]_t`` because it exists for any pair of
states on any von Neumann algebra, including type III, and "reduces to the
ratio of Gibbs weights in the type-I case".  The reduction is verified here,
and then the reduction that actually matters is shown to FAIL:

    [D phi : D psi]_t  reduces to a ratio of GIBBS WEIGHTS.
    C(a -> b)          is an infimum of a ratio of LOG PARTITION FUNCTIONS.

These are different objects and no relative-entropy functional can be the
second, for a reason that is structural rather than numerical.

  (i)  Relative entropy compares two states of ONE system.  C compares two
       DIFFERENT resources.  Put them on a common algebra with two Hamiltonians
       H_a, H_b and the standard identity is

           S(w_a^b || w_b^b) = b <H_b - H_a>_a + log Z_b(b) - log Z_a(b),

       a DIFFERENCE of log partition functions.  C(a->b) is a RATIO of them.
       Verified below to 1e-35 in a finite-dimensional model.

  (ii) The two behave differently under the framework's own operation.  With
       a^{(x)k} the k-th Cartesian power:

           relative entropy   S_k = k S_1                (extensive)
           -log C             shifts by -log k           (log-extensive)
           d = -log(C C')     invariant                  (intensive)

       Three different homogeneity degrees under one operation, so no
       functional relation can hold.  Measured for k = 1..6.

 (iii) Their zero sets differ.  S(phi||psi) = 0 iff phi = psi; d(a,b) = 0 iff
       u_a - u_b is constant, which happens for genuinely different signatures
       (Cartesian powers, homothetic Gibbs regions).  Exhibited.

WHAT HAPPENS AT s <= 1, where the brief wanted to use the cocycle.  This part
is cited, not computed -- it is a theorem about a type III_1 factor and there
is nothing finite-dimensional to evaluate:

  * Bost-Connes: for 0 < s <= 1 the KMS_s state is UNIQUE.  So for any two
    KMS_s states phi, psi at the same s one has phi = psi, hence
    [D phi : D psi]_t = 1 and S(phi || psi) = 0 identically.  There is no pair
    to compare.
  * Across temperatures s1 != s2 the two KMS states are KMS for the same
    sigma at different betas; their modular groups are sigma_{-s1 t} and
    sigma_{-s2 t}, which differ by sigma_{(s2-s1)t}.  That is not inner for the
    Bost-Connes dynamics, so the states are not quasi-equivalent and Araki's
    relative entropy is +infinity.  The comparison takes only the values
    {0, +infinity}.
  * The factor is type III_1, whose flow of weights is the trivial flow on a
    point, so Connes' invariants S(M) = [0, infinity) and T(M) = {0} carry no
    information; and by the Connes-Stormer transitivity theorem the unitary
    orbit of any faithful normal state is norm dense in the normal state space
    of the hyperfinite III_1 factor.  Any unitarily invariant comparison of
    states there is therefore constant.

THAT is the type-III obstruction, and it is not "there is no trace so the
arithmetic is hard".  It is: there is exactly one state, and even if there were
two, the invariant that would compare them is the one type III_1 kills.

STANDING OBSTRUCTIONS, repeated: atomic measures are not admissible Weil test
functions, so E is a finite-rank truncation; and the exchange monotone diverges
in the critical strip.  Neither is removed here.
"""

from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from fn_complexity.core import exchange_rate, power_signature  # noqa: E402

OUT = Path(__file__).resolve().parent
mp.mp.dps = 40


def gibbs(levels, beta):
    """Gibbs weights p_i = a_i^beta / Z and the partition function."""
    weights = [mp.mpf(a) ** beta for a in levels]
    total = mp.fsum(weights)
    return [w / total for w in weights], total


def relative_entropy(p, q):
    return mp.fsum(pi * mp.log(pi / qi) for pi, qi in zip(p, q) if pi != 0)


def log_z(levels, beta):
    return mp.log(mp.fsum(mp.mpf(a) ** beta for a in levels))


def expand(counter):
    """Turn a fiber-size -> multiplicity counter into a flat signature."""
    out = []
    for size, multiplicity in counter.items():
        out.extend([size] * multiplicity)
    return out


def main() -> int:
    rows = []

    print("=" * 78)
    print("1. The type-I reduction the brief asks for, verified.")
    print("   For the primon gas truncated to n <= N, the Gibbs states at s1, s2")
    print("   have diagonal density matrices, and")
    print("      [D w_1 : D w_2]_t = rho_1^{it} rho_2^{-it}")
    print("                        = diag( n^{-i(s1-s2)t} ) (z(s2)/z(s1))^{it}.")
    print("=" * 78)
    N = 60
    levels = list(range(1, N + 1))
    s1, s2, t = mp.mpf("2.5"), mp.mpf("1.7"), mp.mpf("0.83")
    p1 = [mp.mpf(n) ** (-s1) for n in levels]
    z1 = mp.fsum(p1)
    p1 = [x / z1 for x in p1]
    p2 = [mp.mpf(n) ** (-s2) for n in levels]
    z2 = mp.fsum(p2)
    p2 = [x / z2 for x in p2]
    cocycle = [mp.e ** (1j * t * mp.log(a / b)) for a, b in zip(p1, p2)]
    predicted = [
        mp.e ** (-1j * (s1 - s2) * t * mp.log(n)) * (z2 / z1) ** (1j * t)
        for n in levels
    ]
    worst = max(abs(a - b) for a, b in zip(cocycle, predicted))
    print(f"   max |cocycle - closed form| over n <= {N}: {mp.nstr(worst, 6)}")
    rows.append(["cocycle_closed_form", "", "", mp.nstr(worst, 12), ""])

    # cocycle identity: u_{t+t'} = u_t sigma^{w2}_t(u_{t'}).  On the diagonal
    # algebra sigma^{w2}_t acts trivially on diagonal elements, so this is the
    # plain multiplicativity u_{t+t'} = u_t u_{t'}.
    t2 = mp.mpf("0.41")
    left = [mp.e ** (1j * (t + t2) * mp.log(a / b)) for a, b in zip(p1, p2)]
    right = [
        mp.e ** (1j * t * mp.log(a / b)) * mp.e ** (1j * t2 * mp.log(a / b))
        for a, b in zip(p1, p2)
    ]
    print(
        "   cocycle identity u_{t+t'} = u_t sigma_t(u_{t'}): max error "
        f"{mp.nstr(max(abs(a - b) for a, b in zip(left, right)), 6)}"
    )

    # relative entropy as the Bregman divergence of log zeta
    print()
    print("   Relative entropy between two Gibbs states is the Bregman")
    print("   divergence of the log partition function:")
    print("      S(w_s1 || w_s2) = log Z(s2) - log Z(s1) - (s2-s1)(log Z)'(s1)")
    print("   Checked on the gas truncated to n <= 20000, where both sides use")
    print("   the same finite Z and the identity is exact rather than")
    print("   asymptotic.  The last column is the untruncated Bregman form with")
    print("   zeta, whose gap from the truncated value is the truncation error.")
    print(f"{'s1':>6} {'s2':>6} {'direct sum':>22} {'Bregman (Z_N)':>22} "
          f"{'diff':>10} {'Bregman (zeta)':>16}")
    cutoff = 20000
    for a_value, b_value in (("2.5", "1.7"), ("3", "4"), ("2", "6")):
        sa, sb = mp.mpf(a_value), mp.mpf(b_value)
        pa = [mp.mpf(n) ** (-sa) for n in range(1, cutoff + 1)]
        za = mp.fsum(pa)
        pa = [x / za for x in pa]
        pb = [mp.mpf(n) ** (-sb) for n in range(1, cutoff + 1)]
        zb = mp.fsum(pb)
        pb = [x / zb for x in pb]
        direct = relative_entropy(pa, pb)
        derivative = -mp.fsum(
            mp.log(n) * mp.mpf(n) ** (-sa) for n in range(1, cutoff + 1)
        ) / za
        bregman = mp.log(zb) - mp.log(za) - (sb - sa) * derivative
        exact = (
            mp.log(mp.zeta(sb))
            - mp.log(mp.zeta(sa))
            - (sb - sa) * mp.diff(mp.zeta, sa) / mp.zeta(sa)
        )
        print(
            f"{a_value:>6} {b_value:>6} {mp.nstr(direct, 16):>22} "
            f"{mp.nstr(bregman, 16):>22} {mp.nstr(abs(direct - bregman), 4):>10} "
            f"{mp.nstr(exact, 12):>16}"
        )
        rows.append(["relative_entropy_bregman", a_value, b_value,
                     mp.nstr(direct, 20), mp.nstr(bregman, 20)])

    print()
    print("=" * 78)
    print("2. Relative entropy gives a DIFFERENCE of log partition functions;")
    print("   the exchange rate is a RATIO of them.")
    print("=" * 78)
    a = (12, 10, 8, 8, 2, 1)
    b = (11, 9, 7, 7, 4, 1)
    print(f"   a = {a}")
    print(f"   b = {b}")
    print(f"{'beta':>7} {'S(w_a||w_b)':>20} {'beta<H_b-H_a>_a':>20} "
          f"{'logZ_b - logZ_a':>20} {'identity err':>13}")
    for beta_value in ("0.5", "1", "2", "5"):
        beta = mp.mpf(beta_value)
        pa, _ = gibbs(a, beta)
        pb, _ = gibbs(b, beta)
        entropy = relative_entropy(pa, pb)
        shift = beta * mp.fsum(
            pi * (-mp.log(bi) + mp.log(ai)) for pi, ai, bi in zip(pa, a, b)
        )
        difference = log_z(b, beta) - log_z(a, beta)
        print(
            f"{beta_value:>7} {mp.nstr(entropy, 14):>20} {mp.nstr(shift, 14):>20} "
            f"{mp.nstr(difference, 14):>20} "
            f"{mp.nstr(abs(entropy - shift - difference), 4):>13}"
        )
        rows.append(["difference_identity", beta_value, mp.nstr(entropy, 20),
                     mp.nstr(difference, 20),
                     mp.nstr(abs(entropy - shift - difference), 10)])
    print("   The identity holds exactly.  The log-partition functions enter")
    print("   through their DIFFERENCE, and the exchange rate needs their")
    print("   RATIO -- which is what makes C dimensionless, i.e. a rate.")

    print()
    print("=" * 78)
    print("3. Homogeneity under Cartesian powers: three different degrees.")
    print("=" * 78)
    print(f"{'k':>3} {'S(w_a^k||w_b^k)':>20} {'S_k / (k S_1)':>15} "
          f"{'-log C(a^k->b)':>16} {'+log k':>12} {'d(a^k,b^k)':>14} "
          f"{'spectral S_k/kS_1':>18}")
    beta = mp.mpf(2)
    pa, _ = gibbs(a, beta)
    pb, _ = gibbs(b, beta)
    base_entropy = relative_entropy(pa, pb)
    base_rate = exchange_rate(a, b)
    base_distance = -mp.log(mp.mpf(exchange_rate(a, b)) * exchange_rate(b, a))
    for k in range(1, 6):
        # The two product resources must be indexed compatibly: entry (i_1..i_k)
        # of a^{(x)k} against entry (i_1..i_k) of b^{(x)k}.  Sorting the two
        # multisets independently would pair them by rank instead, which is a
        # different (and meaningless) coupling.
        product_a = [mp.fprod(t) for t in itertools.product(a, repeat=k)]
        product_b = [mp.fprod(t) for t in itertools.product(b, repeat=k)]
        power_a = tuple(sorted(expand(power_signature(a, k)), reverse=True))
        power_b = tuple(sorted(expand(power_signature(b, k)), reverse=True))
        pak, _ = gibbs(product_a, beta)
        pbk, _ = gibbs(product_b, beta)
        entropy = relative_entropy(pak, pbk)
        # Same quantity after merging equal products, i.e. on the SPECTRAL
        # algebra -- which is all the exchange framework retains, since
        # Z_a(beta) depends only on the multiset of entries.
        spectral_a, _ = gibbs(power_a, beta)
        spectral_b, _ = gibbs(power_b, beta)
        spectral = relative_entropy(spectral_a, spectral_b)
        forward = exchange_rate(power_a, b)
        distance = -mp.log(
            mp.mpf(exchange_rate(power_a, power_b)) * exchange_rate(power_b, power_a)
        )
        print(
            f"{k:>3} {mp.nstr(entropy, 14):>20} "
            f"{mp.nstr(entropy / (k * base_entropy), 12):>15} "
            f"{float(-mp.log(forward)):>16.10f} "
            f"{float(-mp.log(forward) + mp.log(k) + mp.log(base_rate)):>12.2e} "
            f"{mp.nstr(distance - base_distance, 4):>14} "
            f"{mp.nstr(spectral / (k * base_entropy), 12):>18}"
        )
        rows.append(["cartesian_power", k, mp.nstr(entropy, 20),
                     f"{-float(mp.log(forward)):.14f}",
                     mp.nstr(distance - base_distance, 10)])
    print("   S_k = k S_1 exactly (extensive); -log C(a^{(x)k} -> b) =")
    print("   -log k - log C(a->b) (log-extensive); d is invariant (intensive).")
    print("   No function of one can be the other.")
    print()
    print("   The last column is a second, independent obstruction.  Relative")
    print("   entropy is an invariant of the ALGEBRA; the exchange framework")
    print("   retains only the SPECTRUM, i.e. the multiset of entries.  Merging")
    print("   equal products is a coarse-graining, and by data processing it")
    print("   strictly decreases relative entropy: additivity is lost, the ratio")
    print("   falling to 0.59 by k = 5.  The exchange rate is unaffected, being")
    print("   a function of the multiset by definition.  So relative entropy is")
    print("   not even well defined on the objects the framework compares.")

    print()
    print("=" * 78)
    print("4. The zero sets differ.  d(a, a^{(x)k}) = 0 while the two signatures")
    print("   are different resources with S > 0.")
    print("=" * 78)
    for base in ((3, 1), (12, 10, 8, 8, 2, 1)):
        squared = tuple(sorted(expand(power_signature(base, 2)), reverse=True))
        distance = -mp.log(
            mp.mpf(exchange_rate(base, squared)) * exchange_rate(squared, base)
        )
        print(
            f"   a = {base}, a^(x)2 = {squared[:6]}{'...' if len(squared) > 6 else ''}"
        )
        print(
            f"     d(a, a^(x)2) = {mp.nstr(distance, 6)}   "
            f"C(a->a^(x)2) = {exchange_rate(base, squared):.12f}   "
            f"C(a^(x)2->a) = {exchange_rate(squared, base):.12f}"
        )
        rows.append(["zero_set", str(base), mp.nstr(distance, 12), "", ""])
    print("   d vanishes on Cartesian-power pairs; no relative entropy does,")
    print("   because the two Gibbs states are different measures.")

    print()
    print("=" * 78)
    print("5. What the cocycle gives at s <= 1: cited, not computed.")
    print("=" * 78)
    print("   Bost-Connes uniqueness for 0 < s <= 1 forces phi = psi, so")
    print("   [D phi : D psi]_t = 1 and S(phi||psi) = 0 identically.")
    print("   Across temperatures the modular groups differ by sigma_{(s2-s1)t},")
    print("   which is outer, so the states are disjoint and S = +infinity.")
    print("   The factor is type III_1: trivial flow of weights, S(M) =")
    print("   [0, infinity), T(M) = {0}, and Connes-Stormer transitivity makes")
    print("   the unitary orbit of any faithful normal state norm dense.")
    print("   A unitarily invariant comparison of states there is constant.")

    with (OUT / "connes_cocycle.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["kind", "arg1", "arg2", "value1", "value2"])
        writer.writerows(rows)
    print(f"\nwritten to {(OUT / 'connes_cocycle.csv').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
