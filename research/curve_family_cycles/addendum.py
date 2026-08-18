#!/usr/bin/env python3
"""The checks addendum 1 to brief B asks for, and what they say.

The addendum derives, for a pair of fibrations of ``A^2`` over ``F_q``,

    C(f->g) = 1 + min(-0.08578644 * dm2, dM) / (q log q),      dm2 > 0,

with ``dM = M_f - M_g``, ``M = max_c(-a_c)``, and ``dm2 = m2(f) - m2(g)``, and
concludes in its section 3 that the comparison is therefore governed by the
scalar

    phi~(f) = M_f - 0.08578644 * m2(f),

"still a total order, so it still forbids cycles".  Its section 6 asks for that
to be verified numerically before it is trusted, and its section 4(c) names the
crossover region ``beta ~ sqrt(q)`` as the one place a scalar prediction cannot
exist, with the instruction: histogram the argmin ``beta`` of every computed
rate, and "if the histogram has mass anywhere other than {0} u [0,1] u {inf},
that is the finding".

Five parts, in the addendum's own order.

    A.  phi~ against the certified F_11 cycle          (section 3)
    B.  phi~ as an order on large pools                (section 6, second check)
    C.  the rate formula against measured rates        (section 6, first check)
    D.  the histogram of contact temperatures          (section 4(c))
    E.  the multiplicity term                          (section 4(b))
    F.  is phi~ right on any subpopulation?            (it is not, provably)

    python research/curve_family_cycles/addendum.py
"""

from __future__ import annotations

import collections
import csv
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from common import Engine, beta_grid, build_pool, superelliptic_counts  # noqa: E402
from search import exhaustive_pool, three_cycles  # noqa: E402
from fn_complexity import exchange_rate_result  # noqa: E402

PSI_MIN = (3.0 - 2.0 * math.sqrt(2.0)) / 2.0  # 0.08578644, |psi_2| at beta* = sqrt2 - 1
TOLERANCE = 1e-9

A = (18, 16, 15, 15, 14, 12, 9, 6, 6, 5, 5)
B = (18, 18, 14, 13, 12, 9, 9, 9, 8, 7, 4)
C = (19, 14, 12, 11, 11, 10, 10, 10, 9, 9, 6)


def largest_deviation(signature: tuple[int, ...], q: int) -> int:
    return max(signature) - q


def second_moment(signature: tuple[int, ...], q: int) -> float:
    return sum((q - n) ** 2 for n in signature) / q**2


def phi_tilde(signature: tuple[int, ...], q: int) -> float:
    return largest_deviation(signature, q) - PSI_MIN * second_moment(signature, q)


def part_a() -> None:
    print("\n" + "=" * 78)
    print("A --- phi~ against the certified F_11 cycle  (addendum section 3)")
    print("=" * 78)
    q = 11
    names = {A: "A", B: "B", C: "C"}
    print(f"\n  {'':4}{'M':>4}{'m2':>12}{'phi~':>14}")
    for signature in (A, B, C):
        print(
            f"  {names[signature]:4}{largest_deviation(signature, q):>4}"
            f"{second_moment(signature, q):>12.6f}{phi_tilde(signature, q):>14.8f}"
        )
    order = sorted((A, B, C), key=lambda s: phi_tilde(s, q))
    print(f"\n  phi~ orders them {' < '.join(names[s] for s in order)}, a chain.")
    print("  The computed comparison is a cycle:\n")
    failures = 0
    for first, second in ((A, B), (B, C), (C, A)):
        forward = exchange_rate_result(implementer=first, implemented=second, grid_size=16384).rate
        backward = exchange_rate_result(implementer=second, implemented=first, grid_size=16384).rate
        predicted = phi_tilde(first, q) < phi_tilde(second, q)
        agrees = predicted == (forward < backward)
        failures += not agrees
        print(
            f"    {names[first]} < {names[second]}  computed "
            f"(margin {backward - forward:+.4e});  phi~ says "
            f"{names[first] if predicted else names[second]} first  "
            f"[{'agrees' if agrees else 'DISAGREES'}]"
        )
    print(
        f"\n  {failures} of the three edges contradict phi~, so phi~ is not the\n"
        "  comparison.  It fails on exactly the edge phi fails on, and by the same\n"
        "  mechanism: an interior minimum at beta = 3.83, which is neither beta* =\n"
        "  0.414 nor infinity, and which neither term of phi~ describes."
    )


def part_b(rows: list) -> None:
    print("\n" + "=" * 78)
    print("B --- phi~ as an order on large pools  (addendum section 6)")
    print("=" * 78)
    print(
        "\n  phi~ differs from phi only in breaking phi's ties, and phi's ties are\n"
        "  exactly the pairs of equal largest fiber, on which phi~ reduces to\n"
        "  'the larger m2 precedes'.  So the sharp test is that subpopulation ---\n"
        "  which is also where every cycle of the census lives.\n"
    )
    print(
        f"  {'q':>5}{'strict pairs':>14}{'phi wrong':>11}{'phi blind':>11}"
        f"{'phi~ wrong':>12}{'tie pairs':>11}{'phi~ wrong there':>18}"
    )
    for q, signatures in pools():
        rates, _ = Engine(signatures, beta_grid(q, points=20_000)).rate_matrix(chunk=64)
        strict = (rates - rates.T) < -TOLERANCE
        identity = np.eye(len(signatures), dtype=bool)
        phi = np.array([math.log(len(s)) * math.log(max(s)) for s in signatures])
        tilde = np.array([phi_tilde(s, q) for s in signatures])
        tops = np.array([max(s) for s in signatures])
        moments = np.array([second_moment(s, q) for s in signatures])
        gap_phi = phi[:, None] - phi[None, :]
        gap_tilde = tilde[:, None] - tilde[None, :]
        # pairs phi cannot separate but phi~ does: equal largest fiber, distinct m2
        decided = (
            (tops[:, None] == tops[None, :])
            & (np.abs(moments[:, None] - moments[None, :]) > 1e-12)
            & ~identity
        )
        total = int(strict.sum())
        wrong_phi = int((strict & (gap_phi > 1e-12)).sum())
        blind_phi = int((strict & (np.abs(gap_phi) <= 1e-12)).sum())
        wrong_tilde = int((strict & (gap_tilde > 1e-12)).sum())
        on_ties = int((strict & decided).sum())
        wrong_on_ties = int((strict & decided & (gap_tilde > 1e-12)).sum())
        print(
            f"  {q:>5}{total:>14}{wrong_phi:>11}{blind_phi:>11}{wrong_tilde:>12}"
            f"{on_ties:>11}{100 * wrong_on_ties / max(on_ties, 1):>17.1f}%"
        )
        rows.append([q, len(signatures), total, wrong_phi, wrong_tilde, on_ties, wrong_on_ties])
    print(
        "\n  On the pairs it is meant to decide, phi~ is not merely inexact: it is\n"
        "  anti-predictive, wrong on about nine in ten.  The rule it encodes is\n"
        "  'larger m2 precedes'; the truth is the opposite about that often, because\n"
        "  a signature whose Psi rises faster at small tau must flatten before the\n"
        "  common endpoint, and the later excursion is usually the deeper one.\n"
        "  So phi~ is a worse order than phi overall, despite being one term deeper:\n"
        "  it converts phi's honest ties into confident errors."
    )


def pools():
    """Signature pools: exhaustive genus-two at small q, sampled at larger q."""

    for q in (11, 13):
        pool, _ = exhaustive_pool(q)
        yield q, sorted(pool)
    for q in (101,):
        pool = build_pool(
            q,
            budget=dict(hyperelliptic=1500, superelliptic=600, twist=600, additive=600, dense=600),
        )
        by_max: dict[int, list] = collections.defaultdict(list)
        for entry in pool:
            by_max[entry.max_fiber].append(entry)
        yield q, [entry.signature for members in by_max.values() for entry in members[:12]]


def part_c() -> None:
    print("\n" + "=" * 78)
    print("C --- the predicted rate formula against measured rates  (section 6)")
    print("=" * 78)
    print(
        "\n  D(f,g) = (C(f->g) - 1) * q log q, measured against the addendum's\n"
        "  min(-0.0858 dm2, dM), over pairs of 40 genus-2 pencils with dm2 > 0,\n"
        "  split by whether the largest fibers differ.\n"
    )
    print(
        f"  {'q':>6}{'population':>17}{'pairs':>7}{'median |D|':>13}"
        f"{'median rel. error':>19}{'1/sqrt(q)':>11}{'median beta*/sqrt(q)':>22}"
    )
    for q in (211, 1009):
        rng = np.random.default_rng(11)
        signatures = []
        while len(signatures) < 40:
            coefficients = [int(v) for v in rng.integers(0, q, size=5)] + [1]
            counts = superelliptic_counts(q, 2, coefficients)
            if counts.min() == 0:
                continue
            signatures.append(tuple(sorted((int(v) for v in counts), reverse=True)))
        # The package solver costs O(grid * q) Python operations per pair, which is
        # prohibitive at q = 1009; the grid engine agrees with it to 1e-13 here.
        rates, contacts = Engine(signatures, beta_grid(q, points=30_000)).rate_matrix(chunk=32)
        exact = exchange_rate_result(
            implementer=signatures[0], implemented=signatures[1], grid_size=16384
        ).rate
        tops = np.array([max(s) - q for s in signatures], dtype=float)
        moments = np.array([second_moment(s, q) for s in signatures])
        deviation = (rates - 1.0) * q * math.log(q)
        gap_top = tops[:, None] - tops[None, :]
        gap_m2 = moments[:, None] - moments[None, :]
        predicted = np.minimum(-PSI_MIN * gap_m2, gap_top)
        base = (gap_m2 > 0) & ~np.eye(len(signatures), dtype=bool)
        for label, mask in (
            ("largest differ", base & (gap_top != 0)),
            ("largest tie", base & (gap_top == 0)),
        ):
            if not mask.any():
                continue
            error = np.abs(deviation[mask] - predicted[mask]) / np.abs(deviation[mask])
            finite = contacts[mask]
            finite = finite[np.isfinite(finite) & (finite > 0)]
            print(
                f"  {q:>6}{label:>17}{int(mask.sum()):>7}"
                f"{np.median(np.abs(deviation[mask])):>13.4f}"
                f"{np.median(error):>18.1%}{1 / math.sqrt(q):>11.4f}"
                f"{float(np.median(finite)) / math.sqrt(q):>22.3f}"
            )
        print(
            f"         (grid engine against the package solver on one pair: "
            f"{abs(exact - rates[0, 1]):.1e})"
        )
    print(
        "\n  Where the largest fibers differ the formula is a fair description and its\n"
        "  error is roughly flat in q at about 10%.  Where they tie --- the only\n"
        "  pairs that can carry a cycle --- the error GROWS with q, 15% to 95%, and\n"
        "  the contact walks out from beta* to beta of order sqrt(q).  The addendum\n"
        "  predicts an O(1/sqrt(q)) relative error in both cases; on the population\n"
        "  that matters the error moves the other way, because the term it drops\n"
        "  grows like sqrt(q) against the terms it keeps."
    )


def part_d() -> None:
    print("\n" + "=" * 78)
    print("D --- histogram of the contact temperature  (section 4(c))")
    print("=" * 78)
    print(
        "\n  Contact temperature of every strict rate inside the largest class of\n"
        "  equal largest fiber, binned by the addendum's own partition.\n"
    )
    print(
        f"  {'q':>6}{'rates':>8}{'beta=0':>9}{'0<b<=1':>9}{'1<b<sqrt q':>12}"
        f"{'~sqrt q':>10}{'>3 sqrt q':>11}{'beta=inf':>10}{'median b/sqrt q':>17}"
    )
    for q in (31, 101, 211, 401):
        pool = build_pool(
            q,
            budget=dict(hyperelliptic=1500, superelliptic=600, twist=600, additive=600, dense=600),
        )
        classes: dict[int, list] = collections.defaultdict(list)
        for entry in pool:
            classes[entry.max_fiber].append(entry)
        members = max(classes.values(), key=len)[:120]
        signatures = [entry.signature for entry in members]
        rates, contacts = Engine(signatures, beta_grid(q, points=30_000)).rate_matrix(chunk=64)
        mask = ~np.eye(len(signatures), dtype=bool)
        beta = contacts[mask]
        root = math.sqrt(q)
        bins = [
            int((beta == 0.0).sum()),
            int(((beta > 0) & (beta <= 1)).sum()),
            int(((beta > 1) & (beta < root / 3)).sum()),
            int(((beta >= root / 3) & (beta <= 3 * root)).sum()),
            int(((beta > 3 * root) & np.isfinite(beta)).sum()),
            int(np.isinf(beta).sum()),
        ]
        finite = beta[np.isfinite(beta) & (beta > 0)]
        widths = (9, 9, 12, 10, 11, 10)
        cells = "".join(f"{value:>{width}}" for value, width in zip(bins, widths))
        print(
            f"  {q:>6}{beta.size:>8}{cells}"
            f"{float(np.median(finite)) / root:>17.3f}"
        )
    print(
        "\n  The mass is in the crossover column, at beta between sqrt(q)/3 and\n"
        "  3 sqrt(q), and essentially nowhere in {0} u (0,1] u {inf}.  By the\n"
        "  addendum's own criterion, that is the finding: the infimum of a pair of\n"
        "  curve families is attained where neither of its expansions is valid."
    )


def part_e() -> None:
    print("\n" + "=" * 78)
    print("E --- the multiplicity term  (section 4(b))")
    print("=" * 78)
    print(
        "\n  Two signatures identical except that the largest fiber is attained\n"
        "  twice instead of once, so log(max) ties and the multiplicity is the\n"
        "  leading endpoint datum.  Prediction: (R(beta) - 1) * beta * log(max)\n"
        "  tends to log(mu_f) - log(mu_g) = log 2 = 0.693147.\n"
    )
    q = 11
    single = (18, 16, 15, 15, 14, 12, 9, 6, 6, 5, 5)
    doubled = (18, 18, 15, 15, 14, 12, 9, 6, 6, 5, 3)  # 16 -> 18, 5 -> 3, same sum
    assert sum(single) == sum(doubled) == q * q
    print(f"    single: {single}   mu = {single.count(18)}")
    print(f"    double: {doubled}   mu = {doubled.count(18)}\n")

    def log_z(signature, beta):
        logs = np.log(np.asarray(signature, dtype=float))
        top = logs.max()
        return beta * top + math.log(float(np.exp(beta * (logs - top)).sum()))

    print(f"    {'beta':>10}{'(R-1) beta log(max)':>24}")
    top = math.log(max(single))
    for beta in (10.0, 100.0, 1000.0, 10_000.0, 100_000.0):
        first, second = log_z(doubled, beta), log_z(single, beta)
        print(f"    {beta:>10.0f}{(first / second - 1.0) * beta * top:>24.6f}")
    print(f"    {'limit':>10}{math.log(2.0):>24.6f}")
    print(
        "\n  Confirmed, and it decays like 1/beta, so at the crossover beta ~ sqrt(q)\n"
        "  it contributes ~ log 2 / sqrt(q) to (R - 1) log Z --- the same\n"
        "  1/sqrt(q) order as the crossover itself, and sqrt(q) times the moment\n"
        "  terms.  Section 4(b) is right that this is real; parts C and D say it is\n"
        "  one contribution to the crossover rather than a separate regime."
    )


def part_f() -> None:
    print("\n" + "=" * 78)
    print("F --- is there a subpopulation on which phi~ works?")
    print("=" * 78)
    print(
        "\n  Stratifying the pairs phi~ decides that phi cannot --- equal largest\n"
        "  fiber, distinct m2 --- three ways.  No stratum is good, and the third\n"
        "  stratification shows the failure is systematic rather than noisy.\n"
    )
    for q in (11, 13):
        pool, _ = exhaustive_pool(q)
        signatures = sorted(pool)
        rates, _ = Engine(signatures, beta_grid(q, points=30_000)).rate_matrix(chunk=64)
        strict = (rates - rates.T) < -TOLERANCE
        tops = np.array([max(s) - q for s in signatures])
        moments = np.array([second_moment(s, q) for s in signatures])
        mult = np.array([s.count(max(s)) for s in signatures])
        gap = moments[:, None] - moments[None, :]
        decided = (
            (tops[:, None] == tops[None, :])
            & (np.abs(gap) > 1e-12)
            & ~np.eye(len(signatures), dtype=bool)
        )
        selected = strict & decided
        correct = selected & (gap > 0)  # phi~ says the larger m2 precedes
        print(
            f"  q = {q}: {int(selected.sum())} such pairs, phi~ correct on "
            f"{100 * correct.sum() / selected.sum():.1f}%"
        )
        print(f"    {'by largest fiber':>22}{'pairs':>9}{'phi~ correct':>15}")
        for value in sorted(set(tops.tolist())):
            block = (tops[:, None] == value) & (tops[None, :] == value)
            count = int((selected & block).sum())
            if count < 20:
                continue
            print(
                f"    {'N_max = ' + str(q + value):>22}{count:>9}"
                f"{100 * (correct & block).sum() / count:>14.1f}%"
            )
        print(f"    {'by |dm2| quintile':>22}{'pairs':>9}{'phi~ correct':>15}")
        values = np.abs(gap[selected])
        edges = [0.0] + list(np.percentile(values, [20, 40, 60, 80])) + [values.max() + 1]
        for low, high in zip(edges[:-1], edges[1:]):
            block = selected & (np.abs(gap) >= low) & (np.abs(gap) < high)
            count = int(block.sum())
            if count < 20:
                continue
            print(
                f"    {f'[{low:.3f}, {high:.3f})':>22}{count:>9}"
                f"{100 * (correct & block).sum() / count:>14.1f}%"
            )
        print(f"    {'by multiplicity':>22}{'pairs':>9}{'phi~ correct':>15}")
        for label, condition in (
            ("mu equal", mult[:, None] == mult[None, :]),
            ("mu differs", mult[:, None] != mult[None, :]),
        ):
            block = selected & condition
            count = int(block.sum())
            if count < 20:
                continue
            print(
                f"    {label:>22}{count:>9}"
                f"{100 * (correct & block).sum() / count:>14.1f}%"
            )
        print()
    print(
        "  The accuracy falls monotonically as |dm2| grows, to 0.1% in the top\n"
        "  quintile.  That is the signature of a reversed sign, not of a missing\n"
        "  correction: the term phi~ keeps dominates exactly where phi~ is most\n"
        "  wrong.  The correct leading rule on these pairs is the opposite one,\n"
        "  smaller m2 precedes.\n"
    )
    print(
        "  But no rule of that shape can be right either, and this is a proof and\n"
        "  not a statistic: there are cycles inside a single largest-fiber class\n"
        "  whose three m2 are distinct.  On such a triple every scalar F(M, m2) is\n"
        "  a function of m2 alone, hence a total order, hence not the comparison.\n"
    )
    for q in (11, 13):
        pool, _ = exhaustive_pool(q)
        signatures = sorted(pool)
        rates, _ = Engine(signatures, beta_grid(q, points=30_000)).rate_matrix(chunk=64)
        tops = [max(s) - q for s in signatures]
        moments = [second_moment(s, q) for s in signatures]
        cycles = [
            entry
            for entry in three_cycles(rates)
            if len({tops[t] for t in entry[1:]}) == 1
            and len({round(moments[t], 12) for t in entry[1:]}) == 3
        ]
        print(f"  q = {q}: {len(cycles)} such cycles.  The widest:")
        margin, *triple = cycles[0]
        for position, index in enumerate(triple):
            print(
                f"      {'ABC'[position]}:  M = {tops[index]}, m2 = {moments[index]:.6f},"
                f"  sigma = {signatures[index]}"
            )
        for position in range(3):
            first, second = triple[position], triple[(position + 1) % 3]
            forward = exchange_rate_result(
                implementer=signatures[first], implemented=signatures[second], grid_size=16384
            ).rate
            backward = exchange_rate_result(
                implementer=signatures[second], implemented=signatures[first], grid_size=16384
            ).rate
            print(
                f"      {'ABC'[position]} < {'ABC'[(position + 1) % 3]}:"
                f"  margin {backward - forward:+.4e}  (package solver)"
            )
        print()


def main() -> int:
    rows: list = []
    part_a()
    part_b(rows)
    part_c()
    part_d()
    part_e()
    part_f()
    with (HERE / "addendum.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "q",
                "signatures",
                "strict_pairs",
                "phi_errors",
                "phi_tilde_errors",
                "largest_fiber_tie_pairs",
                "phi_tilde_errors_on_ties",
            ]
        )
        writer.writerows(rows)
    print("\nwritten: addendum.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
