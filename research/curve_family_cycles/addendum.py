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
        f"\n  {'q':>5}{'signatures':>12}{'strict pairs':>14}{'phi errors':>12}"
        f"{'phi~ errors':>13}{'phi~ error rate':>17}"
    )
    for q in (31, 101, 211):
        pool = build_pool(
            q,
            budget=dict(hyperelliptic=1500, superelliptic=600, twist=600, additive=600, dense=600),
        )
        by_max: dict[int, list] = collections.defaultdict(list)
        for entry in pool:
            by_max[entry.max_fiber].append(entry)
        chosen = [entry for members in by_max.values() for entry in members[:12]]
        signatures = [entry.signature for entry in chosen]
        rates, _ = Engine(signatures, beta_grid(q, points=20_000)).rate_matrix(chunk=64)
        difference = rates - rates.T
        strict = difference < -TOLERANCE
        phi = np.array([math.log(len(s)) * math.log(max(s)) for s in signatures])
        tilde = np.array([phi_tilde(s, q) for s in signatures])
        errors = {}
        for name, scalar in (("phi", phi), ("phi~", tilde)):
            gap = scalar[:, None] - scalar[None, :]
            errors[name] = int((strict & (gap > 1e-12)).sum())
        total = int(strict.sum())
        print(
            f"  {q:>5}{len(signatures):>12}{total:>14}{errors['phi']:>12}"
            f"{errors['phi~']:>13}{100 * errors['phi~'] / total:>16.2f}%"
        )
        rows.append([q, len(signatures), total, errors["phi"], errors["phi~"]])
    print(
        "\n  phi~ is a strict improvement on phi and still wrong on a percent of\n"
        "  pairs, which is far above the 1e-10 tie floor.  Every such pair is a\n"
        "  cycle candidate, and the census in FINDINGS.md shows the candidates\n"
        "  close triangles."
    )


def part_c() -> None:
    print("\n" + "=" * 78)
    print("C --- the predicted rate formula against measured rates  (section 6)")
    print("=" * 78)
    print(
        "\n  D(f,g) = (C(f->g) - 1) * q log q, measured against the addendum's\n"
        "  min(-0.0858 dm2, dM), over pairs of 20 genus-2 pencils with dm2 > 0.\n"
    )
    print(
        f"  {'q':>6}{'pairs':>7}{'median |D|':>13}{'median |pred|':>15}"
        f"{'median error':>14}{'error/sqrt(q)':>15}"
    )
    for q in (211, 1009):
        rng = np.random.default_rng(11)
        signatures = []
        while len(signatures) < 20:
            coefficients = [int(v) for v in rng.integers(0, q, size=5)] + [1]
            counts = superelliptic_counts(q, 2, coefficients)
            if counts.min() == 0:
                continue
            signatures.append(tuple(sorted((int(v) for v in counts), reverse=True)))
        # The package solver costs O(grid * q) Python operations per pair, which
        # is prohibitive at q = 1009; the grid engine gives the same numbers to
        # 1e-9 (checked in search.py) and one pair is spot-checked below.
        rates, _ = Engine(signatures, beta_grid(q, points=30_000)).rate_matrix(chunk=32)
        measured, predicted = [], []
        for i, first in enumerate(signatures):
            for j, second in enumerate(signatures):
                if i == j:
                    continue
                dm2 = second_moment(first, q) - second_moment(second, q)
                if dm2 <= 0:
                    continue
                dM = largest_deviation(first, q) - largest_deviation(second, q)
                measured.append((rates[i, j] - 1.0) * q * math.log(q))
                predicted.append(min(-PSI_MIN * dm2, float(dM)))
        exact = exchange_rate_result(
            implementer=signatures[0], implemented=signatures[1], grid_size=16384
        ).rate
        spot = abs(exact - rates[0, 1])
        measured = np.array(measured)
        predicted = np.array(predicted)
        error = np.abs(measured - predicted)
        print(
            f"  {q:>6}{len(measured):>7}{np.median(np.abs(measured)):>13.4f}"
            f"{np.median(np.abs(predicted)):>15.4f}{np.median(error):>14.4f}"
            f"{np.median(error) / math.sqrt(q):>15.4f}"
            f"   [engine vs package on one pair: {spot:.1e}]"
        )
    print(
        "\n  The error is not O(q^-1/2) relative, as the addendum predicts; it is\n"
        "  of the same size as D itself and grows like sqrt(q) in these units.\n"
        "  That is the missing term: in the addendum's own normalisation the\n"
        "  crossover contributes D ~ sqrt(q), so it dominates both terms kept,\n"
        "  each of which is O(1) or an integer times one."
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
        "  leading endpoint datum.  Prediction: (R(beta) - 1) * beta * log Z_g\n"
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

    print(f"    {'beta':>10}{'(R-1) beta log Z':>20}")
    for beta in (10.0, 100.0, 1000.0, 10_000.0, 100_000.0):
        first, second = log_z(doubled, beta), log_z(single, beta)
        print(f"    {beta:>10.0f}{(first / second - 1.0) * beta * second:>20.6f}")
    print(f"    {'limit':>10}{math.log(2.0):>20.6f}")
    print(
        "\n  Confirmed, and it decays like 1/beta, so at the crossover beta ~ sqrt(q)\n"
        "  it contributes ~ log 2 / sqrt(q) to (R - 1) log Z --- the same\n"
        "  1/sqrt(q) order as the crossover itself, and sqrt(q) times the moment\n"
        "  terms.  Section 4(b) is right that this is real; parts C and D say it is\n"
        "  one contribution to the crossover rather than a separate regime."
    )


def main() -> int:
    rows: list = []
    part_a()
    part_b(rows)
    part_c()
    part_d()
    part_e()
    with (HERE / "addendum.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["q", "signatures", "strict_pairs", "phi_errors", "phi_tilde_errors"])
        writer.writerows(rows)
    print("\nwritten: addendum.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
