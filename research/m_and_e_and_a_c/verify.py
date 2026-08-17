#!/usr/bin/env python3
"""Independent verification of the certificates in T1_1_negative_type_minimality.md.

Every distance is computed three (optionally four) ways:

* ``exchange_rate`` from ``fn_complexity`` (golden-section refinement of the
  infimum that defines the rate);
* from scratch on a dense uniform beta grid, 2,000,001 points on [0, 600] plus
  the ``beta = infinity`` endpoint ``log max a / log max b``, using
  ``d(a,b) = osc_beta(log log Z_a - log log Z_b)``;
* the same on [0, 2000], to rule out truncation of the large-beta tail;
* with ``mpmath`` at 40 working digits (a coarse high-precision scan followed by
  golden-section refinement of every local basin).

Then the certificates are evaluated: ``x^T D x`` for the negative-type witnesses
(``sum(x) = 0``, including an exactly rational witness) and
``sum_{i<j} b_i b_j d_ij`` for the hypermetric witnesses.
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from exchange_geometry import (
    dense_distance_matrix,
    exact_distance_matrix,
    negative_type_defect,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from fn_complexity import exchange_rate_result  # noqa: E402

# --- certificates ---------------------------------------------------------

# Headline: five signatures, the smallest possible number (see the note).
NEGATIVE_TYPE_FIVE = [
    (12, 10, 8, 8, 2, 1),
    (11, 9, 7, 7, 4, 1),
    (12, 12, 6, 5, 4, 4),
    (12, 10, 7, 4, 3, 3),
    (11, 11, 7, 7, 4, 3),
]
FIVE_WITNESS = [
    0.302625832,
    -0.510642368,
    -0.576418122,
    0.330027051,
    0.454407607,
]

# An independent second five-point family, found from a different seed.
NEGATIVE_TYPE_FIVE_ALT = [
    (12, 10, 4, 4, 4, 4),
    (12, 9, 8, 8, 1, 1),
    (12, 12, 5, 5, 4, 3),
    (11, 9, 7, 6, 3, 1),
    (11, 11, 7, 5, 5, 2),
]

NEGATIVE_TYPE_SIX = [
    (10, 7, 5, 5, 4),
    (12, 9, 7, 5, 2, 1),
    (12, 12, 5, 5, 5, 2),
    (12, 9, 7, 7, 6, 5),
    (11, 8, 8, 5, 3, 1),
    (11, 11, 7, 5, 4, 3),
]

PENTAGONAL_FIVE = [
    (10, 6),
    (8, 8, 1, 1, 1, 1),
    (10, 10, 6, 5, 4, 4),
    (9, 9, 2),
    (10, 5, 5, 3, 3, 1),
]
PENTAGONAL_WEIGHTS = (1, 1, 1, -1, -1)


# --- high precision -------------------------------------------------------


def mp_rate(implementer, implemented, digits=40):
    """``C(implementer -> implemented) = inf_beta log Z_implementer / log Z_implemented``."""

    from mpmath import mp

    mp.dps = digits
    g = [mp.mpf(v) for v in implementer]
    f = [mp.mpf(v) for v in implemented]

    def ratio(beta):
        return mp.log(mp.fsum(v**beta for v in g)) / mp.log(mp.fsum(v**beta for v in f))

    best = min(
        mp.log(len(g)) / mp.log(len(f)),  # beta = 0
        mp.log(max(g)) / mp.log(max(f)),  # beta = infinity
    )
    grid = [mp.mpf(10) ** (mp.mpf(-8) + mp.mpf(13) * k / 3000) for k in range(3001)]
    values = [ratio(beta) for beta in grid]
    inverse_phi = (mp.sqrt(5) - 1) / 2
    for index in range(1, len(grid) - 1):
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]:
            left, right = grid[index - 1], grid[index + 1]
            x1 = right - inverse_phi * (right - left)
            x2 = left + inverse_phi * (right - left)
            y1, y2 = ratio(x1), ratio(x2)
            for _ in range(200):
                if y1 <= y2:
                    right, x2, y2 = x2, x1, y1
                    x1 = right - inverse_phi * (right - left)
                    y1 = ratio(x1)
                else:
                    left, x1, y1 = x1, x2, y2
                    x2 = left + inverse_phi * (right - left)
                    y2 = ratio(x2)
            best = min(best, y1, y2)
    best = min(best, min(values))
    return best


def mp_distance_matrix(family, digits=40):
    from mpmath import mp

    mp.dps = digits
    size = len(family)
    matrix = [[mp.mpf(0) for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(i + 1, size):
            forward = mp_rate(family[i], family[j], digits)
            backward = mp_rate(family[j], family[i], digits)
            value = -mp.log(forward * backward)
            matrix[i][j] = matrix[j][i] = value
    return matrix


def mp_quadratic(matrix, weights, digits=40):
    from mpmath import mp

    mp.dps = digits
    x = [mp.mpf(w) if not isinstance(w, Fraction) else mp.mpf(w.numerator) / w.denominator
         for w in weights]
    total = mp.mpf(0)
    for i in range(len(x)):
        for j in range(len(x)):
            total += x[i] * x[j] * matrix[i][j]
    return total


# --- reporting ------------------------------------------------------------


def report(name, family):
    print(f"\n=== {name}: {len(family)} signatures ===")
    for signature in family:
        print("   ", signature)
    exact = exact_distance_matrix(family)
    dense600 = dense_distance_matrix(family, points=2_000_001, beta_max=600.0)
    dense2000 = dense_distance_matrix(family, points=2_000_001, beta_max=2000.0)
    print(f"  max |exact - dense[0,600]|  = {np.abs(exact - dense600).max():.3e}")
    print(f"  max |exact - dense[0,2000]| = {np.abs(exact - dense2000).max():.3e}")
    print("  distance matrix (exact solver):")
    for row in exact:
        print("     ", "  ".join(f"{value:.9f}" for value in row))
    return exact, dense600, dense2000


def optimal_betas(family):
    print("  beta attaining each infimum (inf = the largest-fibre endpoint):")
    for i, a in enumerate(family):
        for j, b in enumerate(family):
            if i >= j:
                continue
            forward = exchange_rate_result(a, b)
            backward = exchange_rate_result(b, a)
            print(
                f"     {i}->{j}: rate {forward.rate:.12f} at beta {forward.beta:.6g}"
                f"    {j}->{i}: rate {backward.rate:.12f} at beta {backward.beta:.6g}"
            )


def triangle_slack(distances):
    slack = math.inf
    size = distances.shape[0]
    for i in range(size):
        for j in range(size):
            for k in range(size):
                if len({i, j, k}) < 3:
                    continue
                slack = min(slack, distances[i, k] + distances[k, j] - distances[i, j])
    return slack


def rationalize(vector, digits=6):
    """Round the witness to a decimal vector whose entries sum to exactly zero."""

    scaled = [round(value * 10**digits) for value in vector]
    scaled[-1] -= sum(scaled)
    return [Fraction(value, 10**digits) for value in scaled]


def negative_type_report(name, family, witness=None, high_precision=False):
    exact, dense600, dense2000 = report(name, family)
    for label, matrix in (
        ("exact      ", exact),
        ("dense[0,600] ", dense600),
        ("dense[0,2000]", dense2000),
    ):
        value, vector = negative_type_defect(matrix)
        print(f"  [{label}] max x^T D x over centred unit x = {value:+.12e}")
        if label.strip() == "exact":
            optimal = vector
    if witness is None:
        witness = optimal
    witness = np.asarray(witness, dtype=float)
    print("  witness x =", np.round(witness, 9).tolist())
    print(f"  sum(x) = {witness.sum():+.3e},  |x| = {np.linalg.norm(witness):.9f}")
    for label, matrix in (
        ("exact      ", exact),
        ("dense[0,600] ", dense600),
        ("dense[0,2000]", dense2000),
    ):
        print(f"  [{label}] x^T D x at the stated witness = {witness @ matrix @ witness:+.12e}")
    rational = rationalize(witness / np.linalg.norm(witness))
    print("  exactly rational witness (sum = 0):", [str(v) for v in rational])
    numeric = np.array([float(v) for v in rational])
    print(f"    sum = {sum(rational)}   x^T D x (exact solver) = {numeric @ exact @ numeric:+.12e}")
    print(f"  triangle slack (min over distinct triples) = {triangle_slack(exact):+.6e}")
    if high_precision:
        try:
            matrix = mp_distance_matrix(family)
            print("  mpmath (40 digits) distance matrix:")
            for row in matrix:
                print("     ", "  ".join(str(value)[:22] for value in row))
            print(f"    max |mpmath - exact| = "
                  f"{max(abs(float(matrix[i][j]) - exact[i, j]) for i in range(len(family)) for j in range(len(family))):.3e}")
            print(f"    x^T D x at the rational witness = {mp_quadratic(matrix, rational)}")
        except ImportError:
            print("  mpmath not available")
    return exact


def main() -> int:
    exact = negative_type_report(
        "NEGATIVE TYPE, FIVE POINTS (headline)",
        NEGATIVE_TYPE_FIVE,
        FIVE_WITNESS,
        high_precision=True,
    )
    optimal_betas(NEGATIVE_TYPE_FIVE)
    print("  defect of every 4-point subfamily (all must be <= 0: MET_4 = CUT_4):")
    for drop in range(5):
        keep = [i for i in range(5) if i != drop]
        value, _ = negative_type_defect(exact[np.ix_(keep, keep)])
        print(f"     drop {NEGATIVE_TYPE_FIVE[drop]}: {value:+.6e}")
    best, pattern = -math.inf, None
    for negatives in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
        b = np.ones(5)
        b[list(negatives)] = -1.0
        value = 0.5 * float(b @ exact @ b)
        if value > best:
            best, pattern = value, b.astype(int).tolist()
    print(f"  best pentagonal value on the same family = {best:+.6e} at b = {pattern}")

    negative_type_report("negative type, five points (second family)", NEGATIVE_TYPE_FIVE_ALT)
    negative_type_report("negative type, six points (larger margin)", NEGATIVE_TYPE_SIX)

    exact, dense600, dense2000 = report("PENTAGONAL, FIVE POINTS", PENTAGONAL_FIVE)
    b = np.array(PENTAGONAL_WEIGHTS, dtype=float)
    for label, matrix in (
        ("exact      ", exact),
        ("dense[0,600] ", dense600),
        ("dense[0,2000]", dense2000),
    ):
        print(f"  [{label}] sum_(i<j) b_i b_j d_ij = {0.5 * b @ matrix @ b:+.12e}")
    print("  b =", PENTAGONAL_WEIGHTS, "(sum = 1)")
    print(f"  negative-type defect of the same family = {negative_type_defect(exact)[0]:+.6e}")
    print(f"  triangle slack (min over distinct triples) = {triangle_slack(exact):+.6e}")
    try:
        matrix = mp_distance_matrix(PENTAGONAL_FIVE)
        print(f"  mpmath (40 digits) pentagonal value = {mp_quadratic(matrix, [Fraction(v) for v in PENTAGONAL_WEIGHTS]) / 2}")
        print(f"    max |mpmath - exact| = "
              f"{max(abs(float(matrix[i][j]) - exact[i, j]) for i in range(5) for j in range(5)):.3e}")
    except ImportError:
        print("  mpmath not available")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
