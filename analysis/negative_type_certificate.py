#!/usr/bin/env python3
"""The minimal certificate that the exchange metric is not of negative type.

The exchange pseudometric is ``d(a,b) = -log(C(a->b) C(b->a))``.  It obeys the
triangle inequality, but it is not of NEGATIVE TYPE: there is a family and a
vector ``x`` with ``sum(x) = 0`` and ``x^T D x > 0``.  Equivalently, by
Schoenberg, ``exp(-t d)`` fails to be positive semidefinite for small ``t`` and
``sqrt(d)`` admits no isometric embedding into a Hilbert space.

The witness below has **five** signatures, and five is the floor:

* every cut semimetric satisfies ``x^T delta_S x = -2 (sum_{i in S} x_i)^2 <= 0``,
  so ``CUT_n`` is contained in ``NEG_n``;
* ``MET_4 = CUT_4`` (Deza--Laurent), so no four-point metric can violate.

Since ``CUT_5 = HYP_5``, the same family also shows that l1-embeddability breaks
at exactly five points, and the pentagonal inequality -- the first hypermetric
inequality beyond the triangle -- is violated as well.  So the exchange metric
lies in ``MET \\ HYP`` and leaves the hierarchy at the first opportunity.

The exploratory search that found the family is in ``research/m_and_e_and_a_c/``;
this script only verifies, and verifies twice: once through the package solver
and once on an independent dense beta-grid, since the violation is invisible to
any grid truncated below ``beta ~ 500``.
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate, exchange_rate_result  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "negative_type_certificate.csv"

CERTIFICATE = [
    (12, 10, 8, 8, 2, 1),
    (11, 9, 7, 7, 4, 1),
    (12, 12, 6, 5, 4, 4),
    (12, 10, 7, 4, 3, 3),
    (11, 11, 7, 7, 4, 3),
]
WITNESS = np.array([302626, -510642, -576418, 330027, 454407], dtype=float) / 1e6

# The pentagonal inequality  sum_{i<j} b_i b_j d_ij <= 0  with b = (1,1,1,-1,-1).
PENTAGON = [
    (10, 6),
    (8, 8, 1, 1, 1, 1),
    (10, 10, 6, 5, 4, 4),
    (9, 9, 2),
    (10, 5, 5, 3, 3, 1),
]
PENTAGON_WEIGHTS = np.array([1, 1, 1, -1, -1], dtype=float)

GRID_POINTS = 2_000_001
GRID_MAX = 600.0


def distance(a: tuple[int, ...], b: tuple[int, ...]) -> float:
    return -math.log(exchange_rate(a, b) * exchange_rate(b, a))


def distance_matrix(family: list[tuple[int, ...]]) -> np.ndarray:
    size = len(family)
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            matrix[i, j] = matrix[j, i] = distance(family[i], family[j])
    return matrix


def grid_distance_matrix(family: list[tuple[int, ...]]) -> np.ndarray:
    """Independent recomputation, summing from the largest entry to stay stable."""

    betas = np.concatenate([[0.0], np.linspace(1e-9, GRID_MAX, GRID_POINTS)])

    def log_partition(signature: tuple[int, ...]) -> np.ndarray:
        entries = np.asarray(signature, dtype=float)
        top = entries.max()
        return betas * math.log(top) + np.log(
            np.sum((entries[:, None] / top) ** betas[None, :], axis=0)
        )

    profiles = {s: log_partition(s) for s in family}

    def rate(source: tuple[int, ...], target: tuple[int, ...]) -> float:
        ratio = profiles[source] / profiles[target]
        ratio[0] = math.log(len(source)) / math.log(len(target))
        return min(ratio.min(), math.log(max(source)) / math.log(max(target)))

    size = len(family)
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            value = -(math.log(rate(family[i], family[j]))
                      + math.log(rate(family[j], family[i])))
            matrix[i, j] = matrix[j, i] = value
    return matrix


def label(signature: tuple[int, ...]) -> str:
    return "{" + ",".join(map(str, signature)) + "}"


def main() -> int:
    matrix = distance_matrix(CERTIFICATE)
    quadratic = float(WITNESS @ matrix @ WITNESS)

    print("the five-signature certificate")
    for signature, weight in zip(CERTIFICATE, WITNESS):
        print(f"   {label(signature):<24} x = {weight:+.6f}")
    print(f"\n   sum of x        {WITNESS.sum():+.3e}   (zero)")
    print(f"   x^T D x         {quadratic:+.15e}   (positive: NOT of negative type)")

    centering = np.eye(5) - np.ones((5, 5)) / 5
    smallest = np.linalg.eigvalsh(-0.5 * centering @ matrix @ centering).min()
    print(f"   min eig -1/2JDJ {smallest:+.6e}")

    slack = min(
        matrix[i, k] + matrix[k, j] - matrix[i, j]
        for i, j, k in itertools.permutations(range(5), 3)
    )
    print(f"   triangle slack  {slack:+.6e}   (non-negative: a genuine metric)")

    contacts = []
    for i, j in itertools.combinations(range(5), 2):
        beta = exchange_rate_result(
            implemented=CERTIFICATE[j], implementer=CERTIFICATE[i]
        ).beta
        contacts.append("0" if beta == 0 else "inf" if math.isinf(beta) else f"{beta:.2f}")
    print(f"   contacts        {', '.join(contacts)}")

    print("\nindependent recomputation on a dense beta-grid")
    grid = grid_distance_matrix(CERTIFICATE)
    print(f"   grid            {GRID_POINTS} points on [0, {GRID_MAX}] plus beta = infinity")
    print(f"   max |D - D_grid| {abs(matrix - grid).max():.3e}")
    print(f"   x^T D_grid x     {float(WITNESS @ grid @ WITNESS):+.15e}")

    print("\nminimality: MET_4 = CUT_4 and CUT_n is contained in NEG_n,")
    print("so no four-point metric can violate. Our own 4-subsets:")
    worst = 0.0
    for indices in itertools.combinations(range(5), 4):
        sub = matrix[np.ix_(indices, indices)]
        centre = np.eye(4) - np.ones((4, 4)) / 4
        worst = min(worst, np.linalg.eigvalsh(-0.5 * centre @ sub @ centre).min())
    print(f"   worst 4-subset eigenvalue {worst:+.3e}   (non-positive, as forced)")

    print("\npentagonal inequality, b = (1,1,1,-1,-1)")
    pentagon = distance_matrix(PENTAGON)
    value = 0.0
    for i, j in itertools.combinations(range(5), 2):
        value += PENTAGON_WEIGHTS[i] * PENTAGON_WEIGHTS[j] * pentagon[i, j]
    for signature, weight in zip(PENTAGON, PENTAGON_WEIGHTS):
        print(f"   {label(signature):<24} b = {weight:+.0f}")
    print(f"   sum b_i b_j d_ij  {value:+.12e}   (positive: hypermetric fails)")
    print("   CUT_5 = HYP_5, so l1-embeddability breaks at five points too.")

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["family", "i", "j", "signature_i", "signature_j", "distance"])
        for name, fam, mat in (("certificate", CERTIFICATE, matrix),
                               ("pentagon", PENTAGON, pentagon)):
            for i, j in itertools.combinations(range(5), 2):
                writer.writerow(
                    [name, i, j, label(fam[i]), label(fam[j]), f"{mat[i, j]:.15f}"]
                )
    print(f"\nwritten to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
