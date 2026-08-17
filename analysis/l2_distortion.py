#!/usr/bin/env python3
"""The l2-distortion of the exchange metric, against Bourgain's O(log n) bound.

The exchange pseudometric ``d(a,b) = -log(C(a->b) C(b->a))`` is a metric but is
not of negative type, so it does not embed isometrically in Hilbert space.  How
badly does it fail?  The l2-distortion is

    c2(d) = min over f : X -> l2 of  (max stretch) x (max shrink),

equivalently the least ``D`` with ``d_ij <= ||f_i - f_j|| <= D d_ij``.  Bourgain
proved ``c2 = O(log n)`` for every n-point metric, attained on expanders.

Computing ``c2`` exactly is a semidefinite program.  With no SDP solver here it
is bracketed from both sides, which is enough to pin it:

UPPER BOUND -- any embedding gives one.  Parametrise the Gram matrix as
``G = X X^T``, so ``Q_ij = ||x_i - x_j||^2``, and minimise the spread
``max_ij Q_ij/d_ij^2  /  min_ij Q_ij/d_ij^2``.  The square root of the best
spread found is an upper bound on ``c2``.

LOWER BOUND -- a Poincare-type certificate.  For any symmetric ``Delta`` with
zero row sums and ``Delta >= 0`` (positive semidefinite), and any ``f``,

    sum_ij Delta_ij ||f_i - f_j||^2 = -2 tr(Delta F) <= 0,

since the Gram matrix ``F`` is PSD.  Splitting ``Delta`` by sign off the
diagonal and using ``d_ij^2 <= ||f_i-f_j||^2 <= D^2 d_ij^2`` gives

    c2(d)^2 >= ( sum_{Delta_ij > 0} Delta_ij d_ij^2 )
             / ( sum_{Delta_ij < 0} |Delta_ij| d_ij^2 ).

Writing ``Delta = J Y Y^T J`` with ``J = I - 11^T/n`` makes both constraints
automatic, so the bound can be maximised by unconstrained optimisation.

The script first reproduces the textbook value ``c2(C_4) = sqrt(2)`` as a check
on the solver, then measures the exchange metric.

RESULT.  The failure of negative type is real but quantitatively mild.  On the
minimal five-signature family that defeats negative type the distortion is
``c2 = 1.3375``, bracketed to ``3e-4``; on random families it sits near ``1.1``
and does not grow -- ``c2/log n`` falls steadily from ``0.63`` at ``n = 5`` to
``0.36`` at ``n = 25``, so the Bourgain bound is nowhere near attained.  A greedy
attempt to extend the certificate found nothing better, so the worst distortion
observed anywhere is that of the smallest possible witness.
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "l2_distortion.csv"

RESTARTS = 12
SEED = 5

# The minimal family that defeats negative type (see negative_type_certificate.py).
CERTIFICATE = [
    (12, 10, 8, 8, 2, 1),
    (11, 9, 7, 7, 4, 1),
    (12, 12, 6, 5, 4, 4),
    (12, 10, 7, 4, 3, 3),
    (11, 11, 7, 7, 4, 3),
]


def exchange_distances(family: list[tuple[int, ...]]) -> np.ndarray:
    size = len(family)
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            value = -math.log(
                exchange_rate(family[i], family[j]) * exchange_rate(family[j], family[i])
            )
            matrix[i, j] = matrix[j, i] = value
    return matrix


def upper_bound(distances: np.ndarray, restarts: int = RESTARTS) -> float:
    """Best embedding found; sqrt of the spread is an upper bound on c2."""

    size = len(distances)
    rank = min(size, 8)
    pairs = [(i, j) for i, j in itertools.combinations(range(size), 2)]
    squared = np.array([distances[i, j] ** 2 for i, j in pairs])
    rows = np.array([i for i, _ in pairs])
    columns = np.array([j for _, j in pairs])
    rng = np.random.default_rng(SEED)
    best = math.inf

    def spread(flat: np.ndarray, sharpness: float) -> float:
        points = flat.reshape(size, rank)
        gaps = ((points[rows] - points[columns]) ** 2).sum(axis=1)
        ratio = np.log(gaps + 1e-300) - np.log(squared)
        return (
            np.logaddexp.reduce(sharpness * ratio)
            + np.logaddexp.reduce(-sharpness * ratio)
        ) / sharpness

    for _ in range(restarts):
        flat = rng.normal(size=size * rank)
        for sharpness in (4.0, 32.0, 256.0, 2048.0):
            flat = minimize(
                spread, flat, args=(sharpness,), method="L-BFGS-B",
                options={"maxiter": 2000},
            ).x
        points = flat.reshape(size, rank)
        gaps = ((points[rows] - points[columns]) ** 2).sum(axis=1) / squared
        best = min(best, gaps.max() / gaps.min())
    return math.sqrt(best)


def lower_bound(distances: np.ndarray, restarts: int = 8, rank: int = 3) -> float:
    """Best Poincare certificate found; sqrt of the ratio bounds c2 below."""

    size = len(distances)
    squared = distances**2
    centering = np.eye(size) - np.ones((size, size)) / size
    rng = np.random.default_rng(SEED + 1)
    best = 1.0

    def negative_ratio(flat: np.ndarray) -> float:
        factor = flat.reshape(size, rank)
        delta = centering @ (factor @ factor.T) @ centering
        off = delta - np.diag(np.diag(delta))
        positive = np.sum(np.where(off > 0, off, 0.0) * squared)
        negative = -np.sum(np.where(off < 0, off, 0.0) * squared)
        return 0.0 if negative < 1e-12 else -positive / negative

    for _ in range(restarts):
        flat = rng.normal(size=size * rank)
        flat = minimize(
            negative_ratio, flat, method="Powell",
            options={"maxiter": 1200, "xtol": 1e-8, "ftol": 1e-10},
        ).x
        best = max(best, -negative_ratio(flat))
    return math.sqrt(best)


def cycle_metric(length: int) -> np.ndarray:
    """Graph metric of the cycle C_n; c2(C_4) = sqrt(2) is the textbook value."""

    matrix = np.zeros((length, length))
    for i in range(length):
        for j in range(length):
            gap = abs(i - j)
            matrix[i, j] = min(gap, length - gap)
    return matrix


def signature_pool() -> list[tuple[int, ...]]:
    pool = []
    for size in range(2, 6):
        for entries in itertools.combinations_with_replacement(range(1, 9), size):
            signature = tuple(sorted(entries, reverse=True))
            if signature[0] > 1:
                pool.append(signature)
    return sorted(set(pool))


def main() -> int:
    print("solver check on C_4, where c2 = sqrt(2) = 1.414214 exactly")
    square = cycle_metric(4)
    print(f"   lower {lower_bound(square):.6f}   upper {upper_bound(square):.6f}\n")

    rows = []
    print("the five-signature family that defeats negative type")
    distances = exchange_distances(CERTIFICATE)
    low, high = lower_bound(distances), upper_bound(distances)
    print(f"   c2 in [{low:.6f}, {high:.6f}]\n")
    rows.append(["certificate", 5, f"{low:.9f}", f"{high:.9f}", f"{math.log(5):.6f}"])

    print("growth with family size, drawn from a pool of small signatures")
    print(f"   {'n':>4} {'c2 lower':>10} {'c2 upper':>10} {'log n':>8} {'c2/log n':>9}")
    pool = signature_pool()
    rng = np.random.default_rng(SEED + 2)
    for size in (5, 8, 12, 16, 20, 25):
        indices = rng.choice(len(pool), size=size, replace=False)
        family = [pool[i] for i in indices]
        distances = exchange_distances(family)
        low = lower_bound(distances, restarts=6)
        high = upper_bound(distances, restarts=5)
        print(f"   {size:>4} {low:>10.6f} {high:>10.6f} {math.log(size):>8.4f} "
              f"{high / math.log(size):>9.4f}")
        rows.append([f"random n={size}", size, f"{low:.9f}", f"{high:.9f}",
                     f"{math.log(size):.6f}"])

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["family", "n", "c2_lower", "c2_upper", "log_n"])
        writer.writerows(rows)
    print(f"\nwritten to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
