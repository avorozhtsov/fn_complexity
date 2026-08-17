#!/usr/bin/env python3
"""Positivity experiments behind the Weil/exchange-matrix note.

Three computations:

1. The isometry identity.  With ``u_a(beta) = log log Z_a(beta)``,

       -log C(a -> b) = sup_beta (u_b - u_a),
       d(a, b) := -log (C(a -> b) C(b -> a)) = osc_beta (u_a - u_b).

   The script checks the second identity by direct sampling.

2. Failure of negative type.  ``d`` is a pseudometric (that is the triangle
   inequality, equivalently supermultiplicativity of the rates), but it is not
   of negative type: an explicit thirteen-signature family carries a vector
   ``x`` with ``sum(x) = 0`` and ``x^T D x > 0``.  By Schoenberg's theorem the
   exchange metric therefore admits no Hilbert-space embedding, and
   ``exp(-t D)`` fails to be positive semidefinite for small ``t``.

3. Conics over finite fields.  The split and anisotropic binary quadratic
   forms over ``F_q`` have signatures built from the point counts ``2q-1``,
   ``q-1`` and ``q+1``, so their exchange rates are functions of Weil numbers:

       C(anisotropic -> split) = log(q+1) / log(2q-1)              (exact),
       C(split -> anisotropic) = 1 - kappa / (q log q) + O(q^-2),

   with ``kappa = max_{beta>0} (2 beta - 2^beta) / (beta + 1)``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CERTIFICATE_CSV = OUTPUT_DIRECTORY / "exchange_negative_type_certificate.csv"
CONIC_CSV = OUTPUT_DIRECTORY / "finite_field_conic_rates.csv"

# The family found by greedy reduction from a sixty-signature pool.  Dropping
# any single member destroys the violation, and no subfamily of size at most
# twelve violates negative type.
CERTIFICATE_FAMILY = [
    (5, 5, 5, 1),
    (5, 5, 4, 2),
    (7, 7, 6, 1),
    (6, 6, 3),
    (6, 6, 2, 2),
    (5, 4, 4, 1),
    (7, 6, 6, 1),
    (7, 5, 5, 3),
    (6, 3, 2),
    (5, 4, 4),
    (6, 5, 5, 5),
    (6, 4, 2, 1),
    (5, 5, 5),
]

ISOMETRY_PAIRS = [
    ((2, 2), (3, 1)),
    ((5, 3), (3, 1, 1)),
    ((6, 5, 2, 1), (6, 4, 3, 2)),
    ((3, 3, 3), (7, 5)),
    ((9, 1), (3, 2, 2, 2)),
]

CONIC_PRIMES = [3, 5, 7, 11, 13, 17, 23, 31, 41, 61, 101, 151, 251, 509]


def label(signature: tuple[int, ...]) -> str:
    return "{" + ",".join(map(str, signature)) + "}"


def log_partition(signature: tuple[int, ...], betas: np.ndarray) -> np.ndarray:
    entries = np.asarray(signature, dtype=float)
    return np.log(np.sum(entries[:, None] ** betas[None, :], axis=0))


def oscillation(a: tuple[int, ...], b: tuple[int, ...], beta_max: float = 60.0) -> float:
    """max - min of ``u_a - u_b`` over ``[0, inf]``."""

    betas = np.linspace(1e-9, beta_max, 200_001)
    values = np.log(log_partition(a, betas)) - np.log(log_partition(b, betas))
    at_infinity = math.log(math.log(max(a))) - math.log(math.log(max(b)))
    return max(values.max(), at_infinity) - min(values.min(), at_infinity)


def irreversibility(a: tuple[int, ...], b: tuple[int, ...]) -> float:
    return -math.log(exchange_rate(a, b) * exchange_rate(b, a))


def distance_matrix(family: list[tuple[int, ...]]) -> np.ndarray:
    size = len(family)
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            matrix[i, j] = matrix[j, i] = irreversibility(family[i], family[j])
    return matrix


def negative_type_eigenvalue(matrix: np.ndarray) -> tuple[float, np.ndarray]:
    """Smallest eigenvalue of the classical multidimensional-scaling matrix.

    ``d`` is of negative type exactly when this matrix is positive
    semidefinite; the eigenvector of a negative eigenvalue is a certificate.
    """

    size = matrix.shape[0]
    centering = np.eye(size) - np.ones((size, size)) / size
    scaled = -0.5 * centering @ matrix @ centering
    values, vectors = np.linalg.eigh(scaled)
    return float(values.min()), centering @ vectors[:, 0]


def split_conic(q: int) -> tuple[int, ...]:
    """Fiber signature of ``xy`` on ``F_q^2``: one fiber ``2q-1``, then ``q-1``."""

    return tuple(sorted([2 * q - 1] + [q - 1] * (q - 1), reverse=True))


def anisotropic_conic(q: int) -> tuple[int, ...]:
    """Fiber signature of an anisotropic form: ``q-1`` fibers ``q+1``, then ``1``."""

    return tuple(sorted([q + 1] * (q - 1) + [1], reverse=True))


def conic_constant() -> tuple[float, float]:
    """``max_{beta>0} (2 beta - 2^beta)/(beta+1)`` and its maximiser."""

    objective = lambda beta: (2 * beta - 2.0**beta) / (beta + 1)  # noqa: E731
    low, high = 0.5, 3.0
    for _ in range(400):
        first = low + (high - low) / 3
        second = high - (high - low) / 3
        if objective(first) < objective(second):
            low = first
        else:
            high = second
    beta = (low + high) / 2
    return objective(beta), beta


def report_isometry() -> None:
    print("1. d(a,b) = oscillation of u_a - u_b")
    for a, b in ISOMETRY_PAIRS:
        print(
            f"   {label(a):>12} {label(b):>12}"
            f"   d = {irreversibility(a, b):.9f}"
            f"   osc = {oscillation(a, b):.9f}"
        )


def report_certificate() -> np.ndarray:
    print("\n2. failure of negative type")
    matrix = distance_matrix(CERTIFICATE_FAMILY)
    eigenvalue, certificate = negative_type_eigenvalue(matrix)
    quadratic = float(certificate @ matrix @ certificate)
    print(f"   family size            {len(CERTIFICATE_FAMILY)}")
    print(f"   min eigenvalue         {eigenvalue:.6e}   (negative type needs >= 0)")
    print(f"   x^T D x                {quadratic:.6e}   (negative type needs <= 0)")
    print(f"   sum of x               {float(certificate.sum()):.3e}")
    for t in (0.01, 0.05, 0.1, 0.5, 1.0):
        smallest = float(np.linalg.eigvalsh(np.exp(-t * matrix)).min())
        print(f"   min eig exp(-{t:>4} D)  {smallest:+.6e}")

    with CERTIFICATE_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["signature", "certificate_weight"])
        for signature, weight in zip(CERTIFICATE_FAMILY, certificate):
            writer.writerow([label(signature), f"{weight:.15f}"])
    print(f"   certificate written to {CERTIFICATE_CSV.name}")
    return matrix


def report_conics() -> None:
    kappa, beta_star = conic_constant()
    print("\n3. conics over finite fields")
    print(f"   kappa = {kappa:.12f} at beta* = {beta_star:.12f}")
    rows = []
    for q in CONIC_PRIMES:
        split = split_conic(q)
        anisotropic = anisotropic_conic(q)
        forward = exchange_rate(split, anisotropic)
        backward = exchange_rate(anisotropic, split)
        predicted = math.log(q + 1) / math.log(2 * q - 1)
        scaled = (1 - forward) * q * math.log(q)
        rows.append((q, forward, backward, predicted, scaled))
        print(
            f"   q={q:>4}  C(split->aniso)={forward:.10f}"
            f"  (1-C) q log q = {scaled:.6f}"
            f"   C(aniso->split)={backward:.10f}"
            f"  log(q+1)/log(2q-1)={predicted:.10f}"
        )
    with CONIC_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "q",
                "C_split_to_anisotropic",
                "C_anisotropic_to_split",
                "log_q_plus_1_over_log_2q_minus_1",
                "one_minus_C_times_q_log_q",
            ]
        )
        for q, forward, backward, predicted, scaled in rows:
            writer.writerow(
                [q, f"{forward:.15f}", f"{backward:.15f}", f"{predicted:.15f}", f"{scaled:.15f}"]
            )
    print(f"   rates written to {CONIC_CSV.name}")


def main() -> int:
    report_isometry()
    report_certificate()
    report_conics()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
