#!/usr/bin/env python3
"""The Weil pairing of two signatures, and how it differs from the exchange rate.

A signature ``a = (a_1, ..., a_r)`` *is* a Weil test measure: the atomic measure
``sum_i delta_{a_i}`` on the positive reals has Mellin transform
``Z_a(s) = sum_i a_i^s``, which is exactly the partition function.  So the two
matrices of the project live on the same objects:

    M_ab = C(a -> b) = inf_beta log Z_a(beta) / log Z_b(beta)      (real spectrum)
    E_ab = sum over zeta zeros rho of Z_a(rho) conj(Z_b(rho))      (complex spectrum)

``E`` restricted to the first ``N`` zeros is a genuine Gram matrix of the vectors
``v_a = (Z_a(rho_n))_n``, hence positive semidefinite, and it is the finite-rank
Weil matrix for this family of test measures.

THE POINT.  Writing ``rho = 1/2 + i gamma`` gives
``Z_a(rho) = sum_i sqrt(a_i) exp(i gamma log a_i)``, so

    E_ab = sum_{i,j} sqrt(a_i b_j) * S(a_i / b_j),
    S(x) = sum_{0 < gamma <= T} x^{i gamma}.

Landau's theorem evaluates ``S``: it is ``N`` at ``x = 1`` and
``-(T/2pi) Lambda(y)/sqrt(y) + O(log T)`` otherwise, where ``y = max(x, 1/x)``
and ``Lambda`` is von Mangoldt.  Since ``sqrt(a_i b_j)/sqrt(y) = min(a_i, b_j)``,
this collapses to

    E_ab = N * O(a,b) - (T/2pi) * A(a,b) + O(r s log T),
    O(a,b) = sum over pairs with a_i = b_j of a_i,                (multiset overlap)
    A(a,b) = sum over pairs with a_i != b_j of
             min(a_i,b_j) * Lambda(max(a_i,b_j)/min(a_i,b_j)).    (prime-power ratios)

So ``E`` sees only *multiplicative coincidences* between the fiber sizes --
equality, and ratios that are prime powers -- while ``M`` sees the sizes through
their logarithms and is finite for every pair.  Two signatures sharing no entry
and no prime-power ratio are Weil-orthogonal while being perfectly comparable in
the exchange sense.  This script verifies the formula and quantifies the
divergence between the two geometries.
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path
import sys

import numpy as np

RESEARCH_DIRECTORY = Path(__file__).resolve().parent
PROJECT_ROOT = RESEARCH_DIRECTORY.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate  # noqa: E402

ZEROS_PATH = RESEARCH_DIRECTORY / "zeta_zeros_1200.npy"
CSV_PATH = RESEARCH_DIRECTORY / "weil_pairing_on_signatures.csv"

FAMILY = [
    (2, 2), (3, 1), (4, 2), (5, 3), (6, 1), (3, 1, 1), (8, 4), (9, 3),
    (4, 4, 1), (7, 5), (5, 5), (7, 1), (9, 1), (6, 3), (10, 5), (11, 7),
]


def von_mangoldt(x: float, tolerance: float = 1e-9) -> float:
    """``Lambda(x)``: ``log p`` when ``x`` is exactly a prime power, else 0."""

    if abs(x - round(x)) > tolerance:
        return 0.0
    n = int(round(x))
    if n < 2:
        return 0.0
    remainder, divisor = n, 2
    while divisor * divisor <= remainder:
        if remainder % divisor == 0:
            while remainder % divisor == 0:
                remainder //= divisor
            return math.log(divisor) if remainder == 1 else 0.0
        divisor += 1
    return math.log(n)


def overlap(a: tuple[int, ...], b: tuple[int, ...]) -> float:
    return float(sum(x for x in a for y in b if x == y))


def arithmetic_term(a: tuple[int, ...], b: tuple[int, ...]) -> float:
    total = 0.0
    for x in a:
        for y in b:
            if x != y:
                high, low = max(x, y), min(x, y)
                total += low * von_mangoldt(high / low)
    return total


def mellin_vectors(family, gammas: np.ndarray) -> dict:
    vectors = {}
    for signature in family:
        entries = np.asarray(signature, dtype=float)
        phases = np.exp(1j * np.outer(np.log(entries), gammas))
        vectors[signature] = (np.sqrt(entries)[:, None] * phases).sum(axis=0)
    return vectors


def main() -> int:
    gammas = np.load(ZEROS_PATH)
    count = len(gammas)
    height = gammas[-1]
    scale = height / (2 * math.pi)
    print(f"{count} zeta zeros, T = {height:.3f}, T/2pi = {scale:.3f}\n")

    vectors = mellin_vectors(FAMILY, gammas)
    gram = np.array(
        [[np.vdot(vectors[b], vectors[a]).real for b in FAMILY] for a in FAMILY]
    )
    print("E is a Gram matrix, so positive semidefinite by construction:")
    print(f"   min eigenvalue = {np.linalg.eigvalsh(gram).min():.6e}\n")

    print("Landau prediction  E_ab = N*O(a,b) - (T/2pi)*A(a,b):")
    print(f"   {'a':>10} {'b':>10} {'E (actual)':>13} {'predicted':>12} "
          f"{'O':>5} {'A':>8} {'abs err':>9}")
    rows = []
    errors = []
    orthogonal = []
    for index_a, a in enumerate(FAMILY):
        for index_b, b in enumerate(FAMILY):
            if index_b < index_a:
                continue
            actual = gram[index_a, index_b]
            common, arithmetic = overlap(a, b), arithmetic_term(a, b)
            predicted = count * common - scale * arithmetic
            errors.append(abs(actual - predicted))
            rows.append(
                [str(a), str(b), f"{actual:.6f}", f"{predicted:.6f}",
                 f"{common:.0f}", f"{arithmetic:.6f}"]
            )
            if common == 0 and arithmetic == 0 and index_a != index_b:
                orthogonal.append((a, b, actual))
            if index_a < 6 and index_b < 8:
                print(f"   {str(a):>10} {str(b):>10} {actual:>13.2f} "
                      f"{predicted:>12.2f} {common:>5.0f} {arithmetic:>8.4f} "
                      f"{abs(actual - predicted):>9.2f}")
    typical = count * max(overlap(a, a) for a in FAMILY)
    print(f"\n   max abs error {max(errors):.1f}, mean {np.mean(errors):.1f}, "
          f"against a diagonal scale of ~{typical:.0f}")
    print(f"   predicted error size O(r s log T) ~ {math.log(height) * 9:.0f}\n")

    print("Weil-orthogonal pairs (no shared entry, no prime-power ratio):")
    for a, b, value in orthogonal[:10]:
        distance = -math.log(exchange_rate(a, b) * exchange_rate(b, a))
        print(
            f"   {str(a):>10} vs {str(b):>10}   E = {value:>9.2f}  "
            f"(scale {count * 8:.0f})   exchange distance d = {distance:.6f}"
        )

    print("\nthe two geometries, compared:")
    diagonal = np.sqrt(np.diag(gram))
    correlation = gram / np.outer(diagonal, diagonal)
    pairs = []
    for index_a, a in enumerate(FAMILY):
        for index_b, b in enumerate(FAMILY):
            if index_b <= index_a:
                continue
            distance = -math.log(exchange_rate(a, b) * exchange_rate(b, a))
            angle = math.acos(max(-1.0, min(1.0, correlation[index_a, index_b])))
            pairs.append((distance, angle, a, b))
    distances = np.array([p[0] for p in pairs])
    angles = np.array([p[1] for p in pairs])
    print(f"   correlation between exchange distance d and Weil angle: "
          f"{np.corrcoef(distances, angles)[0, 1]:+.4f}")
    pairs.sort(key=lambda p: p[0])
    print(f"   {'d (exchange)':>13} {'angle (Weil)':>13}   pair")
    for distance, angle, a, b in pairs[:4] + pairs[-4:]:
        print(f"   {distance:>13.6f} {angle:>13.6f}   {a} vs {b}")

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["a", "b", "E_actual", "E_predicted", "overlap", "arithmetic"])
        writer.writerows(rows)
    print(f"\nwritten to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
