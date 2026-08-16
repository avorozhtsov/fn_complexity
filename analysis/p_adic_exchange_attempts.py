#!/usr/bin/env python3
"""Reproducible attempts at 2-adic quadratic-map exchange rates.

The affine rate itself is defined using exact affine processors over Q_2.  This
script computes two deliberately weaker diagnostics:

1. finite-map signature rates after reduction modulo 2^m;
2. an exhaustive obstruction modulo 4 for integral 3-source/2-target linear
   processor identities in the standard lattices.

Neither diagnostic is silently promoted to the exact Q_2-affine rate.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations_with_replacement, product
import math
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
ANISOTROPIC_D = (1, 2, -2, 5, -5, 10, -10)
MONOMIAL_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


@dataclass(frozen=True)
class ResidueRate:
    """A finite-signature exchange rate and its minimizing exponent."""

    rate: float
    beta: float


@dataclass(frozen=True)
class ModFourSearch:
    """Summary of one exhaustive integral 3-to-2 search modulo four."""

    source_d: int
    target_d: int
    pullback_forms: int
    contribution_terms: int
    pair_sums: int
    found: bool


def residue_signature(m: int, d: int) -> tuple[int, ...]:
    """Return the nonzero fiber sizes of x^2+d*y^2 modulo 2^m.

    Circular convolution of the two square-residue distributions replaces an
    O(2^(2m)) enumeration.  The integer result is checked after FFT rounding.
    """

    if m < 1 or m > 20:
        raise ValueError("m must lie between 1 and 20")
    modulus = 1 << m
    residues = np.arange(modulus, dtype=np.int64)
    squares = np.bincount((residues * residues) % modulus, minlength=modulus)
    scaled_squares = np.zeros(modulus, dtype=np.int64)
    np.add.at(scaled_squares, (d * residues) % modulus, squares)
    counts = np.rint(
        np.fft.ifft(np.fft.fft(squares) * np.fft.fft(scaled_squares)).real
    ).astype(np.int64)
    if counts.min() < 0 or int(counts.sum()) != modulus * modulus:
        raise AssertionError("FFT convolution did not recover an exact fiber partition")
    return tuple(sorted((int(value) for value in counts if value), reverse=True))


def _compress(signature: tuple[int, ...]) -> tuple[tuple[float, float], ...]:
    values, multiplicities = np.unique(signature, return_counts=True)
    return tuple(
        (math.log(int(value)), math.log(int(multiplicity)))
        for value, multiplicity in zip(values, multiplicities)
    )


def _log_partition(compressed: tuple[tuple[float, float], ...], beta: float) -> float:
    terms = tuple(log_multiplicity + beta * log_value for log_value, log_multiplicity in compressed)
    largest = max(terms)
    return largest + math.log(sum(math.exp(term - largest) for term in terms))


def _golden_minimum(function, left: float, right: float) -> tuple[float, float]:
    inverse_phi = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = right - inverse_phi * (right - left)
    x2 = left + inverse_phi * (right - left)
    y1, y2 = function(x1), function(x2)
    for _ in range(160):
        if right - left <= 1e-12 * (1.0 + abs(x1) + abs(x2)):
            break
        if y1 <= y2:
            right, x2, y2 = x2, x1, y1
            x1 = right - inverse_phi * (right - left)
            y1 = function(x1)
        else:
            left, x1, y1 = x1, x2, y2
            x2 = left + inverse_phi * (right - left)
            y2 = function(x2)
    return (y1, x1) if y1 <= y2 else (y2, x2)


def residue_exchange_rate(
    implementer: tuple[int, ...], implemented: tuple[int, ...]
) -> ResidueRate:
    """Compute C(implementer -> implemented) from compressed fiber data."""

    source = _compress(implementer)
    target = _compress(implemented)
    gaps: list[float] = []
    for compressed in (source, target):
        log_values = [entry[0] for entry in compressed]
        if len(log_values) > 1:
            gaps.extend(log_values[-1] - value for value in log_values[:-1])
    horizon = min(36.0 / min(gaps), 1.0e5) if gaps else 64.0
    grid = (0.0,) + tuple(np.geomspace(1e-10, horizon, 400))

    def ratio(beta: float) -> float:
        return _log_partition(source, beta) / _log_partition(target, beta)

    values = tuple(ratio(beta) for beta in grid)
    candidates = [
        (values[0], 0.0),
        (
            math.log(implementer[0]) / math.log(implemented[0]),
            math.inf,
        ),
    ]
    for index in range(1, len(grid) - 1):
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]:
            candidates.append(
                _golden_minimum(ratio, grid[index - 1], grid[index + 1])
            )
    rate, beta = min(candidates)
    return ResidueRate(rate, beta)


def _quadratic_coefficients(
    first_row: tuple[int, ...],
    second_row: tuple[int, ...],
    d: int,
    modulus: int,
) -> tuple[int, ...]:
    diagonal = tuple(
        (first_row[index] ** 2 + d * second_row[index] ** 2) % modulus
        for index in range(4)
    )
    cross = tuple(
        (
            2
            * (
                first_row[left] * first_row[right]
                + d * second_row[left] * second_row[right]
            )
        )
        % modulus
        for left, right in MONOMIAL_PAIRS
    )
    return diagonal + cross


@lru_cache(maxsize=None)
def pullback_forms_mod_four(d: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate all q_d pullbacks from four variables over Z/4Z."""

    rows = tuple(product(range(4), repeat=4))
    return tuple(
        sorted(
            {
                _quadratic_coefficients(first, second, d, 4)
                for first in rows
                for second in rows
            }
        )
    )


_LOW_BITS_40 = sum(1 << (2 * index) for index in range(20))
_HIGH_BITS_40 = _LOW_BITS_40 << 1


def _encode_base_four(values: tuple[int, ...]) -> int:
    return sum((value & 3) << (2 * index) for index, value in enumerate(values))


def _add_base_four(left: int, right: int) -> int:
    low = (left ^ right) & _LOW_BITS_40
    high = ((left ^ right) & _HIGH_BITS_40) ^ ((left & right & _LOW_BITS_40) << 1)
    return low | high


def _negate_base_four(value: int) -> int:
    return value ^ ((value & _LOW_BITS_40) << 1)


def integral_three_to_two_mod_four(source_d: int, target_d: int) -> ModFourSearch:
    """Exhaust the standard-lattice integral 3-to-2 identity modulo four.

    Each source copy contributes ``b_i q_source(A_i z)`` to two target
    outputs, with arbitrary A_i in Mat(2,4,Z/4) and b_i in (Z/4)^2.  A
    meet-in-the-middle three-sum decides whether the two target quadrics occur.
    """

    forms = pullback_forms_mod_four(source_d)
    contributions = {
        _encode_base_four(
            tuple((first * value) % 4 for value in form)
            + tuple((second * value) % 4 for value in form)
        )
        for first, second in product(range(4), repeat=2)
        for form in forms
    }
    terms = tuple(sorted(contributions))
    pair_sums = {
        _add_base_four(terms[left], terms[right])
        for left, right in combinations_with_replacement(range(len(terms)), 2)
    }
    first_target = (1, target_d % 4, 0, 0, 0, 0, 0, 0, 0, 0)
    second_target = (0, 0, 1, target_d % 4, 0, 0, 0, 0, 0, 0)
    target = _encode_base_four(first_target + second_target)
    found = any(
        _add_base_four(target, _negate_base_four(term)) in pair_sums
        for term in terms
    )
    return ModFourSearch(
        source_d,
        target_d,
        len(forms),
        len(terms),
        len(pair_sums),
        found,
    )


def write_pair_convergence(max_m: int, output: Path) -> None:
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "m",
                "q1_outputs",
                "q2_outputs",
                "C_mod_2m_q1_to_q2",
                "beta_q1_to_q2",
                "C_mod_2m_q2_to_q1",
                "beta_q2_to_q1",
            )
        )
        for m in range(2, max_m + 1):
            q1 = residue_signature(m, 1)
            q2 = residue_signature(m, 2)
            forward = residue_exchange_rate(q1, q2)
            reverse = residue_exchange_rate(q2, q1)
            writer.writerow(
                (
                    m,
                    len(q1),
                    len(q2),
                    f"{forward.rate:.15f}",
                    "inf" if math.isinf(forward.beta) else f"{forward.beta:.15f}",
                    f"{reverse.rate:.15f}",
                    "inf" if math.isinf(reverse.beta) else f"{reverse.beta:.15f}",
                )
            )


def write_rate_matrix(m: int, output: Path) -> None:
    signatures = {d: residue_signature(m, d) for d in ANISOTROPIC_D}
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("implementer_d",) + ANISOTROPIC_D)
        for implementer_d in ANISOTROPIC_D:
            writer.writerow(
                (implementer_d,)
                + tuple(
                    f"{residue_exchange_rate(signatures[implementer_d], signatures[implemented_d]).rate:.15f}"
                    for implemented_d in ANISOTROPIC_D
                )
            )


def write_mod_four_searches(output: Path) -> tuple[ModFourSearch, ModFourSearch]:
    searches = (
        integral_three_to_two_mod_four(1, 2),
        integral_three_to_two_mod_four(2, 1),
    )
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "source_d",
                "target_d",
                "pullback_forms",
                "contribution_terms",
                "pair_sums",
                "solution_found",
            )
        )
        for search in searches:
            writer.writerow(
                (
                    search.source_d,
                    search.target_d,
                    search.pullback_forms,
                    search.contribution_terms,
                    search.pair_sums,
                    str(search.found).lower(),
                )
            )
    return searches


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-m", type=int, default=14)
    parser.add_argument("--matrix-m", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, default=ANALYSIS_DIR)
    parser.add_argument("--skip-mod-four", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pair_path = args.output_dir / "p_adic_residue_rates_q1_q2.csv"
    matrix_path = args.output_dir / f"p_adic_residue_rate_matrix_m{args.matrix_m}.csv"
    write_pair_convergence(args.max_m, pair_path)
    write_rate_matrix(args.matrix_m, matrix_path)
    print(pair_path)
    print(matrix_path)
    if not args.skip_mod_four:
        search_path = args.output_dir / "p_adic_integral_3_to_2_mod4.csv"
        searches = write_mod_four_searches(search_path)
        print(search_path)
        for search in searches:
            print(
                f"q_{search.source_d} -> q_{search.target_d}: "
                f"mod-4 solution found = {search.found}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
