#!/usr/bin/env python3
"""Reproduce exchange-matrix examples used by the polynomial-map paper.

The first matrix is the finite-map *signature* exchange matrix attached to the
six nonzero nodes of the homogeneous quadratic tensor poset

    Sym^2(F_3^2) -> F_3^2.

It is not promoted to the operational affine-processor rate.  The second
calculation reads the already generated modulo-2^10 residue-signature matrix
for the seven anisotropic Q_2 classes.  For each matrix we enumerate simple
directed cycles and record the largest geometric-mean return at every length.
"""

from __future__ import annotations

import argparse
import csv
from itertools import permutations
import math
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate_result  # noqa: E402


Q3_NAMES = (
    "<x^2>",
    "<xy>",
    "<x^2+y^2>",
    "<xy;x^2>",
    "<y^2;x^2>",
    "<-x^2+y^2;xy>",
)
Q3_SIGNATURES = (
    (6, 3),
    (5, 2, 2),
    (4, 4, 1),
    (3, 2, 2, 2),
    (4, 2, 2, 1),
    (2, 2, 2, 2, 1),
)


def signature_matrix() -> list[list[float]]:
    return [
        [
            exchange_rate_result(
                implementer,
                implemented,
                grid_size=4096,
                tolerance=1e-14,
            ).rate
            for implemented in Q3_SIGNATURES
        ]
        for implementer in Q3_SIGNATURES
    ]


def read_matrix(path: Path) -> tuple[tuple[str, ...], list[list[float]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    names = tuple(rows[0][1:])
    row_names = tuple(row[0] for row in rows[1:])
    if row_names != names:
        raise ValueError(f"matrix row and column labels differ in {path}")
    matrix = [[float(value) for value in row[1:]] for row in rows[1:]]
    return names, matrix


def best_simple_cycles(
    names: tuple[str, ...], matrix: list[list[float]]
) -> list[tuple[int, float, float, tuple[int, ...]]]:
    """Return the best oriented simple cycle of each possible length.

    Rotations are removed by requiring the smallest vertex index to occur
    first.  Reverse orientations remain distinct, as exchange rates are
    directed.
    """

    size = len(names)
    result: list[tuple[int, float, float, tuple[int, ...]]] = []
    for length in range(2, size + 1):
        best: tuple[float, float, tuple[int, ...]] | None = None
        for cycle in permutations(range(size), length):
            if cycle[0] != min(cycle):
                continue
            product = math.prod(
                matrix[cycle[index]][cycle[(index + 1) % length]]
                for index in range(length)
            )
            mean = product ** (1.0 / length)
            candidate = (mean, product, cycle)
            if best is None or candidate > best:
                best = candidate
        if best is None:
            raise AssertionError("every requested cycle length must exist")
        result.append((length, best[0], best[1], best[2]))
    return result


def maximum_composition_residual(matrix: list[list[float]]) -> float:
    """Return max(M_ij M_jk - M_ik); it should be at most roundoff."""

    size = len(matrix)
    return max(
        matrix[first][middle] * matrix[middle][last] - matrix[first][last]
        for first in range(size)
        for middle in range(size)
        for last in range(size)
    )


def write_matrix(
    path: Path, names: tuple[str, ...], matrix: list[list[float]]
) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("implementer / implemented",) + names)
        for name, row in zip(names, matrix):
            writer.writerow((name,) + tuple(f"{value:.15f}" for value in row))


def write_cycles(
    path: Path,
    names: tuple[str, ...],
    cycles: list[tuple[int, float, float, tuple[int, ...]]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("length", "geometric_mean", "product", "cycle"))
        for length, mean, product, cycle in cycles:
            labels = tuple(names[index] for index in cycle)
            writer.writerow(
                (
                    length,
                    f"{mean:.15f}",
                    f"{product:.15f}",
                    " -> ".join(labels + labels[:1]),
                )
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "paper_finite_fields_maps" / "anc",
    )
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    q3_matrix = signature_matrix()
    q3_cycles = best_simple_cycles(Q3_NAMES, q3_matrix)
    q3_matrix_path = arguments.output_dir / "q3_quadratic_p1_signature_matrix.csv"
    q3_cycles_path = arguments.output_dir / "q3_quadratic_p1_cycle_means.csv"
    write_matrix(q3_matrix_path, Q3_NAMES, q3_matrix)
    write_cycles(q3_cycles_path, Q3_NAMES, q3_cycles)

    q2_source = PROJECT_ROOT / "analysis" / "p_adic_residue_rate_matrix_m10.csv"
    q2_names, q2_matrix = read_matrix(q2_source)
    q2_cycles = best_simple_cycles(q2_names, q2_matrix)
    q2_cycles_path = arguments.output_dir / "q2_anisotropic_residue_cycle_means_m10.csv"
    write_cycles(q2_cycles_path, q2_names, q2_cycles)

    print(f"q3 maximum composition residual: {maximum_composition_residual(q3_matrix):.3e}")
    print(f"q2 maximum composition residual: {maximum_composition_residual(q2_matrix):.3e}")
    print(q3_matrix_path)
    print(q3_cycles_path)
    print(q2_cycles_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
