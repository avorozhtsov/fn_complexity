#!/usr/bin/env python3
"""T1.2 part 3: structure of the PSD set ``S(d) = {t > 0 : exp(-t d) is PSD}``.

Three things:

A. A solvable family.  For the complete-bipartite graph metric ``K_{m,n}``
   (distance 1 across the parts, 2 inside a part),

       exp(-t d) is PSD  <=>  t >= (1/2) log((m-1)(n-1)),

   so ``S`` is exactly a closed ray and ``t*`` is unbounded at fixed diameter.
   Proof: writing ``p = exp(-2t)``, ``q = exp(-t)``, the form on
   ``x = a 1_A + b 1_B + (centred)`` is minimised by the centred part zero and
   reduces to the 2x2 form ``(1-p+pm)A + (1-p+pn)B - 2 q sqrt(mn AB)``
   in ``A = m a^2``, ``B = n b^2``; its determinant condition simplifies to
   ``p (m-1)(n-1) <= 1``.  The script checks the formula numerically.

   Note the critical direction is NOT centred: PSD of ``K_t`` at a fixed ``t``
   is strictly stronger than negative type of the linearisation.

B. Is ``S`` ever disconnected?  ``S`` is a *closed additive sub-semigroup* of
   ``(0, inf)`` -- if ``K_s`` and ``K_u`` are PSD then so is their Schur product
   ``K_s o K_u = K_{s+u}`` (Schur product theorem) -- and it contains a ray,
   because ``K_t -> I``.  That leaves a gap logically possible.  The script
   searches random metrics and adversarial two-scale metrics for one.

C. Which direction is critical at ``t*`` for the exchange families -- centred
   (a negative-type certificate) or not.

Writes ``t1_2_part3_bipartite.csv``, ``t1_2_part3_disconnection_search.csv``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from t1_2_common import (  # noqa: E402
    ZERO_TOLERANCE,
    build_families,
    distance_matrix,
    kernel_min_eigenvalue,
    negative_type_defect,
    psd_threshold,
)

HERE = Path(__file__).resolve().parent
BIPARTITE_CSV = HERE / "t1_2_part3_bipartite.csv"
SEARCH_CSV = HERE / "t1_2_part3_disconnection_search.csv"

FINE = np.exp(np.linspace(np.log(1e-4), np.log(3e4), 1201))


def bipartite_metric(m: int, n: int) -> np.ndarray:
    size = m + n
    matrix = np.ones((size, size))
    matrix[:m, :m] = 2.0
    matrix[m:, m:] = 2.0
    np.fill_diagonal(matrix, 0.0)
    return matrix


def count_sign_changes(distances: np.ndarray, grid: np.ndarray = FINE) -> int:
    """Number of sign flips of ``lambda_min(K_t)`` along the grid.

    One flip means the PSD set is the ray ``[t*, inf)``; more would mean a gap.
    The kernels are stacked and diagonalised in one batched call -- per-``t``
    calls are dominated by numpy overhead at these matrix sizes.
    """

    scale = float(distances.max())
    unit = distances / scale
    stack = np.exp(-grid[:, None, None] * unit[None, :, :])
    values = np.linalg.eigvalsh(stack)[:, 0]
    signs = np.where(values < -ZERO_TOLERANCE, -1, 1)
    return int(np.count_nonzero(signs[1:] != signs[:-1]))


def report_bipartite() -> None:
    print("A. complete bipartite graph metrics: t* = (1/2) log((m-1)(n-1))")
    rows = []
    print(f"   {'m':>3} {'n':>3} {'t* numeric':>14} {'t* formula':>14} {'rel err':>10} {'flips':>6}")
    for m in range(2, 9):
        for n in range(m, 13):
            distances = bipartite_metric(m, n)
            star = psd_threshold(distances)
            product = (m - 1) * (n - 1)
            formula = 0.5 * math.log(product) if product > 1 else None
            flips = count_sign_changes(distances)
            numeric = float("nan") if star is None else star
            error = (
                float("nan")
                if star is None or formula in (None, 0.0)
                else abs(star - formula) / formula
            )
            rows.append(
                [m, n, "-" if star is None else f"{star:.12f}",
                 "-" if formula is None else f"{formula:.12f}",
                 f"{error:.3e}", flips]
            )
            print(
                f"   {m:>3} {n:>3} {numeric:>14.10f}"
                f" {float('nan') if formula is None else formula:>14.10f}"
                f" {error:>10.2e} {flips:>6}"
            )
    with BIPARTITE_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["m", "n", "t_star_numeric", "t_star_formula", "relative_error",
                         "sign_changes"])
        writer.writerows(rows)


def random_metric(size: int, rng: np.random.Generator) -> np.ndarray:
    """A random metric: shortest paths of a random weighted complete graph."""

    weights = rng.uniform(0.05, 1.0, size=(size, size))
    weights = np.triu(weights, 1)
    weights = weights + weights.T
    matrix = weights.copy()
    for k in range(size):  # Floyd-Warshall
        matrix = np.minimum(matrix, matrix[:, [k]] + matrix[[k], :])
    np.fill_diagonal(matrix, 0.0)
    return matrix


def two_scale_metric(rng: np.random.Generator) -> np.ndarray:
    """Two bipartite violators glued at very different scales.

    If a gap in ``S`` can be produced by superposing thresholds, this is where
    it should show up: the blocks want ``exp(-t d)`` PSD at two separated ``t``.
    """

    m1, n1 = int(rng.integers(2, 6)), int(rng.integers(2, 8))
    m2, n2 = int(rng.integers(2, 6)), int(rng.integers(2, 8))
    first = bipartite_metric(m1, n1)
    second = bipartite_metric(m2, n2) * float(rng.uniform(0.02, 0.5))
    size = first.shape[0] + second.shape[0]
    gap = float(rng.uniform(1.0, 40.0)) + max(first.max(), second.max())
    matrix = np.full((size, size), gap)
    matrix[: first.shape[0], : first.shape[0]] = first
    matrix[first.shape[0] :, first.shape[0] :] = second
    np.fill_diagonal(matrix, 0.0)
    return matrix


def report_search() -> None:
    print("\nB. search for a disconnected PSD set")
    rng = np.random.default_rng(31337)
    rows = []
    for kind, maker, trials in (
        ("random_metric_5", lambda r: random_metric(5, r), 4000),
        ("random_metric_6", lambda r: random_metric(6, r), 4000),
        ("random_metric_8", lambda r: random_metric(8, r), 3000),
        ("random_metric_12", lambda r: random_metric(12, r), 2000),
        ("two_scale_bipartite", two_scale_metric, 4000),
    ):
        violating = 0
        multi = 0
        worst = 0
        for _ in range(trials):
            distances = maker(rng)
            if negative_type_defect(distances / distances.max()) <= 1e-12:
                continue
            violating += 1
            flips = count_sign_changes(distances)
            worst = max(worst, flips)
            if flips > 1:
                multi += 1
        rows.append([kind, trials, violating, multi, worst])
        print(
            f"   {kind:<22} trials={trials:>5} non-negative-type={violating:>5}"
            f"  with >1 sign change={multi:>3}  max sign changes={worst}"
        )
    with SEARCH_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["kind", "trials", "non_negative_type", "with_multiple_sign_changes",
                         "max_sign_changes"])
        writer.writerows(rows)


def report_critical_direction() -> None:
    print("\nC. is the critical direction at t* centred (a negative-type certificate)?")
    families = build_families()
    for name, family in families.items():
        distances = distance_matrix(family)
        star = psd_threshold(distances)
        if star is None:
            continue
        kernel = np.exp(-star * distances)
        values, vectors = np.linalg.eigh(kernel)
        null = vectors[:, int(np.argmin(values))]
        size = len(family)
        share = abs(float(null.sum())) / math.sqrt(size)
        print(
            f"   {name:<12} t* = {star:>12.6f}   lambda_min = {values.min():+.2e}"
            f"   |<v,1>|/sqrt(N) = {share:.6f}"
            f"   ({'centred' if share < 1e-3 else 'NOT centred'})"
        )
    for m, n in ((2, 3), (3, 3), (3, 7), (5, 5)):
        distances = bipartite_metric(m, n)
        star = psd_threshold(distances)
        kernel = np.exp(-star * distances)
        values, vectors = np.linalg.eigh(kernel)
        null = vectors[:, int(np.argmin(values))]
        share = abs(float(null.sum())) / math.sqrt(m + n)
        print(
            f"   K_{m},{n:<9} t* = {star:>12.6f}   |<v,1>|/sqrt(N) = {share:.6f}"
            f"   ({'centred' if share < 1e-3 else 'NOT centred'})"
        )


def main() -> int:
    report_bipartite()
    report_search()
    report_critical_direction()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
