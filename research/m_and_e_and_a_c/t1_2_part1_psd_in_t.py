#!/usr/bin/env python3
"""T1.2 part 1: for which ``t > 0`` is ``K_t = exp(-t d)`` positive semidefinite?

Scans ``lambda_min(K_t)`` on a log grid ``t in [0.001, 50]`` for the published
13-signature certificate family, for pseudo-random families of size 5..40 and
for greedily-optimised violating families, bisects the critical ``t*``, checks
whether the PSD set is a single ray, and brackets ``t*`` between two provable
bounds (a Gershgorin upper bound and a certificate lower bound).

Writes ``t1_2_part1_lambda_min.csv``, ``t1_2_part1_thresholds.csv``,
``t1_2_distance_matrices.json``.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from t1_2_common import (  # noqa: E402
    ZERO_TOLERANCE,
    build_families,
    distance_matrix,
    distance_matrix_grid,
    kernel_min_eigenvalue,
    label,
    negative_type_defect,
    negative_type_eigenvalue,
    psd_threshold,
    second_order_estimate,
)

HERE = Path(__file__).resolve().parent
SCAN_CSV = HERE / "t1_2_part1_lambda_min.csv"
THRESHOLD_CSV = HERE / "t1_2_part1_thresholds.csv"
DISTANCES_JSON = HERE / "t1_2_distance_matrices.json"

GRID = np.exp(np.linspace(np.log(1e-3), np.log(50.0), 601))


def sign_changes(values: np.ndarray, grid: np.ndarray) -> list[tuple[float, float]]:
    """Bracketing intervals where the sign of ``lambda_min`` flips."""

    signs = np.where(values < -ZERO_TOLERANCE, -1, 1)
    flips = np.nonzero(signs[1:] != signs[:-1])[0]
    return [(float(grid[i]), float(grid[i + 1])) for i in flips]


def gershgorin_bound(distances: np.ndarray) -> float:
    """Smallest ``t`` with ``max_i sum_{j!=i} exp(-t d_ij) <= 1``.

    ``K_t`` is diagonally dominant there, hence PSD; a proved upper bound on ``t*``.
    """

    off = distances + np.diag(np.full(distances.shape[0], np.inf))

    def excess(t: float) -> float:
        return float(np.exp(-t * off).sum(axis=1).max()) - 1.0

    low, high = 1e-6, 1.0
    while excess(high) > 0:
        high *= 2
        if high > 1e12:
            raise RuntimeError("Gershgorin bracket overflow")
    for _ in range(200):
        middle = 0.5 * (low + high)
        if excess(middle) > 0:
            low = middle
        else:
            high = middle
    return high


def certificate_bound(distances: np.ndarray) -> float:
    """Largest ``t`` at which the negative-type certificate still proves failure.

    With ``x`` centred, ``|x| = 1`` and ``q = x^T d x > 0``,
    ``x^T K_t x = -t q + sum_{k>=2} (-t)^k x^T d^{o k} x / k!`` and
    ``|x^T d^{o k} x| <= N diam^k``, so ``K_t`` is not PSD whenever
    ``t q > N (exp(t diam) - 1 - t diam)``.  Hence ``t*`` exceeds the root.
    """

    size = distances.shape[0]
    centring = np.eye(size) - np.ones((size, size)) / size
    values, vectors = np.linalg.eigh(-0.5 * centring @ distances @ centring)
    if values.min() >= 0:
        return 0.0
    x = centring @ vectors[:, int(np.argmin(values))]
    x = x / np.linalg.norm(x)
    q = float(x @ distances @ x)
    if q <= 0:
        return 0.0
    diameter = float(distances.max())

    def slack(t: float) -> float:
        return t * q - size * (math.expm1(t * diameter) - t * diameter)

    low, high = 1e-12, 1.0
    while slack(high) > 0:
        high *= 2
    for _ in range(200):
        middle = 0.5 * (low + high)
        if slack(middle) > 0:
            low = middle
        else:
            high = middle
    return low


def main() -> int:
    families = build_families()
    scan_rows: list[list[object]] = []
    threshold_rows: list[list[object]] = []
    stored: dict[str, list[list[float]]] = {}

    header = (
        f"{'family':<12} {'N':>3} {'defect':>11} {'t*':>13} {'flips':>5}"
        f" {'diam':>9} {'t2nd':>12} {'t2nd/t*':>8} {'lower':>10} {'upper':>10}"
    )
    print(header)
    print("-" * len(header))
    for name, family in families.items():
        distances = distance_matrix(family)
        stored[name] = distances.tolist()
        values = np.array([kernel_min_eigenvalue(distances, t) for t in GRID])
        for t, value in zip(GRID, values):
            scan_rows.append([name, f"{t:.8g}", f"{value:.12e}"])

        flips = sign_changes(values, GRID)
        star = psd_threshold(distances)
        defect = negative_type_defect(distances)
        off = distances[~np.eye(len(family), dtype=bool)]
        diameter = float(distances.max())
        lower = certificate_bound(distances)
        upper = gershgorin_bound(distances)
        second = second_order_estimate(distances)
        threshold_rows.append(
            [
                name,
                len(family),
                f"{negative_type_eigenvalue(distances):.6e}",
                f"{defect:.6e}",
                "-" if star is None else f"{star:.9f}",
                len(flips),
                f"{diameter:.6f}",
                f"{off.min():.6f}",
                f"{off.mean():.6f}",
                "-" if star is None else f"{star * diameter:.6f}",
                f"{second:.9f}",
                "-" if star is None else f"{second / star:.6f}",
                f"{lower:.6e}",
                f"{upper:.4f}",
            ]
        )
        star_text = "     negtype" if star is None else f"{star:>13.6f}"
        ratio_text = "       -" if star is None else f"{second / star:>8.4f}"
        print(
            f"{name:<12} {len(family):>3} {defect:>11.3e} {star_text}"
            f" {len(flips):>5} {diameter:>9.4f} {second:>12.6f}"
            f" {ratio_text} {lower:>10.2e} {upper:>10.2f}"
        )

    with SCAN_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["family", "t", "lambda_min_exp_minus_t_d"])
        writer.writerows(scan_rows)
    with THRESHOLD_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "family",
                "size",
                "negative_type_min_eigenvalue",
                "negative_type_defect",
                "t_star",
                "grid_sign_changes",
                "diameter",
                "min_distance",
                "mean_distance",
                "t_star_times_diameter",
                "second_order_estimate",
                "second_order_over_t_star",
                "provable_lower_bound",
                "gershgorin_upper_bound",
            ]
        )
        writer.writerows(threshold_rows)
    DISTANCES_JSON.write_text(
        json.dumps(
            {
                name: {
                    "signatures": [label(s) for s in families[name]],
                    "distances": stored[name],
                }
                for name in families
            }
        ),
        encoding="utf-8",
    )

    # -------------------------------------------------- semigroup consistency
    print("\nSchur/semigroup check: t, s PSD => t+s PSD (200 random pairs/family)")
    rng = np.random.default_rng(2024)
    for name in families:
        distances = np.array(stored[name])
        star = psd_threshold(distances)
        if star is None:
            continue
        pairs = np.exp(rng.uniform(np.log(star), np.log(50.0), size=(200, 2)))
        worst = min(kernel_min_eigenvalue(distances, float(t + s)) for t, s in pairs)
        print(f"   {name:<12} min lambda_min(K_(t+s)) = {worst:+.3e}")

    # ------------------------------------------- fine scan just below/above t*
    print("\nFine scan around t* (is the transition a single crossing?)")
    for name in families:
        distances = np.array(stored[name])
        star = psd_threshold(distances)
        if star is None:
            continue
        fine = star * np.exp(np.linspace(np.log(0.5), np.log(2.0), 4001))
        values = np.array([kernel_min_eigenvalue(distances, float(t)) for t in fine])
        flips = sign_changes(values, fine)
        print(f"   {name:<12} sign changes in [t*/2, 2t*]: {len(flips)}")

    # -------------------------------------------------- independent check
    print("\nIndependent grid recomputation of d (2e6 points on [0,60] + beta=inf)")
    for name, family in families.items():
        solver = distance_matrix(family)
        star_solver = psd_threshold(solver)
        if star_solver is None:
            continue
        grid = distance_matrix_grid(family)
        star_grid = psd_threshold(grid)
        print(
            f"   {name:<12} max|d_solver - d_grid| = {np.abs(solver - grid).max():.3e}"
            f"   defect_grid = {negative_type_defect(grid):.6e}"
            f"   t*_solver = {star_solver:.9f}   t*_grid = {star_grid:.9f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
