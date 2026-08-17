#!/usr/bin/env python3
"""T1.2 part 2: the ``L^t`` bridge from the Gram (``L^2``) end to the tropical
(``L^infinity``) end.

For a probability measure ``mu`` on the spectrum, with ``phi_a = log Z_a`` and
``u_a = log phi_a``:

    G_t(i,j) = int (phi_i phi_j)^(t/2) dmu          Gram of phi^(t/2): PSD, all t
    P_t(i,j) = ( int (phi_i/phi_j)^t dmu )^(1/t)    normalised price
    D_t(i,j) = log P_t(i,j) + log P_t(j,i)          soft irreversibility

``D_t(i,i) = 0`` and ``D_t >= 0`` by Jensen.  Writing ``g = u_i - u_j``,
``log P_t(g) = sum_k t^(k-1) kappa_k(g) / k!`` (cumulants under ``mu``) and
``kappa_k(-g) = (-1)^k kappa_k(g)``, so the odd cumulants cancel in ``D_t``:

    D_t(i,j) = t kappa_2(g) + t^3 kappa_4(g)/12 + t^5 kappa_6(g)/360 + ...

Hence two exact limits:

    D_t / t -> Var_mu(u_i - u_j)      as t -> 0    (squared L^2(mu) distance)
    D_t     -> osc_supp(u_i - u_j)    as t -> inf  (the exchange metric d)

The first is a squared Hilbert-space distance, so it is of negative type and
``exp(-s D_t)`` is PSD for every ``s``; the second is not of negative type.
The phase boundary is the curve ``s*(t) = inf{s : exp(-s D_t) is PSD}``.

Writes ``t1_2_part2_price_limit.csv``, ``t1_2_part2_phase.csv``,
``t1_2_part2_l2_spectrum.csv``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from t1_2_common import (  # noqa: E402
    CERTIFICATE_FAMILY,
    build_families,
    distance_matrix,
    exchange_rate,
    label,
    min_eigenvalue,
    negative_type_defect,
    psd_threshold,
    u_values,
)

HERE = Path(__file__).resolve().parent
PRICE_CSV = HERE / "t1_2_part2_price_limit.csv"
PHASE_CSV = HERE / "t1_2_part2_phase.csv"
SPECTRUM_CSV = HERE / "t1_2_part2_l2_spectrum.csv"

GRID_POINTS = 10_001
T_GRID = np.exp(np.linspace(np.log(0.03), np.log(3.0e4), 41))


def logsumexp(values: np.ndarray, axis: int | None = None):
    if axis is None:
        largest = float(np.max(values))
        return largest + float(np.log(np.exp(values - largest).sum()))
    largest = np.max(values, axis=axis, keepdims=True)
    total = largest + np.log(np.exp(values - largest).sum(axis=axis, keepdims=True))
    return np.squeeze(total, axis=axis)


def measure(name: str, points: int = GRID_POINTS):
    """Nodes and normalised log-weights of ``mu``.

    ``full`` is the push-forward of the uniform law on ``(0,1)`` under
    ``beta = x/(1-x)``: it charges the whole spectrum ``[0, inf]``, so the
    ``t -> inf`` limit of ``D_t`` is the exchange metric ``d`` itself.  The
    others are uniform on the stated window.
    """

    if name == "full":
        x = (np.arange(points) + 0.5) / points
        betas = x / (1.0 - x)
    else:
        low, high = WINDOWS[name]
        betas = np.linspace(low, high, points)
    return betas, np.full(points, -math.log(points))


WINDOWS = {
    "u_0.05_20": (0.05, 20.0),
    "u_0.01_5": (0.01, 5.0),
    "u_0.5_60": (0.5, 60.0),
    "u_0.001_200": (0.001, 200.0),
}
MEASURES = ["full", *WINDOWS]


def u_table(family, betas: np.ndarray) -> np.ndarray:
    return np.array([u_values(signature, betas) for signature in family])


class Kernels:
    """Precomputed pairwise differences and sums of ``u`` on the quadrature grid.

    Every ``D_t`` and ``G_t`` is a logsumexp over the same tensors, so caching
    them turns the whole ``t``-scan into one pass of elementwise work per ``t``.
    """

    def __init__(self, table: np.ndarray, log_weights: np.ndarray) -> None:
        self.size = table.shape[0]
        self.log_weights = log_weights
        self.difference = table[:, None, :] - table[None, :, :]
        self.total = table[:, None, :] + table[None, :, :]

    def soft_distance(self, t: float) -> np.ndarray:
        price = (
            logsumexp(t * self.difference + self.log_weights[None, None, :], axis=2) / t
        )
        matrix = price + price.T
        np.fill_diagonal(matrix, 0.0)
        return matrix

    def gram_correlation(self, t: float) -> np.ndarray:
        log_gram = logsumexp(
            (t / 2) * self.total + self.log_weights[None, None, :], axis=2
        )
        diagonal = np.diag(log_gram)
        return np.exp(log_gram - 0.5 * (diagonal[:, None] + diagonal[None, :]))

    def oscillation(self) -> np.ndarray:
        return self.difference.max(axis=2) - self.difference.min(axis=2)


def soft_distance(table: np.ndarray, log_weights: np.ndarray, t: float) -> np.ndarray:
    return Kernels(table, log_weights).soft_distance(t)


def oscillation_matrix(table: np.ndarray) -> np.ndarray:
    """``osc`` of ``u_i - u_j`` over the support: the ``t -> inf`` limit of ``D_t``."""

    size = table.shape[0]
    out = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            difference = table[i] - table[j]
            out[i, j] = out[j, i] = float(difference.max() - difference.min())
    return out


def l2_squared_distance(table: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """``Var_mu(u_i - u_j)``: the ``t -> 0`` limit of ``D_t/t``."""

    mean = table @ weights
    centred = table - mean[:, None]
    gram = (centred * weights[None, :]) @ centred.T
    diagonal = np.diag(gram)
    return diagonal[:, None] + diagonal[None, :] - 2 * gram


def l2_gram(table: np.ndarray, weights: np.ndarray) -> np.ndarray:
    mean = table @ weights
    centred = table - mean[:, None]
    return (centred * weights[None, :]) @ centred.T


def triangle_defect(matrix: np.ndarray) -> float:
    """``max_ijk (D_ij - D_ik - D_kj)``: positive iff the triangle inequality fails."""

    size = matrix.shape[0]
    worst = -np.inf
    for k in range(size):
        worst = max(worst, float((matrix - matrix[:, [k]] - matrix[[k], :]).max()))
    return worst


# ------------------------------------------------------------------- reports


def report_price_limit() -> None:
    """``P_t(i,j) -> 1/C(f_j -> f_i)``: the Laplace / Varadhan principle."""

    print("1. Laplace limit   P_t(a,b) -> 1 / C(b -> a) = sup_beta phi_a/phi_b")
    pairs = [
        ((5, 5, 5, 1), (6, 3, 2)),
        ((7, 7, 6, 1), (5, 4, 4)),
        ((6, 6, 2, 2), (6, 5, 5, 5)),
        ((9, 1), (3, 2, 2, 2)),
        ((2, 2), (3, 1)),
    ]
    betas, log_weights = measure("full", 400_001)
    rows = []
    print(
        f"   {'a':>12} {'b':>12} {'P_100':>14} {'P_1e4':>14} {'P_1e6':>14}"
        f" {'1/C (solver)':>14} {'rel err':>10}"
    )
    for a, b in pairs:
        table = u_table([a, b], betas)
        difference = table[0] - table[1]
        exact = 1.0 / exchange_rate(b, a)
        values = {}
        for t in (1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6):
            value = math.exp(logsumexp(t * difference + log_weights) / t)
            values[t] = value
            rows.append(
                [
                    label(a),
                    label(b),
                    f"{t:g}",
                    f"{value:.12f}",
                    f"{exact:.12f}",
                    f"{(value - exact) / exact:.6e}",
                    f"{(exact - value) * t / math.log(t) if t > 1 else float('nan'):.6e}",
                ]
            )
        print(
            f"   {label(a):>12} {label(b):>12} {values[100.0]:>14.10f}"
            f" {values[1e4]:>14.10f} {values[1e6]:>14.10f} {exact:>14.10f}"
            f" {abs(values[1e6] - exact) / exact:>10.2e}"
        )
    with PRICE_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["a", "b", "t", "P_t", "one_over_C", "relative_error", "gap_times_t_over_log_t"]
        )
        writer.writerows(rows)


def report_l2_spectrum() -> None:
    """The ``L^2(mu)`` Gram spectrum of ``{u_a}``: how singular is the Hilbert end?"""

    print("\n2. the L^2 end: spectrum of the centred Gram of {u_a} (cert13, mu = full)")
    betas, log_weights = measure("full")
    weights = np.exp(log_weights)
    table = u_table(CERTIFICATE_FAMILY, betas)
    rows = []
    print(f"   {'N':>3} {'lambda_max':>12} {'lambda_min':>12} {'ratio':>11} {'decay':>8}")
    previous = None
    for n in range(2, len(CERTIFICATE_FAMILY) + 1):
        values = np.linalg.eigvalsh(l2_gram(table[:n], weights))
        ratio = float(values.min() / values.max())
        decay = "" if previous is None else f"{previous / ratio:8.1f}"
        print(
            f"   {n:>3} {values.max():>12.4e} {values.min():>12.4e}"
            f" {ratio:>11.3e} {decay:>8}"
        )
        rows.append([n, f"{values.max():.9e}", f"{values.min():.9e}", f"{ratio:.9e}"])
        previous = ratio
    with SPECTRUM_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["N", "lambda_max", "lambda_min", "ratio"])
        writer.writerows(rows)

    print("\n   small-t asymptotics: D_t/t - E should be O(t^2)")
    limit_zero = l2_squared_distance(table, weights)
    for t in (0.01, 0.02, 0.04, 0.08, 0.16, 0.32):
        gap = float(np.abs(soft_distance(table, log_weights, t) / t - limit_zero).max())
        print(f"      t = {t:5.3f}   max|D_t/t - E| = {gap:.6e}   /t^2 = {gap / t**2:.6e}")


def scan(name: str, family, measure_name: str, writer) -> None:
    betas, log_weights = measure(measure_name)
    weights = np.exp(log_weights)
    table = u_table(family, betas)

    kernels = Kernels(table, log_weights)
    limit_zero = l2_squared_distance(table, weights)
    limit_infinity = kernels.oscillation()
    solver = distance_matrix(family)
    star = psd_threshold(solver)

    print(f"\n   family {name}   mu = {measure_name}")
    print(
        f"      t -> 0   limit E = Var(u_i - u_j):   defect = "
        f"{negative_type_defect(limit_zero):.3e}  (0 = negative type, as proved)"
    )
    print(
        f"      t -> inf limit osc:  max|osc - d| = "
        f"{np.abs(limit_infinity - solver).max():.3e}   defect = "
        f"{negative_type_defect(limit_infinity):.3e}"
        f"   t*(osc) = {psd_threshold(limit_infinity)}"
    )
    print(f"      reference t*(d) from part 1 = {star}")
    print(
        f"      {'t':>10} {'defect(D_t)':>13} {'s*(t)':>13} {'max D_t':>10}"
        f" {'s* max D_t':>11} {'min eig G_t':>12} {'triangle':>11}"
    )
    for index, t in enumerate(T_GRID):
        matrix = kernels.soft_distance(float(t))
        soft_star = psd_threshold(matrix)
        defect = negative_type_defect(matrix)
        # G_t is a Gram matrix, hence PSD for every t; spot-check it rather than
        # paying for the (expensive) exponentials at all 41 nodes.
        gram_min = (
            min_eigenvalue(kernels.gram_correlation(float(t)))
            if index % 5 == 0
            else float("nan")
        )
        triangle = triangle_defect(matrix)
        writer.writerow(
            [
                name,
                measure_name,
                f"{t:.6g}",
                f"{defect:.9e}",
                "-" if soft_star is None else f"{soft_star:.9f}",
                f"{matrix.max():.9f}",
                "-" if soft_star is None else f"{soft_star * matrix.max():.9f}",
                f"{gram_min:.6e}",
                f"{triangle:.6e}",
                "-" if star is None else f"{star:.9f}",
            ]
        )
        star_text = "     negtype" if soft_star is None else f"{soft_star:>13.6f}"
        scaled = (
            "          -"
            if soft_star is None
            else f"{soft_star * matrix.max():>11.6f}"
        )
        print(
            f"      {t:>10.4g} {defect:>13.4e} {star_text} {matrix.max():>10.6f}"
            f" {scaled} {gram_min:>12.3e} {triangle:>11.2e}"
        )


def main() -> int:
    report_price_limit()
    report_l2_spectrum()
    families = build_families()
    print("\n3. the D_t bridge:  s*(t) = inf{s : exp(-s D_t) is PSD}")
    with PHASE_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "family",
                "measure",
                "t",
                "negative_type_defect_D_t",
                "s_star",
                "max_D_t",
                "s_star_times_max_D_t",
                "min_eig_G_t_rescaled",
                "triangle_defect_D_t",
                "t_star_of_d",
            ]
        )
        for name in ("cert13", "greedy25", "rand30_61"):
            scan(name, families[name], "full", writer)
        for measure_name in WINDOWS:
            scan("cert13", families["cert13"], measure_name, writer)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
