#!/usr/bin/env python3
"""Fast vectorized geometry of the exchange pseudometric.

The exchange rate between signatures ``g`` and ``f`` is

    C(g -> f) = inf_{beta in [0, inf]} log Z_g(beta) / log Z_f(beta),
    Z_a(beta) = sum_i a_i^beta,

with the ``beta = inf`` value read as ``log max(g) / log max(f)``.  Setting
``u_a(beta) = log log Z_a(beta)`` one has ``-log C(g -> f) = sup_beta (u_f - u_g)``
and therefore

    d(a, b) = -log( C(a -> b) C(b -> a) ) = osc_beta ( u_a - u_b ),

the oscillation (max minus min) of a single smooth function on ``[0, inf]``.
That identity makes ``d`` computable by evaluating one matrix ``U`` of shape
``(#signatures, #betas)`` once and taking row differences, which is orders of
magnitude faster than calling the exact solver pairwise.  The exact solver is
still used to certify every headline number.
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate  # noqa: E402


# ---------------------------------------------------------------- beta grids


def beta_grid(count: int = 1200, beta_max: float = 400.0, beta_min: float = 1e-7) -> np.ndarray:
    """Geometric grid on ``(0, beta_max]`` prefixed by ``beta = 0``."""

    grid = np.geomspace(beta_min, beta_max, count)
    return np.concatenate(([0.0], grid))


def log_log_partition(signatures, betas: np.ndarray) -> np.ndarray:
    """Matrix ``U[i, j] = log log Z_{a_i}(beta_j)`` with an appended beta = inf column."""

    rows = np.empty((len(signatures), betas.size + 1))
    for index, signature in enumerate(signatures):
        entries = np.asarray(signature, dtype=float)
        log_entries = np.log(entries)
        # log-sum-exp keeps large beta stable
        scaled = log_entries[:, None] * betas[None, :]
        peak = scaled.max(axis=0)
        log_z = peak + np.log(np.exp(scaled - peak).sum(axis=0))
        rows[index, : betas.size] = np.log(log_z)
        rows[index, betas.size] = math.log(log_entries.max())
    return rows


def grid_distance_matrix(rows: np.ndarray) -> np.ndarray:
    """Pairwise oscillation distances from a ``U`` matrix (chunked for memory)."""

    size = rows.shape[0]
    matrix = np.zeros((size, size))
    chunk = max(1, int(6e7 // max(1, size * rows.shape[1])))
    for start in range(0, size, chunk):
        stop = min(size, start + chunk)
        block = rows[start:stop, None, :] - rows[None, :, :]
        matrix[start:stop] = block.max(axis=2) - block.min(axis=2)
    np.fill_diagonal(matrix, 0.0)
    return matrix


def grid_distance(rows: np.ndarray, indices) -> np.ndarray:
    """Distance matrix of a small subset, given global ``U`` rows."""

    block = rows[np.asarray(indices)]
    difference = block[:, None, :] - block[None, :, :]
    matrix = difference.max(axis=2) - difference.min(axis=2)
    np.fill_diagonal(matrix, 0.0)
    return matrix


# ------------------------------------------------------------ exact distances


def exact_distance(a, b) -> float:
    """``d(a, b)`` from the repository's exact solver."""

    return -math.log(exchange_rate(a, b) * exchange_rate(b, a))


def exact_distance_matrix(family) -> np.ndarray:
    size = len(family)
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            matrix[i, j] = matrix[j, i] = exact_distance(family[i], family[j])
    return matrix


def dense_distance_matrix(family, points: int = 2_000_001, beta_max: float = 60.0) -> np.ndarray:
    """Independent recomputation of ``d`` on a dense uniform beta grid."""

    betas = np.linspace(0.0, beta_max, points)
    rows = log_log_partition(family, betas)
    return grid_distance_matrix(rows)


# ------------------------------------------------------------- cone functions


def negative_type_matrix(distances: np.ndarray) -> np.ndarray:
    """``-1/2 J D J``; negative type holds exactly when this is PSD."""

    size = distances.shape[0]
    centering = np.eye(size) - np.ones((size, size)) / size
    return -0.5 * centering @ distances @ centering


def helmert(size: int) -> np.ndarray:
    """Orthonormal basis of ``{x : sum(x) = 0}`` as a ``(size-1, size)`` matrix."""

    basis = np.zeros((size - 1, size))
    for row in range(size - 1):
        count = row + 1
        basis[row, :count] = 1.0
        basis[row, count] = -float(count)
        basis[row] /= math.sqrt(count * (count + 1.0))
    return basis


def negative_type_defect(distances: np.ndarray) -> tuple[float, np.ndarray]:
    """Return ``(max_x x^T D x, x)`` over unit vectors with ``sum(x) = 0``.

    Positive value = negative type fails.  Working in an orthonormal basis of the
    centred subspace is essential: ``-1/2 J D J`` always has the constant vector
    in its kernel, so its smallest eigenvalue is capped at zero and carries no
    information about how comfortably negative type holds.
    """

    size = distances.shape[0]
    basis = helmert(size)
    values, vectors = np.linalg.eigh(basis @ distances @ basis.T)
    vector = basis.T @ vectors[:, -1]
    return float(values[-1]), vector


def hypermetric_value(distances: np.ndarray, weights) -> float:
    """``sum_{i<j} b_i b_j d_ij`` for an integer vector ``b`` with ``sum(b) = 1``."""

    b = np.asarray(weights, dtype=float)
    return float(0.5 * b @ distances @ b)


def triangle_defect(distances: np.ndarray) -> float:
    """Largest violation of ``d_ij <= d_ik + d_kj`` (positive means failure)."""

    return float((distances[:, None, :] + distances[None, :, :] - distances[:, :, None]).min()) * -1


if __name__ == "__main__":  # small self-check
    family = [(5, 5, 5, 1), (5, 5, 4, 2), (7, 7, 6, 1), (6, 6, 3), (6, 3, 2)]
    betas = beta_grid()
    rows = log_log_partition(family, betas)
    approximate = grid_distance_matrix(rows)
    exact = exact_distance_matrix(family)
    print("max |grid - exact| =", np.abs(approximate - exact).max())
