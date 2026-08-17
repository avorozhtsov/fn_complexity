#!/usr/bin/env python3
"""Shared machinery for T1.2 (the L^t interpolation between E and M).

Objects
-------
``phi_a(beta) = log Z_a(beta)``, ``Z_a(beta) = sum_i a_i**beta``.
``u_a = log phi_a`` (defined for beta > 0; ``phi_a(0) = log(len(a)) > 0`` too,
so ``u_a`` is finite on all of ``[0, inf)`` as soon as ``len(a) >= 2``).

``C(a -> b) = inf_beta phi_a/phi_b`` (repo solver, ~1e-13),
``d(a,b) = -log(C(a->b) C(b->a)) = osc_beta (u_a - u_b)``.

Families
--------
``FAMILIES`` holds the published 13-signature negative-type certificate plus a
set of pseudo-random families of sizes 5..40 drawn from the pool of decreasing
signatures with length 2..6 and entries 1..12.  Members are filtered so that
all pairwise ``d`` exceed 1e-9 (the pool contains genuine distance-zero pairs,
e.g. ``a`` and ``a (x) a``, which would make every ``K_t`` singular).
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate  # noqa: E402

ZERO_TOLERANCE = 1e-10  # eigenvalues smaller than this in modulus are "zero"

# ---------------------------------------------------------------- signatures

CERTIFICATE_FAMILY: list[tuple[int, ...]] = [
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


def label(signature: tuple[int, ...]) -> str:
    return "{" + ",".join(map(str, signature)) + "}"


def signature_pool(
    lengths: tuple[int, ...] = (2, 3, 4, 5, 6),
    max_entry: int = 12,
) -> list[tuple[int, ...]]:
    """All decreasing signatures of the given lengths with entries 1..max_entry.

    Constant signatures (all fibers of size one) are dropped: their ``phi`` is
    constant and every distance to a non-constant signature is infinite.
    """

    pool: list[tuple[int, ...]] = []
    for length in lengths:
        for combo in itertools.combinations_with_replacement(
            range(max_entry, 0, -1), length
        ):
            if combo[0] == 1:
                continue
            pool.append(tuple(combo))
    return pool


def sample_family(
    size: int,
    seed: int,
    pool: list[tuple[int, ...]] | None = None,
) -> list[tuple[int, ...]]:
    """Draw ``size`` signatures at pairwise exchange distance > 1e-9."""

    pool = pool if pool is not None else signature_pool()
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(pool))
    chosen: list[tuple[int, ...]] = []
    for index in order:
        candidate = pool[int(index)]
        if all(irreversibility(candidate, member) > 1e-9 for member in chosen):
            chosen.append(candidate)
        if len(chosen) == size:
            return chosen
    raise RuntimeError(f"pool exhausted before reaching size {size}")


# ------------------------------------------------------------- the metric d


def irreversibility(a: tuple[int, ...], b: tuple[int, ...]) -> float:
    """``d(a,b) = -log(C(a->b) C(b->a))`` from the repo solver."""

    return -math.log(exchange_rate(a, b) * exchange_rate(b, a))


def distance_matrix(family: list[tuple[int, ...]]) -> np.ndarray:
    size = len(family)
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            matrix[i, j] = matrix[j, i] = irreversibility(family[i], family[j])
    return matrix


# ------------------------------------------- independent grid recomputation


def log_partition(signature: tuple[int, ...], betas: np.ndarray) -> np.ndarray:
    """``phi_a(beta)`` evaluated stably on a grid."""

    logs = np.log(np.asarray(signature, dtype=float))[:, None]
    scaled = logs * betas[None, :]
    largest = scaled.max(axis=0)
    return largest + np.log(np.exp(scaled - largest).sum(axis=0))


def u_values(signature: tuple[int, ...], betas: np.ndarray) -> np.ndarray:
    """``u_a = log phi_a`` on a grid."""

    return np.log(log_partition(signature, betas))


def distance_matrix_grid(
    family: list[tuple[int, ...]],
    beta_max: float = 60.0,
    points: int = 2_000_000,
    include_infinity: bool = True,
) -> np.ndarray:
    """``d`` recomputed as ``osc_beta(u_a - u_b)`` on a dense grid.

    Independent of the solver; used to certify the headline claims.
    """

    betas = np.linspace(1e-9, beta_max, points)
    table = np.array([u_values(signature, betas) for signature in family])
    at_infinity = np.array(
        [math.log(math.log(max(signature))) for signature in family]
    )
    size = len(family)
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            difference = table[i] - table[j]
            edge = at_infinity[i] - at_infinity[j]
            high = max(float(difference.max()), edge)
            low = min(float(difference.min()), edge)
            matrix[i, j] = matrix[j, i] = high - low
    return matrix


# ------------------------------------------------------------- linear algebra


def min_eigenvalue(matrix: np.ndarray) -> float:
    return float(np.linalg.eigvalsh(matrix).min())


def kernel_min_eigenvalue(distances: np.ndarray, t: float) -> float:
    """Smallest eigenvalue of ``K_t = exp(-t d)``."""

    return min_eigenvalue(np.exp(-t * distances))


def _centred(distances: np.ndarray) -> np.ndarray:
    """``-1/2 J d J`` with ``J`` the centring projector."""

    size = distances.shape[0]
    centring = np.eye(size) - np.ones((size, size)) / size
    return -0.5 * centring @ distances @ centring


def negative_type_eigenvalue(distances: np.ndarray) -> float:
    """Smallest eigenvalue of ``-1/2 J d J``.

    ``d`` is of negative type exactly when this is >= 0.  The value is
    1-homogeneous in ``d``, so it must be normalised before being compared
    across matrices of different scale -- see ``negative_type_defect``.
    """

    return min_eigenvalue(_centred(distances))


def negative_type_ratio(distances: np.ndarray) -> float:
    """Signed scale-free ``lambda_min / lambda_max`` of ``-1/2 J d J``.

    Non-negative exactly when ``d`` is of negative type.  Unlike the defect it
    keeps varying inside the negative-type region, so a search can descend on it.
    """

    values = np.linalg.eigvalsh(_centred(distances))
    top = float(values.max())
    return float(values.min()) / top if top > 0 else float("-inf")


def negative_type_defect(distances: np.ndarray) -> float:
    """Scale-free violation ``-lambda_min / lambda_max`` of ``-1/2 J d J``.

    Zero (up to roundoff) exactly when ``d`` is of negative type.  Invariant
    under ``d -> c d``, so it is the right quantity to compare across the
    ``L^t`` family, whose scale varies by orders of magnitude.
    """

    values = np.linalg.eigvalsh(_centred(distances))
    top = float(values.max())
    if top <= 0:
        return float("inf")
    return max(0.0, -float(values.min()) / top)


def psd_threshold(
    distances: np.ndarray,
    low: float = 1e-6,
    high: float = 1e6,
    scan_points: int = 481,
    defect_tolerance: float = 1e-12,
) -> float | None:
    """``t* = sup {t : exp(-t d) is NOT PSD}``.

    ``d`` is rescaled to unit diameter first (``t*(c d) = t*(d)/c`` exactly), so
    the bracket means the same thing at every input scale.  The threshold is
    located by scanning a log grid for the *largest* ``t`` at which ``K_t`` has
    a resolvable negative eigenvalue and then bisecting; scanning from above is
    essential because ``lambda_min(K_t) -> 0`` as ``t -> 0`` (``K_t -> J``), so a
    small-``t`` bracket cannot be certified numerically.

    Returns ``None`` when ``d`` is of negative type, i.e. when ``exp(-t d)`` is
    PSD for every ``t > 0`` (Schoenberg).
    """

    scale = float(distances.max())
    if scale <= 0:
        return None
    unit = distances / scale
    if negative_type_defect(unit) <= defect_tolerance:
        return None
    grid = np.exp(np.linspace(math.log(low), math.log(high), scan_points))
    negative = [
        index
        for index, t in enumerate(grid)
        if kernel_min_eigenvalue(unit, float(t)) < -ZERO_TOLERANCE
    ]
    if not negative:
        return None
    last = negative[-1]
    if last == len(grid) - 1:
        raise RuntimeError("no PSD point found below the upper bracket")
    left, right = float(grid[last]), float(grid[last + 1])
    for _ in range(200):
        middle = math.sqrt(left * right)
        if kernel_min_eigenvalue(unit, middle) < -ZERO_TOLERANCE:
            left = middle
        else:
            right = middle
        if right - left < 1e-13 * right:
            break
    return math.sqrt(left * right) / scale


def second_order_estimate(distances: np.ndarray) -> float:
    """``max 2 (x^T d x) / (x^T d^{o2} x)`` over centred ``x``: a closed-form scale
    for ``t*``.

    Truncating ``x^T exp(-t d) x = sum_{k>=1} (-t)^k x^T d^{ok} x / k!`` (valid
    for centred ``x``) after two terms, the certificate ``x`` shows failure for
    ``t < 2 (x^T d x) / (x^T d^{o2} x)``.  The truncation is not controlled, so
    this is a heuristic scale, not a bound.
    """

    size = distances.shape[0]
    centring = np.eye(size) - np.ones((size, size)) / size
    values, vectors = np.linalg.eigh(_centred(distances))
    best = 0.0
    squared = distances**2
    for index in range(size):
        x = centring @ vectors[:, index]
        norm = np.linalg.norm(x)
        if norm < 1e-12:
            continue
        x = x / norm
        numerator = float(x @ distances @ x)
        denominator = float(x @ squared @ x)
        if numerator > 0 and denominator > 0:
            best = max(best, 2 * numerator / denominator)
    return best


# ------------------------------------------------- cached working pool + search

WORKING_POOL_NPZ = HERE / "t1_2_working_pool.npz"


def working_pool(size: int = 400, seed: int = 7) -> tuple[list[tuple[int, ...]], np.ndarray]:
    """A cached pool of signatures together with its full distance matrix."""

    if WORKING_POOL_NPZ.exists():
        blob = np.load(WORKING_POOL_NPZ, allow_pickle=True)
        if int(blob["size"]) == size and int(blob["seed"]) == seed:
            members = [tuple(int(v) for v in row if v > 0) for row in blob["members"]]
            return members, blob["distances"]
    members = sample_family(size, seed)
    distances = distance_matrix(members)
    width = max(len(m) for m in members)
    padded = np.zeros((len(members), width), dtype=int)
    for index, member in enumerate(members):
        padded[index, : len(member)] = member
    np.savez(
        WORKING_POOL_NPZ,
        members=padded,
        distances=distances,
        size=size,
        seed=seed,
    )
    return members, distances


def greedy_violating_subset(
    distances: np.ndarray,
    target_size: int,
    seed: int,
    restarts: int = 40,
    sweeps: int = 12,
) -> np.ndarray:
    """Hill-climb a subset of the given size minimising ``negative_type_ratio``.

    The scale-free ratio is the right objective: the absolute eigenvalue is
    1-homogeneous in ``d`` and a search on it simply collapses onto whichever
    sub-family has the smallest diameter.
    """

    rng = np.random.default_rng(seed)
    population = distances.shape[0]
    best_indices: np.ndarray | None = None
    best_value = np.inf
    for _ in range(restarts):
        indices = rng.choice(population, size=target_size, replace=False)
        value = negative_type_ratio(distances[np.ix_(indices, indices)])
        for _ in range(sweeps):
            improved = False
            for position in range(target_size):
                outside = rng.choice(population, size=24, replace=False)
                for candidate in outside:
                    if candidate in indices:
                        continue
                    trial = indices.copy()
                    trial[position] = candidate
                    trial_value = negative_type_ratio(
                        distances[np.ix_(trial, trial)]
                    )
                    if trial_value < value - 1e-15:
                        indices, value, improved = trial, trial_value, True
            if not improved:
                break
        if value < best_value:
            best_value, best_indices = value, indices
    assert best_indices is not None
    return best_indices


def build_families() -> dict[str, list[tuple[int, ...]]]:
    pool = signature_pool()
    families: dict[str, list[tuple[int, ...]]] = {"cert13": CERTIFICATE_FAMILY}
    for size, seed in [
        (5, 11),
        (5, 12),
        (8, 21),
        (8, 22),
        (10, 31),
        (13, 41),
        (13, 42),
        (20, 51),
        (20, 52),
        (25, 55),
        (30, 61),
        (30, 62),
        (40, 71),
        (40, 72),
    ]:
        families[f"rand{size}_{seed}"] = sample_family(size, seed, pool)
    members, distances = working_pool()
    for size, seed in [
        (5, 100),
        (6, 101),
        (8, 102),
        (10, 103),
        (12, 107),
        (13, 104),
        (16, 108),
        (20, 105),
        (25, 109),
        (30, 106),
        (40, 110),
    ]:
        indices = greedy_violating_subset(distances, size, seed)
        families[f"greedy{size}"] = [members[int(i)] for i in indices]
    return families
