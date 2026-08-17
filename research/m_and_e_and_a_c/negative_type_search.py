#!/usr/bin/env python3
"""Search for small families of signatures violating negative type / hypermetricity.

The exchange pseudometric is ``d(a,b) = -log(C(a->b) C(b->a))``.  Negative type
asks for ``x^T D x <= 0`` whenever ``sum(x) = 0``; equivalently ``-1/2 J D J``
is positive semidefinite.  This module searches for finite families where that
fails, and for violations of the hypermetric inequalities
``sum_{i<j} b_i b_j d_ij <= 0`` (``b`` integral, ``sum(b) = 1``).

Search uses the fast grid evaluation of ``d`` (error below 1e-6, always an
underestimate); every reported certificate is recomputed with the repository's
exact solver and with an independent dense-grid evaluation.

Modes::

    python negative_type_search.py pool
    python negative_type_search.py reduce   --core 600 --rounds 20
    python negative_type_search.py exhaust5 --core 70
    python negative_type_search.py anneal   --k 5 --steps 40000
    python negative_type_search.py hyper    --k 5   (pentagonal) / --k 7
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from exchange_geometry import (  # noqa: E402
    beta_grid,
    exact_distance_matrix,
    grid_distance_matrix,
    helmert,
    log_log_partition,
    negative_type_defect,
)

SEARCH_BETAS = beta_grid(count=4000, beta_max=400.0, beta_min=1e-3)


# ------------------------------------------------------------------- pools


def enumerate_signatures(max_length: int, max_entry: int, min_length: int = 2):
    """All decreasing tuples with entries in ``1..max_entry`` and largest entry >= 2."""

    out = []
    for length in range(min_length, max_length + 1):
        for combo in itertools.combinations_with_replacement(range(1, max_entry + 1), length):
            signature = tuple(sorted(combo, reverse=True))
            if signature[0] >= 2:
                out.append(signature)
    return out


def deduplicate(signatures, rows, digits: int = 9):
    """Drop signatures at distance zero from an earlier one.

    ``d(a,b) = 0`` means ``u_a - u_b`` is constant, so the row modulo an additive
    constant is a complete invariant.  Collisions really occur: tensor powers
    satisfy ``log Z_{a^{(k)}} = k log Z_a``, hence ``d(a, a^{(k)}) = 0``.
    """

    seen = set()
    keep = []
    columns = np.linspace(0, rows.shape[1] - 1, 24).astype(int)
    for index in range(len(signatures)):
        key = tuple(np.round(rows[index, columns] - rows[index, 0], digits))
        if key in seen:
            continue
        seen.add(key)
        keep.append(index)
    return [signatures[i] for i in keep], rows[np.array(keep)]


def sample_signatures(count, max_length, max_entry, min_length=2, seed=0):
    """Random decreasing tuples, for pools too large to enumerate."""

    rng = random.Random(seed)
    out = set()
    while len(out) < count:
        length = rng.randint(min_length, max_length)
        signature = tuple(
            sorted((rng.randint(1, max_entry) for _ in range(length)), reverse=True)
        )
        if signature[0] >= 2:
            out.add(signature)
    return sorted(out, reverse=True)


def build_pool(max_length=6, max_entry=12, betas=None, sample=0, min_length=2, seed=0):
    if sample:
        signatures = sample_signatures(sample, max_length, max_entry, min_length, seed)
    else:
        signatures = enumerate_signatures(max_length, max_entry, min_length)
    rows = log_log_partition(signatures, SEARCH_BETAS if betas is None else betas)
    return deduplicate(signatures, rows)


# --------------------------------------------------------------- objectives


def defect(distances: np.ndarray) -> float:
    """``max_{sum x = 0, |x| = 1} x^T D x``; positive means negative type fails.

    Computed as the largest eigenvalue of ``Q D Q^T`` for an orthonormal basis
    ``Q`` of the centred hyperplane, which keeps the trivial constant direction
    (always an eigenvector of ``-1/2 J D J`` with eigenvalue 0) out of the way.
    """

    basis = helmert(distances.shape[0])
    return float(np.linalg.eigvalsh(basis @ distances @ basis.T)[-1])


def batched_defect(blocks: np.ndarray) -> np.ndarray:
    basis = helmert(blocks.shape[-1])
    return np.linalg.eigvalsh(np.einsum("ai,...ij,bj->...ab", basis, blocks, basis))[..., -1]


def _det2(m):
    return m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0]


def _det3(m):
    return (
        m[..., 0, 0] * (m[..., 1, 1] * m[..., 2, 2] - m[..., 1, 2] * m[..., 2, 1])
        - m[..., 0, 1] * (m[..., 1, 0] * m[..., 2, 2] - m[..., 1, 2] * m[..., 2, 0])
        + m[..., 0, 2] * (m[..., 1, 0] * m[..., 2, 1] - m[..., 1, 1] * m[..., 2, 0])
    )


def _det4(m):
    total = 0.0
    for column in range(4):
        keep = [c for c in range(4) if c != column]
        minor = m[..., 1:, :][..., keep]
        sign = -1.0 if column % 2 else 1.0
        total = total + sign * m[..., 0, column] * _det3(minor)
    return total


def batched_violation_mask(blocks: np.ndarray) -> np.ndarray:
    """Cheap necessary test for ``lambda_max(Q D Q^T) > 0`` on stacks of blocks.

    ``-Q D Q^T`` is positive definite exactly when all its leading principal
    minors are positive; checking those with closed-form determinants runs two
    orders of magnitude faster than a batched eigenvalue solve, and the (rare)
    survivors are then handed to the exact routine.
    """

    basis = helmert(blocks.shape[-1])
    matrix = -np.einsum("ai,...ij,bj->...ab", basis, blocks, basis)
    order = matrix.shape[-1]
    mask = matrix[..., 0, 0] <= 0
    if order >= 2:
        mask |= _det2(matrix[..., :2, :2]) <= 0
    if order >= 3:
        mask |= _det3(matrix[..., :3, :3]) <= 0
    if order >= 4:
        mask |= _det4(matrix[..., :4, :4]) <= 0
    if order >= 5:  # fall back to the eigen test for larger blocks
        mask |= np.linalg.eigvalsh(-matrix)[..., -1] > 0
    return mask


def hypermetric_patterns(k: int, positives: int):
    """All +-1 patterns of length ``k`` with ``positives`` ones (sum = 1 needs 2p = k+1)."""

    patterns = []
    for chosen in itertools.combinations(range(k), positives):
        b = -np.ones(k)
        b[list(chosen)] = 1.0
        patterns.append(b)
    return np.array(patterns)


def batched_hypermetric(blocks: np.ndarray, patterns: np.ndarray):
    """Best ``sum_{i<j} b_i b_j d_ij`` over the given patterns, per block."""

    values = 0.5 * np.einsum("pi,bij,pj->bp", patterns, blocks, patterns)
    return values.max(axis=1), values.argmax(axis=1)


def subset(distances, indices):
    idx = np.asarray(indices)
    return distances[np.ix_(idx, idx)]


# ----------------------------------------------------------------- searches


VIOLATION_TOLERANCE = 1e-9


def eigen_shrink(distances, indices, floor=5, tolerance=VIOLATION_TOLERANCE):
    """Shrink a violating family using the most-negative eigenvector's support.

    Keeps the coordinates carrying the largest weight, halving the size while the
    violation survives.  Much faster than one-at-a-time greedy on large cores.
    """

    current = list(indices)
    while len(current) > floor:
        block = subset(distances, current)
        size = len(current)
        basis = helmert(size)
        values, vectors = np.linalg.eigh(basis @ block @ basis.T)
        if values[-1] <= tolerance:
            return current
        weight = np.abs(basis.T @ vectors[:, -1])
        order = np.argsort(-weight)
        target = max(floor, size // 2)
        found = False
        while target < size:
            trial = [current[i] for i in sorted(order[:target])]
            if defect(subset(distances, trial)) > tolerance:
                current = trial
                found = True
                break
            target = target + max(1, (size - target) // 2)
        if not found:
            break
    return current


def greedy_reduce(distances, indices, floor=5, tolerance=VIOLATION_TOLERANCE):
    """Drop points one at a time while the family still violates negative type."""

    current = list(indices)
    while len(current) > floor:
        best = None
        for position in range(len(current)):
            trial = current[:position] + current[position + 1 :]
            value = defect(subset(distances, trial))
            if value > tolerance and (best is None or value > best[0]):
                best = (value, trial)
        if best is None:
            break
        current = best[1]
    return current


def mode_reduce(signatures, rows, core_size, rounds, seed, floor=5):
    rng = random.Random(seed)
    best_overall = None
    histogram = {}
    for round_index in range(rounds):
        sample = rng.sample(range(len(signatures)), core_size)
        distances = grid_distance_matrix(rows[np.array(sample)])
        value = defect(distances)
        if value <= VIOLATION_TOLERANCE:
            print(f"  round {round_index:>3}: core defect {value:.3e} (no violation)")
            continue
        shrunk = eigen_shrink(distances, list(range(core_size)), floor)
        reduced_local = greedy_reduce(distances, shrunk, floor)
        reduced = [sample[i] for i in reduced_local]
        final = defect(subset(distances, reduced_local))
        histogram[len(reduced)] = histogram.get(len(reduced), 0) + 1
        print(
            f"  round {round_index:>3}: core defect {value:.3e} -> {len(reduced)} points, "
            f"defect {final:.3e}",
            flush=True,
        )
        if best_overall is None or (len(reduced), -final) < (
            len(best_overall[0]),
            -best_overall[1],
        ):
            best_overall = (reduced, final)
    print("  reduced-size histogram:", dict(sorted(histogram.items())))
    return best_overall


def mode_exhaust(signatures, rows, core, k=5, chunk=100_000, patterns=None):
    """Exhaustive sweep of all k-subsets of ``core``."""

    core = list(core)
    distances = grid_distance_matrix(rows[np.array(core)])
    n = len(core)
    total = math.comb(n, k)
    print(f"  sweeping {total} subsets of size {k} from a core of {n}")
    best_defect = (-math.inf, None)
    best_hyper = (-math.inf, None, None)
    combos = itertools.combinations(range(n), k)
    processed = 0
    while True:
        batch = list(itertools.islice(combos, chunk))
        if not batch:
            break
        idx = np.array(batch)
        blocks = distances[idx[:, :, None], idx[:, None, :]]
        defects = batched_defect(blocks)
        pick = int(defects.argmax())
        if defects[pick] > best_defect[0]:
            best_defect = (float(defects[pick]), [core[i] for i in idx[pick]])
        if patterns is not None:
            values, which = batched_hypermetric(blocks, patterns)
            pick = int(values.argmax())
            if values[pick] > best_hyper[0]:
                best_hyper = (
                    float(values[pick]),
                    [core[i] for i in idx[pick]],
                    patterns[which[pick]].astype(int).tolist(),
                )
        processed += len(batch)
        if processed % (chunk * 20) == 0:
            print(f"    {processed}/{total}", flush=True)
    return best_defect, best_hyper


def anneal(distances_getter, n, k, steps, restarts, seed, score):
    """Generic simulated annealing over k-subsets of ``range(n)``."""

    rng = random.Random(seed)
    overall = (-math.inf, None)
    for restart in range(restarts):
        current = rng.sample(range(n), k)
        value = score(distances_getter(current))
        best = (value, list(current))
        scale = max(abs(value), 1e-3)
        for step in range(steps):
            temperature = scale * (1 - step / steps) ** 3 + 1e-14
            trial = list(current)
            trial[rng.randrange(k)] = rng.randrange(n)
            if len(set(trial)) < k:
                continue
            candidate = score(distances_getter(trial))
            if candidate > value or rng.random() < math.exp(
                min(0.0, (candidate - value) / temperature)
            ):
                current, value = trial, candidate
                if value > best[0]:
                    best = (value, list(current))
        print(f"  restart {restart}: best {best[0]:.6e}", flush=True)
        if best[0] > overall[0]:
            overall = best
    return overall


# ------------------------------------------------------------------- driver


def report_family(signatures, indices, patterns=None):
    family = [signatures[i] for i in indices]
    exact = exact_distance_matrix(family)
    value, vector = negative_type_defect(exact)
    print("  family:", family)
    print(f"  exact negative-type defect: {value:.6e}")
    print("  x =", np.round(vector, 8).tolist())
    return family, exact, value, vector


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode")
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--steps", type=int, default=40_000)
    parser.add_argument("--restarts", type=int, default=6)
    parser.add_argument("--rounds", type=int, default=20)
    parser.add_argument("--core", type=int, default=600)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    parser.add_argument("--floor", type=int, default=5)
    args = parser.parse_args()

    signatures, rows = build_pool(args.max_length, args.max_entry)
    print(
        f"pool: {len(signatures)} distinct signatures "
        f"(length <= {args.max_length}, entries <= {args.max_entry}), "
        f"{rows.shape[1]} beta samples"
    )

    if args.mode == "pool":
        return 0

    if args.mode == "reduce":
        best = mode_reduce(signatures, rows, args.core, args.rounds, args.seed, args.floor)
        if best:
            print(f"\nsmallest family found: {len(best[0])} points, grid defect {best[1]:.6e}")
            report_family(signatures, best[0])
        return 0

    if args.mode in {"exhaust5", "exhaust"}:
        rng = random.Random(args.seed)
        core = rng.sample(range(len(signatures)), args.core)
        patterns = hypermetric_patterns(args.k, (args.k + 1) // 2) if args.k % 2 == 1 else None
        best_defect, best_hyper = mode_exhaust(signatures, rows, core, args.k, patterns=patterns)
        print(f"\nbest negative-type defect over {args.k}-subsets: {best_defect[0]:.6e}")
        if best_defect[1]:
            print("   ", [signatures[i] for i in best_defect[1]])
        if patterns is not None:
            print(f"best hypermetric value: {best_hyper[0]:.6e}  b = {best_hyper[2]}")
            if best_hyper[1]:
                print("   ", [signatures[i] for i in best_hyper[1]])
        return 0

    if args.mode == "anneal":
        distances_getter = lambda idx: grid_distance_matrix(rows[np.array(idx)])  # noqa: E731
        value, indices = anneal(
            distances_getter, len(signatures), args.k, args.steps, args.restarts, args.seed, defect
        )
        print(f"\nbest grid defect for k = {args.k}: {value:.6e}")
        if indices:
            report_family(signatures, indices)
        return 0

    if args.mode == "hyper":
        patterns = hypermetric_patterns(args.k, (args.k + 1) // 2)
        score = lambda D: float(  # noqa: E731
            (0.5 * np.einsum("pi,ij,pj->p", patterns, D, patterns)).max()
        )
        distances_getter = lambda idx: grid_distance_matrix(rows[np.array(idx)])  # noqa: E731
        value, indices = anneal(
            distances_getter, len(signatures), args.k, args.steps, args.restarts, args.seed, score
        )
        print(f"\nbest hypermetric value for k = {args.k}: {value:.6e}")
        if indices:
            family = [signatures[i] for i in indices]
            exact = exact_distance_matrix(family)
            values = 0.5 * np.einsum("pi,ij,pj->p", patterns, exact, patterns)
            best = int(values.argmax())
            print("  family:", family)
            print(f"  exact hypermetric value: {values[best]:.6e}")
            print("  b =", patterns[best].astype(int).tolist())
        return 0

    raise SystemExit(f"unknown mode {args.mode}")


if __name__ == "__main__":
    raise SystemExit(main())
