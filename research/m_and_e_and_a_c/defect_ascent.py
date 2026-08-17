#!/usr/bin/env python3
"""Alternating maximization of the negative-type defect over k-point families.

The defect of a family ``F`` is ``max{ x^T D_F x : sum(x) = 0, |x| = 1 }``, so
the quantity to maximize over both the family and the vector is

    max_F max_x  x^T D_F x .

With ``x`` held fixed the objective is *linear* in the row of ``D`` belonging to
any one point, so the best replacement for that point is a single argmax over
the entire core; with ``F`` held fixed the best ``x`` is the top eigenvector of
``Q D_F Q^T`` in an orthonormal basis of the centred hyperplane.  Alternating
the two is far stronger than simulated annealing on the raw defect, which drifts
towards families of tiny diameter where every value is near zero.

Usage::

    python defect_ascent.py --k 5 --restarts 2000 --core 2000 --seed 0
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from campaign import load_core
from exchange_geometry import exact_distance_matrix, helmert, negative_type_defect


def top_direction(block):
    basis = helmert(block.shape[0])
    values, vectors = np.linalg.eigh(basis @ block @ basis.T)
    return float(values[-1]), basis.T @ vectors[:, -1]


def seeded_ascent(distances, index, sweeps=40, normalize=False):
    """Run the alternating ascent from a given family; returns (defect, family)."""

    index = list(index)
    k = len(index)
    block = distances[np.ix_(index, index)]
    value, x = top_direction(block)
    for _ in range(sweeps):
        moved = False
        for position in range(k):
            others = [i for i in range(k) if i != position]
            columns = distances[:, [index[i] for i in others]]
            rest = block[np.ix_(others, others)]
            form = float(x[others] @ rest @ x[others]) + 2.0 * x[position] * (
                columns @ x[others]
            )
            if normalize:
                total = float(rest.sum()) + 2.0 * columns.sum(axis=1)
                gains = form / np.maximum(total / (k * (k - 1)), 1e-15)
            else:
                gains = form
            gains[index] = -np.inf
            candidate = int(gains.argmax())
            if candidate != index[position]:
                index[position] = candidate
                moved = True
        block = distances[np.ix_(index, index)]
        value, x = top_direction(block)
        if not moved:
            break
    return value, index


def ascent(distances, k, rng, restarts, sweeps=25, normalize=True):
    """Return ``(best relative defect, best absolute defect, family)``."""

    n = distances.shape[0]
    best = (-np.inf, -np.inf, None)
    for _ in range(restarts):
        index = rng.sample(range(n), k)
        block = distances[np.ix_(index, index)]
        value, x = top_direction(block)
        for _ in range(sweeps):
            moved = False
            for position in range(k):
                others = [i for i in range(k) if i != position]
                columns = distances[:, [index[i] for i in others]]
                rest = block[np.ix_(others, others)]
                fixed_form = float(x[others] @ rest @ x[others])
                fixed_sum = float(rest.sum())
                # x^T D x and the mean distance are both affine in the row of the
                # candidate that replaces position p, so the exact objective can
                # be evaluated for every candidate at once.
                form = fixed_form + 2.0 * x[position] * (columns @ x[others])
                if normalize:
                    total = fixed_sum + 2.0 * columns.sum(axis=1)
                    gains = form / np.maximum(total / (k * (k - 1)), 1e-15)
                else:
                    gains = form
                gains[index] = -np.inf
                candidate = int(gains.argmax())
                if candidate != index[position]:
                    index[position] = candidate
                    moved = True
            block = distances[np.ix_(index, index)]
            value, x = top_direction(block)
            if not moved:
                break

        block = distances[np.ix_(index, index)]
        absolute, _ = top_direction(block)
        scale = block.sum() / (k * (k - 1))
        rel = absolute / scale if scale > 0 else -np.inf
        if rel > best[0]:
            best = (rel, absolute, list(index))
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--core", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--restarts", type=int, default=500)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    parser.add_argument("--absolute", action="store_true")
    args = parser.parse_args()

    core, distances = load_core(args.seed, args.core, args.max_length, args.max_entry)
    print(f"core {len(core)} signatures (seed {args.seed}); k = {args.k}")
    rng = random.Random(args.seed + 101)
    rel, absolute, index = ascent(
        distances, args.k, rng, args.restarts, normalize=not args.absolute
    )
    family = [core[i] for i in index]
    exact = exact_distance_matrix(family)
    exact_value, vector = negative_type_defect(exact)
    print(f"  best relative defect {rel:+.6e}, absolute {absolute:+.6e}")
    print("  family:", family)
    print(f"  exact defect {exact_value:+.6e}")
    if exact_value > 0:
        print("  x =", np.round(vector, 9).tolist())
        print("  VIOLATION")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
