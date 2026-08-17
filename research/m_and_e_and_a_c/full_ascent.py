#!/usr/bin/env python3
"""Alternating ascent for the negative-type defect over the *whole* pool.

``defect_ascent`` works on a cached core of a couple of thousand signatures.
This module runs the same alternating maximization

    max over families F of size k,  max over x with sum(x) = 0, |x| = 1,  x^T D_F x

but lets each coordinate range over all ~18.5k pool signatures.  The distance
column from every pool member to one fixed signature is a single vectorized
oscillation over the ``log log Z`` rows; those columns are cached, so a sweep
costs a handful of them.

Usage::

    python full_ascent.py --k 5 --restarts 40 --seed 0
    python full_ascent.py --k 5 --seeds-from six_point_families.json
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from exchange_geometry import exact_distance_matrix, helmert, negative_type_defect
from negative_type_search import build_pool


class ColumnCache:
    def __init__(self, rows):
        self.rows = rows
        self.cache = {}

    def column(self, index):
        column = self.cache.get(index)
        if column is None:
            difference = self.rows - self.rows[index]
            column = difference.max(axis=1) - difference.min(axis=1)
            self.cache[index] = column
        return column

    def block(self, indices):
        size = len(indices)
        out = np.zeros((size, size))
        for a in range(size):
            column = self.column(indices[a])
            for b in range(size):
                out[a, b] = column[indices[b]]
        np.fill_diagonal(out, 0.0)
        return out


def top_direction(block):
    basis = helmert(block.shape[0])
    values, vectors = np.linalg.eigh(basis @ block @ basis.T)
    return float(values[-1]), basis.T @ vectors[:, -1]


def ascend(cache, index, sweeps=30):
    index = list(index)
    k = len(index)
    block = cache.block(index)
    value, x = top_direction(block)
    for _ in range(sweeps):
        moved = False
        for position in range(k):
            others = [i for i in range(k) if i != position]
            columns = np.stack([cache.column(index[i]) for i in others], axis=1)
            rest = block[np.ix_(others, others)]
            form = float(x[others] @ rest @ x[others]) + 2.0 * x[position] * (
                columns @ x[others]
            )
            form[index] = -np.inf
            candidate = int(form.argmax())
            if candidate != index[position]:
                index[position] = candidate
                moved = True
                block = cache.block(index)
                value, x = top_direction(block)
        if not moved:
            break
    return value, index


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--restarts", type=int, default=30)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    parser.add_argument("--seeds-from", default="")
    args = parser.parse_args()

    signatures, rows = build_pool(args.max_length, args.max_entry)
    lookup = {signature: i for i, signature in enumerate(signatures)}
    print(f"pool: {len(signatures)} signatures, k = {args.k}", flush=True)
    cache = ColumnCache(rows)
    rng = random.Random(args.seed + 31)

    starts = []
    if args.seeds_from:
        families = json.loads(Path(args.seeds_from).read_text())
        for family in families:
            index = [lookup[tuple(s)] for s in family]
            if len(index) == args.k:
                starts.append(index)
            elif len(index) == args.k + 1:  # every k-subset of a bigger family
                for drop in range(len(index)):
                    starts.append([index[i] for i in range(len(index)) if i != drop])
    starts.extend(rng.sample(range(len(signatures)), args.k) for _ in range(args.restarts))

    best = (-np.inf, None)
    for number, start in enumerate(starts):
        value, index = ascend(cache, start)
        if value > best[0]:
            best = (value, list(index))
        print(
            f"  start {number:>3}: defect {value:+.6e}   (best {best[0]:+.6e}, "
            f"{len(cache.cache)} columns cached)",
            flush=True,
        )

    family = [signatures[i] for i in best[1]]
    exact = exact_distance_matrix(family)
    exact_value, vector = negative_type_defect(exact)
    print(f"\nbest defect for k = {args.k}: grid {best[0]:+.6e}, exact {exact_value:+.6e}")
    print("family:", family)
    if exact_value > 0:
        print("x =", np.round(vector, 9).tolist())
        print("VIOLATION")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
