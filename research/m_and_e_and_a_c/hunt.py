#!/usr/bin/env python3
"""The protocol that actually finds small violations, run on a given pool.

Each round:

1. draw a random core and test it for negative type;
2. if it violates, shrink it to an inclusion-minimal violating family;
3. use every k-subset of that family as a seed for the alternating ascent
   (family and witness updated in turn) over the whole pool.

Reporting the same protocol on restricted pools makes the comparison fair:
purely random ascent starts almost never reach a violation even on the full
pool, so a "no violation found" from random starts alone means very little.

Usage::  python hunt.py --k 5 --rounds 8 --max-length 6 --max-entry 12
"""

from __future__ import annotations

import argparse
import itertools
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from exchange_geometry import exact_distance_matrix, grid_distance_matrix, negative_type_defect
from full_ascent import ColumnCache, ascend
from negative_type_search import build_pool, defect, eigen_shrink, greedy_reduce, subset


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--rounds", type=int, default=8)
    parser.add_argument("--core", type=int, default=300)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    parser.add_argument("--max-seeds", type=int, default=12)
    parser.add_argument("--sample", type=int, default=0, help="sample a pool instead of enumerating")
    parser.add_argument("--min-length", type=int, default=2)
    args = parser.parse_args()

    signatures, rows = build_pool(
        args.max_length,
        args.max_entry,
        sample=args.sample,
        min_length=args.min_length,
        seed=args.seed,
    )
    print(
        f"pool: {len(signatures)} signatures "
        f"(length <= {args.max_length}, entries <= {args.max_entry}), k = {args.k}",
        flush=True,
    )
    cache = ColumnCache(rows)
    rng = random.Random(args.seed + 4211)
    best = (-np.inf, None)
    successes = 0
    attempts = 0

    for round_index in range(args.rounds):
        sample = rng.sample(range(len(signatures)), min(args.core, len(signatures)))
        distances = grid_distance_matrix(rows[np.array(sample)])
        if defect(distances) <= 1e-9:
            print(f"  round {round_index}: core of {len(sample)} is of negative type", flush=True)
            continue
        shrunk = eigen_shrink(distances, list(range(len(sample))), floor=args.k)
        reduced = greedy_reduce(distances, shrunk, floor=args.k)
        minimal = [sample[i] for i in reduced]
        print(
            f"  round {round_index}: inclusion-minimal family of {len(minimal)}, "
            f"defect {defect(subset(distances, reduced)):+.4e}",
            flush=True,
        )
        seeds = list(itertools.combinations(minimal, args.k))
        rng.shuffle(seeds)
        for seed in seeds[: args.max_seeds]:
            attempts += 1
            value, index = ascend(cache, list(seed))
            if value > 1e-9:
                successes += 1
            if value > best[0]:
                best = (value, list(index))
            print(f"      seed -> defect {value:+.6e}", flush=True)

    print(f"\nascent runs {attempts}, of which violating: {successes}")
    if best[1] is None:
        print("no family found")
        return 0
    family = [signatures[i] for i in best[1]]
    exact = exact_distance_matrix(family)
    exact_value, vector = negative_type_defect(exact)
    print(f"best defect: grid {best[0]:+.6e}, exact {exact_value:+.6e}")
    print("family:", family)
    if exact_value > 0:
        print("x =", np.round(vector, 9).tolist())
        print("VIOLATION")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
