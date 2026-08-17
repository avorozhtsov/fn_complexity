#!/usr/bin/env python3
"""How few signatures can violate negative type?

Pipeline per round:

1. draw a random core from the pool and compute its grid distance matrix;
2. shrink the violating core with the most-negative eigenvector, then drop
   points one at a time while the violation survives (inclusion-minimal family);
3. for every inclusion-minimal family of size six, try hard to get to five:
   each of its six 5-subsets seeds a coordinate ascent that replaces one member
   at a time by the best of thousands of candidates.

Everything reported is recomputed with the exact solver at the end.

Usage::

    python minimality.py --rounds 40 --core 400 --refine-pool 2000
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from campaign import load_core
from defect_ascent import seeded_ascent
from exchange_geometry import exact_distance_matrix, grid_distance_matrix, negative_type_defect
from negative_type_search import build_pool, defect, eigen_shrink, greedy_reduce, subset

OUTPUT = Path(__file__).resolve().parent / "minimality_results.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rounds", type=int, default=40)
    parser.add_argument("--core", type=int, default=400)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--refine-pool", type=int, default=2000)
    parser.add_argument("--refine-rounds", type=int, default=4)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    args = parser.parse_args()

    signatures, rows = build_pool(args.max_length, args.max_entry)
    print(f"pool: {len(signatures)} signatures", flush=True)
    refine_core, refine_distances = load_core(0, args.refine_pool, args.max_length, args.max_entry)
    refine_index = {signature: i for i, signature in enumerate(refine_core)}

    rng = random.Random(args.seed)
    histogram = {}
    best_by_size = {}
    five_point_best = (-np.inf, None)
    five_subset_best = (-np.inf, None)

    for round_index in range(args.rounds):
        sample = rng.sample(range(len(signatures)), args.core)
        distances = grid_distance_matrix(rows[np.array(sample)])
        if defect(distances) <= 1e-9:
            print(f"round {round_index}: core of {args.core} is of negative type", flush=True)
            continue
        shrunk = eigen_shrink(distances, list(range(args.core)), floor=5)
        reduced = greedy_reduce(distances, shrunk, floor=5)
        family = [signatures[sample[i]] for i in reduced]
        size = len(family)
        value = defect(subset(distances, reduced))
        histogram[size] = histogram.get(size, 0) + 1
        if size not in best_by_size or value > best_by_size[size][0]:
            best_by_size[size] = (value, family)
        print(f"round {round_index}: minimal size {size}, defect {value:.4e}", flush=True)

        if size == 6:
            for drop in range(6):
                sub = [reduced[i] for i in range(6) if i != drop]
                plain = defect(subset(distances, sub))
                if plain > five_subset_best[0]:
                    five_subset_best = (plain, [signatures[sample[i]] for i in sub])
                seed_index = [
                    refine_index[signatures[sample[i]]]
                    for i in sub
                    if signatures[sample[i]] in refine_index
                ]
                if len(seed_index) < 5:
                    continue
                refined_value, refined = seeded_ascent(refine_distances, seed_index)
                if refined_value > five_point_best[0]:
                    five_point_best = (refined_value, [refine_core[i] for i in refined])

    print("\nhistogram of inclusion-minimal sizes:", dict(sorted(histogram.items())))
    results = {"histogram": {str(k): v for k, v in sorted(histogram.items())}}

    for size in sorted(best_by_size):
        value, family = best_by_size[size]
        exact = exact_distance_matrix(family)
        exact_value, vector = negative_type_defect(exact)
        print(f"\nbest {size}-point family (grid defect {value:.6e})")
        print("  family:", family)
        print(f"  exact defect: {exact_value:+.6e}")
        print("  x =", np.round(vector, 9).tolist())
        results[f"best_{size}"] = {
            "family": [list(s) for s in family],
            "exact_defect": exact_value,
            "x": vector.tolist(),
        }

    print("\nbest negative-type defect over 5-point families")
    print(f"  from 5-subsets of minimal 6-families: {five_subset_best[0]:+.6e}")
    print("   ", five_subset_best[1])
    print(f"  after seeded alternating ascent:      {five_point_best[0]:+.6e}")
    print("   ", five_point_best[1])
    if five_point_best[1]:
        exact = exact_distance_matrix(five_point_best[1])
        print(f"  exact defect of that 5-family: {negative_type_defect(exact)[0]:+.6e}")
        results["best_five"] = {
            "family": [list(s) for s in five_point_best[1]],
            "exact_defect": negative_type_defect(exact)[0],
        }
    OUTPUT.write_text(json.dumps(results, indent=2))
    print(f"\nwritten to {OUTPUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
