#!/usr/bin/env python3
"""Where does the exchange metric sit in the Deza-Laurent hierarchy?

    CUT (l1) subset HYPERMETRIC subset NEGATIVE TYPE subset METRIC

The hypermetric inequalities are ``sum_{i<j} b_i b_j d_ij <= 0`` for integer
vectors ``b`` with ``sum(b) = 1``; ``b = (1,1,-1)`` is the triangle inequality
and ``b = (1,1,1,-1,-1)`` the pentagonal inequality.  This script tests a list
of small ``b`` patterns against the exchange metric and reports, for each, the
largest value it can find together with the witnessing family.

The search is a coordinate ascent that is exact in each coordinate: with ``b``
fixed, the objective is *linear* in the row of ``D`` belonging to one point, so
the best replacement for that point is a single argmax over the whole core.

Usage::

    python hierarchy.py --core 2000 --restarts 400
"""

from __future__ import annotations

import argparse
import itertools
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from campaign import load_core
from exchange_geometry import exact_distance_matrix, negative_type_defect

# (name, b) with sum(b) = 1.  Ordering inside b is irrelevant: the search
# assigns signatures to positions.
PATTERNS = [
    ("triangle            b=(1,1,-1)", (1, 1, -1)),
    ("4-point             b=(1,1,1,-2)", (1, 1, 1, -2)),
    ("4-point             b=(2,1,-1,-1)", (2, 1, -1, -1)),
    ("pentagonal          b=(1,1,1,-1,-1)", (1, 1, 1, -1, -1)),
    ("5-point             b=(2,1,-1,-1,0)", (2, 1, -1, -1)),
    ("5-point             b=(2,1,1,-1,-2)", (2, 1, 1, -1, -2)),
    ("5-point             b=(1,1,1,1,-3)", (1, 1, 1, 1, -3)),
    ("6-point             b=(2,1,1,-1,-1,-1)", (2, 1, 1, -1, -1, -1)),
    ("6-point             b=(1,1,1,1,-1,-2)", (1, 1, 1, 1, -1, -2)),
    ("7-point             b=(1,1,1,1,-1,-1,-1)", (1, 1, 1, 1, -1, -1, -1)),
    ("9-point             b=(1,1,1,1,1,-1,-1,-1,-1)", (1, 1, 1, 1, 1, -1, -1, -1, -1)),
]


def coordinate_ascent(distances, b, rng, restarts, sweeps=12):
    """Maximize ``0.5 b^T D b`` over injective assignments of ``b`` to the core."""

    n = distances.shape[0]
    weights = np.asarray(b, dtype=float)
    k = len(weights)
    best = (-np.inf, None)
    for _ in range(restarts):
        index = rng.sample(range(n), k)
        for _ in range(sweeps):
            improved = False
            for position in range(k):
                others = [i for i in range(k) if i != position]
                # objective as a function of the point at `position`
                scores = weights[position] * (
                    distances[:, [index[i] for i in others]] @ weights[others]
                )
                scores[index] = -np.inf  # keep the assignment injective
                candidate = int(scores.argmax())
                if candidate != index[position]:
                    index[position] = candidate
                    improved = True
            if not improved:
                break
        block = distances[np.ix_(index, index)]
        value = 0.5 * float(weights @ block @ weights)
        if value > best[0]:
            best = (value, list(index))
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--restarts", type=int, default=300)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    args = parser.parse_args()

    core, distances = load_core(args.seed, args.core, args.max_length, args.max_entry)
    scale = distances.sum() / (len(core) * (len(core) - 1))
    print(f"core {len(core)} signatures, mean distance {scale:.4f}\n")
    rng = random.Random(args.seed + 17)

    for name, b in PATTERNS:
        value, index = coordinate_ascent(distances, b, rng, args.restarts)
        family = [core[i] for i in index]
        exact = exact_distance_matrix(family)
        weights = np.asarray(b, dtype=float)
        exact_value = 0.5 * float(weights @ exact @ weights)
        status = "VIOLATED" if exact_value > 1e-9 else "holds   "
        print(f"{name:<44} {status}  grid {value:+.6e}  exact {exact_value:+.6e}")
        if exact_value > 1e-9:
            print("    family:", family)
            print(f"    negative-type defect of the same family: {negative_type_defect(exact)[0]:+.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
