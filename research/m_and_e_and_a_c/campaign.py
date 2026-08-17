#!/usr/bin/env python3
"""Search campaign: how small can a negative-type / hypermetric violation be?

Works on cached cores: a core of a few thousand signatures is drawn from the
pool, its grid distance matrix is computed once, and all subset searches then
run on that matrix (microseconds per evaluation).

Usage::

    python campaign.py core   --core 2000 --seed 0        # build + cache a core
    python campaign.py k5     --seed 0                    # hunt 5-point violations
    python campaign.py k6     --seed 0                    # optimize 6-point violations
    python campaign.py hyper  --seed 0 --k 5|7|9          # hypermetric violations
    python campaign.py sizes  --seed 0                    # minimal-size statistics
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

from exchange_geometry import exact_distance_matrix, grid_distance_matrix, negative_type_defect
from negative_type_search import (
    batched_defect,
    batched_hypermetric,
    build_pool,
    defect,
    eigen_shrink,
    greedy_reduce,
    hypermetric_patterns,
    subset,
)

CACHE = Path(__file__).resolve().parent / "cache"
CACHE.mkdir(exist_ok=True)


def relative(score):
    """Make a quadratic score scale-free by dividing by the mean distance.

    Both the negative-type defect and the hypermetric forms are homogeneous of
    degree one in ``D``, so an unnormalized search drifts towards families of
    tiny diameter (where every value is near zero) instead of towards families
    whose violation is structurally large.
    """

    def wrapped(block):
        size = block.shape[0]
        scale = block.sum() / max(1.0, size * (size - 1))
        if scale <= 0:
            return -math.inf
        return score(block) / scale

    return wrapped


class LazyDistances:
    """Pairwise grid distances over the whole pool, cached on demand."""

    def __init__(self, rows):
        self.rows = rows
        self.shape = (rows.shape[0], rows.shape[0])
        self._cache = {}

    def pair(self, i, j):
        if i == j:
            return 0.0
        key = (i, j) if i < j else (j, i)
        value = self._cache.get(key)
        if value is None:
            difference = self.rows[key[0]] - self.rows[key[1]]
            value = float(difference.max() - difference.min())
            self._cache[key] = value
        return value

    def block(self, indices):
        size = len(indices)
        out = np.zeros((size, size))
        for a in range(size):
            for b in range(a + 1, size):
                out[a, b] = out[b, a] = self.pair(indices[a], indices[b])
        return out


def core_path(seed, size, max_length, max_entry):
    return CACHE / f"core_s{seed}_n{size}_L{max_length}_E{max_entry}.npz"


def load_core(seed, size, max_length, max_entry, structured=False):
    path = core_path(seed, size, max_length, max_entry)
    if path.exists():
        data = np.load(path, allow_pickle=True)
        return [tuple(s) for s in data["signatures"]], data["distances"]
    signatures, rows = build_pool(max_length, max_entry)
    rng = random.Random(seed)
    if structured:
        chosen = structured_core(signatures, size, rng)
    else:
        chosen = rng.sample(range(len(signatures)), min(size, len(signatures)))
    core = [signatures[i] for i in chosen]
    distances = grid_distance_matrix(rows[np.array(chosen)])
    np.savez_compressed(
        path, signatures=np.array(core, dtype=object), distances=distances
    )
    return core, distances


def structured_core(signatures, size, rng):
    """Prefer near-flat signatures and ones with a single dominant entry."""

    scored = []
    for index, signature in enumerate(signatures):
        spread = signature[0] / signature[-1]
        flat = abs(spread - 1.0) < 0.35
        spike = signature[0] >= 3 * max(signature[1], 1)
        if flat or spike:
            scored.append(index)
    if len(scored) >= size:
        return rng.sample(scored, size)
    rest = [i for i in range(len(signatures)) if i not in set(scored)]
    return scored + rng.sample(rest, size - len(scored))


# ------------------------------------------------------------------ searches


def block_of(distances, indices):
    if isinstance(distances, LazyDistances):
        return distances.block(list(indices))
    return subset(distances, indices)


def anneal_subsets(distances, k, steps, restarts, seed, score, verbose=True):
    rng = random.Random(seed)
    n = distances.shape[0]
    overall = (-math.inf, None)
    for restart in range(restarts):
        current = rng.sample(range(n), k)
        value = score(block_of(distances, current))
        best = (value, list(current))
        scale = max(abs(value), 1e-4)
        for step in range(steps):
            temperature = scale * (1 - step / steps) ** 3 + 1e-15
            trial = list(current)
            trial[rng.randrange(k)] = rng.randrange(n)
            if len(set(trial)) < k:
                continue
            candidate = score(block_of(distances, trial))
            if candidate > value or rng.random() < math.exp(
                min(0.0, (candidate - value) / temperature)
            ):
                current, value = trial, candidate
                if value > best[0]:
                    best = (value, list(current))
        if verbose:
            print(f"    restart {restart}: {best[0]:+.6e}", flush=True)
        if best[0] > overall[0]:
            overall = best
    return overall


def exhaustive(distances, indices, k, score_batch, chunk=200_000):
    """Sweep every k-subset of ``indices``; returns (best value, best subset)."""

    idx_all = list(indices)
    n = len(idx_all)
    block_matrix = subset(distances, idx_all)
    best = (-math.inf, None)
    combos = itertools.combinations(range(n), k)
    total = math.comb(n, k)
    processed = 0
    while True:
        batch = list(itertools.islice(combos, chunk))
        if not batch:
            break
        idx = np.array(batch)
        blocks = block_matrix[idx[:, :, None], idx[:, None, :]]
        values = score_batch(blocks)
        pick = int(values.argmax())
        if values[pick] > best[0]:
            best = (float(values[pick]), [idx_all[i] for i in idx[pick]])
        processed += len(batch)
        if processed % (chunk * 25) == 0:
            print(f"      {processed}/{total}", flush=True)
    return best


def local_pool_refine(distances, family, score, rounds=8):
    """Coordinate ascent: replace one member at a time by the best pool element."""

    current = list(family)
    value = score(block_of(distances, current))
    n = distances.shape[0]
    for _ in range(rounds):
        improved = False
        for position in range(len(current)):
            best_swap = (value, current[position])
            for candidate in range(n):
                if candidate in current:
                    continue
                trial = list(current)
                trial[position] = candidate
                v = score(block_of(distances, trial))
                if v > best_swap[0]:
                    best_swap = (v, candidate)
            if best_swap[1] != current[position]:
                current[position] = best_swap[1]
                value = best_swap[0]
                improved = True
        if not improved:
            break
    return value, current


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode")
    parser.add_argument("--core", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--steps", type=int, default=60_000)
    parser.add_argument("--restarts", type=int, default=10)
    parser.add_argument("--rounds", type=int, default=30)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    parser.add_argument("--structured", action="store_true")
    parser.add_argument("--sweep-core", type=int, default=90)
    parser.add_argument("--full", action="store_true", help="search the whole pool lazily")
    parser.add_argument("--absolute", action="store_true", help="do not normalize by scale")
    parser.add_argument("--start", default="", help="python literal list of starting signatures")
    args = parser.parse_args()

    if args.full:
        core, rows = build_pool(args.max_length, args.max_entry)
        distances = LazyDistances(rows)
        print(f"pool: {len(core)} signatures (lazy distances)")
    else:
        core, distances = load_core(
            args.seed, args.core, args.max_length, args.max_entry, args.structured
        )
        print(f"core: {len(core)} signatures, seed {args.seed}, structured={args.structured}")

    wrap = (lambda f: f) if args.absolute else relative

    if args.mode == "core":
        triple = math.inf
        for start in range(0, len(core), 50):  # chunked: the full triple tensor is huge
            stop = min(len(core), start + 50)
            block = (
                distances[start:stop, None, :]
                + distances[None, :, :]
                - distances[start:stop, :, None]
            )
            triple = min(triple, float(block.min()))
        print(f"  min over triples of d_ik + d_kj - d_ij = {triple:.3e} (>= 0 is the triangle law)")
        print(f"  full-core negative-type defect = {defect(distances):+.6e}")
        return 0

    if args.mode.startswith("k") and args.mode[1:].isdigit():
        k = int(args.mode[1:])
        score = wrap(defect)
        print(f"  annealing negative-type defect over {k}-subsets")
        if args.start:
            start = [core.index(tuple(s)) for s in eval(args.start)]  # noqa: S307
            value, family = local_pool_refine(distances, start, score)
            print(f"  refined from the given start: {value:+.6e}")
        else:
            value, family = anneal_subsets(
                distances, k, args.steps, args.restarts, args.seed, score
            )
            value, family = local_pool_refine(distances, family, score)
        signatures = [core[i] for i in family]
        raw = defect(block_of(distances, family))
        print(f"  best grid defect for k={k}: {raw:+.6e}   (relative {value:+.6e})")
        print("  family:", signatures)
        exact = exact_distance_matrix(signatures)
        exact_value, vector = negative_type_defect(exact)
        print(f"  exact defect: {exact_value:+.6e}")
        if exact_value > 0:
            print("  x =", np.round(vector, 9).tolist())
        return 0

    if args.mode == "sweep":
        # exhaustive k-subset sweep over a sub-core chosen by pool degree
        rng = random.Random(args.seed + 7)
        sub = rng.sample(range(len(core)), args.sweep_core)
        print(f"  exhaustive sweep k={args.k} over {args.sweep_core} signatures")
        best = exhaustive(distances, sub, args.k, batched_defect)
        print(f"  best defect: {best[0]:+.6e}")
        print("  family:", [core[i] for i in best[1]])
        return 0

    if args.mode == "hyper":
        patterns = hypermetric_patterns(args.k, (args.k + 1) // 2)
        raw_score = lambda D: float(  # noqa: E731
            (0.5 * np.einsum("pi,ij,pj->p", patterns, D, patterns)).max()
        )
        score = wrap(raw_score)
        print(f"  annealing hypermetric value over {args.k}-subsets ({len(patterns)} patterns)")
        if args.start:
            start = [core.index(tuple(s)) for s in eval(args.start)]  # noqa: S307
            value, family = local_pool_refine(distances, start, score)
        else:
            value, family = anneal_subsets(
                distances, args.k, args.steps, args.restarts, args.seed, score
            )
            value, family = local_pool_refine(distances, family, score)
        print(f"  best grid hypermetric value (relative): {value:+.6e}")
        signatures = [core[i] for i in family]
        print("  family:", signatures)
        exact = exact_distance_matrix(signatures)
        values = 0.5 * np.einsum("pi,ij,pj->p", patterns, exact, patterns)
        pick = int(values.argmax())
        print(f"  exact hypermetric value: {values[pick]:+.6e}")
        print("  b =", patterns[pick].astype(int).tolist())
        print(f"  exact negative-type defect of the same family: {negative_type_defect(exact)[0]:+.6e}")
        return 0

    if args.mode == "sizes":
        rng = random.Random(args.seed + 11)
        histogram = {}
        best = None
        for round_index in range(args.rounds):
            sample = rng.sample(range(len(core)), 300)
            block = subset(distances, sample)
            if defect(block) <= 0:
                continue
            shrunk = eigen_shrink(block, list(range(len(sample))), floor=5)
            reduced = greedy_reduce(block, shrunk, floor=5)
            size = len(reduced)
            histogram[size] = histogram.get(size, 0) + 1
            value = defect(subset(block, reduced))
            family = [core[sample[i]] for i in reduced]
            if best is None or (size, -value) < (len(best[1]), -best[0]):
                best = (value, family)
            print(f"    round {round_index}: {size} points, defect {value:.3e}", flush=True)
        print("  histogram of inclusion-minimal sizes:", dict(sorted(histogram.items())))
        print(f"  best: {len(best[1])} points, defect {best[0]:.6e}")
        print("  family:", best[1])
        return 0

    raise SystemExit(f"unknown mode {args.mode}")


if __name__ == "__main__":
    raise SystemExit(main())
