#!/usr/bin/env python3
"""Collect distinct five-signature families that violate negative type.

Runs the full-pool alternating ascent from many random starts and records every
run that lands on a violating family, verifying each with the exact solver.
Writes ``five_point_violations.csv``.

Usage::  python five_point_census.py --restarts 200 --seed 0
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from exchange_geometry import exact_distance_matrix, negative_type_defect
from full_ascent import ColumnCache, ascend
from negative_type_search import build_pool

OUTPUT = Path(__file__).resolve().parent / "five_point_violations.csv"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--restarts", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-entry", type=int, default=12)
    parser.add_argument("--max-length", type=int, default=6)
    parser.add_argument(
        "--seeds-from",
        default="",
        help="json list of families; every k-subset of each becomes an ascent start",
    )
    args = parser.parse_args()

    signatures, rows = build_pool(args.max_length, args.max_entry)
    cache = ColumnCache(rows)
    rng = random.Random(args.seed + 977)
    lookup = {signature: i for i, signature in enumerate(signatures)}
    starts = []
    if args.seeds_from:
        import itertools
        import json

        for family in json.loads(Path(args.seeds_from).read_text()):
            index = [lookup[tuple(s)] for s in family]
            starts.extend(list(combo) for combo in itertools.combinations(index, args.k))
    starts.extend(rng.sample(range(len(signatures)), args.k) for _ in range(args.restarts))

    found = {}
    for number, start in enumerate(starts):
        value, index = ascend(cache, start)
        if value > 1e-9:
            family = tuple(sorted(signatures[i] for i in index))
            if family not in found:
                exact = exact_distance_matrix(list(family))
                defect, witness = negative_type_defect(exact)
                found[family] = (defect, witness)
                print(f"  run {number:>4}: NEW violating family, exact defect {defect:+.6e}")
                for signature in family:
                    print("      ", signature)
        if (number + 1) % 25 == 0:
            print(f"  {number + 1}/{len(starts)} runs, {len(found)} distinct families", flush=True)

    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["family", "exact_defect", "witness"])
        for family, (defect, witness) in sorted(found.items(), key=lambda item: -item[1][0]):
            writer.writerow(
                [
                    " | ".join("{" + ",".join(map(str, s)) + "}" for s in family),
                    f"{defect:.15e}",
                    " ".join(f"{value:.12f}" for value in witness),
                ]
            )
    print(f"\n{len(found)} distinct violating {args.k}-point families written to {OUTPUT.name}")
    if found:
        best = max(found.items(), key=lambda item: item[1][0])
        print(f"best exact defect {best[1][0]:+.9e}")
        for signature in best[0]:
            print("   ", signature)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
