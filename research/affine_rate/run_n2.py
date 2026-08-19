#!/usr/bin/env python3
"""Exact ``N_2(g -> f)`` for the quadratic F_3 classes, up to a search limit.

``N_2`` is the smallest number of g-atoms on ``F_3^4`` whose span, together
with the constants, contains both ``f(x_1)`` and ``f(x_2)``.  Since
``k_{g->f}(r) = max{k : N_k <= r}``, the pair ``(N_1, N_2)`` already gives
``k(1), ..., k(N_2)`` and the Fekete lower bound ``C_aff >= 2 / N_2``.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from affine_atoms import (  # noqa: E402
    atom_tables, enumerate_classes, minimal_atom_count, target_tables,
)

OUTPUT = Path(__file__).resolve().parent / "n2_matrix.json"
QUADRATIC = ("linear", "rank1", "split", "anisotropic", "parabolic")
# Plain subset search is only feasible for the small atom sets; for the others
# use run_n2_small.py (the exhaustive N_2 <= 3 decision) together with N_1.
FEASIBLE = ("linear", "rank1")


def main(limit: int = 4) -> None:
    classes = {c.key: c for c in enumerate_classes()}
    atoms = {}
    for key in FEASIBLE:
        start = time.time()
        atoms[key] = atom_tables(classes[key].table, 2)
        print(f"  atoms(k=2) {key:<12} {len(atoms[key]):>6}"
              f"   {time.time() - start:5.1f}s", flush=True)

    result: dict[str, dict[str, str]] = {}
    for g in FEASIBLE:
        row: dict[str, str] = {}
        for f in QUADRATIC:
            start = time.time()
            targets = target_tables(classes[f].table, 2)
            value = minimal_atom_count(atoms[g], targets, limit)
            row[f] = str(value) if value is not None else f">{limit}"
            print(f"  N_2({g:<12}-> {f:<12}) = {row[f]:>4}"
                  f"   {time.time() - start:6.1f}s", flush=True)
        result[g] = row

    OUTPUT.write_text(json.dumps({"limit": limit, "n2": result}, indent=1))
    print()
    print(f"N_2(g -> f), searched to {limit}:")
    print(f"{'':<14}" + "".join(f"{k:>14}" for k in QUADRATIC))
    for g in FEASIBLE:
        print(f"{g:<14}" + "".join(f"{result[g][f]:>14}" for f in QUADRATIC))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
