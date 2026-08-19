#!/usr/bin/env python3
"""Exact ``N_1(g -> f)`` for all 14 x 14 ordered pairs of F_3 map classes."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from affine_atoms import (  # noqa: E402
    atom_tables, enumerate_classes, minimal_atom_count_k1,
)

OUTPUT = Path(__file__).resolve().parent / "n1_matrix.json"


def main() -> None:
    classes = enumerate_classes()
    keys = [c.key for c in classes]
    atoms = {c.key: atom_tables(c.table, 1) for c in classes}

    start = time.time()
    table: dict[str, dict[str, int | None]] = {}
    for g in classes:
        row: dict[str, int | None] = {}
        for f in classes:
            row[f.key] = minimal_atom_count_k1(atoms[g.key],
                                               np.array(f.table, dtype=np.int8))
        table[g.key] = row
        print(f"  {g.key:<13} done  {time.time() - start:6.1f}s", flush=True)

    OUTPUT.write_text(json.dumps({
        "classes": [{"key": c.key, "representative": c.representative,
                     "signature": list(c.signature), "degree": c.degree,
                     "size": c.size, "atoms_k1": int(len(atoms[c.key]))}
                    for c in classes],
        "n1": table,
    }, indent=1))

    width = max(len(k) for k in keys) + 2
    print()
    print("N_1(g -> f):  rows g (resource), columns f (target)")
    print(" " * width + "".join(f"{k[:12]:>13}" for k in keys))
    for g in keys:
        cells = []
        for f in keys:
            value = table[g][f]
            cells.append("inf" if value is None else str(value))
        print(f"{g:<{width}}" + "".join(f"{c:>13}" for c in cells))


if __name__ == "__main__":
    main()
