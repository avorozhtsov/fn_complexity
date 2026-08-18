#!/usr/bin/env python3
"""Is ``N_2(x^2 -> x^2+y) = 5`` or ``6``?

``N_1 = 3`` (computed) and ``N_k >= 2k+1`` (proved in FINDINGS.md), so
``N_2`` is 5 or 6; block-diagonal composition gives the upper bound 6.
The x^2-atoms on F_3^4 are just the 122 squares of affine forms, so the
exhaustive subset search at depth 5 is feasible.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from affine_atoms import (  # noqa: E402
    atom_tables, enumerate_classes, minimal_atom_count, target_tables,
)


def main() -> None:
    classes = {c.key: c for c in enumerate_classes()}
    atoms = atom_tables(classes["rank1"].table, 2)
    targets = target_tables(classes["parabolic"].table, 2)
    print(f"atoms(k=2, g = x^2): {len(atoms)}", flush=True)
    start = time.time()
    value = minimal_atom_count(atoms, targets, 5)
    print(f"N_2(x^2 -> x^2+y) {'=' if value else '>'} "
          f"{value if value else 5}   ({time.time() - start:.0f}s)")


if __name__ == "__main__":
    main()
