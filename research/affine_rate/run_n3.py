#!/usr/bin/env python3
"""Decide ``N_3(g -> f) <= 4`` exhaustively for the F_3 quadratic pairs.

``N_3 = 3`` would force ``C_aff = 1``; ``N_3 = 4`` gives the Fekete lower bound
``C_aff >= 3/4``.  Both are decided here, because a solution with four atoms has
span ``U`` of dimension at most 4 containing the three-dimensional hull
``H = span{J(f(x_1)), J(f(x_2)), J(f(x_3))}``, so ``U = H`` or ``U = H + <a>``.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from isotropy_atoms import (  # noqa: E402
    BINARY, JetSpace, atom_codes, n_k_le_hull_plus_one, target_jets,
)

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "n3_le4.json"
KEYS = ("linear", "rank1", "split", "sum2", "parabolic")


def main() -> None:
    q, k = 3, 3
    space = JetSpace.make(q, 2 * k)
    print(f"jet space: q={q}, n={2*k}, dim={space.dim}")
    atoms = {}
    for g in KEYS:
        start = time.time()
        atoms[g] = atom_codes(space, BINARY[g])
        print(f"  atoms(k={k}) {g:<10} {len(atoms[g]):>9}   {time.time()-start:6.1f}s",
              flush=True)

    result: dict[str, dict[str, str]] = {}
    for g in KEYS:
        row: dict[str, str] = {}
        for f in KEYS:
            start = time.time()
            value = n_k_le_hull_plus_one(space, atoms[g], target_jets(space, BINARY[f], k))
            row[f] = str(value) if value is not None else ">4"
            print(f"  N_3({g:<10}-> {f:<10}) = {row[f]:>3}   {time.time()-start:7.1f}s",
                  flush=True)
        result[g] = row

    OUTPUT.write_text(json.dumps(result, indent=1))
    print()
    print("N_3(g -> f) decided up to 4:")
    print(f"{'':<12}" + "".join(f"{k2:>12}" for k2 in KEYS))
    for g in KEYS:
        print(f"{g:<12}" + "".join(f"{result[g][f]:>12}" for f in KEYS))


if __name__ == "__main__":
    main()
