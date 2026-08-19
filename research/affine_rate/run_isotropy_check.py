#!/usr/bin/env python3
"""Validate the jet representation against the committed F_3 value-table results.

``affine_atoms.py`` computes ``N_1`` and ``N_2`` on ``F_3`` by exhaustive search
over value tables of functions ``F_3^{2k} -> F_3``.  ``isotropy_atoms.py``
computes the same quantity from jets (quadratic part, linear part).  The two
must agree on the whole quadratic block; this script checks that.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from isotropy_atoms import (  # noqa: E402
    BINARY, JetSpace, atom_codes, n_k_le_hull_plus_one, n_k_small, target_jets,
)

HERE = Path(__file__).resolve().parent
KEYS = ("linear", "rank1", "split", "sum2", "parabolic")
OLD_NAME = {"linear": "linear", "rank1": "rank1", "split": "split",
            "sum2": "anisotropic", "parabolic": "parabolic"}


def main() -> None:
    q = 3
    n1_old = json.loads((HERE / "n1_matrix.json").read_text())["n1"]
    n2_old = json.loads((HERE / "n2_le3.json").read_text())

    print(f"q = {q}: N_1 from jets versus the committed value-table N_1")
    space1 = JetSpace.make(q, 2)
    atoms1 = {g: atom_codes(space1, BINARY[g]) for g in KEYS}
    ok = True
    for g in KEYS:
        for f in KEYS:
            got = n_k_small(space1, atoms1[g], target_jets(space1, BINARY[f], 1), 6)
            want = n1_old[OLD_NAME[g]][OLD_NAME[f]]
            flag = "ok" if got == want else "MISMATCH"
            ok &= got == want
            print(f"  N_1({g:<10}-> {f:<10}) = {str(got):>4}   committed {str(want):>4}  {flag}")

    print()
    print(f"q = {q}: N_2 from jets versus the committed exhaustive 'N_2 <= 3' decision")
    space2 = JetSpace.make(q, 4)
    atoms2 = {g: atom_codes(space2, BINARY[g]) for g in KEYS}
    for g in KEYS:
        for f in KEYS:
            got = n_k_le_hull_plus_one(space2, atoms2[g],
                                       target_jets(space2, BINARY[f], 2))
            want = n2_old[OLD_NAME[g]][OLD_NAME[f]]
            got_text = str(got) if got is not None else ">3"
            flag = "ok" if got_text == want else "MISMATCH"
            ok &= got_text == want
            print(f"  N_2({g:<10}-> {f:<10}) = {got_text:>4}   committed {want:>4}  {flag}")

    print()
    print("all agree" if ok else "DISAGREEMENT -- do not trust the jet code")


if __name__ == "__main__":
    main()
