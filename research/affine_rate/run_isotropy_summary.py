#!/usr/bin/env python3
"""The `C_aff` bracket table for the F_3 quadratic pool, after ISOTROPY.md.

Same shape as ``run_summary.py``, with the upper bounds proved in
``ISOTROPY.md`` folded in:

* ``C_aff <= 2/3``  whenever the resource is anisotropic and the target admits an
  affine line on which it is a bijective affine function (Theorem 2); over F_3
  the targets with such a line are ``x``, ``xy`` and ``x^2+y``;
* ``C_aff = 1/rho`` for the parabolic resource ``x^2+y``, ``rho`` the rank of the
  target's quadratic part (Theorem 3);
* ``C_aff = 1/rho`` for the rank-one resource ``x^2`` against a purely quadratic
  target, and ``C_aff = 1/2`` for ``x^2 -> x`` (section 4).

Lower bounds are the Fekete ratios ``k/N_k`` from the committed exhaustive
``N_1`` and ``N_2``.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

HERE = Path(__file__).resolve().parent
QUADRATIC = ("linear", "rank1", "split", "anisotropic", "parabolic")
REPRESENTATIVE = {"linear": "x", "rank1": "x^2", "split": "xy",
                  "anisotropic": "x^2+y^2", "parabolic": "x^2+y"}

# exhaustive; see run_n2.py, run_n2_small.py, run_rank1_parabolic.py
N2 = {
    ("linear", "linear"): 2,
    ("rank1", "linear"): 4, ("rank1", "rank1"): 2, ("rank1", "split"): 4,
    ("rank1", "anisotropic"): 4, ("rank1", "parabolic"): 6,
    ("split", "linear"): 2, ("split", "rank1"): 2, ("split", "split"): 2,
    ("split", "anisotropic"): 4, ("split", "parabolic"): 4,
    ("anisotropic", "linear"): 3, ("anisotropic", "rank1"): 2,
    ("anisotropic", "split"): 3, ("anisotropic", "anisotropic"): 2,
    ("anisotropic", "parabolic"): 4,
    ("parabolic", "linear"): 2, ("parabolic", "rank1"): 2,
    ("parabolic", "split"): 4, ("parabolic", "anisotropic"): 4,
    ("parabolic", "parabolic"): 2,
}

# proved upper bounds, with the section of ISOTROPY.md that proves them
PROVED_UPPER = {
    ("rank1", "linear"): (Fraction(1, 2), "s4"),
    ("rank1", "split"): (Fraction(1, 2), "s4"),
    ("rank1", "anisotropic"): (Fraction(1, 2), "s4"),
    ("rank1", "parabolic"): (Fraction(1, 2), "s4"),
    ("anisotropic", "linear"): (Fraction(2, 3), "Thm 1"),
    ("anisotropic", "split"): (Fraction(2, 3), "Thm 2"),
    ("anisotropic", "parabolic"): (Fraction(2, 3), "Thm 2"),
    ("parabolic", "split"): (Fraction(1, 2), "Thm 3"),
    ("parabolic", "anisotropic"): (Fraction(1, 2), "Thm 3"),
}


def main() -> None:
    n1 = json.loads((HERE / "n1_matrix.json").read_text())["n1"]
    csig = json.loads((HERE / "csig_matrix.json").read_text())["csig"]

    print(f"{'g':<10}{'f':<10}{'N_1':>4}{'N_2':>4}"
          f"{'C_aff lower':>13}{'C_aff upper':>13}{'by':>8}{'C_sig':>12}   verdict")
    below = above = equal = undetermined = 0
    for g in QUADRATIC:
        for f in QUADRATIC:
            value = n1[g][f]
            source = ""
            if value is None:
                low = high = Fraction(0)
            elif value == 1:
                low = high = Fraction(1)
            else:
                low = max(Fraction(1, value), Fraction(2, N2[(g, f)]))
                high, source = PROVED_UPPER.get((g, f), (Fraction(1), ""))
            sig = float(csig[g][f])
            if high < sig:
                verdict, _ = "C_aff < C_sig", (below := below + 1)
            elif low > sig:
                verdict, _ = "C_aff > C_sig", (above := above + 1)
            elif low == high == 1 and abs(sig - 1) < 1e-30:
                verdict, _ = "equal (both 1)", (equal := equal + 1)
            else:
                verdict, _ = "undetermined", (undetermined := undetermined + 1)
            print(f"{REPRESENTATIVE[g]:<10}{REPRESENTATIVE[f]:<10}"
                  f"{str(value):>4}{str(N2.get((g, f))):>4}"
                  f"{str(low):>13}{str(high):>13}{source:>8}{sig:>12.7f}   {verdict}")
    print()
    print(f"C_aff < C_sig: {below}   C_aff > C_sig: {above}   "
          f"both 1: {equal}   undetermined: {undetermined}")


if __name__ == "__main__":
    main()
