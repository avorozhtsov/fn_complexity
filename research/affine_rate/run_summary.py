#!/usr/bin/env python3
"""Assemble the C_aff bracket table and compare it with C_sig.

Lower bounds come from the computed ``N_k`` (Fekete: ``C_aff >= k / N_k``).
Upper bounds are the proved ones recorded in ``FINDINGS.md``:

* ``C_aff <= 1``                        for every pair (output dimensions equal);
* ``C_aff = 0``                         when ``N_1`` is infinite;
* ``C_aff <= 1/2``   for ``g = x^2``    and ``f`` in {x, xy, x^2+y^2, x^2+y};
* ``C_aff <= 2/3``   for ``g = x^2+y^2`` and ``f = x``.
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

# N_2 values: exhaustive.  "4" entries are  N_1 = 2  (so N_2 <= 4 by
# block-diagonal composition) together with the exhaustive refutation of
# N_2 <= 3 in run_n2_small.py.  "6" is N_1 = 3 with N_2 > 5 from
# run_rank1_parabolic.py.
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

PROVED_UPPER = {("rank1", f): Fraction(1, 2)
                for f in ("linear", "split", "anisotropic", "parabolic")}
PROVED_UPPER[("anisotropic", "linear")] = Fraction(2, 3)


def main() -> None:
    n1 = json.loads((HERE / "n1_matrix.json").read_text())["n1"]
    csig = json.loads((HERE / "csig_matrix.json").read_text())["csig"]

    rows = []
    for g in QUADRATIC:
        for f in QUADRATIC:
            value = n1[g][f]
            if value is None:
                low, high = Fraction(0), Fraction(0)
            elif value == 1:
                low, high = Fraction(1), Fraction(1)
            else:
                low = max(Fraction(1, value), Fraction(2, N2[(g, f)]))
                high = PROVED_UPPER.get((g, f), Fraction(1))
            rows.append((g, f, value, N2.get((g, f)), low, high,
                         float(csig[g][f])))

    print(f"{'g':<12}{'f':<12}{'N_1':>5}{'N_2':>5}"
          f"{'C_aff lower':>14}{'C_aff upper':>14}{'C_sig':>12}   verdict")
    strictly_below = strictly_above = undetermined = 0
    for g, f, one, two, low, high, sig in rows:
        if high < sig:
            verdict = "C_aff < C_sig"
            strictly_below += 1
        elif low > sig:
            verdict = "C_aff > C_sig"
            strictly_above += 1
        elif low == high == 1 and abs(sig - 1) < 1e-30:
            verdict = "equal (both 1)"
        else:
            verdict = "undetermined"
            undetermined += 1
        print(f"{g:<12}{f:<12}{str(one):>5}{str(two):>5}"
              f"{str(low):>14}{str(high):>14}{sig:>12.6f}   {verdict}")
    print()
    print(f"strictly below: {strictly_below}   strictly above: {strictly_above}   "
          f"undetermined: {undetermined}")


if __name__ == "__main__":
    main()
