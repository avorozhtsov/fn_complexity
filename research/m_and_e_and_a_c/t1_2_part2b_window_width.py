#!/usr/bin/env python3
"""T1.2 part 2b: how wide must the spectral window be to keep the violation?

``D_t -> osc_window(u_i - u_j)`` as ``t -> inf``.  If the window is bounded,
that limit is only an approximation of the exchange metric ``d``, and the
negative-type violation of ``d`` is small enough that a truncation error of
1e-3 already destroys it.  This script sweeps the upper cut ``beta_max`` and
records the truncation error, the negative-type defect and ``t*``.

Writes ``t1_2_part2b_window_width.csv``.
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from t1_2_common import (  # noqa: E402
    CERTIFICATE_FAMILY,
    build_families,
    distance_matrix,
    negative_type_defect,
    psd_threshold,
    u_values,
)

HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "t1_2_part2b_window_width.csv"

BETA_MAX = [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 20000, 100000]
BETA_MIN = 1e-4
POINTS = 200_001


def oscillation(family, betas: np.ndarray) -> np.ndarray:
    """``osc_[beta_min, beta_max](u_i - u_j)`` pair by pair.

    Only the support of the grid matters for an oscillation, so the nodes are
    log-spaced: that resolves both ends of the spectrum at moderate cost.
    """

    table = np.array([u_values(s, betas) for s in family])
    size = len(family)
    out = np.zeros((size, size))
    for i in range(size):
        for j in range(i + 1, size):
            difference = table[i] - table[j]
            out[i, j] = out[j, i] = float(difference.max() - difference.min())
    return out


def main() -> int:
    families = build_families()
    selected = {
        "cert13": CERTIFICATE_FAMILY,
        "rand30_61": families["rand30_61"],
        "greedy25": families["greedy25"],
    }
    rows = []
    for name, family in selected.items():
        reference = distance_matrix(family)
        print(
            f"\n{name}: d has defect {negative_type_defect(reference):.4e},"
            f" t* = {psd_threshold(reference)}"
        )
        print(f"   {'beta_max':>9} {'max|osc-d|':>12} {'defect':>12} {'t*(osc)':>14}")
        for cut in BETA_MAX:
            betas = np.exp(np.linspace(np.log(BETA_MIN), np.log(float(cut)), POINTS))
            osc = oscillation(family, betas)
            error = float(np.abs(osc - reference).max())
            defect = negative_type_defect(osc)
            star = psd_threshold(osc)
            rows.append(
                [name, cut, f"{error:.9e}", f"{defect:.9e}",
                 "-" if star is None else f"{star:.9f}"]
            )
            print(
                f"   {cut:>9} {error:>12.3e} {defect:>12.3e}"
                f" {'  negative type' if star is None else f'{star:>14.6f}'}"
            )
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["family", "beta_max", "max_abs_osc_minus_d",
                         "negative_type_defect", "t_star"])
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
