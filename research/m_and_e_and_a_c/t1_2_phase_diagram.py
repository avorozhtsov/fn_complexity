#!/usr/bin/env python3
"""The PSD threshold as a critical inverse temperature.

T1.2 established that ``S(d) = {t : exp(-t d) >= 0}`` was a single closed ray
``[t*, inf)`` in every one of 6886 cases, with exactly one sign change of
``lambda_min``.  Read as quantum mechanics that is a phase diagram, not a
positivity scan.  The matrix

    rho_t = exp(-t d) / tr exp(-t d)

is a Gibbs state whose Hamiltonian is the exchange metric and whose inverse
temperature is ``t``.  By Schoenberg it is a legitimate state for every ``t > 0``
iff ``d`` is of negative type; Theorem 2 says it is not.  So there is a critical
``t*`` below which ``rho_t`` has a negative eigenvalue -- a negative probability
-- and the family admits no Hilbert-space description at all.  ``lambda_min(t)``
is the order parameter and the single sign change is the transition.

Two points on the ``t`` axis are distinguished:

  * ``t = 1/2``, where ``exp(-d/2) = sqrt(M o M^T)`` is the discriminant matrix
    of Szegedy's quantum walk (see `szegedy_walk.py`);
  * ``t -> inf``, the tropical limit, where ``exp(-t d) -> I`` and positivity is
    trivial.  Positivity is EASY in the classical limit and hard under strong
    quantisation, which is the opposite of the naive expectation and is worth
    stating explicitly.

Writes `t1_2_phase_diagram.csv` and, if matplotlib is available,
`figures/psd_phase_diagram.png`.

    python research/m_and_e_and_a_c/t1_2_phase_diagram.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))

import t1_2_common as T  # noqa: E402
from gauge_decomposition import CERT5  # noqa: E402

CSV_PATH = HERE / "t1_2_phase_diagram.csv"
FIG_PATH = HERE.parents[1] / "figures" / "psd_phase_diagram.png"


def curve(distances: np.ndarray, ts: np.ndarray) -> np.ndarray:
    return np.array([float(np.linalg.eigvalsh(np.exp(-t * distances))[0]) for t in ts])


def main() -> None:
    families = T.build_families()
    families["cert5"] = CERT5
    chosen = [n for n in ("cert5", "cert13", "greedy25", "greedy40",
                          "rand25_55", "rand40_71") if n in families]

    ts = np.geomspace(1e-3, 1e3, 400)
    rows, curves = [], {}
    for name in chosen:
        d = T.distance_matrix(families[name])
        lam = curve(d, ts)
        star = T.psd_threshold(d)
        curves[name] = (lam, star)
        signs = int(np.sum(np.diff(np.sign(lam)) != 0))
        half = float(np.linalg.eigvalsh(np.exp(-0.5 * d))[0])
        print(f"{name:>10} n={len(families[name]):>3}  "
              f"t* = {'-' if star is None else f'{star:.6f}':>12}   "
              f"sign changes on the grid = {signs}   "
              f"lambda_min(t=1/2) = {half: .4e}")
        for t, v in zip(ts, lam):
            rows.append({"family": name, "size": len(families[name]),
                         "t": f"{t:.6e}", "lambda_min": f"{v:.9e}",
                         "t_star": "" if star is None else f"{star:.9f}"})

    with CSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["family", "size", "t",
                                                    "lambda_min", "t_star"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH.relative_to(HERE.parents[1])} ({len(rows)} rows)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping the figure")
        return

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for name, (lam, star) in curves.items():
        scale = np.abs(lam).max()
        ax.plot(ts, lam / scale, label=f"{name}"
                + ("" if star is None else f"  ($t^*={star:.3g}$)"))
        if star is not None:
            ax.axvline(star, color="0.75", lw=0.6, zorder=0)
    ax.axhline(0.0, color="0.2", lw=0.8)
    ax.axvline(0.5, color="C3", lw=1.0, ls=":", label="Szegedy point $t=1/2$")
    ax.set_xscale("log")
    ax.set_xlabel(r"inverse temperature $t$")
    ax.set_ylabel(r"$\lambda_{\min}(e^{-td})$, normalised")
    ax.set_title("PSD threshold as a phase transition: below $t^*$ the Gibbs\n"
                 "state of the exchange metric has a negative eigenvalue")
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    FIG_PATH.parent.mkdir(exist_ok=True)
    fig.savefig(FIG_PATH, dpi=160)
    print(f"wrote {FIG_PATH.relative_to(HERE.parents[1])}")


if __name__ == "__main__":
    main()
