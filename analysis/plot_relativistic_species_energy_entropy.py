#!/usr/bin/env python3
"""Plot the exactly homothetic E-S curves of relativistic radiation sectors."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT_ROOT / "paper" / "figures" / "relativistic-radiation-energy-entropy.pdf"


def radiation_curve(temperatures: np.ndarray, degrees_of_freedom: float):
    """Return entropy and energy densities in natural units."""
    energy = (np.pi**2 / 30.0) * degrees_of_freedom * temperatures**4
    entropy = (2.0 * np.pi**2 / 45.0) * degrees_of_freedom * temperatures**3
    return entropy, energy


def main() -> None:
    temperatures = np.linspace(0.0, 1.5, 500)
    photon_dof = 2.0
    electron_positron_dof = (7.0 / 8.0) * 4.0
    scale = electron_positron_dof / photon_dof

    photon_entropy, photon_energy = radiation_curve(temperatures, photon_dof)
    pair_entropy, pair_energy = radiation_curve(
        temperatures, electron_positron_dof
    )

    figure, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    axis.plot(
        photon_energy,
        photon_entropy,
        color="#1769aa",
        linewidth=2.4,
        label=r"photons: $g_\gamma=2$",
    )
    axis.plot(
        pair_energy,
        pair_entropy,
        color="#d95f02",
        linewidth=3.0,
        label=r"$e^\pm$ (ultrarelativistic): $g_{e^\pm}=7/2$",
    )
    axis.plot(
        scale * photon_energy,
        scale * photon_entropy,
        color="#1769aa",
        linewidth=1.8,
        linestyle=(0, (2, 2)),
        label=r"homothety $(E,S)\mapsto\frac{7}{4}(E,S)$",
    )

    reference_temperature = 1.0
    photon_point = radiation_curve(
        np.array([reference_temperature]), photon_dof
    )
    pair_point = radiation_curve(
        np.array([reference_temperature]), electron_positron_dof
    )
    axis.plot(
        [0.0, pair_point[1][0]],
        [0.0, pair_point[0][0]],
        color="#555555",
        linewidth=0.9,
        linestyle="--",
        zorder=1,
    )
    axis.scatter(
        [photon_point[1][0], pair_point[1][0]],
        [photon_point[0][0], pair_point[0][0]],
        s=42,
        facecolors="white",
        edgecolors=["#1769aa", "#d95f02"],
        linewidths=1.8,
        zorder=5,
    )
    axis.annotate(
        r"same $T$",
        xy=(pair_point[1][0], pair_point[0][0]),
        xytext=(12, 3),
        textcoords="offset points",
        fontsize=9,
        color="#444444",
    )

    axis.set_title("Energy-entropy curves of early-universe radiation sectors")
    axis.set_xlabel(r"energy density $\rho$")
    axis.set_ylabel(r"entropy density $s$")
    axis.set_xlim(left=0.0)
    axis.set_ylim(bottom=0.0)
    axis.grid(True, color="#d8dde6", linewidth=0.7, alpha=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(loc="upper left", frameon=True, facecolor="white", fontsize=9)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, format="pdf", facecolor="white")
    plt.close(figure)
    print(OUTPUT)


if __name__ == "__main__":
    main()
