#!/usr/bin/env python3
"""Generate the zeta Gibbs entropy-energy curve used by the paper."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "paper_finite_fields_maps"
    / "images"
    / "zeta-entropy-energy.svg"
)


def normalize_svg(output: Path) -> None:
    """Remove Matplotlib's line-ending spaces from generated SVG markup."""

    text = output.read_text(encoding="utf-8")
    output.write_text(
        "\n".join(line.rstrip() for line in text.splitlines()) + "\n",
        encoding="utf-8",
    )


def thermodynamic_point(beta: float) -> tuple[float, float]:
    """Return mean energy U and entropy S for beta > 1."""

    zeta = mp.zeta(beta)
    mean_energy = -mp.diff(mp.zeta, beta) / zeta
    entropy = mp.log(zeta) + beta * mean_energy
    return float(mean_energy), float(entropy)


def render(output: Path, *, show_title: bool) -> None:
    mp.mp.dps = 40
    beta = 1.0 + np.geomspace(0.03, 12.0, 600)
    points = np.array([thermodynamic_point(float(value)) for value in beta])

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "svg.fonttype": "none",
        }
    )
    figure, (curve_axis, domain_axis) = plt.subplots(
        1,
        2,
        figsize=(7.4, 3.25),
        gridspec_kw={"width_ratios": (1.55, 1.15)},
    )
    if show_title:
        figure.suptitle(r"Zeta Gibbs entropy-energy curve")

    curve_axis.plot(
        points[:, 0],
        points[:, 1],
        color="#1f5f8b",
        linewidth=2.0,
    )
    marker_betas = (1.05, 1.1, 1.2, 1.5, 2.0)
    marker_offsets = {
        1.05: (-34, -2),
        1.1: (6, -2),
        1.2: (6, -2),
        1.5: (6, -2),
        2.0: (8, 4),
    }
    for value in marker_betas:
        energy, entropy = thermodynamic_point(value)
        curve_axis.scatter(
            [energy],
            [entropy],
            color="#1f5f8b",
            edgecolor="white",
            linewidth=0.7,
            s=24,
            zorder=3,
        )
        curve_axis.annotate(
            rf"$\beta={value:g}$",
            (energy, entropy),
            xytext=marker_offsets[value],
            textcoords="offset points",
            va="center",
        )
    curve_axis.annotate(
        r"$\beta\downarrow1$: $U,S\to\infty$",
        xy=(19.4, 23.4),
        xytext=(11.0, 21.0),
        arrowprops={"arrowstyle": "->", "color": "#444444", "linewidth": 0.8},
        ha="center",
    )
    curve_axis.annotate(
        r"$\beta\to\infty$",
        xy=(0.03, 0.17),
        xytext=(4.0, 0.9),
        arrowprops={"arrowstyle": "->", "color": "#444444", "linewidth": 0.8},
    )
    curve_axis.set_xlim(0.0, 21.0)
    curve_axis.set_ylim(0.0, 25.0)
    curve_axis.set_xlabel(r"mean energy $U=-\partial_\beta\log\zeta(\beta)$")
    curve_axis.set_ylabel(r"entropy $S=\log\zeta(\beta)+\beta U$")
    curve_axis.grid(color="#d8d8d8", linewidth=0.6)
    curve_axis.text(
        0.02,
        0.95,
        r"canonical branch: $\beta>1$",
        transform=curve_axis.transAxes,
        va="top",
    )

    domain_axis.axvspan(-4.0, 0.0, color="#c05a5a", alpha=0.20)
    domain_axis.axvspan(0.0, 1.0, color="#c05a5a", alpha=0.10)
    domain_axis.axvspan(1.0, 5.0, color="#4b8b63", alpha=0.18)
    domain_axis.axvline(0.0, color="#666666", linewidth=0.8)
    domain_axis.axvline(1.0, color="#333333", linewidth=1.1)
    domain_axis.hlines(0.48, -4.0, 1.0, color="#9b3f3f", linewidth=2.0)
    domain_axis.hlines(0.48, 1.0, 5.0, color="#39764f", linewidth=2.0)
    domain_axis.text(
        -1.5,
        0.76,
        "no Gibbs state\n" + r"$Z(\beta)=\infty$",
        ha="center",
        va="center",
        color="#7e2f2f",
    )
    domain_axis.text(
        -2.0,
        0.25,
        "negative T\n" + r"$\beta<0$",
        ha="center",
        va="center",
    )
    domain_axis.text(
        0.5,
        0.25,
        r"$0\leq\beta\leq1$",
        ha="center",
        va="center",
    )
    domain_axis.text(
        3.0,
        0.76,
        "Gibbs state\n" + r"$Z(\beta)=\zeta(\beta)<\infty$",
        ha="center",
        va="center",
        color="#2d6340",
    )
    domain_axis.annotate(
        r"$\beta=1$ threshold",
        xy=(1.0, 0.49),
        xytext=(2.6, 0.22),
        arrowprops={"arrowstyle": "->", "color": "#444444", "linewidth": 0.8},
        ha="center",
    )
    domain_axis.set_xlim(-4.0, 5.0)
    domain_axis.set_ylim(0.0, 1.08)
    domain_axis.set_xlabel(r"inverse temperature $\beta=1/T$")
    domain_axis.set_yticks([])
    domain_axis.spines["left"].set_visible(False)
    domain_axis.spines["right"].set_visible(False)
    domain_axis.spines["top"].set_visible(False)
    domain_axis.set_xticks((-4, -2, 0, 1, 3, 5))

    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, format="svg", bbox_inches="tight")
    plt.close(figure)
    normalize_svg(output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    title_group = parser.add_mutually_exclusive_group()
    title_group.add_argument("--titles", dest="show_title", action="store_true")
    title_group.add_argument("--no-titles", dest="show_title", action="store_false")
    parser.set_defaults(show_title=True)
    arguments = parser.parse_args()
    render(arguments.output, show_title=arguments.show_title)
    print(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
