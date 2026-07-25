#!/usr/bin/env python3
"""Plot a pair whose two directed rates have finite-temperature tangencies."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate_result, gibbs_point  # noqa: E402

SIGNATURES = ((6, 5, 2, 1), (6, 4, 3, 2))
COLORS = ("#0072B2", "#D55E00")
OUTPUT_STEM = (
    PROJECT_ROOT / "images" / "exchange-homotheties_6-5-2-1_6-4-3-2"
)
FILL_ALPHA = 0.15


def temperatures() -> tuple[float, ...]:
    finite = tuple(
        10.0 ** (-4.0 + 8.0 * index / 799.0)
        for index in range(800)
    )
    return (0.0,) + finite + (math.inf,)


def curve(signature: tuple[int, ...]) -> tuple[list[float], list[float]]:
    points = tuple(
        gibbs_point(signature, temperature)
        for temperature in temperatures()
    )
    return (
        [point.entropy for point in points],
        [point.energy for point in points],
    )


def signature_text(signature: tuple[int, ...]) -> str:
    return r"\{" + ",".join(map(str, signature)) + r"\}"


def fill_region(
    axis,
    entropies: list[float],
    energies: list[float],
    color: str,
    *,
    alpha: float = FILL_ALPHA,
) -> None:
    """Fill the closed Gibbs region horizontally from its boundary to E=0."""

    axis.fill(
        energies + [0.0, 0.0, energies[0]],
        entropies + [entropies[-1], 0.0, 0.0],
        color=color,
        alpha=alpha,
        linewidth=0.0,
        zorder=0,
    )


def endpoint_lines(
    axis,
    entropies: list[float],
    energies: list[float],
    color: str,
    *,
    alpha: float,
    infinity_linestyle: str | tuple = "solid",
    zero_temperature_linestyle: str | tuple = "solid",
) -> None:
    axis.hlines(
        entropies[-1],
        0.0,
        energies[-1],
        color=color,
        linewidth=1.0,
        linestyle=infinity_linestyle,
        alpha=alpha,
        zorder=1,
    )
    if entropies[0] > 1.0e-12:
        axis.vlines(
            energies[0],
            0.0,
            entropies[0],
            color=color,
            linewidth=1.0,
            linestyle=zero_temperature_linestyle,
            alpha=alpha,
            zorder=1,
        )


def main() -> int:
    curves = {
        signature: curve(signature)
        for signature in SIGNATURES
    }

    plt.style.use("default")
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12.0, 5.4),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    figure.patch.set_facecolor("white")

    for target_index, (axis, target) in enumerate(zip(axes, SIGNATURES)):
        axis.set_facecolor("white")
        for source_index, (signature, color) in enumerate(
            zip(SIGNATURES, COLORS)
        ):
            entropies, energies = curves[signature]
            strong = source_index == target_index
            if strong:
                fill_region(axis, entropies, energies, color)
            axis.plot(
                energies,
                entropies,
                color=color,
                linewidth=3.0 if strong else 1.8,
                alpha=1.0 if strong else 0.45,
                zorder=3 if strong else 2,
            )
            endpoint_lines(
                axis,
                entropies,
                energies,
                color,
                alpha=0.68 if strong else 0.28,
            )

        source_index = 1 - target_index
        source = SIGNATURES[source_index]
        result = exchange_rate_result(
            source,
            target,
            grid_size=8192,
            tolerance=1e-14,
        )
        scale = 1.0 / result.rate
        temperature = result.temperature

        entropies, energies = curves[source]
        scaled_entropies = [scale * value for value in entropies]
        scaled_energies = [scale * value for value in energies]
        fill_region(
            axis,
            scaled_entropies,
            scaled_energies,
            COLORS[source_index],
        )
        axis.plot(
            scaled_energies,
            scaled_entropies,
            color=COLORS[source_index],
            linewidth=2.6,
            linestyle=(0, (2, 3)),
            zorder=5,
        )
        endpoint_lines(
            axis,
            scaled_entropies,
            scaled_energies,
            COLORS[source_index],
            alpha=0.78,
            infinity_linestyle=(0, (2, 3)),
            zero_temperature_linestyle=(0, (2, 3)),
        )

        target_contact = gibbs_point(target, temperature)
        source_contact = gibbs_point(source, temperature)
        scaled_source_contact = (
            scale * source_contact.energy,
            scale * source_contact.entropy,
        )
        axis.scatter(
            [target_contact.energy],
            [target_contact.entropy],
            facecolor="white",
            edgecolor=COLORS[target_index],
            linewidth=2.0,
            s=72,
            zorder=7,
        )
        axis.scatter(
            [scaled_source_contact[0]],
            [scaled_source_contact[1]],
            color=COLORS[source_index],
            marker="x",
            linewidth=2.0,
            s=58,
            zorder=8,
        )
        axis.annotate(
            rf"$T_*={temperature:.6f}$",
            xy=(target_contact.energy, target_contact.entropy),
            xytext=(12, 18) if target_index == 0 else (-112, -42),
            textcoords="offset points",
            color="#111827",
            fontsize=9,
            arrowprops={
                "arrowstyle": "->",
                "color": "#6b7280",
                "linewidth": 0.9,
            },
            zorder=9,
        )

        target_text = signature_text(target)
        source_text = signature_text(source)
        axis.set_title(
            rf"Target ${target_text}$",
            loc="left",
            fontsize=16,
            color="#111827",
            pad=12,
        )
        axis.legend(
            handles=[
                Line2D(
                    [0],
                    [0],
                    color=COLORS[source_index],
                    linewidth=2.6,
                    linestyle=(0, (2, 3)),
                    label=rf"${scale:.6f}\,{source_text}$",
                ),
                Line2D(
                    [0],
                    [0],
                    color="#111827",
                    marker="x",
                    linestyle="None",
                    markersize=7,
                    label="finite-temperature contact",
                ),
            ],
            title="tight homothety",
            loc="lower right",
            fontsize=8.5,
            title_fontsize=10,
            frameon=True,
            facecolor="white",
            edgecolor="#9ca3af",
            labelcolor="#111827",
        )
        axis.axvline(0.0, color="#6b7280", linewidth=1.0)
        axis.grid(color="#d1d5db", alpha=0.72, linewidth=0.8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.spines[["bottom", "left"]].set_color("#6b7280")
        axis.tick_params(colors="#374151")
        axis.set_xlabel(r"energy $E$", labelpad=10)
        axis.set_xlim(-1.93, 0.07)
        axis.set_ylim(-0.03, 1.53)

    axes[0].set_ylabel(r"entropy $H$", labelpad=10)
    figure.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=color,
                linewidth=2.6,
                label="{" + ",".join(map(str, signature)) + "}",
            )
            for signature, color in zip(SIGNATURES, COLORS)
        ],
        loc="outside upper center",
        ncols=2,
        title="solid original Gibbs curves",
        frameon=False,
        labelcolor="#111827",
    )

    OUTPUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        OUTPUT_STEM.with_suffix(".pdf"),
        facecolor=figure.get_facecolor(),
    )
    plt.close(figure)
    print(OUTPUT_STEM.with_suffix(".pdf"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
