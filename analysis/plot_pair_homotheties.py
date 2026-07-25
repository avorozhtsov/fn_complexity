#!/usr/bin/env python3
"""Plot the two directed homotheties for {2,2} and {3,1}."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import ExchangeRateCache, gibbs_point  # noqa: E402

SIGNATURES = ((2, 2), (3, 1))
COLORS = ("#0072B2", "#D55E00")
OUTPUT_STEM = PROJECT_ROOT / "images" / "exchange-homotheties_2-2_3-1"
FILL_ALPHA = 0.15


def temperatures() -> tuple[float, ...]:
    finite = tuple(
        10.0 ** (-4.0 + 8.0 * index / 599.0)
        for index in range(600)
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


def main() -> int:
    cache = ExchangeRateCache()
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
            axis.hlines(
                entropies[-1],
                0.0,
                energies[-1],
                color=color,
                linewidth=1.0,
                alpha=0.68 if strong else 0.28,
            )
            if entropies[0] > 1.0e-12:
                axis.vlines(
                    energies[0],
                    0.0,
                    entropies[0],
                    color=color,
                    linewidth=1.0,
                    alpha=0.68 if strong else 0.28,
                )
            axis.scatter(
                [energies[-1]],
                [entropies[-1]],
                color=color,
                s=35,
                zorder=4,
            )

        source_index = 1 - target_index
        source = SIGNATURES[source_index]
        scale = 1.0 / cache.get(source, target)
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
        axis.hlines(
            scaled_entropies[-1],
            0.0,
            scaled_energies[-1],
            color=COLORS[source_index],
            linewidth=1.0,
            linestyle=(0, (2, 3)),
            alpha=0.78,
            zorder=4,
        )
        if scaled_entropies[0] > 1.0e-12:
            axis.vlines(
                scaled_energies[0],
                0.0,
                scaled_entropies[0],
                color=COLORS[source_index],
                linewidth=1.0,
                linestyle=(0, (2, 3)),
                alpha=0.78,
                zorder=4,
            )
        axis.scatter(
            [scaled_energies[-1]],
            [scaled_entropies[-1]],
            color=COLORS[source_index],
            marker="D",
            s=38,
            zorder=6,
        )

        target_text = r"\{" + ",".join(map(str, target)) + r"\}"
        source_text = r"\{" + ",".join(map(str, source)) + r"\}"
        axis.set_title(
            rf"Target ${target_text}$",
            loc="left",
            fontsize=16,
            color="#111827",
            pad=12,
        )
        dotted_handle = Line2D(
            [0],
            [0],
            color=COLORS[source_index],
            linewidth=2.6,
            linestyle=(0, (2, 3)),
            label=rf"${scale:.6f}\,{source_text}$",
        )
        axis.legend(
            handles=[dotted_handle],
            title="tight homothety",
            loc="lower right",
            fontsize=9,
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
        axis.set_xlim(-1.82, 0.07)
        axis.set_ylim(-0.03, 1.18)

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
