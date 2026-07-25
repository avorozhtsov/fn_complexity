#!/usr/bin/env python3
"""Plot pairwise exchange-rate homotheties for a three-signature cycle."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import ExchangeRateCache, gibbs_point  # noqa: E402

SIGNATURES = ((5, 3), (3, 1, 1), (6, 1))
COLORS = ("#0072B2", "#D55E00", "#7E57C2")
OUTPUT_STEM = PROJECT_ROOT / "images" / "exchange-homotheties_cycle-5-3_3-1-1_6-1"


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


def signature_text(signature: tuple[int, ...]) -> str:
    return r"\{" + ",".join(map(str, signature)) + r"\}"


def main() -> int:
    cache = ExchangeRateCache()
    curves = {
        signature: curve(signature)
        for signature in SIGNATURES
    }

    plt.style.use("default")
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(18.0, 6.4),
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
            axis.plot(
                entropies,
                energies,
                color=color,
                linewidth=3.0 if source_index == target_index else 1.7,
                alpha=1.0 if source_index == target_index else 0.48,
                zorder=3 if source_index == target_index else 2,
            )

        dotted_handles: list[Line2D] = []
        for source_index, (source, color) in enumerate(
            zip(SIGNATURES, COLORS)
        ):
            if source == target:
                continue
            scale = 1.0 / cache.get(source, target)
            entropies, energies = curves[source]
            axis.plot(
                [scale * entropy for entropy in entropies],
                [scale * energy for energy in energies],
                color=color,
                linewidth=2.3,
                linestyle=(0, (2, 3)),
                zorder=4,
            )
            dotted_handles.append(
                Line2D(
                    [0],
                    [0],
                    color=color,
                    linewidth=2.3,
                    linestyle=(0, (2, 3)),
                    label=(
                        rf"${scale:.6f}\,f_{source_index + 1}$"
                        rf"$=f_{source_index + 1}/"
                        rf"C(f_{source_index + 1}\mid f_{target_index + 1})$"
                    ),
                )
            )

        axis.set_title(
            rf"Embrace $f_{target_index + 1}={signature_text(target)}$",
            loc="left",
            fontsize=16,
            color="#111827",
            pad=12,
        )
        axis.axhline(0.0, color="#6b7280", linewidth=1.0)
        axis.grid(color="#d1d5db", alpha=0.72, linewidth=0.8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.spines[["bottom", "left"]].set_color("#6b7280")
        axis.tick_params(colors="#374151")
        axis.legend(
            handles=dotted_handles,
            title="tight dotted homotheties",
            loc="lower right",
            fontsize=9,
            title_fontsize=10,
            frameon=True,
            facecolor="white",
            edgecolor="#9ca3af",
            labelcolor="#111827",
        )

    axes[0].set_ylabel(r"energy $E$", labelpad=10)
    for axis in axes:
        axis.set_xlabel(r"entropy $H$", labelpad=10)
        axis.set_xlim(-0.03, 1.84)
        axis.set_ylim(-3.0, 0.08)

    solid_handles = [
        Line2D(
            [0],
            [0],
            color=color,
            linewidth=2.6,
            label=rf"$f_{index}={signature_text(signature)}$",
        )
        for index, (signature, color) in enumerate(
            zip(SIGNATURES, COLORS), 1
        )
    ]
    figure.legend(
        handles=solid_handles,
        loc="outside upper center",
        ncols=3,
        title="solid original curves",
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
