#!/usr/bin/env python3
"""Plot Gibbs energy as a function of entropy for selected signatures."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import gibbs_point  # noqa: E402

SIGNATURES = ((5, 3), (3, 1, 1), (6, 1))
COLORS = ("#0072B2", "#D55E00", "#7E57C2")
OUTPUT_STEM = (
    PROJECT_ROOT / "images" / "energy-entropy_5-3_3-1-1_6-1_ho"
)
HO_GAMMA = 1.0
HO_E0 = -math.log(6.0)


def temperatures() -> tuple[float, ...]:
    """Resolve both low- and high-temperature ends of each Gibbs curve."""

    finite = tuple(
        10.0 ** (-4.0 + 8.0 * index / 599.0)
        for index in range(600)
    )
    return (0.0,) + finite + (math.inf,)


def signature_label(signature: tuple[int, ...]) -> str:
    return r"$\{" + ",".join(map(str, signature)) + r"\}$"


def harmonic_oscillator_curve() -> tuple[tuple[float, float], ...]:
    """Exact infinite-spectrum HO curve, parameterized by temperature."""

    result = [(0.0, HO_E0)]
    for index in range(600):
        temperature = 10.0 ** (
            -4.0 + math.log10(1.25e4) * index / 599.0
        )
        inverse_temperature_gap = HO_GAMMA / temperature
        if inverse_temperature_gap > 40.0:
            mean_occupation = math.exp(-inverse_temperature_gap)
        else:
            mean_occupation = 1.0 / math.expm1(inverse_temperature_gap)
        entropy = 0.0 if mean_occupation == 0.0 else (
            (mean_occupation + 1.0) * math.log1p(mean_occupation)
            - mean_occupation * math.log(mean_occupation)
        )
        energy = HO_E0 + HO_GAMMA * mean_occupation
        result.append((entropy, energy))
    return tuple(result)


def main() -> int:
    plt.style.use("default")
    figure, axis = plt.subplots(figsize=(10.8, 7.0), constrained_layout=True)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    for signature, color in zip(SIGNATURES, COLORS):
        points = tuple(
            gibbs_point(signature, temperature)
            for temperature in temperatures()
        )
        entropies = [point.entropy for point in points]
        energies = [point.energy for point in points]
        axis.plot(
            energies,
            entropies,
            color=color,
            linewidth=2.6,
            label=signature_label(signature),
        )
        axis.hlines(
            entropies[-1],
            0.0,
            energies[-1],
            color=color,
            linewidth=2.0,
            linestyle=(0, (2, 3)),
            alpha=0.72,
            zorder=1,
        )
        axis.scatter(
            [energies[0], energies[-1]],
            [entropies[0], entropies[-1]],
            color=color,
            edgecolor="white",
            linewidth=0.8,
            s=42,
            zorder=3,
        )

    ho_points = harmonic_oscillator_curve()
    axis.plot(
        [point[1] for point in ho_points],
        [point[0] for point in ho_points],
        color="#009E73",
        linewidth=2.6,
        label=r"$\mathrm{HO}(1,-\log 6)$",
    )
    axis.scatter(
        [HO_E0],
        [0.0],
        color="#009E73",
        marker="s",
        s=47,
        edgecolor="white",
        linewidth=0.8,
        zorder=4,
    )

    axis.annotate(
        r"$T=0$",
        xy=(-math.log(6), 0.0),
        xytext=(-1.72, 0.10),
        color="#374151",
        arrowprops={"arrowstyle": "->", "color": "#6b7280"},
    )
    axis.annotate(
        r"$T\to\infty$",
        xy=(-math.log(3) / 3, math.log(3)),
        xytext=(-0.72, 1.02),
        color="#374151",
        arrowprops={"arrowstyle": "->", "color": "#6b7280"},
    )

    axis.set_title(
        "Gibbs energy–entropy curves",
        loc="left",
        fontsize=19,
        fontweight="semibold",
        color="#111827",
        pad=18,
    )
    axis.text(
        0.0,
        1.015,
        r"Map levels $\epsilon_i=-\log n_i$; oscillator levels "
        r"$\epsilon_n=-\log 6+n$",
        transform=axis.transAxes,
        color="#4b5563",
        fontsize=11,
    )
    axis.set_xlabel(r"energy $E=\sum_i p_i\epsilon_i$", labelpad=12)
    axis.set_ylabel(r"entropy $H=-\sum_i p_i\log p_i$", labelpad=12)
    axis.set_xlim(-1.87, 0.06)
    axis.set_ylim(-0.025, 1.15)
    axis.grid(color="#d1d5db", alpha=0.72, linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["bottom", "left"]].set_color("#6b7280")
    axis.tick_params(colors="#374151")
    axis.legend(
        title="signature",
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor="#d1d5db",
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
