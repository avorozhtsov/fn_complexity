#!/usr/bin/env python3
"""Generate entropy-energy curves for the first prime oscillator modes."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "paper_finite_fields_maps"
    / "images"
    / "prime-mode-entropy-energy-curves.svg"
)
DEFAULT_PRIMES = (2, 3, 5, 7, 11, 13)


def normalize_svg(output: Path) -> None:
    """Remove Matplotlib's line-ending spaces from generated SVG markup."""

    text = output.read_text(encoding="utf-8")
    output.write_text(
        "\n".join(line.rstrip() for line in text.splitlines()) + "\n",
        encoding="utf-8",
    )


def entropy_from_energy(energy: np.ndarray, prime: int) -> np.ndarray:
    """Return the oscillator entropy at mean energy U for one prime mode."""

    occupation = energy / math.log(prime)
    entropy = np.zeros_like(occupation)
    positive = occupation > 0.0
    values = occupation[positive]
    entropy[positive] = (1.0 + values) * np.log1p(values) - values * np.log(values)
    return entropy


def render(output: Path, *, primes: tuple[int, ...], show_title: bool) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "svg.fonttype": "none",
        }
    )
    figure, axis = plt.subplots(figsize=(7.4, 3.25))
    if show_title:
        axis.set_title("Prime-mode entropy-energy curves")

    energy = np.linspace(0.0, 20.0, 700)
    colors = ("#1f5f8b", "#b3532d", "#4b7f52", "#795a9b", "#9a7b20", "#447a78")
    line_styles = ("-", "--", "-.", ":", (0, (5, 1.5)), (0, (3, 1, 1, 1)))
    label_offsets = {11: 0.06, 13: -0.07}

    for index, prime in enumerate(primes):
        entropy = entropy_from_energy(energy, prime)
        color = colors[index % len(colors)]
        axis.plot(
            energy,
            entropy,
            color=color,
            linestyle=line_styles[index % len(line_styles)],
            linewidth=2.0,
        )
        axis.text(
            20.35,
            entropy[-1] + label_offsets.get(prime, 0.0),
            rf"$p={prime}$",
            color=color,
            va="center",
        )

    axis.text(
        0.03,
        0.95,
        r"levels $E_k=k\log p$, $k=0,1,2,\ldots$",
        transform=axis.transAxes,
        va="top",
    )
    axis.annotate(
        r"$\beta\downarrow0$: $U_p,S_p\to\infty$",
        xy=(19.0, entropy_from_energy(np.array([19.0]), primes[0])[0]),
        xytext=(12.2, 4.42),
        arrowprops={"arrowstyle": "->", "color": "#444444", "linewidth": 0.8},
        ha="center",
    )
    axis.annotate(
        r"$\beta\to\infty$",
        xy=(0.03, 0.03),
        xytext=(2.7, 0.55),
        arrowprops={"arrowstyle": "->", "color": "#444444", "linewidth": 0.8},
    )
    axis.set_xlim(0.0, 22.0)
    axis.set_ylim(0.0, 4.7)
    axis.set_xlabel(r"mean energy $U_p=\log p/(p^\beta-1)$")
    axis.set_ylabel(r"entropy $S_p$")
    axis.grid(color="#d8d8d8", linewidth=0.6)
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, format="svg", bbox_inches="tight")
    plt.close(figure)
    normalize_svg(output)


def parse_primes(value: str) -> tuple[int, ...]:
    primes = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not primes or any(prime < 2 for prime in primes):
        raise argparse.ArgumentTypeError("primes must be comma-separated integers >= 2")
    return primes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--primes",
        type=parse_primes,
        default=DEFAULT_PRIMES,
        help="comma-separated primes (default: 2,3,5,7,11,13)",
    )
    title_group = parser.add_mutually_exclusive_group()
    title_group.add_argument("--titles", dest="show_title", action="store_true")
    title_group.add_argument("--no-titles", dest="show_title", action="store_false")
    parser.set_defaults(show_title=True)
    arguments = parser.parse_args()
    render(arguments.output, primes=arguments.primes, show_title=arguments.show_title)
    print(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
