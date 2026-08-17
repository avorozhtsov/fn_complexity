#!/usr/bin/env python3
"""The energy--entropy curve of the completed zeta function.

The Riemann zeta function is a poor partition function: the Euler product
converges only for ``Re s > 1``, the pole at ``s = 1`` is a Hagedorn temperature
at which energy and entropy both diverge, and the region under the curve is
unbounded, so every exchange rate against it degenerates to ``0`` or ``infinity``.

The completion repairs all of it.  With

    xi(s) = (1/2) s (s-1) pi^{-s/2} Gamma(s/2) zeta(s),

Riemann's kernel gives

    xi(1/2 + b) = integral Phi(u) e^{b u} du,
    Phi(u) = sum_n (4 pi^2 n^4 e^{9u/2} - 6 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u}),

with ``Phi > 0``, ``Phi(-u) = Phi(u)``, and doubly exponential decay.  So ``xi``
recentred on the critical line is the two-sided Laplace transform of a positive
even measure -- a bona fide partition function, with

    Z(b) = xi(1/2 + b),   E(b) = -(log Z)'(b),   H(b) = log Z(b) + b E(b).

Three consequences, all checked here: ``log Z`` is finite and convex for every
real ``b`` (no Hagedorn), ``Z(b) = Z(-b)`` so the spectrum is symmetric, and
``E`` is odd while ``H`` is even, so the Gibbs curve is symmetric about the
entropy axis.  That symmetry *is* the functional equation.
"""

from __future__ import annotations

import csv
from pathlib import Path

from mpmath import diff, exp, gamma, inf, log, mp, mpf, nsum, pi, quad, zeta

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
IMAGE_DIRECTORY = Path(__file__).resolve().parents[1] / "paper_finite_fields_maps" / "images"
CSV_PATH = OUTPUT_DIRECTORY / "xi_gibbs_curve.csv"
SVG_PATH = IMAGE_DIRECTORY / "riemann-xi-gibbs-curve.svg"

mp.dps = 15
INK = "#172033"
CURVE = "#7c3aed"
ACCENT = "#b45309"


def xi(s):
    return mpf(1) / 2 * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def polya_kernel(u):
    return nsum(
        lambda n: (4 * pi**2 * n**4 * exp(mpf(9) * u / 2) - 6 * pi * n**2 * exp(mpf(5) * u / 2))
        * exp(-pi * n**2 * exp(2 * u)),
        [1, inf],
    )


def log_partition(beta):
    return log(xi(mpf(1) / 2 + beta))


def check_kernel() -> None:
    print("Riemann kernel: xi(1/2+b) = int Phi(u) e^{bu} du")
    for value in ("0", "1", "-2"):
        beta = mpf(value)
        left = xi(mpf(1) / 2 + beta)
        right = quad(lambda u: polya_kernel(u) * exp(beta * u), [-4, 4])
        print(
            f"   b={value:>5}   xi={mp.nstr(left, 14)}   integral={mp.nstr(right, 14)}"
            f"   rel. diff={mp.nstr(abs(left - right) / left, 3)}"
        )
    print("\nPhi is positive and even:")
    for value in ("0", "0.5", "1", "1.5"):
        u = mpf(value)
        print(
            f"   Phi(±{value:>4}) = {mp.nstr(polya_kernel(u), 10)}"
            f"   /  {mp.nstr(polya_kernel(-u), 10)}"
        )


def curve(betas):
    """Gibbs data, normalised so that H(0) = 0.

    Only non-negative beta is evaluated.  The naive product for ``xi`` hits the
    poles of ``Gamma(s/2)`` at ``s = 0, -2, -4, ...`` -- points where ``xi``
    itself is perfectly regular, the poles being cancelled by the trivial zeros
    of ``zeta`` -- so the negative half is obtained from the functional equation
    instead: ``log Z`` and ``H`` are even in beta and ``E`` is odd.
    """

    reference = log_partition(mpf(0))
    positive = []
    for beta in betas:
        b = mpf(beta)
        if b < 0:
            continue
        if abs(b - mpf(1) / 2) < mpf(10) ** -9:  # s = 1, the pole of zeta
            continue
        value = log_partition(b) - reference
        energy = -diff(log_partition, b)
        entropy = value + b * energy
        curvature = diff(log_partition, b, 2)
        positive.append((float(b), float(value), float(energy), float(entropy), float(curvature)))
    mirrored = [(-b, v, -e, h, c) for (b, v, e, h, c) in reversed(positive) if b > 0]
    return mirrored + positive


def check_reflection() -> None:
    """Confirm the mirrored half against direct evaluation.

    ``xi`` is regular at negative beta away from the Gamma poles, so the
    reflection used in ``curve`` can be checked rather than assumed.
    """

    reference = log_partition(mpf(0))
    print("\nreflection check: negative beta evaluated directly")
    print(f"{'beta':>7} {'log Z':>14} {'E':>14} {'H':>14}")
    for value in ("3", "-3", "2", "-2"):
        b = mpf(value)
        L = log_partition(b) - reference
        E = -diff(log_partition, b)
        print(f"{value:>7} {mp.nstr(L, 10):>14} {mp.nstr(E, 10):>14} {mp.nstr(L + b * E, 10):>14}")


def render(rows) -> Path:
    width, height = 1100, 720
    left, right, top, bottom = 130, width - 300, 150, height - 90
    energies = [row[2] for row in rows]
    entropies = [row[3] for row in rows]
    min_energy, max_energy = min(energies) * 1.08, max(energies) * 1.08
    min_entropy, max_entropy = min(entropies) * 1.08, 0.08

    def project(energy, entropy):
        x = left + (energy - min_energy) / (max_energy - min_energy) * (right - left)
        y = bottom - (entropy - min_entropy) / (max_entropy - min_entropy) * (bottom - top)
        return x, y

    pieces = []
    for index in range(7):
        value = min_energy + (max_energy - min_energy) * index / 6
        x, _ = project(value, 0)
        pieces.append(f'<line class="grid" x1="{x:g}" y1="{top}" x2="{x:g}" y2="{bottom}"/>')
        pieces.append(
            f'<text class="tick" x="{x:g}" y="{bottom + 26}" text-anchor="middle">{value:.2f}</text>'
        )
    for index in range(6):
        value = min_entropy + (max_entropy - min_entropy) * index / 5
        _, y = project(min_energy, value)
        pieces.append(f'<line class="grid" x1="{left}" y1="{y:g}" x2="{right}" y2="{y:g}"/>')
        pieces.append(
            f'<text class="tick" x="{left - 14}" y="{y + 5:g}" text-anchor="end">{value:.2f}</text>'
        )
    zero_x, _ = project(0.0, 0.0)
    pieces.append(f'<line class="axis-zero" x1="{zero_x:g}" y1="{top}" x2="{zero_x:g}" y2="{bottom}"/>')

    points = " ".join(f"{x:g},{y:g}" for x, y in (project(r[2], r[3]) for r in rows))
    pieces.append(f'<polyline class="curve" stroke="{CURVE}" points="{points}"/>')
    peak_x, peak_y = project(0.0, 0.0)
    pieces.append(
        f'<circle cx="{peak_x:g}" cy="{peak_y:g}" r="6" fill="#ffffff" '
        f'stroke="{ACCENT}" stroke-width="3"/>'
    )
    pieces.append(
        f'<text class="note" x="{peak_x + 14:g}" y="{peak_y - 12:g}">'
        "β = 0, the critical line</text>"
    )
    for beta, marker in ((3.0, "β = 3"), (-3.0, "β = −3")):
        row = min(rows, key=lambda r: abs(r[0] - beta))
        x, y = project(row[2], row[3])
        pieces.append(f'<circle cx="{x:g}" cy="{y:g}" r="4.5" fill="{CURVE}"/>')
        anchor = "start" if beta < 0 else "end"
        offset = 12 if beta < 0 else -12
        pieces.append(
            f'<text class="mark" x="{x + offset:g}" y="{y + 18:g}" '
            f'text-anchor="{anchor}">{marker}</text>'
        )

    legend_x = right + 40
    pieces.append(f'<text class="legend" x="{legend_x}" y="{top + 20}">Z(β) = ξ(½ + β)</text>')
    pieces.append(
        f'<text class="legend" x="{legend_x}" y="{top + 46}">= ∫ Φ(u) e^{{βu}} du,  Φ &gt; 0</text>'
    )
    pieces.append(f'<text class="legend" x="{legend_x}" y="{top + 82}">Φ even ⟹ Z(β) = Z(−β)</text>')
    pieces.append(f'<text class="legend" x="{legend_x}" y="{top + 108}">⟹ E odd, H even</text>')
    pieces.append(
        f'<text class="legend" x="{legend_x}" y="{top + 144}">the mirror symmetry is</text>'
    )
    pieces.append(
        f'<text class="legend" x="{legend_x}" y="{top + 170}">the functional equation</text>'
    )
    pieces.append(f'<text class="legend" x="{legend_x}" y="{top + 196}">ξ(s) = ξ(1 − s)</text>')

    body = "\n  ".join(pieces)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"
     viewBox="0 0 {width} {height}" role="img"
     aria-labelledby="title description">
  <title id="title">Energy-entropy curve of the completed zeta function</title>
  <desc id="description">The Gibbs curve of Z(b) = xi(1/2 + b). Because the
  Riemann kernel Phi is positive and even, the curve exists for every real b and
  is symmetric about the entropy axis; that symmetry is the functional
  equation.</desc>
  <defs>
    <style>
      .grid {{ stroke: #e2e8f0; stroke-width: 1; }}
      .axis {{ stroke: #94a3b8; stroke-width: 1.8; }}
      .axis-zero {{ stroke: #cbd5e1; stroke-width: 1.6; stroke-dasharray: 5 4; }}
      .curve {{ fill: none; stroke-width: 2.8; stroke-linejoin: round; }}
      .tick {{ fill: #64748b; font: 13px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .axis-label {{ fill: #475569; font: 15px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .legend {{ fill: #475569; font: 14px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .note {{ fill: {ACCENT}; font: 600 13px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .mark {{ fill: #475569; font: 13px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .title {{ fill: {INK}; font: 700 29px Inter, ui-sans-serif, system-ui, sans-serif; }}
      .subtitle {{ fill: #475569; font: 16px Inter, ui-sans-serif, system-ui, sans-serif; }}
    </style>
  </defs>
  <rect width="{width}" height="{height}" fill="#ffffff"/>
  <text class="title" x="{width / 2:g}" y="58" text-anchor="middle">
    The energy–entropy curve of the completed zeta function
  </text>
  <text class="subtitle" x="{width / 2:g}" y="90" text-anchor="middle">
    ζ has none — its region is unbounded at the pole. ξ does, for every real β.
  </text>
  {body}
  <line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}"/>
  <line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>
  <text class="axis-label" x="{(left + right) / 2:g}" y="{bottom + 58}"
        text-anchor="middle">mean energy E</text>
  <text class="axis-label" x="{left - 78}" y="{(top + bottom) / 2:g}"
        text-anchor="middle" transform="rotate(-90 {left - 78} {(top + bottom) / 2:g})">
    entropy H</text>
</svg>
'''
    SVG_PATH.parent.mkdir(parents=True, exist_ok=True)
    SVG_PATH.write_text(svg, encoding="utf-8")
    return SVG_PATH


def main() -> int:
    check_kernel()
    check_reflection()
    betas = [mpf(k) / 8 for k in range(0, 97)]
    rows = curve(betas)
    print("\nGibbs data, normalised so that H(0) = 0:")
    print(f"{'beta':>7} {'log Z':>12} {'E':>12} {'H':>12} {'(log Z)″':>12}")
    for row in rows:
        if abs(row[0] * 2 % 6) < 1e-9:
            print(f"{row[0]:>7.2f} {row[1]:>12.6f} {row[2]:>12.6f} {row[3]:>12.6f} {row[4]:>12.6f}")
    print(f"\nconvexity: min (log Z)″ = {min(r[4] for r in rows):.6f}  (positive ⟹ log Z convex)")
    print("symmetry: E is odd and H even by construction of the negative half;")
    print("          the functional equation is what licenses that reflection.")

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["beta", "log_Z", "energy", "entropy", "log_Z_second_derivative"])
        for row in rows:
            writer.writerow([f"{value:.15f}" for value in row])
    print(f"\nwritten to {CSV_PATH.name} and {render(rows).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
