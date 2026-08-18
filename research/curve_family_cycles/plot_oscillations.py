#!/usr/bin/env python3
"""The picture of the cycle: three curves whose midranges are all negative.

Put ``u_a(beta) = log log Z_a(beta)``.  Then, exactly,

    C(a->b) = exp(min_beta (u_a - u_b)),   C(b->a) = exp(-max_beta (u_a - u_b)),

so the exchange distance is the oscillation of ``u_a - u_b`` and the comparison
is the sign of its midrange.  The index ``phi`` only sees the two endpoint
values ``u_a(0) - u_b(0)`` and ``u_a(inf) - u_b(inf)``.

The figure draws the three differences of the certified ``F_11`` cycle against
``log beta``.  Each has its deepest excursion at a different temperature, each
midrange is below zero, and the three curves sum to zero pointwise --- which is
the whole mechanism.

    python research/curve_family_cycles/plot_oscillations.py
"""

from __future__ import annotations

import math
from pathlib import Path

Q = 11
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "curve_family_cycle_f11.svg"

A = (18, 16, 15, 15, 14, 12, 9, 6, 6, 5, 5)
B = (18, 18, 14, 13, 12, 9, 9, 9, 8, 7, 4)
C = (19, 14, 12, 11, 11, 10, 10, 10, 9, 9, 6)
PAIRS = (("A", A, "B", B, "#2563eb"), ("B", B, "C", C, "#be123c"), ("C", C, "A", A, "#0f766e"))

WIDTH, HEIGHT = 780, 470
LEFT, RIGHT, TOP, BOTTOM = 86, 660, 40, 350
LOG_LOW, LOG_HIGH = -1.0, 4.0


def u(signature: tuple[int, ...], beta: float) -> float:
    logs = [math.log(v) for v in signature]
    top = max(logs)
    return math.log(beta * top + math.log(sum(math.exp(beta * (l - top)) for l in logs)))


def difference(first, second, beta: float) -> float:
    return u(first, beta) - u(second, beta)


def lam(signature: tuple[int, ...], beta: float) -> float:
    """``Lambda(beta) = log((1/q) sum_c (N_c/q)^beta)``, the exact reduced form."""

    ells = [math.log(v / Q) for v in signature]
    top = max(ells)
    return beta * top + math.log(sum(math.exp(beta * (e - top)) for e in ells) / len(ells))


def psi_value(signature: tuple[int, ...], tau: float) -> float:
    """``Psi(tau) = Lambda_hat(tau)/tau`` for the normalised traces."""

    alphas = [(v - Q) / math.sqrt(Q) for v in signature]
    top = max(alphas)
    total = math.log(sum(math.exp(tau * (a - top)) for a in alphas) / len(alphas))
    return (tau * top + total) / tau


def main() -> int:
    samples = [LOG_LOW + (LOG_HIGH - LOG_LOW) * i / 900 for i in range(901)]
    series = []
    for name_a, first, name_b, second, colour in PAIRS:
        values = [difference(first, second, 10.0**s) for s in samples]
        endpoint = math.log(math.log(max(first))) - math.log(math.log(max(second)))
        low = min(min(values), endpoint, 0.0)
        high = max(max(values), endpoint, 0.0)
        series.append((name_a, name_b, colour, values, endpoint, low, high))

    span = max(max(abs(s[5]), abs(s[6])) for s in series) * 1.15

    def x_of(log_beta: float) -> float:
        return LEFT + (log_beta - LOG_LOW) / (LOG_HIGH - LOG_LOW) * (RIGHT - LEFT)

    def y_of(value: float) -> float:
        middle = (TOP + BOTTOM) / 2
        return middle - value / span * (BOTTOM - TOP) / 2

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" font-family="Georgia, serif">',
        f'<rect width="{WIDTH}" height="{HEIGHT}" fill="#ffffff"/>',
        f'<text x="{LEFT}" y="24" font-size="15" fill="#111827">'
        "Three genus-two pencils over F_11: the differences "
        "u_a - u_b, u_a = log log Z_a</text>",
        f'<line x1="{LEFT}" y1="{y_of(0):.1f}" x2="{RIGHT}" y2="{y_of(0):.1f}" '
        'stroke="#9ca3af" stroke-width="1"/>',
        f'<line x1="{LEFT}" y1="{TOP}" x2="{LEFT}" y2="{BOTTOM}" stroke="#9ca3af" stroke-width="1"/>',
    ]
    for decade in range(int(LOG_LOW), int(LOG_HIGH) + 1):
        x = x_of(decade)
        out.append(
            f'<line x1="{x:.1f}" y1="{BOTTOM}" x2="{x:.1f}" y2="{BOTTOM + 5}" stroke="#9ca3af"/>'
        )
        out.append(
            f'<text x="{x:.1f}" y="{BOTTOM + 20}" font-size="11" fill="#374151" '
            f'text-anchor="middle">10^{decade}</text>'
        )
    out.append(
        f'<text x="{(LEFT + RIGHT) / 2:.0f}" y="{BOTTOM + 40}" font-size="12" '
        'fill="#374151" text-anchor="middle">inverse temperature beta</text>'
    )
    for tick in (-0.02, -0.01, 0.0, 0.01, 0.02):
        y = y_of(tick)
        if not TOP - 5 <= y <= BOTTOM + 5:
            continue
        out.append(f'<line x1="{LEFT - 5}" y1="{y:.1f}" x2="{LEFT}" y2="{y:.1f}" stroke="#9ca3af"/>')
        out.append(
            f'<text x="{LEFT - 9}" y="{y + 4:.1f}" font-size="11" fill="#374151" '
            f'text-anchor="end">{tick:+.2f}</text>'
        )

    legend = TOP + 6
    for order, (name_a, name_b, colour, values, endpoint, low, high) in enumerate(series):
        path = " ".join(
            f"{'M' if index == 0 else 'L'}{x_of(s):.1f},{y_of(v):.1f}"
            for index, (s, v) in enumerate(zip(samples, values))
        )
        out.append(f'<path d="{path}" fill="none" stroke="{colour}" stroke-width="2"/>')
        # the beta = infinity endpoint, drawn as a stub past the right edge
        out.append(
            f'<line x1="{RIGHT}" y1="{y_of(endpoint):.1f}" x2="{RIGHT + 26}" '
            f'y2="{y_of(endpoint):.1f}" stroke="{colour}" stroke-width="2" '
            'stroke-dasharray="4 3"/>'
        )
        midrange = (low + high) / 2
        bar = RIGHT + 34 + 13 * order
        out.append(
            f'<line x1="{bar}" y1="{y_of(midrange):.1f}" x2="{bar + 10}" '
            f'y2="{y_of(midrange):.1f}" stroke="{colour}" stroke-width="3"/>'
        )
        out.append(
            f'<text x="{LEFT + 8}" y="{legend}" font-size="12" fill="{colour}">'
            f"u_{name_a} - u_{name_b}:  midrange {midrange:+.2e}  "
            f"({name_a} precedes {name_b})</text>"
        )
        legend += 17

    out.append(
        f'<text x="{RIGHT + 46}" y="{BOTTOM + 20}" font-size="11" fill="#374151" '
        'text-anchor="middle">inf</text>'
    )
    out.append(
        f'<text x="{LEFT}" y="{HEIGHT - 42}" font-size="12" fill="#374151">'
        "Dashed stubs are the beta = infinity endpoint values, the only part of each"
        " curve the index phi reads;</text>"
    )
    out.append(
        f'<text x="{LEFT}" y="{HEIGHT - 25}" font-size="12" fill="#374151">'
        "thick bars are the midranges, whose signs are the comparison. The three "
        "curves sum to zero pointwise.</text>"
    )
    out.append("</svg>")
    OUTPUT.write_text("\n".join(out), encoding="utf-8")
    print(f"written {OUTPUT.relative_to(HERE.parents[1])}")
    for name_a, name_b, _, values, endpoint, low, high in series:
        print(
            f"  u_{name_a} - u_{name_b}: min {low:+.6f}  max {high:+.6f}  "
            f"midrange {(low + high) / 2:+.3e}  oscillation {high - low:.6f}"
        )

    print(
        "\nThe two leading-order forms of the midrange, against the exact one.\n"
        "Both drop Lambda_v/((1+beta) log q), nominally O(1/(sqrt(q) log q)) = 0.13\n"
        "at q = 11 but in fact reaching 0.23, and both get the third edge wrong:\n"
    )
    print(f"  {'edge':<10}{'exact':>14}{'reduced Lambda':>18}{'alpha-scaled Psi':>20}")
    betas = [10.0**s for s in samples]
    taus = [b / math.sqrt(Q) for b in betas]
    for name_a, first, name_b, second, _ in PAIRS:
        exact = next(
            (low + high) / 2 for a, b, _, _, _, low, high in series if (a, b) == (name_a, name_b)
        )
        reduced = [
            (lam(first, b) - lam(second, b)) / (1.0 + b) for b in betas
        ] + [math.log(max(first) / max(second)), 0.0]
        scaled = [
            psi_value(first, t) - psi_value(second, t) for t in taus
        ] + [(max(first) - max(second)) / math.sqrt(Q), 0.0]
        print(
            f"  {name_a} < {name_b}   {exact:>+13.3e}"
            f"{(min(reduced) + max(reduced)) / (2 * math.log(Q)):>+18.3e}"
            f"{(min(scaled) + max(scaled)) / (2 * math.sqrt(Q) * math.log(Q)):>+20.3e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
