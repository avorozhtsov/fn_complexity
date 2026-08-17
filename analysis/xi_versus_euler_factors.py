#!/usr/bin/env python3
"""Exchange rates between the completed zeta function and truncated Euler factors.

The completion makes ``xi`` a genuine partition function, but it does not make
it comparable to a finite resource.  With the mass-one normalisation
``Z(b) = xi(1/2+b)/xi(1/2)`` the two ends of the temperature range both fail:

    b -> 0     log Z(b) ~ (1/2) (log Z)''(0) b^2 = 0.0231 b^2  ->  0,
    b -> inf   log Z(b) ~ (b/2) log b                          ->  faster than linear,

whereas a truncated Euler factor ``P(p,K) = {1, p, ..., p^K}`` has
``log Z(0) = log(K+1) > 0`` and ``log Z(b) ~ K b log p``, linear.  So the ratio
of the two profiles tends to ``0`` at one end and ``infinity`` at the other, and
both unrestricted rates vanish:

    C(xi -> P) = C(P -> xi) = 0.

This is the first limitation stated in the companion paper -- an unbounded
spectrum with a different growth exponent gives nothing -- and it is why the
comparison has to be restricted to a temperature window, exactly as was done
there for the Standard Model plasma.

On a window ``W = [b1, b2]`` the rates are finite, attained at opposite ends of
the window, and the script reports them together with the irreversibility
``d = -log(C(xi->P) C(P->xi))``.  Minimising ``d`` over ``K`` answers the
currency question: how many units of ``p``-adic currency the completed zeta
function is worth on that window.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

from mpmath import gamma, log, mp, mpf, pi, zeta

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "xi_versus_euler_factors.csv"

mp.dps = 15
WINDOWS = [(0.5, 5.0), (1.0, 10.0), (2.0, 20.0)]
PRIMES = [2, 3, 5, 7]
TRUNCATIONS = [1, 2, 4, 8, 16, 32, 64, 128]
GRID = 240


def xi(s):
    return mpf(1) / 2 * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def xi_profile(betas: list[float]) -> list[float]:
    """``log Z(b)`` for the mass-one normalisation, positive for ``b > 0``."""

    reference = log(xi(mpf(1) / 2))
    values = []
    for beta in betas:
        b = mpf(beta)
        if abs(b - mpf(1) / 2) < mpf(10) ** -6:  # s = 1, the pole of zeta
            b = b + mpf(10) ** -6
        values.append(float(log(xi(mpf(1) / 2 + b)) - reference))
    return values


def euler_profile(p: int, truncation: int, betas: list[float]) -> list[float]:
    """``log Z`` of ``{1, p, ..., p^K}``, summed from the top to avoid overflow."""

    log_p = math.log(p)
    return [
        truncation * beta * log_p
        + math.log(sum(math.exp(-i * beta * log_p) for i in range(truncation + 1)))
        for beta in betas
    ]


def windowed_rates(source: list[float], target: list[float]) -> tuple[float, float]:
    ratios = [s / t for s, t in zip(source, target)]
    return min(ratios), 1.0 / max(ratios)


def main() -> int:
    rows = []
    print("Unrestricted rates between xi and any truncated Euler factor: both 0.")
    print("  b -> 0   : log Z_xi ~ 0.0231 b^2 -> 0 while log Z_P(0) = log(K+1) > 0")
    print("  b -> inf : log Z_xi ~ (b/2) log b grows faster than log Z_P ~ K b log p")
    print("  so inf over [0, inf] of each ratio is 0, in both directions.\n")

    for low, high in WINDOWS:
        betas = [low + (high - low) * index / (GRID - 1) for index in range(GRID)]
        xi_values = xi_profile(betas)
        print(f"window beta in [{low}, {high}]")
        print(
            f"  {'p':>3} {'K':>3} {'C(xi->P)':>11} {'C(P->xi)':>11} "
            f"{'product':>10} {'d':>9}"
        )
        for p in PRIMES:
            for truncation in TRUNCATIONS:
                target = euler_profile(p, truncation, betas)
                forward, backward = windowed_rates(xi_values, target)
                product = forward * backward
                distance = -math.log(product) if product > 0 else math.inf
                rows.append(
                    [
                        low,
                        high,
                        p,
                        truncation,
                        f"{forward:.12f}",
                        f"{backward:.12f}",
                        f"{product:.12f}",
                        f"{distance:.12f}",
                    ]
                )
                if truncation in (1, 8, 128):
                    print(
                        f"  {p:>3} {truncation:>3} {forward:>11.6f} {backward:>11.6f} "
                        f"{product:>10.6f} {distance:>9.6f}"
                    )
        # As K -> infinity, log Z_P -> K b log p, and K and p cancel from the
        # product of the two rates.  What survives is the departure of xi's own
        # profile from a straight line through the origin.
        slopes = [value / beta for value, beta in zip(xi_values, betas)]
        limit_product = min(slopes) / max(slopes)
        print(
            f"  K -> infinity limit, independent of p and K: "
            f"product={limit_product:.6f}, d={-math.log(limit_product):.6f}\n"
        )

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["window_low", "window_high", "p", "K", "C_xi_to_P", "C_P_to_xi", "product", "d"]
        )
        writer.writerows(rows)
    print(f"written to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
