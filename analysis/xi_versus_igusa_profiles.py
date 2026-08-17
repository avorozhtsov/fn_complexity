#!/usr/bin/env python3
"""Exchange rates between the completed zeta function and p-adic Igusa profiles.

An Igusa profile is the cost curve of a p-adic map,

    Z_{f,p}(s) = integral over Z_p^n of |f(x)|_p^s dx,     Phi_{f,p} = -log Z_{f,p},

positive and increasing, with ``Phi(0) = 0`` because the measure is normalised.
Its two ends are honest invariants of the map:

    Phi'(0)      = log p * (average valuation of f),   the profile is LINEAR at 0,
    Phi(infty)   = -log measure{ |f|_p = 1 },          the profile SATURATES.

The completed zeta profile ``log Z_xi(b) = log(xi(1/2+b)/xi(1/2))`` has neither
shape: it is quadratic at the origin and grows like ``(b/2) log b``.  So both
unrestricted rates vanish again, but for a new pair of reasons -- quadratic
against linear at ``s = 0``, divergent against bounded at ``s = infinity``.
Against the truncated Euler factors the failure was quadratic-against-positive
and slower-than-linear; here it is the opposite end that saturates.

On a window the rates are finite.  Taking ``r`` independent copies of a map
multiplies its profile by ``r``, and ``r`` cancels from the product of the two
rates, so the irreversibility is a pure shape mismatch, exactly as the
truncation ``K`` cancelled for the Euler factors.

Profiles are the closed forms of the local-currency note; the script first
reproduces its constant ``C_zeta(x^2 -> x^2-y^2) = 0.9397027875459163...`` as a
check.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

from mpmath import gamma, log, mp, mpf, pi, zeta

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "xi_versus_igusa_profiles.csv"

mp.dps = 15
WINDOWS = [(0.5, 5.0), (1.0, 10.0), (2.0, 20.0)]
GRID = 240


def xi(s):
    return mpf(1) / 2 * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def xi_profile(betas):
    reference = log(xi(mpf(1) / 2))
    values = []
    for beta in betas:
        b = mpf(beta)
        if abs(b - mpf(1) / 2) < mpf(10) ** -6:
            b = b + mpf(10) ** -6
        values.append(float(log(xi(mpf(1) / 2 + b)) - reference))
    return values


def phi_identity(p, s):
    """f(x) = x on Z_p:  Z = (1 - 1/p)/(1 - p^{-s-1})."""

    return -math.log((1 - 1 / p) / (1 - p ** (-s - 1)))


def phi_square_2(s):
    """f(x,y) = x^2 on Z_2^2:  Z = 1/(2 - t^2), t = 2^{-s}."""

    t = 2.0 ** (-s)
    return math.log(2 - t * t)


def phi_split_2(s):
    """f(x,y) = x^2 - y^2 on Z_2^2:  Z = (t^2 - 2t + 2)/(2 - t)^2."""

    t = 2.0 ** (-s)
    return 2 * math.log(2 - t) - math.log(t * t - 2 * t + 2)


PROFILES = {
    "x on Z_2": lambda s: phi_identity(2, s),
    "x on Z_3": lambda s: phi_identity(3, s),
    "x on Z_5": lambda s: phi_identity(5, s),
    "x on Z_7": lambda s: phi_identity(7, s),
    "x^2 on Z_2 (degenerate)": phi_square_2,
    "x^2-y^2 on Z_2 (split, bad prime)": phi_split_2,
    "xy on Z_3 (= x^2-y^2)": lambda s: 2 * phi_identity(3, s),
    "x^3 on Z_2": lambda s: phi_identity(2, 3 * s),
    "x^2 y on Z_2": lambda s: phi_identity(2, 2 * s) + phi_identity(2, s),
}


def check_documented_constant() -> None:
    objective = lambda s: phi_square_2(s) / phi_split_2(s)  # noqa: E731
    low, high = 0.1, 3.0
    for _ in range(300):
        first = low + (high - low) / 3
        second = high - (high - low) / 3
        if objective(first) < objective(second):
            high = second
        else:
            low = first
    argument = (low + high) / 2
    print("check against the local-currency note:")
    print(f"   C_zeta(x^2 -> x^2-y^2) = {objective(argument):.15f}")
    print("   documented             = 0.939702787545916")
    print(f"   minimiser t* = {2.0 ** (-argument):.15f}   documented 0.696541929482172\n")


def endpoints() -> None:
    print("the two ends of each Igusa profile (both finite -- this is the obstruction):")
    print(f"   {'profile':<36} {'Phi(0.001)/0.001':>17} {'Phi(60)':>10}")
    for name, phi in PROFILES.items():
        print(f"   {name:<36} {phi(0.001) / 0.001:>17.6f} {phi(60.0):>10.6f}")
    print("   xi, for contrast: quadratic at 0, and unbounded:")
    reference = float(log(xi(mpf(1) / 2)))
    small = float(log(xi(mpf(1) / 2 + mpf("0.001")))) - reference
    print(f"   {'log Z_xi':<36} {small / 0.001:>17.6f} {'divergent':>10}\n")


def main() -> int:
    check_documented_constant()
    endpoints()

    print("Unrestricted rates against xi: both vanish, for two new reasons.")
    print("   s -> 0    : Phi_f is linear, log Z_xi is quadratic  => C(xi->f) = 0")
    print("   s -> inf  : Phi_f saturates, log Z_xi diverges      => C(f->xi) = 0\n")

    rows = []
    for low, high in WINDOWS:
        betas = [low + (high - low) * index / (GRID - 1) for index in range(GRID)]
        xi_values = xi_profile(betas)
        print(f"window s in [{low}, {high}]")
        print(f"   {'profile':<36} {'C(xi->f)':>11} {'C(f->xi)':>11} {'product':>10} {'d':>9}")
        results = []
        for name, phi in PROFILES.items():
            target = [phi(beta) for beta in betas]
            ratios = [x / t for x, t in zip(xi_values, target)]
            forward, backward = min(ratios), 1.0 / max(ratios)
            product = forward * backward
            distance = -math.log(product)
            results.append((distance, name, forward, backward, product))
            rows.append(
                [
                    low,
                    high,
                    name,
                    f"{forward:.12f}",
                    f"{backward:.12f}",
                    f"{product:.12f}",
                    f"{distance:.12f}",
                ]
            )
            print(
                f"   {name:<36} {forward:>11.6f} {backward:>11.6f} "
                f"{product:>10.6f} {distance:>9.6f}"
            )
        best = min(results)
        print(f"   closest: {best[1]}, d={best[0]:.6f}")
        # r independent copies multiply the profile by r, which cancels from the
        # product; the irreversibility is a pure shape mismatch.
        target = [PROFILES["x on Z_2"](beta) for beta in betas]
        for copies in (1, 5, 50):
            ratios = [x / (copies * t) for x, t in zip(xi_values, target)]
            product = min(ratios) / max(ratios)
            print(f"   r={copies:>3} copies of 'x on Z_2': d={-math.log(product):.6f}")
        print()

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["window_low", "window_high", "profile", "C_xi_to_f", "C_f_to_xi", "product", "d"]
        )
        writer.writerows(rows)
    print(f"written to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
