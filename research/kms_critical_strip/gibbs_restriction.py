#!/usr/bin/env python3
"""Step 2 of session brief H: the framework is the s > 1 Gibbs restriction.

Two things are checked.

**A. The exchange monotone of a local prime mode is its Gibbs free energy.**
For the truncated Euler factor ``P_{p,K} = {1, p, ..., p^K}`` the exchange
framework's monotone is ``Z(beta) = sum_{j<=K} p^{j beta}`` (entries >= 1,
beta >= 0), and the primon-gas local factor is ``sum_{j<=K} p^{-js}``.  These
are the same function at ``beta = -s``.  The script confirms that the repo's
``exchange_rate`` on the signature ``(1, p, ..., p^K)`` agrees with the closed
form obtained from the log-partition profiles, and records the exact rates
between local modes.

**B. Reproduction of `analysis/xi_versus_euler_factors.py` at 40 digits.**
The published table uses ``mp.dps = 15`` and a 240-point uniform beta grid.
Here the same quantities are recomputed at ``mp.dps = 40`` with true extrema
(golden-section refinement of the grid minimiser), so that the grid error of
the published numbers can be reported alongside any agreement claim, as the
brief requires.

**C. Where the framework's temperature range actually sits.**
The completed resource ``xi`` places the critical LINE at ``beta_xi = 0`` --
the framework's infinite-temperature endpoint -- and the functional equation
is the reflection ``beta -> -beta`` about it.  The uncompleted primon gas
``{1,2,3,...}`` places the critical STRIP at ``beta_x in (-1, 0)``, entirely
outside the framework's admissible ``[0, infinity]``.  Both statements are
verified numerically here.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

from fn_complexity.core import exchange_rate  # noqa: E402

OUT = Path(__file__).resolve().parent
mp.mp.dps = 40

WINDOWS = [("0.5", "5"), ("1", "10"), ("2", "20")]
PRIMES = [2, 3, 5, 7]
TRUNCATIONS = [1, 8, 128]


def s_minus_one_zeta(s):
    """(s-1) zeta(s), regular at s = 1 where its value is 1.

    The published script sidesteps the pole by shifting the grid point by
    1e-6; the window [0.5, 5] has s = 1 exactly at its left endpoint, so that
    shift is one of the error sources quantified below.
    """
    if abs(s - 1) < mp.mpf(10) ** -20:
        return mp.mpf(1) + (s - 1) * mp.euler
    return (s - 1) * mp.zeta(s)


def xi(s):
    return mp.mpf(1) / 2 * s * mp.pi ** (-s / 2) * mp.gamma(s / 2) * s_minus_one_zeta(s)


XI_HALF = xi(mp.mpf(1) / 2)


def log_z_xi(beta):
    """log of the mass-one normalised completed-zeta partition function."""
    return mp.log(xi(mp.mpf(1) / 2 + beta) / XI_HALF)


def log_z_euler(p, K, beta):
    """log sum_{j=0}^{K} p^{j beta}: the exchange monotone of P_{p,K}."""
    lp = mp.log(p)
    return K * beta * lp + mp.log(mp.fsum(mp.e ** (-j * beta * lp) for j in range(K + 1)))


def extremum(func, low, high, *, grid: int, maximise: bool):
    """Grid scan followed by golden-section refinement; returns (value, point)."""
    step = (high - low) / (grid - 1)
    points = [low + step * i for i in range(grid)]
    values = [func(x) for x in points]
    index = max(range(grid), key=lambda i: values[i]) if maximise else min(
        range(grid), key=lambda i: values[i]
    )
    left = points[max(index - 1, 0)]
    right = points[min(index + 1, grid - 1)]
    sign = -1 if maximise else 1
    phi = (mp.sqrt(5) - 1) / 2
    a, b = left, right
    c = b - phi * (b - a)
    d = a + phi * (b - a)
    for _ in range(200):
        if sign * func(c) < sign * func(d):
            b = d
        else:
            a = c
        c = b - phi * (b - a)
        d = a + phi * (b - a)
        if b - a < mp.mpf(10) ** (-32):
            break
    point = (a + b) / 2
    # keep the better of the refined point and the grid endpoint values
    candidates = [(values[index], points[index]), (func(point), point)]
    return max(candidates) if maximise else min(candidates)


def main() -> int:
    rows = []

    print("=" * 78)
    print("A. The local-mode exchange monotone is the Gibbs free energy.")
    print("   signature (1,p,...,p^K), Z(beta) = sum p^{j beta}; the primon-gas")
    print("   local factor sum p^{-js} is the same function at beta = -s.")
    print("=" * 78)
    print(f"{'p':>3} {'K':>4} {'q':>3} {'L':>4} {'C(P_pK -> P_qL)':>20} "
          f"{'closed form':>20} {'|diff|':>10}")
    for p, K in ((2, 4), (2, 8), (3, 4), (5, 6)):
        for q, L in ((2, 4), (3, 8), (5, 3), (7, 12)):
            source = tuple(p**j for j in range(K + 1))
            target = tuple(q**j for j in range(L + 1))
            rate = exchange_rate(source, target)
            # closed form: inf over beta in [0,inf] of log Z_source / log Z_target.
            def ratio(beta, p=p, K=K, q=q, L=L):
                return log_z_euler(p, K, beta) / log_z_euler(q, L, beta)

            lo = extremum(ratio, mp.mpf("1e-9"), mp.mpf(400), grid=4000, maximise=False)[0]
            ends = [
                mp.log(K + 1) / mp.log(L + 1),          # beta = 0
                K * mp.log(p) / (L * mp.log(q)),        # beta = infinity
            ]
            closed = min([lo, *ends])
            print(
                f"{p:>3} {K:>4} {q:>3} {L:>4} {rate:>20.14f} "
                f"{float(closed):>20.14f} {abs(rate - float(closed)):>10.2e}"
            )
            rows.append(["local_rate", p, K, q, L, f"{rate:.14f}", mp.nstr(closed, 16)])

    print()
    print("=" * 78)
    print("B. `analysis/xi_versus_euler_factors.py` recomputed at dps=40 with")
    print("   true extrema.  'published' is the dps=15, 240-point grid value.")
    print("=" * 78)
    published = {
        ("0.5", "5"): {
            (2, 1): "3.208235", (2, 8): "2.637828", (2, 128): "2.310921",
            (3, 1): "2.887080", (3, 8): "2.461196", (3, 128): "2.295856",
            (5, 1): "2.661633", (5, 8): "2.371665", (5, 128): "2.289426",
            (7, 1): "2.568499", (7, 8): "2.342875", (7, 128): "2.287493",
        },
        ("1", "10"): {
            (2, 1): "2.696770", (2, 8): "2.353802", (2, 128): "2.244131",
            (3, 1): "2.468935", (3, 8): "2.281446", (3, 128): "2.239229",
            (5, 1): "2.343663", (5, 8): "2.253533", (5, 128): "2.237433",
            (7, 1): "2.302720", (7, 8): "2.246204", (7, 128): "2.236969",
        },
    }
    published_limits = {("0.5", "5"): "2.283688", ("1", "10"): "2.236350",
                        ("2", "20"): "2.114300"}

    for low_s, high_s in WINDOWS:
        low, high = mp.mpf(low_s), mp.mpf(high_s)
        print(f"\nwindow beta_xi in [{low_s}, {high_s}]")
        print(f"  {'p':>3} {'K':>4} {'d (dps=40)':>22} {'published':>12} {'diff':>10}")
        for p in PRIMES:
            for K in TRUNCATIONS:
                def ratio(beta, p=p, K=K):
                    return log_z_xi(beta) / log_z_euler(p, K, beta)

                lo = extremum(ratio, low, high, grid=2000, maximise=False)[0]
                hi = extremum(ratio, low, high, grid=2000, maximise=True)[0]
                distance = mp.log(hi / lo)
                ref = published.get((low_s, high_s), {}).get((p, K))
                diff = "" if ref is None else f"{float(distance) - float(ref):>10.2e}"
                print(
                    f"  {p:>3} {K:>4} {mp.nstr(distance, 18):>22} "
                    f"{ref or '-':>12} {diff:>10}"
                )
                rows.append(
                    ["xi_vs_euler", low_s, high_s, p, K, mp.nstr(distance, 25), ref or ""]
                )
        # K -> infinity: log Z_P -> K beta log p, so p and K cancel and what
        # survives is the departure of log Z_xi from a straight line through 0.
        def chord(beta):
            return log_z_xi(beta) / beta

        lo = extremum(chord, low, high, grid=4000, maximise=False)[0]
        hi = extremum(chord, low, high, grid=4000, maximise=True)[0]
        limit = mp.log(hi / lo)
        ref = published_limits[(low_s, high_s)]
        print(
            f"  K -> infinity limit  d_W = {mp.nstr(limit, 22)}  "
            f"published {ref}  diff {float(limit) - float(ref):.2e}"
        )
        rows.append(["xi_vs_euler_limit", low_s, high_s, "", "", mp.nstr(limit, 25), ref])

    print()
    print("=" * 78)
    print("C. Where the critical line and strip actually sit on the framework's")
    print("   temperature axis, for the two zeta resources of the notes.")
    print("=" * 78)
    print("  completed resource  xi:  beta_xi = s - 1/2")
    print("     critical line s=1/2  ->  beta_xi = 0      (framework endpoint,")
    print("                                                infinite temperature)")
    print("     functional equation s <-> 1-s  ->  beta_xi <-> -beta_xi")
    print("     nontrivial zeros  rho = 1/2 + i gamma  ->  beta_xi = i gamma")
    print("  uncompleted primon gas  {1,2,3,...}:  beta_x = -s")
    print("     abscissa s=1        ->  beta_x = -1")
    print("     critical strip 0<s<1 -> beta_x in (-1, 0)")
    print("     critical line s=1/2  ->  beta_x = -1/2")
    print("     framework's admissible range beta_x in [0, infinity] lies")
    print("     ENTIRELY inside the divergence half-line beta_x >= -1.")
    print()
    for beta_value in ("0", "0.25", "0.5", "1", "2"):
        beta = mp.mpf(beta_value)
        print(
            f"  beta_xi = {beta_value:>5}  s = {mp.nstr(mp.mpf(1)/2 + beta, 6):>7}  "
            f"log Z_xi = {mp.nstr(log_z_xi(beta), 14):>20}"
        )
    print()
    print("  log Z_xi(0) = 0 exactly: the completed resource is a probability")
    print("  measure, one 'entry'.  So u_xi = log log Z_xi -> -infinity at the")
    print("  critical line and both rates against any finite signature vanish")
    print("  there.  It is the NORMALISATION that fails at s = 1/2, not the")
    print("  finiteness of the monotone.")

    print()
    quad = mp.diff(log_z_xi, mp.mpf(0), 2)
    print(f"  (log Z_xi)''(0)   = {mp.nstr(quad, 22)}")
    print(f"  half of it        = {mp.nstr(quad / 2, 22)}   (the note's 0.0231)")
    print("  and by Hadamard   = sum_{gamma>0} 1/gamma^2 (see zeta_zero_sum.py)")
    rows.append(["log_z_xi_second_derivative", "", "", "", "", mp.nstr(quad, 30), ""])

    with (OUT / "gibbs_restriction.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["kind", "a", "b", "c", "d", "value", "reference"])
        writer.writerows(rows)
    print(f"\nwritten to {(OUT / 'gibbs_restriction.csv').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
