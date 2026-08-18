#!/usr/bin/env python3
"""Brief M item 4, continued -- what replaces (log 2)/2 at finite lambda?

Brief G proves that around any directed cycle of the tropical comparison the
mean edge weight is at most ``(log 2)/2``: *the geometric mean of the rate
asymmetry around any preference cycle is at most 2*.  Here that constant is
measured for ``A_lambda``.

The test bed is the **flat locus** -- signatures with all fibers equal -- for
three reasons: brief G proves the tropical curl is *exactly* zero there, the
``lambda = 0`` curl is exactly zero by the potential identity, and a flat
signature is a two-parameter object ``u(s) = log(R + e^s Lambda)``, so the
search is over three ``sigma = log(R/Lambda)`` and nothing else (an overall
scaling of ``(R, Lambda)`` is a Cartesian power and shifts ``u`` by a constant,
which cancels in every ``A_lambda``).

``logaddexp`` is used throughout: the naive ``log(1 + e^{s-sigma})`` overflows
for the ``sigma`` spreads that matter and silently returns garbage.

Writes ``m4b_cycle_strength.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import common as C

W = 640.0  # prior half-width; the constant has saturated well before this
NGRID = 12001


def uflat(s, sigma):
    """``u`` of a flat signature with ``R = 1``, ``Lambda = e^{-sigma}``."""

    return np.logaddexp(0.0, s - sigma)


def main():
    rows, out = [], []

    def say(t=""):
        print(t)
        out.append(t)

    s, w = C.uniform_prior(-W, W, NGRID)

    def curl(sig3, lam):
        u = [uflat(s, t) for t in sig3]
        return sum(C.soft_mid(u[i] - u[(i + 1) % 3], w, lam) for i in range(3))

    def search(lam, coarse=29):
        V = np.linspace(-1.4 * W, 1.4 * W, coarse)
        cache = {float(t): uflat(s, float(t)) for t in V}
        best, arg = 0.0, None
        for a in V:
            for b in V:
                f1 = C.soft_mid(cache[float(a)] - cache[float(b)], w, lam)
                for c in V:
                    z = abs(f1 + C.soft_mid(cache[float(b)] - cache[float(c)], w, lam)
                            + C.soft_mid(cache[float(c)] - cache[float(a)], w, lam))
                    if z > best:
                        best, arg = z, [float(a), float(b), float(c)]
        step = 0.35 * W
        while step > 1e-6:
            moved = False
            for i in range(3):
                for d in (+step, -step):
                    y = list(arg)
                    y[i] += d
                    z = abs(curl(y, lam))
                    if z > best + 1e-15:
                        best, arg, moved = z, y, True
            if not moved:
                step *= 0.5
        return best, arg

    say("=" * 88)
    say(f"Sharp cycle strength on the FLAT locus, prior = uniform on [-{W:g}, {W:g}]")
    say("  brief G: mean_e |A_inf| around a directed cycle <= (log 2)/2 = 0.3465736,")
    say("           and on the flat locus the tropical curl is exactly 0.")
    say("=" * 88)
    say(f"  {'lambda':>10s} {'max mean|A_lam|':>17s} {'/ (log2)/2':>12s}   sigma triple")
    for lam in (0.03, 0.1, 0.3, 1.0, 2.0, 3.0, 10.0, 30.0, 100.0, 1000.0, None):
        best, arg = search(lam)
        tag = "inf" if lam is None else f"{lam:g}"
        say(f"  {tag:>10s} {best/3:17.9f} {best/3/(math.log(2)/2):12.5f}   "
            + "(" + ", ".join(f"{t:.3f}" for t in arg) + ")")
        rows.append(dict(block="flat_sharp", lam=tag, mean_abs_A=best / 3,
                         over_log2half=best / 3 / (math.log(2) / 2), sigmas=str(arg)))

    say()
    say("  saturation in the prior width at lambda = 1 (the constant density cancels")
    say("  between softmax and softmin, so widening the prior stops mattering):")
    say(f"  {'half-width':>12s} {'max mean|A_1|':>16s}")
    for hw in (10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 1280.0):
        ss, ww = C.uniform_prior(-hw, hw, 8001)

        def curl2(t3):
            u = [np.logaddexp(0.0, ss - t) for t in t3]
            return abs(sum(C.soft_mid(u[i] - u[(i + 1) % 3], ww, 1.0) for i in range(3)))

        V = np.linspace(-1.4 * hw, 1.4 * hw, 29)
        best, arg = 0.0, None
        for a in V:
            for b in V:
                for c in V:
                    z = curl2([a, b, c])
                    if z > best:
                        best, arg = z, [float(a), float(b), float(c)]
        step = 0.35 * hw
        while step > 1e-6:
            moved = False
            for i in range(3):
                for d in (+step, -step):
                    y = list(arg)
                    y[i] += d
                    z = curl2(y)
                    if z > best + 1e-15:
                        best, arg, moved = z, y, True
            if not moved:
                step *= 0.5
        say(f"  {hw:12.0f} {best/3:16.9f}")
        rows.append(dict(block="saturation", half_width=hw, mean_abs_A=best / 3))

    with open("m4b_cycle_strength.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    say("\nwrote m4b_cycle_strength.csv")


if __name__ == "__main__":
    main()
