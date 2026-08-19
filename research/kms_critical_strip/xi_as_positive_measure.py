#!/usr/bin/env python3
"""What the completion actually does to the resource, and why it is not a
semiring operation.

The two-positivities note says the completed zeta "repairs zeta as a partition
function": ``xi(1/2 + b) = int Phi(u) e^{bu} du`` with ``Phi > 0`` even and
doubly exponentially decaying, so ``log Z_xi`` is finite and convex for every
real ``b``.  That is true, and it is verified here at 40 digits.  What the note
does not say is *how* the repair works, and the answer decides whether the KMS
route can reach the critical strip from inside the framework.

Riemann's measure is a DIFFERENCE of two positive measures,

    Phi(u) = A(u) - B(u),
    A(u) = 4 pi^2 e^{9u/2} sum_n n^4 exp(-pi n^2 e^{2u}),
    B(u) = 6 pi   e^{5u/2} sum_n n^2 exp(-pi n^2 e^{2u}),

and their two-sided Laplace transforms are, in closed form,

    A^(b) = 2 pi^{-1/4-b/2} Gamma(9/4 + b/2) zeta(b + 1/2),
    B^(b) = 3 pi^{-1/4-b/2} Gamma(5/4 + b/2) zeta(b + 1/2),

so that with ``s = b + 1/2``

    A^/B^ = (2/3)(5/4 + b/2) = (s + 2)/3,        Phi^ = B^ (s-1)/3 = xi(s).

Consequences, each checked below:

  * A and B are honest positive resources whose partition functions have
    abscissa of convergence exactly ``s = 1`` -- the Hagedorn point of the
    primon gas.  Both diverge on the whole critical strip.
  * Their PRICE ratio ``A^/B^ = (s+2)/3`` equals 1 exactly at ``s = 1``.  The
    completion is the cancellation of two divergent resources that become
    equal in price precisely at the phase transition.
  * On the measure side, multiplying a transform by ``(s-1) = (b - 1/2)``
    is the first-order differential operator ``f -> -f' - f/2``.  It carries a
    minus sign, and the exchange semiring has ``(x)`` and ``(+)`` but no
    subtraction and no differentiation.  So ``xi`` is a legitimate resource,
    but it is NOT reachable from ``(x)_p P_{p,inf}`` by any operation of the
    framework.

Also verified: ``(1/2)(log Z_xi)''(0) = sum_{gamma>0} gamma^{-2}``, the note's
``0.0231``, computed independently from 1200 zeta zeros with a tail estimate.
"""

from __future__ import annotations

import csv
from pathlib import Path

import mpmath as mp
import numpy as np

OUT = Path(__file__).resolve().parent
ZEROS = OUT.parent / "m_and_e_and_a_c" / "zeta_zeros_1200.npy"
mp.mp.dps = 40


def theta_sum(power: int, x):
    """sum_{n>=1} n^{power} exp(-pi n^2 x), power even."""
    total = mp.mpf(0)
    n = 1
    while True:
        term = mp.mpf(n) ** power * mp.e ** (-mp.pi * n * n * x)
        total += term
        if n > 4 and term < mp.mpf(10) ** (-mp.mp.dps - 10) * max(total, mp.mpf(1)):
            break
        n += 1
        if n > 200000:
            break
    return total


def part_a(u):
    x = mp.e ** (2 * u)
    return 4 * mp.pi**2 * mp.e ** (9 * u / 2) * theta_sum(4, x)


def part_b(u):
    x = mp.e ** (2 * u)
    return 6 * mp.pi * mp.e ** (5 * u / 2) * theta_sum(2, x)


def phi(u):
    return part_a(u) - part_b(u)


def s_minus_one_zeta(s):
    if abs(s - 1) < mp.mpf(10) ** -20:
        return mp.mpf(1) + (s - 1) * mp.euler
    return (s - 1) * mp.zeta(s)


def xi(s):
    return mp.mpf(1) / 2 * s * mp.pi ** (-s / 2) * mp.gamma(s / 2) * s_minus_one_zeta(s)


def main() -> int:
    rows = []

    print("=" * 78)
    print("1. Phi is positive and even, and its two-sided Laplace transform is")
    print("   xi(1/2 + b).")
    print("=" * 78)
    print(f"{'u':>7} {'Phi(u)':>26} {'Phi(-u)-Phi(u)':>16}")
    for u_value in ("0", "0.25", "0.5", "0.75", "1", "1.5"):
        u = mp.mpf(u_value)
        value = phi(u)
        print(f"{u_value:>7} {mp.nstr(value, 20):>26} {mp.nstr(phi(-u) - value, 4):>16}")
        rows.append(["phi", u_value, mp.nstr(value, 25), mp.nstr(phi(-u) - value, 6), ""])
    print("   Phi > 0 at every sample point; evenness is Riemann's functional-")
    print("   equation form.  Note that A - B is a catastrophic cancellation for")
    print("   u < -1 -- Phi there is doubly exponentially small -- so the")
    print("   transform below is evaluated as 2 int_0^inf Phi(u) cosh(bu) du,")
    print("   which uses evenness and only touches the stable half-line.")

    print()
    print(f"{'b':>6} {'int Phi e^{bu} du':>28} {'xi(1/2+b)':>28} {'rel.err':>10}")
    for b_value in ("0", "0.5", "1", "2", "5"):
        b = mp.mpf(b_value)
        integral = 2 * mp.quad(
            lambda u: phi(u) * mp.cosh(b * u), [0, mp.mpf("0.5"), 1, 2, 3]
        )
        target = xi(mp.mpf(1) / 2 + b)
        print(
            f"{b_value:>6} {mp.nstr(integral, 22):>28} {mp.nstr(target, 22):>28} "
            f"{mp.nstr(abs(integral - target) / abs(target), 4):>10}"
        )
        rows.append(
            ["transform", b_value, mp.nstr(integral, 25), mp.nstr(target, 25),
             mp.nstr(abs(integral - target) / abs(target), 6)]
        )

    print()
    print("=" * 78)
    print("2. A and B separately: closed-form transforms, and abscissa s = 1.")
    print("   A^(b) = 2 pi^{-1/4-b/2} Gamma(9/4+b/2) zeta(b+1/2)")
    print("   B^(b) = 3 pi^{-1/4-b/2} Gamma(5/4+b/2) zeta(b+1/2)")
    print("=" * 78)

    def a_hat(b):
        return 2 * mp.pi ** (-mp.mpf(1) / 4 - b / 2) * mp.gamma(
            mp.mpf(9) / 4 + b / 2
        ) * mp.zeta(b + mp.mpf(1) / 2)

    def b_hat(b):
        return 3 * mp.pi ** (-mp.mpf(1) / 4 - b / 2) * mp.gamma(
            mp.mpf(5) / 4 + b / 2
        ) * mp.zeta(b + mp.mpf(1) / 2)

    print("   Quadrature is taken over u in [-3, 4] plus the analytic left tail")
    print("   int_{-inf}^{-3} (3/2) e^{(b-1/2)u} du = (3/2) e^{-3(b-1/2)}/(b-1/2),")
    print("   whose error is doubly exponentially small.  The tail formula is")
    print("   itself the statement that the abscissa of convergence is b = 1/2.")
    print(f"{'b':>6} {'s':>6} {'quad A':>24} {'closed A^':>24} {'rel.err':>10}")
    for b_value in ("0.75", "1", "2", "4"):
        b = mp.mpf(b_value)
        integral = mp.quad(
            lambda u: part_a(u) * mp.e ** (b * u), [-3, -1, 0, 1, 4]
        ) + mp.mpf(3) / 2 * mp.e ** (-3 * (b - mp.mpf(1) / 2)) / (b - mp.mpf(1) / 2)
        closed = a_hat(b)
        print(
            f"{b_value:>6} {mp.nstr(b + mp.mpf(1)/2, 4):>6} {mp.nstr(integral, 18):>24} "
            f"{mp.nstr(closed, 18):>24} "
            f"{mp.nstr(abs(integral - closed) / abs(closed), 4):>10}"
        )
        rows.append(["a_hat", b_value, mp.nstr(integral, 25), mp.nstr(closed, 25), ""])

    print()
    print("   Divergence of A and B on the strip: as u -> -infinity both behave")
    print("   like (3/2) e^{-u/2}, so int A e^{bu} du converges iff b > 1/2,")
    print("   i.e. iff s > 1.  Measured A(u) e^{u/2} against the predicted 3/2:")
    print(f"{'u':>8} {'A(u) e^{u/2}':>22} {'B(u) e^{u/2}':>22} "
          f"{'(A-B) e^{u/2}':>22}")
    for u_value in ("-1", "-1.5", "-2", "-2.5"):
        u = mp.mpf(u_value)
        scale = mp.e ** (u / 2)
        va, vb = part_a(u) * scale, part_b(u) * scale
        print(
            f"{u_value:>8} {mp.nstr(va, 16):>22} {mp.nstr(vb, 16):>22} "
            f"{mp.nstr(va - vb, 6):>22}"
        )
        rows.append(["tail", u_value, mp.nstr(va, 25), mp.nstr(vb, 25),
                     mp.nstr(va - vb, 10)])
    print("   Both tend to 3/2 = 1.5 and the difference to 0 doubly")
    print("   exponentially: the leading divergent parts cancel exactly.")

    print()
    print("   Price ratio A^/B^ = (s+2)/3, equal to 1 exactly at s = 1:")
    print(f"{'s':>7} {'A^/B^':>22} {'(s+2)/3':>22} {'diff':>10}")
    for s_value in ("1.05", "1.5", "2", "4"):
        s = mp.mpf(s_value)
        b = s - mp.mpf(1) / 2
        ratio = a_hat(b) / b_hat(b)
        pred = (s + 2) / 3
        print(
            f"{s_value:>7} {mp.nstr(ratio, 18):>22} {mp.nstr(pred, 18):>22} "
            f"{mp.nstr(abs(ratio - pred), 3):>10}"
        )
        rows.append(["price_ratio", s_value, mp.nstr(ratio, 25), mp.nstr(pred, 25), ""])
    print("   At s = 1 the ratio is exactly 1: the two divergent resources have")
    print("   equal price precisely at the Hagedorn temperature, and their")
    print("   difference (s-1)/3 * B^ is entire.")

    print()
    print("=" * 78)
    print("3. (1/2)(log Z_xi)''(0) = sum_{gamma>0} gamma^{-2}, the note's 0.0231.")
    print("=" * 78)
    zeros = np.load(ZEROS)
    partial = sum(mp.mpf(float(g)) ** -2 for g in zeros)
    last = mp.mpf(float(zeros[-1]))
    # Riemann-von Mangoldt density dN = (1/2pi) log(t/2pi) dt gives the tail
    tail = mp.quad(
        lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi * t**2), [last, mp.inf]
    )

    def log_z_xi(b):
        return mp.log(xi(mp.mpf(1) / 2 + b) / xi(mp.mpf(1) / 2))

    second = mp.diff(log_z_xi, mp.mpf(0), 2)
    print(f"  1200 zeros            {mp.nstr(partial, 20)}")
    print(f"  tail estimate         {mp.nstr(tail, 8)}")
    print(f"  total                 {mp.nstr(partial + tail, 20)}")
    print(f"  (1/2)(log Z_xi)''(0)  {mp.nstr(second / 2, 20)}")
    print(f"  difference            {mp.nstr(abs(partial + tail - second / 2), 6)}")
    print(f"  classical sum_rho 1/rho = 1 + gamma_E/2 - (1/2)log(4 pi) = "
          f"{mp.nstr(1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2, 20)}")
    print("  (the two constants differ in the 5th digit and must not be")
    print("   confused: sum 1/gamma^2 versus sum 1/(1/4 + gamma^2).)")
    rows.append(["zero_sum", "", mp.nstr(partial + tail, 25), mp.nstr(second / 2, 25), ""])

    with (OUT / "xi_as_positive_measure.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["kind", "argument", "value1", "value2", "value3"])
        writer.writerows(rows)
    print(f"\nwritten to {(OUT / 'xi_as_positive_measure.csv').name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
