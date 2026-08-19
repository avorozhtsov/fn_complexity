#!/usr/bin/env python3
"""Explicit witnesses for the upper bounds of Theorems 2 and 3, re-verified.

Everything here is checked by evaluating the actual functions on every point of
``F_q^{2k}`` -- no jet algebra, no search.

Theorem 2 (anisotropic resource, hyperbolic target).  ``Q = x^2 - n y^2``.
With ``u_i = x_i + y_i`` and ``v_i = x_i - y_i``,

    h_1 = Q(u_1, u_2),  h_2 = Q(v_1, u_2),  h_3 = Q(u_1, v_2)
    =>  4 x_1 y_1 = h_1 - h_2 ,   -4 n x_2 y_2 = h_1 - h_3 ,

so ``N_2(Q -> xy) <= 3``; and ``x y = (Q(u,0) - Q(v,0)) / 4`` gives
``N_1(Q -> xy) <= 2``.  Together ``N_k <= ceil(3k/2)``.

Theorem 3 (parabolic resource ``g = x^2 + y``).  ``g(l, 0) = l^2`` and
``g(0, m) = m``, so a rank-``rho`` quadratic part costs ``rho`` atoms per copy:

    x y      = ( g(x+y, 0) - g(x-y, 0) ) / 4          (rho = 2)
    x^2+y^2  =   g(x, 0) + g(y, 0)                    (rho = 2)
    x^2+y    =   g(x, y)                              (rho = 1)
    x^2      =   g(x, 0)                              (rho = 1)
    x        =   g(0, x)                              (rho = 0)
"""

from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from isotropy_atoms import least_non_residue  # noqa: E402

PRIMES = (3, 5, 7, 11, 13, 17, 19)


def inverse(value, q):
    return pow(value % q, q - 2, q)


def check_theorem2(q: int) -> tuple[bool, bool]:
    n = least_non_residue(q)

    def Q(u, v):
        return (u * u - n * v * v) % q

    ok1 = True
    for x, y in product(range(q), repeat=2):
        u, v = (x + y) % q, (x - y) % q
        value = (Q(u, 0) - Q(v, 0)) * inverse(4, q) % q
        ok1 &= value == (x * y) % q

    ok2 = True
    for x1, y1, x2, y2 in product(range(q), repeat=4):
        u1, v1 = (x1 + y1) % q, (x1 - y1) % q
        u2, v2 = (x2 + y2) % q, (x2 - y2) % q
        h1, h2, h3 = Q(u1, u2), Q(v1, u2), Q(u1, v2)
        ok2 &= (h1 - h2) * inverse(4, q) % q == (x1 * y1) % q
        ok2 &= (h1 - h3) * inverse((-4 * n) % q, q) % q == (x2 * y2) % q
    return ok1, ok2


def check_theorem3(q: int) -> dict[str, bool]:
    def g(u, v):
        return (u * u + v) % q

    results = {}
    ok = True
    for x, y in product(range(q), repeat=2):
        ok &= (g((x + y) % q, 0) - g((x - y) % q, 0)) * inverse(4, q) % q == x * y % q
    results["xy from 2 atoms"] = ok
    ok = True
    for x, y in product(range(q), repeat=2):
        ok &= (g(x, 0) + g(y, 0)) % q == (x * x + y * y) % q
    results["x^2+y^2 from 2 atoms"] = ok
    ok = True
    for x, y in product(range(q), repeat=2):
        ok &= g(x, y) == (x * x + y) % q
    results["x^2+y from 1 atom"] = ok
    ok = True
    for x, y in product(range(q), repeat=2):
        ok &= g(0, x) == x % q
    results["x from 1 atom"] = ok
    return results


def main() -> None:
    print("Theorem 2 witnesses (anisotropic Q = x^2 - n y^2, target x y):")
    for q in PRIMES:
        ok1, ok2 = check_theorem2(q)
        print(f"  q={q:>3}  n={least_non_residue(q):>2}  "
              f"N_1 <= 2 verified on q^2 points: {ok1}   "
              f"N_2 <= 3 verified on q^4 points: {ok2}")
    print()
    print("Theorem 3 witnesses (parabolic resource g = x^2 + y):")
    for q in PRIMES:
        results = check_theorem3(q)
        summary = "  ".join(f"{name}: {value}" for name, value in results.items())
        print(f"  q={q:>3}  {summary}")


if __name__ == "__main__":
    main()
