#!/usr/bin/env python3
"""Independent checks on the Sato--Tate MGF library of ``st_lib.py``.

Three checks, none of which reuses the determinant machinery it is testing.

1.  **Brute-force Weyl integration.**  The ``g``-dimensional integral
    ``int prod_{i<j}(cos t_i - cos t_j)^2 prod w(t_j) e^{tau tr} dt`` is
    evaluated by tensor-product Gauss--Legendre quadrature in the angles, for
    ``g = 1, 2, 3`` and several ``tau``, and compared with the determinant.

2.  **Exactly known moments.**  ``E[tr^{2k}]`` over ``USp(2g)`` is the dimension
    of the space of invariants in ``V^{otimes 2k}``, which is the number of
    perfect matchings of ``2k`` points -- ``(2k-1)!! = 1, 3, 15, 105`` -- once
    ``2g >= 2k``, with Brauer-algebra relations cutting it down below that.  The
    exact table is

        g = 1:  1, 2,  5,  14        (Catalan)
        g = 2:  1, 3, 14,  84
        g = 3:  1, 3, 15, 104
        g = 4:  1, 3, 15, 105

    for ``E[tr^2], E[tr^4], E[tr^6], E[tr^8]``.  Odd moments vanish.  For
    ``U(1)`` the moments are ``binom(2k,k) = 2, 6, 20, 70``, and for ``SO(n)``
    the second moment is ``1`` (the standard representation is irreducible).
    These are checked by finite differences of the MGF at 40 digits.

3.  **Edge exponents.**  ``t`` in ``M(tau) ~ A e^{tau alpha_max} tau^{-t}`` is
    read off numerically from ``tau (alpha_max - Psi(tau))`` at large ``tau``
    and compared with the closed form ``N^2 + aN``.

    python research/sato_tate_limit/validate_library.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from mpmath import mp, mpf, taylor

import st_lib as S

HERE = Path(__file__).resolve().parent


# ------------------------------------------------------- 1. brute force Weyl


def weyl_density(kind: str, theta: np.ndarray) -> np.ndarray:
    if kind == "sp":
        return np.sin(theta) ** 2
    if kind == "so_even":
        return np.ones_like(theta)
    if kind == "so_odd":
        return np.sin(theta / 2) ** 2
    raise ValueError(kind)


def brute_mgf(kind: str, rank: int, epsilon: int, tau: float, nodes: int) -> float:
    x, w = np.polynomial.legendre.leggauss(nodes)
    theta = 0.5 * math.pi * (x + 1.0)
    weight = 0.5 * math.pi * w * weyl_density(kind, theta)
    c = 2.0 * np.cos(theta)
    if rank == 1:
        num = float(np.sum(weight * np.exp(tau * c)))
        den = float(np.sum(weight))
    elif rank == 2:
        C1, C2 = np.meshgrid(c, c, indexing="ij")
        W = np.outer(weight, weight)
        vd = (C1 - C2) ** 2
        num = float(np.sum(W * vd * np.exp(tau * (C1 + C2))))
        den = float(np.sum(W * vd))
    elif rank == 3:
        C1 = c[:, None, None]
        C2 = c[None, :, None]
        C3 = c[None, None, :]
        W = weight[:, None, None] * weight[None, :, None] * weight[None, None, :]
        vd = ((C1 - C2) * (C1 - C3) * (C2 - C3)) ** 2
        num = float(np.sum(W * vd * np.exp(tau * (C1 + C2 + C3))))
        den = float(np.sum(W * vd))
    else:
        raise ValueError(rank)
    return math.exp(tau * epsilon) * num / den


# ------------------------------------------------------------- 2. moments


EXACT_MOMENTS = {
    ("USp2", 2): 1, ("USp2", 4): 2, ("USp2", 6): 5, ("USp2", 8): 14,
    ("USp4", 2): 1, ("USp4", 4): 3, ("USp4", 6): 14, ("USp4", 8): 84,
    ("USp6", 2): 1, ("USp6", 4): 3, ("USp6", 6): 15, ("USp6", 8): 104,
    ("USp8", 2): 1, ("USp8", 4): 3, ("USp8", 6): 15, ("USp8", 8): 105,
    ("U1", 2): 2, ("U1", 4): 6, ("U1", 6): 20, ("U1", 8): 70,
    ("SO3", 2): 1, ("SO4", 2): 1, ("SO5", 2): 1, ("SO6", 2): 1, ("SO7", 2): 1,
}


def moments_of(name: str, order: int) -> list[float]:
    """Taylor coefficients of the MGF at 0, times ``n!`` -- the raw moments."""

    grp = S.GROUPS[name]
    mp.dps = 90

    def f(t):
        if grp.kind == "sp" and grp.rank == 1:
            from mpmath import besseli
            return besseli(1, 2 * t) / t
        if grp.kind == "so_even" and grp.rank == 1:
            from mpmath import besseli
            return besseli(0, 2 * t)
        return S.mgf_classical(grp.kind, grp.rank, grp.epsilon, t)

    coeffs = taylor(f, mpf("1e-30"), order, method="step", h=mpf("1e-6"))
    out = [float(coeffs[k] * mp.factorial(k)) for k in range(order + 1)]
    mp.dps = 60
    return out


def main() -> int:
    rows: list[list] = []

    print("=" * 78)
    print("1.  determinant MGF against brute-force Weyl quadrature")
    print("=" * 78)
    print(f"{'group':>7} {'tau':>8} {'determinant':>22} {'quadrature':>22} {'rel':>10}")
    worst = 0.0
    for name in ("USp2", "USp4", "USp6", "SO3", "SO4", "SO5", "SO6", "U1"):
        grp = S.GROUPS[name]
        if grp.rank > 3:
            continue
        nodes = {1: 6000, 2: 3000, 3: 400}[grp.rank]
        for tau in (0.25, 1.0, 3.0, 6.0):
            det = float(S.mgf_classical(grp.kind, grp.rank, grp.epsilon, tau))
            bru = brute_mgf(grp.kind, grp.rank, grp.epsilon, tau, nodes)
            rel = abs(det - bru) / det
            worst = max(worst, rel)
            print(f"{name:>7} {tau:>8.2f} {det:>22.14g} {bru:>22.14g} {rel:>10.2e}")
            rows.append(["quadrature", name, tau, f"{det:.16g}", f"{bru:.16g}",
                         f"{rel:.3e}"])
    print(f"\n  worst relative deviation: {worst:.3e}\n")

    print("=" * 78)
    print("2.  moments against their exact values")
    print("=" * 78)
    print(f"{'group':>7} {'k':>3} {'computed E[tr^k]':>22} {'exact':>8} {'abs err':>10}")
    worst_m = 0.0
    for name in ("USp2", "USp4", "USp6", "USp8", "U1", "SO3", "SO4", "SO5", "SO6"):
        want = {k: v for (g, k), v in EXACT_MOMENTS.items() if g == name}
        if not want:
            continue
        order = max(want)
        got = moments_of(name, order)
        for k in sorted(want):
            err = abs(got[k] - want[k])
            worst_m = max(worst_m, err)
            print(f"{name:>7} {k:>3} {got[k]:>22.14g} {want[k]:>8} {err:>10.2e}")
            rows.append(["moment", name, k, f"{got[k]:.16g}", want[k], f"{err:.3e}"])
    print(f"\n  worst absolute deviation: {worst_m:.3e}\n")

    print("=" * 78)
    print("2b. the mean-zero constraint  E[alpha] = 0")
    print("=" * 78)
    print("  Every limit measure of a fibration of A^2 has mean zero exactly,")
    print("  because sum_c a_c = 0 identically.  E[tr] is the dimension of the")
    print("  invariants of the standard representation, hence 0 for every")
    print("  classical group in the library -- checked by a central difference.")
    mp.dps = 60
    h = mpf("1e-10")
    print(f"\n  {'group':>7} {'E[tr] (central difference)':>30}")
    for name in ("USp2", "USp4", "USp6", "USp8", "U1", "SO3", "SO4", "SO5",
                 "SO6", "SO7"):
        grp = S.GROUPS[name]

        def mgf(t):
            if grp.kind == "sp" and grp.rank == 1:
                from mpmath import besseli
                return besseli(1, 2 * t) / t
            if grp.kind == "so_even" and grp.rank == 1:
                from mpmath import besseli
                return besseli(0, 2 * t)
            return S.mgf_classical(grp.kind, grp.rank, grp.epsilon, t)

        mean = (mgf(h) - mgf(-h)) / (2 * h)
        print(f"  {name:>7} {mp.nstr(mean, 8):>30}")
        rows.append(["mean", name, "", mp.nstr(mean, 12), 0, ""])
    print()

    print("=" * 78)
    print("3.  edge exponent t from tau (alpha_max - Psi(tau))  ~  t log tau - log A")
    print("=" * 78)
    print(f"{'group':>7} {'closed form t':>14} {'fitted t':>12} {'log A':>12} {'err':>10}")
    taus = np.array([2.0e3, 4.0e3, 8.0e3, 1.6e4])
    for name in ("USp2", "USp4", "USp6", "USp8", "U1", "SO3", "SO4", "SO5", "SO6"):
        grp = S.GROUPS[name]
        k = S.group_K(name, taus)
        r = grp.alpha_max * taus - k          # = t log tau - log A + o(1)
        slope, intercept = np.polyfit(np.log(taus), r, 1)
        err = abs(slope - grp.tail)
        print(f"{name:>7} {grp.tail:>14.4f} {slope:>12.6f} {-intercept:>12.6f} "
              f"{err:>10.2e}")
        rows.append(["edge", name, grp.tail, f"{slope:.9f}", f"{-intercept:.9f}",
                     f"{err:.3e}"])

    print()
    print("=" * 78)
    print("4.  Psi endpoints:  Psi(0+) = 0 and Psi(inf) = alpha_max")
    print("=" * 78)
    small = np.array([1e-6, 1e-5, 1e-4])
    for name in ("USp2", "USp4", "USp6", "U1", "SO3", "SO4"):
        grp = S.GROUPS[name]
        psi0 = S.group_K(name, small) / small
        big = np.array([1e5])
        psi_inf = float(S.group_K(name, big) / big)
        print(f"{name:>7}  Psi(1e-6) = {psi0[0]:>12.3e}   "
              f"Psi(1e5) = {psi_inf:.6f}  (alpha_max = {grp.alpha_max})")
        rows.append(["endpoint", name, f"{psi0[0]:.6e}", f"{psi_inf:.9f}",
                     grp.alpha_max, ""])

    path = HERE / "validation.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["check", "group", "arg", "computed", "reference", "error"])
        writer.writerows(rows)
    print(f"\nwritten: {path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
