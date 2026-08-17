#!/usr/bin/env python3
"""Shared helpers for T2.3 (symmetry type and exact arithmetic detection).

Adds to ``ffmaps`` the fibrations that are *not* of the form ``y^2 = P(x) + c``
(twist families, the degenerate constant map), the Dickson branch maps, and the
bookkeeping used by the three T2.3 scripts.

Exact identities used throughout (proved in the note, re-checked in
``t2_3_cm_families.py``):

*   ``sum_c a_c = 0``            for every ``f : A^2 -> A^1`` with full image;
*   ``sum_c a_c^2 = Z_f(2) - q^3``   for every such ``f``;
*   ``sum_c a_c^2 = q (nu(P) - q)``  for ``f = y^2 - P(x)``, with
    ``nu(P) = #{(x,x') : P(x) = P(x')}``  (also found independently in T2.2,
    where it is written ``K_P``).
"""

from __future__ import annotations

import math

import numpy as np

import ffmaps as F


def primes_upto(hi: int, lo: int = 3) -> list[int]:
    sieve = np.ones(hi + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(hi ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p:: p] = False
    return [int(p) for p in np.nonzero(sieve)[0] if p >= lo]


def _grid(q: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(q, dtype=np.int64)
    return x[:, None], x[None, :]


# ------------------------------------------------------------- twist families


def quartic_twist_map(q: int) -> np.ndarray:
    """``f(x,y) = x y^2 - x^2``.

    Over ``c != 0`` the substitution ``Y = x y`` is a bijection of ``{x != 0}``
    onto itself and turns the fiber into ``{Y^2 = x^3 + c x, x != 0}``: the
    quartic twist ``E_c : Y^2 = x^3 + c x`` minus its rational 2-torsion point.
    The fiber over ``0`` is the degenerate one, of size ``2q - 1``.
    """
    X, Y = _grid(q)
    return F.counts_of((X * Y * Y - X * X) % q, q)


def sextic_twist_map(q: int) -> np.ndarray:
    """``f(x,y) = y^2 - x^3``; fibers are the sextic twists ``y^2 = x^3 + c``."""
    return F.hyperelliptic(q, [0, 0, 0, 1])


def quadratic_twist_map(q: int, coeffs) -> np.ndarray:
    """``f(x,y) = P(x) y^2`` for a Weierstrass cubic ``P`` of a fixed curve ``E``.

    For ``c != 0`` the fiber is in bijection (via ``Y = 1/y``) with the affine
    quadratic twist ``E^{(c)} : c Y^2 = P(x)`` minus the points with
    ``P(x) = 0``.  Exactly

        N_0 = q + z (q - 1),    N_c = q - z - chi(c) a_E   (c != 0),

    with ``z = #{x : P(x) = 0}`` and ``a_E = q - #{(x,y) : y^2 = P(x)}``.
    """
    X, Y = _grid(q)
    return F.counts_of((F.poly_eval(coeffs, X, q) * Y * Y) % q, q)


def constant_map(q: int) -> np.ndarray:
    """``f(x,y) = 0``: a single fiber of size ``q^2``."""
    out = np.zeros(q, dtype=np.int64)
    out[0] = q * q
    return out


# --------------------------------------------------------------- branch maps


def dickson(n: int) -> list[int]:
    """Dickson polynomial ``D_n(x, 1)``, coefficients low-to-high.

    ``D_n(u + u^{-1}) = u^n + u^{-n}``; ``D_n`` permutes ``F_q`` exactly when
    ``gcd(n, q^2 - 1) = 1``.
    """
    a: list[int] = [2]
    b: list[int] = [0, 1]
    if n == 0:
        return a
    if n == 1:
        return b
    for _ in range(2, n + 1):
        c = [0] + b
        m = max(len(c), len(a))
        c = c + [0] * (m - len(c))
        aa = a + [0] * (m - len(a))
        c = [c[i] - aa[i] for i in range(m)]
        a, b = b, c
    return b


def nu(coeffs, q: int) -> int:
    """``nu(P) = #{(x, x') in F_q^2 : P(x) = P(x')} = sum_u n_P(u)^2``."""
    x = np.arange(q, dtype=np.int64)
    n = np.bincount(F.poly_eval(coeffs, x, q), minlength=q).astype(np.int64)
    return int((n * n).sum())


def is_permutation_polynomial(coeffs, q: int) -> bool:
    return nu(coeffs, q) == q


# ------------------------------------------------------------------- moments


def trace_moments(counts: np.ndarray, q: int, kmax: int = 4) -> dict[str, float]:
    a = (q - counts).astype(np.float64)
    out = {f"m{k}": float((a ** k).sum()) / (q * q ** (k / 2.0)) for k in range(1, kmax + 1)}
    out["sum_a"] = float(a.sum())
    out["sum_a2"] = float((a * a).sum())
    out["absmax"] = float(np.abs(a).max())
    out["max_fiber"] = int(counts.max())
    out["n_image"] = int((counts > 0).sum())
    return out


def sig_key(values: np.ndarray, mults: np.ndarray) -> tuple:
    return (tuple(int(v) for v in values), tuple(int(m) for m in mults))


def residue_classes(qs: list[int], flag: dict[int, bool], max_mod: int = 60):
    """Smallest modulus ``m`` making ``{q : flag[q]}`` a union of classes mod m.

    Returns ``(m, sorted_classes)`` or ``None`` when no modulus up to
    ``max_mod`` is consistent with the observed pattern.
    """
    for m in range(2, max_mod + 1):
        by_class: dict[int, set[bool]] = {}
        for q in qs:
            by_class.setdefault(q % m, set()).add(bool(flag[q]))
        if all(len(v) == 1 for v in by_class.values()):
            good = sorted(r for r, v in by_class.items() if True in v)
            if 0 < len(good) < len(by_class):
                return m, good
    return None


def weil_scale(rate_value: float, q: int) -> float:
    """``(1 - C) sqrt(q) log q``: the genus scale carried by ``C(L -> f)``."""
    return (1.0 - rate_value) * math.sqrt(q) * math.log(q)


def m2_from_rate(rate_value: float, q: int, kappa: float) -> float:
    """Invert ``C(f -> L) = 1 - kappa m2 / (2 q log q)`` for ``m2``."""
    return 2.0 * q * math.log(q) * (1.0 - rate_value) / kappa
