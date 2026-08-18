#!/usr/bin/env python3
"""Arithmetic of the fibers of a genus-two pencil ``y^2 = P(x) + c`` over ``F_p``.

Everything here is elementary point counting; nothing needs a computer algebra
system.  For each fiber the smooth projective model is counted over ``F_p`` and
``F_{p^2}``, which determines the L-polynomial of a genus-two curve completely:

    L(T) = 1 - s1 T + ((s1^2 - s2)/2) T^2 - p s1 T^3 + p^2 T^4,
    s1 = p + 1 - #C(F_p),   s2 = p^2 + 1 - #C(F_{p^2}),

and hence the isogeny class of the Jacobian (Tate).  The real Weil polynomial is
``h(x) = x^2 - s1 x + (e2 - 2p)``; the Jacobian is isogenous over ``F_p`` to a
product of two elliptic curves exactly when ``h`` has integer roots.  The
``p``-rank is read off the Newton polygon of ``L``.

Points at infinity on the smooth model of ``y^2 = f(x)``:  one when ``deg f`` is
odd, two when ``deg f`` is even with square leading coefficient (all our ``f``
are monic).  So ``#C(F_p) = N_affine + 1`` or ``+ 2``, and the "trace"
``a_c = p - N_affine`` of ``research/m_and_e_and_a_c/FINDINGS.md`` equals the
genuine Frobenius trace ``s1`` for degree 5 and ``s1 + 1`` for degree 6.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def nonresidue(p: int) -> int:
    squares = {(x * x) % p for x in range(1, p)}
    return next(r for r in range(2, p) if r not in squares)


class F2:
    """``F_{p^2} = F_p[t]/(t^2 - r)`` on numpy pairs ``(u, v)`` for ``u + v t``."""

    def __init__(self, p: int):
        self.p = p
        self.r = nonresidue(p)

    def mul(self, x, y):
        (u1, v1), (u2, v2) = x, y
        return ((u1 * u2 + self.r * v1 * v2) % self.p, (u1 * v2 + u2 * v1) % self.p)

    def pow(self, x, e: int):
        result = (np.ones_like(x[0]), np.zeros_like(x[1]))
        base = x
        while e:
            if e & 1:
                result = self.mul(result, base)
            base = self.mul(base, base)
            e >>= 1
        return result

    def elements(self):
        u, v = np.meshgrid(np.arange(self.p), np.arange(self.p), indexing="ij")
        return (u.ravel().astype(np.int64), v.ravel().astype(np.int64))

    def chi(self, x):
        """``+1`` on nonzero squares, ``-1`` on nonsquares, ``0`` on zero."""

        u, v = x
        e = (self.p * self.p - 1) // 2
        pu, pv = self.pow(x, e)
        out = np.where(pu == 1, 1, -1)
        return np.where((u == 0) & (v == 0), 0, out)


def poly_values_fp(coeffs: dict[int, int], p: int, xs: np.ndarray) -> np.ndarray:
    out = np.zeros_like(xs)
    for k, c in coeffs.items():
        out = (out + c * pow_mod(xs, k, p)) % p
    return out


def pow_mod(xs: np.ndarray, k: int, p: int) -> np.ndarray:
    out = np.ones_like(xs)
    base = xs % p
    while k:
        if k & 1:
            out = (out * base) % p
        base = (base * base) % p
        k >>= 1
    return out


@dataclass
class Fiber:
    c: int
    n_affine: int
    points_fp: int
    points_fp2: int
    s1: int
    e2: int
    smooth: bool

    @property
    def L(self) -> tuple[int, int]:
        return (self.s1, self.e2)

    def real_weil_discriminant(self, p: int) -> int:
        return self.s1 * self.s1 - 4 * (self.e2 - 2 * p)

    def splits(self, p: int) -> bool:
        """Jacobian isogenous over ``F_p`` to a product of elliptic curves."""

        d = self.real_weil_discriminant(p)
        if d < 0:
            return False
        root = math.isqrt(d)
        return root * root == d

    def p_rank(self, p: int) -> int:
        if self.s1 % p != 0:
            return 2
        return 1 if self.e2 % p != 0 else 0


def squarefree(coeffs: dict[int, int], p: int, degree: int) -> bool:
    """Is ``P(x)`` squarefree over ``F_p``?  (gcd with the derivative.)"""

    f = [0] * (degree + 1)
    for k, c in coeffs.items():
        f[k] = c % p
    f = _trim(f)
    df = _trim([(k * f[k]) % p for k in range(1, len(f))])
    return _degree(_gcd(f, df, p)) == 0


def _trim(f):
    while len(f) > 1 and f[-1] == 0:
        f = f[:-1]
    return f


def _degree(f):
    f = _trim(f)
    return -1 if f == [0] else len(f) - 1


def _gcd(a, b, p):
    a, b = _trim(a[:]), _trim(b[:])
    while _degree(b) >= 0:
        a, b = b, _trim(_rem(a, b, p))
    return a


def _rem(a, b, p):
    a = a[:]
    db = _degree(b)
    inv = pow(b[db], p - 2, p)
    while _degree(a) >= db and _degree(a) >= 0:
        da = _degree(a)
        factor = (a[da] * inv) % p
        for i in range(db + 1):
            a[da - db + i] = (a[da - db + i] - factor * b[i]) % p
        a = _trim(a)
        if all(v == 0 for v in a):
            return [0]
    return a


def pencil_fibers(coeffs: dict[int, int], p: int, degree: int) -> list[Fiber]:
    """The eleven (``p``) fibers of ``y^2 = P(x) + c``, fully described."""

    at_infinity = 1 if degree % 2 else 2
    xs = np.arange(p, dtype=np.int64)
    Pv = poly_values_fp(coeffs, p, xs)
    squares = {(x * x) % p for x in range(1, p)}
    chi = np.array([0] + [1 if x in squares else -1 for x in range(1, p)])

    field = F2(p)
    X2 = field.elements()
    P2 = (np.zeros_like(X2[0]), np.zeros_like(X2[1]))
    for k, c in coeffs.items():
        term = field.pow(X2, k)
        P2 = ((P2[0] + c * term[0]) % p, (P2[1] + c * term[1]) % p)

    out = []
    for c in range(p):
        n_aff = p + int(chi[(Pv + c) % p].sum())
        shifted = ((P2[0] + c) % p, P2[1])
        n_aff2 = p * p + int(field.chi(shifted).sum())
        n1 = n_aff + at_infinity
        n2 = n_aff2 + at_infinity
        s1 = p + 1 - n1
        s2 = p * p + 1 - n2
        e2 = (s1 * s1 - s2) // 2
        shifted_coeffs = dict(coeffs)
        shifted_coeffs[0] = (shifted_coeffs.get(0, 0) + c) % p
        out.append(
            Fiber(c, n_aff, n1, n2, s1, e2, squarefree(shifted_coeffs, p, degree))
        )
    return out


def value_multiplicities(coeffs: dict[int, int], p: int) -> np.ndarray:
    xs = np.arange(p, dtype=np.int64)
    return np.bincount(poly_values_fp(coeffs, p, xs), minlength=p)


def factorisation_type(coeffs: dict[int, int], p: int, degree: int, u: int):
    """Degrees of the irreducible factors of ``P(x) - u`` over ``F_p``,
    with multiplicity, as a sorted tuple of ``(degree, multiplicity)``."""

    f = [0] * (degree + 1)
    for k, c in coeffs.items():
        f[k] = c % p
    f[0] = (f[0] - u) % p
    f = _trim(f)
    return tuple(sorted(_factor_degrees(f, p)))


def _factor_degrees(f, p):
    """Distinct-degree factorisation, enough for the degree multiset."""

    out = []
    # squarefree part first, by repeated gcd with the derivative
    remaining = f
    mult = 1
    while _degree(remaining) > 0:
        df = _trim([(k * remaining[k]) % p for k in range(1, len(remaining))])
        if _degree(df) < 0:  # p-th power
            root = [remaining[i] for i in range(0, len(remaining), p)]
            remaining = _trim(root)
            mult *= p
            continue
        g = _gcd(remaining, df, p)
        sf = _quo(remaining, g, p)
        for d in _distinct_degree(sf, p):
            out.append((d, mult))
        if _degree(g) <= 0:
            break
        remaining = g
        mult += 1
    return out


def _quo(a, b, p):
    a = a[:]
    db = _degree(b)
    inv = pow(b[db], p - 2, p)
    quotient = [0] * (max(_degree(a) - db, 0) + 1)
    while _degree(a) >= db and _degree(a) >= 0:
        da = _degree(a)
        factor = (a[da] * inv) % p
        quotient[da - db] = factor
        for i in range(db + 1):
            a[da - db + i] = (a[da - db + i] - factor * b[i]) % p
        a = _trim(a)
        if all(v == 0 for v in a):
            break
    return _trim(quotient)


def _distinct_degree(f, p):
    """Degrees of the irreducible factors of a squarefree ``f``."""

    degrees = []
    x_power = [0, 1]
    d = 0
    while _degree(f) > 0:
        d += 1
        x_power = _powmod_poly(x_power, p, f, p)
        g = _gcd(_trim([(x_power[i] if i < len(x_power) else 0) - (1 if i == 1 else 0)
                        for i in range(max(len(x_power), 2))]), f, p)
        deg_g = _degree(g)
        if deg_g > 0:
            degrees.extend([d] * (deg_g // d))
            f = _quo(f, g, p)
        if d > _degree(f) + 1 and _degree(f) > 0:
            degrees.append(_degree(f))
            break
    return degrees


def derivative_type(coeffs: dict[int, int], p: int, degree: int):
    """Degrees of the irreducible factors of ``P'(x)`` -- the critical points."""

    f = [0] * (degree + 1)
    for k, c in coeffs.items():
        f[k] = c % p
    df = _trim([(k * f[k]) % p for k in range(1, len(f))])
    return tuple(sorted(_factor_degrees(df, p)))


def rational_critical_values(coeffs: dict[int, int], p: int, degree: int) -> int:
    """How many ``u in F_p`` have ``P(x) - u`` non-squarefree.

    These are the ``c = -u`` for which the fiber ``y^2 = P(x) + c`` is singular.
    """

    count = 0
    for u in range(p):
        shifted = dict(coeffs)
        shifted[0] = (shifted.get(0, 0) - u) % p
        if not squarefree(shifted, p, degree):
            count += 1
    return count


def _mulmod_poly(a, b, m, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] = (out[i + j] + ai * bj) % p
    return _trim(_rem(_trim(out), m, p))


def _powmod_poly(base, e, m, p):
    result = [1]
    while e:
        if e & 1:
            result = _mulmod_poly(result, base, m, p)
        base = _mulmod_poly(base, base, m, p)
        e >>= 1
    return result
