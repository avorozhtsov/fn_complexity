#!/usr/bin/env python3
"""Exact point counting for one-parameter families of hyperelliptic curves.

Shared by ``repeated_factor.py`` and ``witness_search.py``.  Everything is over
a prime field ``F_q``; a family is given as a callable ``c -> coefficient
vector`` and the engine returns the Frobenius trace
``a = q + 1 - #C(F_q)`` of the smooth projective model of ``y^2 = f(x)`` for
every parameter, together with the moments of ``alpha = -a/sqrt(q)``.

The evaluation is vectorised over a ``(chunk of parameters) x F_q`` grid, so
the cost is ``O(q^2)`` field operations done inside numpy.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "chi_table", "is_prime", "family_traces", "moments", "squarefree",
    "polymul", "even_palindromic", "compose_square",
]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def chi_table(q: int) -> np.ndarray:
    chi = -np.ones(q, dtype=np.int64)
    chi[0] = 0
    chi[(np.arange(1, q, dtype=np.int64) ** 2) % q] = 1
    return chi


def polymul(a, b, q: int) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai % q == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % q
    return out


def _polymod(a, b, q):
    a = list(a)
    inv = pow(int(b[0]) % q, q - 2, q)
    while len(a) >= len(b):
        f = (a[0] * inv) % q
        for i in range(len(b)):
            a[i] = (a[i] - f * b[i]) % q
        a = a[1:]
        while a and a[0] == 0:
            a = a[1:]
    return a


def squarefree(coeffs, q: int) -> bool:
    """``gcd(f, f') = 1`` over ``F_q``; coefficients highest degree first."""
    a = [int(c) % q for c in coeffs]
    while a and a[0] == 0:
        a = a[1:]
    if len(a) <= 1:
        return False
    n = len(a) - 1
    d = [((n - i) * a[i]) % q for i in range(n)]
    while d and d[0] == 0:
        d = d[1:]
    if not d:
        return False          # inseparable
    while d:
        a, d = d, _polymod(a, d, q)
        while d and d[0] == 0:
            d = d[1:]
    return len(a) == 1


def compose_square(g):
    """Coefficients of ``g(x^2)`` from those of ``g(u)`` (highest first)."""
    out = []
    for c in g:
        out.append(c)
        out.append(0)
    return out[:-1]


def even_palindromic(m: int, params):
    """``f(x) = sum_j c_j x^{2j}`` with ``c_j = c_{m-j}``, ``c_0 = c_m = 1``.

    ``params`` gives ``c_1, ..., c_{floor(m/2)}``; the returned list is the
    coefficient vector of ``f``, degree ``2m``, highest first.
    """
    half = [1] + list(params)
    g = [half[min(j, m - j)] for j in range(m + 1)]     # c_0 .. c_m
    return compose_square(list(reversed(g)))


def family_traces(q: int, coeff_rows: np.ndarray, deg: int,
                  chunk: int = 4096) -> np.ndarray:
    """Frobenius traces of ``y^2 = f_i(x)`` for a stack of coefficient rows.

    ``coeff_rows`` has shape ``(n, deg+1)``, highest degree first, entries in
    ``[0, q)``.  Leading coefficients must be non-zero.
    """
    chi = chi_table(q)
    x = np.arange(q, dtype=np.int64)
    n = coeff_rows.shape[0]
    out = np.empty(n, dtype=np.int64)
    for lo in range(0, n, chunk):
        hi = min(n, lo + chunk)
        blk = coeff_rows[lo:hi]
        acc = np.zeros((hi - lo, q), dtype=np.int64)
        for j in range(deg + 1):
            acc = (acc * x[None, :] + blk[:, j][:, None]) % q
        pts = np.sum(1 + chi[acc], axis=1)
        if deg % 2 == 1:
            pts = pts + 1
        else:
            pts = pts + 1 + chi[blk[:, 0] % q]
        out[lo:hi] = q + 1 - pts
    return out


def moments(a: np.ndarray, q: int, upto: int = 6) -> dict[int, float]:
    alpha = -a.astype(np.float64) / math.sqrt(q)
    return {j: float(np.mean(alpha ** j)) for j in range(1, upto + 1)}
