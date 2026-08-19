#!/usr/bin/env python3
"""The map ``pencil -> signature`` and its fibers.

The enumeration of brief B and of :mod:`pools` runs over ``P`` monic of degree
``d in {5, 6}`` with ``P(0) = 0``, i.e. over ``q^(d-1)`` polynomials.  Different
polynomials can define *the same pencil up to isomorphism*: substituting
``x -> a x + b`` and renormalising gives

    Q(x) = a^{-d} (P(a x + b) - P(b)),      Q monic, Q(0) = 0,

and the fiber of the ``Q``-pencil over ``c`` is the quadratic twist by ``a^d``
of the fiber of the ``P``-pencil over ``a^{-d}(P(b) + c)``.  So the whole affine
group of order ``q (q-1)`` acts on the enumeration, and the subgroup

    G0 = {(a, b) : a^d is a square in F_q}

acts by *isomorphisms of the fibration*, hence preserves the signature.  On the
complement the signature is sent to its twist ``N_c -> 2q - N_c``.

Counting fibers of ``pencil -> signature`` in ``G0``-orbits rather than in raw
polynomials is the honest measurement: an orbit is one pencil.
"""

from __future__ import annotations

import itertools
import math

import numpy as np


def coefficient_grid(q: int, degree: int) -> np.ndarray:
    """Rows ``(c_1, ..., c_{degree-1})`` in the order :mod:`pools` enumerates."""

    return np.array(list(itertools.product(range(q), repeat=degree - 1)), dtype=np.int64)


def _transform(coeffs: np.ndarray, a: int, b: int, q: int, degree: int) -> np.ndarray:
    """``Q_j = a^{j-d} sum_{k=j}^{d} c_k binom(k, j) b^{k-j}``, ``c_d = 1``."""

    m = coeffs.shape[0]
    out = np.zeros_like(coeffs)
    inv_a = pow(a, q - 2, q)
    for j in range(1, degree):
        acc = np.zeros(m, dtype=np.int64)
        for k in range(j, degree + 1):
            ck = coeffs[:, k - 1] if k < degree else 1
            acc = (acc + ck * math.comb(k, j) * pow(b, k - j, q)) % q
        out[:, j - 1] = (acc * pow(a, j, q) * pow(inv_a, degree, q)) % q
    return out


def canonical(q: int, degree: int) -> np.ndarray:
    """For every ``P``, the least index in its ``G0``-orbit."""

    coeffs = coefficient_grid(q, degree)
    # The weights must invert :func:`coefficient_grid`, so that the value stored
    # here is the *row index* of the canonical representative and the polynomial
    # can be read back.  ``itertools.product`` varies the first coordinate
    # slowest, hence the descending powers.
    weights = q ** np.arange(degree - 2, -1, -1, dtype=np.int64)
    squares = {(x * x) % q for x in range(1, q)}
    best = coeffs @ weights
    for a in range(1, q):
        if pow(a, degree, q) not in squares:
            continue
        for b in range(q):
            if a == 1 and b == 0:
                continue
            best = np.minimum(best, _transform(coeffs, a, b, q, degree) @ weights)
    return best


def group_order(q: int, degree: int) -> int:
    squares = {(x * x) % q for x in range(1, q)}
    return q * sum(1 for a in range(1, q) if pow(a, degree, q) in squares)
