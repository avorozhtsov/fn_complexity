#!/usr/bin/env python3
"""Fiber-count vectors of maps ``f : A^2 -> A^1`` over ``F_q``, ``q`` prime.

Everything is a length-``q`` numpy array ``counts[c] = #f^{-1}(c)``, computed by
``bincount`` over all ``q^2`` input pairs.  ``sum(counts) == q^2`` always.
"""

from __future__ import annotations

import numpy as np


def _grid(q: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(q, dtype=np.int64)
    return x[:, None], x[None, :]


def counts_of(values: np.ndarray, q: int) -> np.ndarray:
    return np.bincount((values % q).ravel(), minlength=q).astype(np.int64)


def poly_eval(coeffs, t: np.ndarray, q: int) -> np.ndarray:
    """Horner evaluation of ``sum coeffs[i] t^i`` mod q (coeffs low-to-high)."""
    out = np.zeros_like(t)
    for c in reversed(coeffs):
        out = (out * t + int(c)) % q
    return out


# ------------------------------------------------------------ named families


def hyperelliptic(q: int, coeffs) -> np.ndarray:
    """``f(x,y) = y^2 - P(x)``; fiber over c is the affine curve ``y^2 = P(x)+c``."""
    X, Y = _grid(q)
    return counts_of((Y * Y - poly_eval(coeffs, X, q)) % q, q)


def superelliptic(q: int, r: int, coeffs) -> np.ndarray:
    """``f(x,y) = y^r - P(x)``."""
    X, Y = _grid(q)
    return counts_of((pow_mod(Y, r, q) - poly_eval(coeffs, X, q)) % q, q)


def pow_mod(t: np.ndarray, r: int, q: int) -> np.ndarray:
    out = np.ones_like(t)
    base = t % q
    e = r
    while e:
        if e & 1:
            out = (out * base) % q
        base = (base * base) % q
        e >>= 1
    return out


def additive(q: int, pc, qc) -> np.ndarray:
    """``f(x,y) = P(x) + Q(y)``."""
    X, Y = _grid(q)
    return counts_of((poly_eval(pc, X, q) + poly_eval(qc, Y, q)) % q, q)


def bilinear_family(q: int, coeffs2d: np.ndarray) -> np.ndarray:
    """``f(x,y) = sum_{i,j} A[i,j] x^i y^j`` for a small coefficient matrix."""
    X, Y = _grid(q)
    dx, dy = coeffs2d.shape
    xp = [np.ones_like(X)]
    for _ in range(dx - 1):
        xp.append((xp[-1] * X) % q)
    yp = [np.ones_like(Y)]
    for _ in range(dy - 1):
        yp.append((yp[-1] * Y) % q)
    out = np.zeros((q, q), dtype=np.int64)
    for i in range(dx):
        for j in range(dy):
            a = int(coeffs2d[i, j]) % q
            if a:
                out = (out + a * xp[i] * yp[j]) % q
    return counts_of(out, q)


def flat_map(q: int) -> np.ndarray:
    """``f(x,y) = x``."""
    return np.full(q, q, dtype=np.int64)


def split_conic(q: int) -> np.ndarray:
    """``f(x,y) = x*y``."""
    X, Y = _grid(q)
    return counts_of(X * Y, q)


def aniso_conic(q: int) -> np.ndarray:
    """``f(x,y) = x^2 - d y^2`` with d a non-residue."""
    d = non_residue(q)
    X, Y = _grid(q)
    return counts_of(X * X - d * Y * Y, q)


def non_residue(q: int) -> int:
    sq = {(i * i) % q for i in range(q)}
    for d in range(2, q):
        if d not in sq:
            return d
    raise ValueError("no non-residue")


# ----------------------------------------------------------------- samplers


def random_maps(q: int, count: int, dx: int, dy: int, rng: np.random.Generator):
    """Random dense bivariate polynomials of bidegree < (dx, dy)."""
    for _ in range(count):
        A = rng.integers(0, q, size=(dx, dy))
        yield A, bilinear_family(q, A)


def random_hyperelliptic(q: int, count: int, deg: int, rng: np.random.Generator):
    """``y^2 = P(x) + c`` with P monic of degree ``deg`` and random lower terms."""
    for _ in range(count):
        c = list(rng.integers(0, q, size=deg)) + [1]
        yield tuple(int(v) for v in c), hyperelliptic(q, c)
