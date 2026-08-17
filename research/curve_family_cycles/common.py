#!/usr/bin/env python3
"""Shared machinery for the cycle search among families of curves.

Two things live here.

**A pool of fiber signatures.**  Every map is an ``f : A^2 -> A^1`` over a prime
field ``F_q``: hyperelliptic and superelliptic pencils ``y^r = P(x) + c``,
quadratic twist families ``P(x) y^2``, additive maps ``P(x) + Q(y)``, and dense
bivariate polynomials.  A signature is the sorted multiset of fiber cardinalities
``N_c = #f^{-1}(c)``, always summing to ``q^2``.

**A vectorised exchange-rate engine.**  ``exchange_rate`` in the package is exact
but costs ``O(grid * q)`` Python operations per pair, which is too slow for the
tens of thousands of pairs a cycle search needs.  The engine here evaluates
``log Z`` for a whole pool on one shared beta-grid with numpy, then minimises
ratios.  Everything it reports as a candidate is re-verified against the package
and against mpmath by ``verify.py``.

Scaling used to direct the search (derived in NOTES.md):  writing
``a_c = q - N_c`` and ``alpha_c = -a_c / sqrt(q)``, and ``beta = tau * sqrt(q)``,

    log Z_f(beta) = (1 + beta) log q + Lambda_f(tau) + O(1/q),
    Lambda_f(tau) = log( (1/q) sum_c exp(tau * alpha_c) ),

so with ``Psi_f(tau) = Lambda_f(tau) / tau``,

    C(u -> v) = 1 + inf_tau (Psi_u - Psi_v) / (sqrt(q) log q) + ...

``Psi_f`` increases from ``0`` (the first moment vanishes identically) to
``alpha_max = (max_c N_c - q)/sqrt(q)`` (the endpoint the index ``phi`` sees).
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


# --------------------------------------------------------------- signatures


@dataclass(frozen=True)
class Entry:
    """One pool member: a signature with a witness map that realises it."""

    signature: tuple[int, ...]
    family: str
    witness: str

    @property
    def fibers(self) -> int:
        return len(self.signature)

    @property
    def max_fiber(self) -> int:
        return self.signature[0]

    def phi(self) -> float:
        return math.log(self.fibers) * math.log(self.max_fiber)


def _grid(q: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(q, dtype=np.int64)
    return x[:, None], x[None, :]


def _counts(values: np.ndarray, q: int) -> np.ndarray:
    return np.bincount((values % q).ravel(), minlength=q).astype(np.int64)


def _poly(coeffs, t: np.ndarray, q: int) -> np.ndarray:
    out = np.zeros_like(t)
    for c in reversed(coeffs):
        out = (out * t + int(c)) % q
    return out


def _powmod(t: np.ndarray, r: int, q: int) -> np.ndarray:
    out = np.ones_like(t)
    base = t % q
    while r:
        if r & 1:
            out = (out * base) % q
        base = (base * base) % q
        r >>= 1
    return out


def _label(coeffs) -> str:
    terms = [
        (f"{c}x^{i}" if i > 1 else f"{c}x" if i == 1 else f"{c}")
        for i, c in enumerate(coeffs)
        if c
    ]
    return " + ".join(reversed(terms)) or "0"


def superelliptic_counts(q: int, r: int, coeffs) -> np.ndarray:
    """``f(x, y) = y^r - P(x)``; the fiber over ``c`` is ``y^r = P(x) + c``."""

    X, Y = _grid(q)
    return _counts((_powmod(Y, r, q) - _poly(coeffs, X, q)) % q, q)


def twist_counts(q: int, coeffs) -> np.ndarray:
    """``f(x, y) = P(x) y^2``; the fibers are the quadratic twist pencil."""

    X, Y = _grid(q)
    return _counts((_poly(coeffs, X, q) * Y * Y) % q, q)


def additive_counts(q: int, pc, qc) -> np.ndarray:
    """``f(x, y) = P(x) + Q(y)``."""

    X, Y = _grid(q)
    return _counts((_poly(pc, X, q) + _poly(qc, Y, q)) % q, q)


def dense_counts(q: int, matrix: np.ndarray) -> np.ndarray:
    """``f(x, y) = sum_{ij} A_ij x^i y^j``."""

    X, Y = _grid(q)
    dx, dy = matrix.shape
    xp = [np.ones_like(X)]
    for _ in range(dx - 1):
        xp.append((xp[-1] * X) % q)
    yp = [np.ones_like(Y)]
    for _ in range(dy - 1):
        yp.append((yp[-1] * Y) % q)
    out = np.zeros((q, q), dtype=np.int64)
    for i in range(dx):
        for j in range(dy):
            a = int(matrix[i, j]) % q
            if a:
                out = (out + a * xp[i] * yp[j]) % q
    return _counts(out, q)


def build_pool(q: int, *, seed: int = 20260817, budget: dict | None = None) -> list[Entry]:
    """Sample many maps ``A^2 -> A^1`` over ``F_q`` and keep distinct signatures.

    Only surjective maps are kept, so every entry has exactly ``q`` fibers and
    the endpoint ``beta = 0`` is a tie across the whole pool: the index ``phi``
    is then a strictly increasing function of the largest fiber alone.
    """

    rng = np.random.default_rng(seed)
    budget = budget or {}
    hyper = budget.get("hyperelliptic", 3000)
    supers = budget.get("superelliptic", 1200)
    twists = budget.get("twist", 1200)
    adds = budget.get("additive", 1200)
    dense = budget.get("dense", 1200)

    seen: dict[tuple[int, ...], Entry] = {}

    def offer(counts: np.ndarray, family: str, witness: str) -> None:
        if counts.min() == 0:
            return
        signature = tuple(sorted((int(v) for v in counts), reverse=True))
        seen.setdefault(signature, Entry(signature, family, witness))

    # y^2 = P(x) + c, deg P = 5..11: genus 2 through 5.
    for degree in (5, 6, 7, 8, 9, 11):
        for _ in range(hyper // 6):
            coeffs = [int(v) for v in rng.integers(0, q, size=degree)] + [1]
            offer(
                superelliptic_counts(q, 2, coeffs),
                f"hyperelliptic deg {degree}",
                f"y^2 - ({_label(coeffs)})",
            )

    # y^r = P(x) + c for r | q - 1 and r not dividing out to a permutation.
    for r in (3, 4, 5, 6):
        if math.gcd(r, q - 1) == 1:
            continue
        for degree in (3, 4, 5, 7):
            for _ in range(supers // 16):
                coeffs = [int(v) for v in rng.integers(0, q, size=degree)] + [1]
                offer(
                    superelliptic_counts(q, r, coeffs),
                    f"superelliptic r={r} deg {degree}",
                    f"y^{r} - ({_label(coeffs)})",
                )

    # P(x) y^2: the quadratic twist pencil.
    for degree in (3, 4, 5, 6, 7):
        for _ in range(twists // 5):
            coeffs = [int(v) for v in rng.integers(0, q, size=degree)] + [1]
            offer(
                twist_counts(q, coeffs),
                f"quadratic twist deg {degree}",
                f"({_label(coeffs)}) y^2",
            )

    # P(x) + Q(y): genus given by the two degrees.
    for _ in range(adds):
        dp = int(rng.integers(2, 8))
        dq = int(rng.integers(2, 8))
        pc = [int(v) for v in rng.integers(0, q, size=dp)] + [1]
        qc = [int(v) for v in rng.integers(0, q, size=dq)] + [1]
        offer(
            additive_counts(q, pc, qc),
            f"additive ({dp},{dq})",
            f"({_label(pc)}) + ({_label(qc)})[y]",
        )

    # Dense bivariate polynomials of small bidegree.
    for _ in range(dense):
        dx = int(rng.integers(3, 6))
        dy = int(rng.integers(3, 6))
        matrix = rng.integers(0, q, size=(dx, dy))
        offer(dense_counts(q, matrix), f"dense ({dx},{dy})", f"A={matrix.tolist()}")

    return list(seen.values())


# ------------------------------------------------------------- rate engine


def beta_grid(q: int, *, points: int = 24_000, low: float = 1e-3, high: float | None = None) -> np.ndarray:
    """Geometric grid.  The default reaches ``36 q``, the package's horizon.

    ``exp(-36)`` is below ``3e-16``, and adjacent fiber sizes differ by at least
    one, so ``beta = 36 q`` already isolates the largest fiber to double
    precision.  The grid is taken an order of magnitude beyond that.
    """

    high = high if high is not None else 360.0 * q
    return np.geomspace(low, high, points)


class Engine:
    """Vectorised ``log Z`` on a shared beta-grid, plus pairwise rates."""

    def __init__(self, signatures: list[tuple[int, ...]], betas: np.ndarray):
        self.signatures = signatures
        self.betas = betas
        self.log_z = np.empty((len(signatures), betas.size))
        for index, signature in enumerate(signatures):
            values, mults = np.unique(np.asarray(signature, dtype=np.int64), return_counts=True)
            logs = np.log(values.astype(float))
            top = logs.max()
            weights = np.exp(np.outer(logs - top, betas))
            self.log_z[index] = betas * top + np.log(mults.astype(float) @ weights)
        self.log_max = np.array([math.log(max(s)) for s in signatures])
        self.log_count = np.array([math.log(len(s)) for s in signatures])

    def rate_matrix(self, chunk: int = 512) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(C, beta_star)`` with ``C[u, v] = C(u -> v)`` on the grid.

        Endpoints are included exactly: ``beta = 0`` contributes
        ``log(#fibers_u)/log(#fibers_v)`` and ``beta = inf`` contributes
        ``log(max_u)/log(max_v)``.
        """

        n = len(self.signatures)
        best = np.full((n, n), np.inf)
        argbest = np.zeros((n, n))
        for start in range(0, self.betas.size, chunk):
            block = self.log_z[:, start : start + chunk]
            ratios = block[:, None, :] / block[None, :, :]
            local = ratios.argmin(axis=2)
            value = np.take_along_axis(ratios, local[:, :, None], axis=2)[:, :, 0]
            improved = value < best
            best = np.where(improved, value, best)
            argbest = np.where(improved, self.betas[start + local], argbest)
        zero = self.log_count[:, None] / self.log_count[None, :]
        infty = self.log_max[:, None] / self.log_max[None, :]
        for endpoint, marker in ((zero, 0.0), (infty, np.inf)):
            improved = endpoint < best
            argbest = np.where(improved, marker, argbest)
            best = np.where(improved, endpoint, best)
        np.fill_diagonal(best, 1.0)
        return best, argbest


def psi(signature: tuple[int, ...], q: int, taus: np.ndarray) -> np.ndarray:
    """``Psi(tau) = Lambda(tau)/tau`` with ``Lambda`` the CGF of ``alpha_c``."""

    alpha = (np.asarray(signature, dtype=float) - q) / math.sqrt(q)
    top = alpha.max()
    out = np.empty_like(taus)
    for index, tau in enumerate(taus):
        shifted = tau * (alpha - top)
        out[index] = (tau * top + math.log(np.exp(shifted).sum() / len(alpha))) / tau
    return out
