#!/usr/bin/env python3
"""The flux ``A``, its Hodge split, and the verification of both.

Conventions are brief D Part 0's, *not* the seed script's (which carried the
opposite sign).  With ``u_a(beta) = log log Z_a(beta)`` and ``g = u_b - u_a``,

    L(a, b) = -log C(a -> b) = max g,        S = (L_ab + L_ba)/2 = osc(g)/2
    A(a, b) = (L_ab - L_ba)/2 = mid(g),      a < b  <=>  A(a, b) > 0

so ``A[i, j] = mid_beta(u_j - u_i)``.  The two endpoints are exact and are
included by hand: ``beta = 0`` gives ``u_j - u_i = 0`` for every pair of
signatures with the same number of fibers, and ``beta = infinity`` gives
``log log max_j - log log max_i``.

HodgeRank: the least-squares potential of an antisymmetric ``A`` is
``psi = -rowmean(A)`` and the gradient part is ``grad[i, j] = psi[j] - psi[i]``.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

TIE = 1e-10  # the package's floor; anything below this is a tie


def beta_grid(q: int, coarse: int = 600, fine: int = 3400, horizon: float = 360.0) -> np.ndarray:
    """Grid on ``(0, horizon*q]``.  Brief B's horizon is ``360 q``."""

    return np.concatenate(
        [np.linspace(0.0, 2.0, coarse)[1:], np.geomspace(2.0, horizon * q, fine)]
    )


def u_grid(signatures, betas: np.ndarray):
    """``u = log log Z`` on the grid, and the two endpoint values.

    ``signatures`` may be a rectangular integer array (all pools here have
    exactly ``q`` fibers) or a ragged list of tuples.  ``u(0) = log log #fibers``
    and ``u(beta) - log beta -> log log max`` as ``beta -> infinity``.
    """

    rows = [np.asarray(s, dtype=float) for s in signatures]
    n = len(rows)
    U = np.empty((n, betas.size))
    for i, row in enumerate(rows):
        z = np.outer(betas, np.log(row))
        m = z.max(axis=1)
        U[i] = np.log(m + np.log(np.exp(z - m[:, None]).sum(axis=1)))
    u_zero = np.array([math.log(math.log(len(r))) for r in rows])
    u_inf = np.array([math.log(math.log(r.max())) for r in rows])
    return U, u_zero, u_inf


def flux_matrix(signatures, betas: np.ndarray) -> np.ndarray:
    """``A[i, j] = mid_beta(u_j - u_i)`` including both endpoints."""

    U, u_zero, u_inf = u_grid(signatures, betas)
    n = U.shape[0]
    A = np.empty((n, n))
    for i in range(n):
        D = U - U[i][None, :]  # (n, grid): u_j - u_i
        e0 = u_zero - u_zero[i]
        ei = u_inf - u_inf[i]
        lo = np.minimum(D.min(axis=1), np.minimum(e0, ei))
        hi = np.maximum(D.max(axis=1), np.maximum(e0, ei))
        A[i] = 0.5 * (lo + hi)
    A = 0.5 * (A - A.T)
    np.fill_diagonal(A, 0.0)
    return A


def hodge(A: np.ndarray) -> dict:
    """Least-squares gradient/curl split of an antisymmetric ``A``."""

    n = A.shape[0]
    psi = -A.mean(axis=1)
    psi -= psi.mean()
    G = psi[None, :] - psi[:, None]
    R = A - G
    nA = float(np.linalg.norm(A))
    off = ~np.eye(n, dtype=bool)
    strict = np.abs(A) > TIE
    agree = int((np.sign(A[off & strict]) == np.sign(G[off & strict])).sum())
    total = int((off & strict).sum())
    return {
        "n": n,
        "normA": nA,
        "grad_frac": float(np.linalg.norm(G)) / nA,
        "curl_frac": float(np.linalg.norm(R)) / nA,
        "grad_energy": float(np.linalg.norm(G)) ** 2 / nA**2,
        "curl_energy": float(np.linalg.norm(R)) ** 2 / nA**2,
        "psi": psi,
        "order_agreement": agree / total if total else float("nan"),
        "ties": int((off & ~strict).sum()) // 2,
    }


def cycle_count(A: np.ndarray, tie: float = TIE) -> tuple[int, int]:
    """Number of strict 3-cycles, and the number of triangles inspected."""

    n = A.shape[0]
    P = A > tie
    total = math.comb(n, 3)
    count = 0
    for i in range(n):
        Pi = P[i]
        for j in range(i + 1, n):
            if P[i, j]:
                # i < j : need j < k and k < i
                count += int(np.count_nonzero(P[j, j + 1 :] & (~Pi[j + 1 :]) & (np.abs(A[i, j + 1 :]) > tie)))
            elif P[j, i]:
                count += int(np.count_nonzero(Pi[j + 1 :] & (~P[j, j + 1 :]) & (np.abs(A[j, j + 1 :]) > tie)))
    return count, total


def cycle_count_reference(A: np.ndarray, tie: float = TIE) -> int:
    """Brute force, for checking :func:`cycle_count` on small pools."""

    n = A.shape[0]
    P = A > tie
    return sum(
        1
        for i, j, k in itertools.combinations(range(n), 3)
        if (P[i, j] and P[j, k] and P[k, i]) or (P[j, i] and P[k, j] and P[i, k])
    )


def verify_against_package(signatures: np.ndarray, A: np.ndarray, pairs, rng) -> float:
    """Max deviation of grid ``A`` from the package's exchange rates."""

    import sys

    sys.path.insert(0, "../../src")
    from fn_complexity import exchange_rate

    worst = 0.0
    for i, j in pairs:
        a = tuple(int(v) for v in signatures[i])
        b = tuple(int(v) for v in signatures[j])
        L_ab = -math.log(exchange_rate(a, b))
        L_ba = -math.log(exchange_rate(b, a))
        worst = max(worst, abs((L_ab - L_ba) / 2 - A[i, j]))
    return worst
