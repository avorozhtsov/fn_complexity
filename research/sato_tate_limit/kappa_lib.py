#!/usr/bin/env python3
"""``K_{USp(2g)}(tau)`` on a grid, fast and to arbitrary precision.

``st_lib.mgf_classical`` rebuilds the whole ``g x g`` moment matrix entry by
entry, which costs ``g^2`` Chebyshev sums per ``tau``; for the genus sweep of
``TRANSITIVITY.md`` we need ``g`` up to 14 on grids of a thousand points, so the
moments are computed once per ``tau`` and shared between all ranks.

Everything is the same Weyl/Andreief determinant as ``st_lib``:

    E_{USp(2g)}[e^{tau tr}] = det( L_tau[x^{i+j}] )_{0<=i,j<g} / (same at tau=0)
    L_tau[T_k] = I_k(2 tau) - (I_{k-2}(2 tau) + I_{k+2}(2 tau))/2

so one call to ``besseli`` per order per ``tau`` suffices.  Agreement with
``st_lib.mgf_classical`` is checked in ``transitivity_dominance.py``.

The Hankel determinants ``H_g(tau) = det(L_tau[x^{i+j}])_{g x g}`` are returned
as well, because the quantity that decides the whole genus question is

    M_{g+1} M_{g-1} / M_g^2  =  ( H_{g+1} H_{g-1} / H_g^2 )(tau)
                                / ( H_{g+1} H_{g-1} / H_g^2 )(0)
                             =  4 b_g(tau)^2,

``b_g`` the off-diagonal Jacobi coefficient of the tilted Chebyshev weight
``e^{2 tau x} sqrt(1-x^2)``, whose free value is ``b_g(0) = 1/2``.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
from mpmath import besseli, matrix, mp, mpf

HERE = Path(__file__).resolve().parent
CACHE = HERE / ".kcache"


def _cheb(n: int) -> list[tuple[int, float]]:
    poly = np.zeros(n + 1)
    poly[n] = 1.0
    c = np.polynomial.chebyshev.poly2cheb(poly)
    return [(k, float(v)) for k, v in enumerate(c) if v != 0.0]


_CHEB: dict[int, list[tuple[int, float]]] = {}


def cheb(n: int) -> list[tuple[int, float]]:
    if n not in _CHEB:
        _CHEB[n] = _cheb(n)
    return _CHEB[n]


def sp_moments(tau, nmax: int) -> list:
    """``L_tau[x^n]`` for ``n = 0 .. nmax``, weight ``(2/pi) sqrt(1-x^2)``."""

    two = 2 * tau
    bes = [besseli(k, two) for k in range(nmax + 3)]

    def T(k: int):
        return bes[k] - (bes[abs(k - 2)] + bes[k + 2]) / 2

    tk = [T(k) for k in range(nmax + 1)]
    return [sum(mpf(c) * tk[k] for k, c in cheb(n)) for n in range(nmax + 1)]


def hankels(tau, gmax: int, dps: int) -> list:
    """``[H_0, H_1, ..., H_gmax]`` with ``H_0 = 1``."""

    saved = mp.dps
    try:
        mp.dps = dps
        mom = sp_moments(mpf(tau), 2 * gmax - 2 if gmax >= 1 else 0)
        out = [mpf(1)]
        for g in range(1, gmax + 1):
            A = matrix(g, g)
            for i in range(g):
                for j in range(g):
                    A[i, j] = mom[i + j]
            out.append(mp.det(A))
        return out
    finally:
        mp.dps = saved


def working_dps(gmax: int, tau: float) -> int:
    cancel = gmax * (gmax - 1) * math.log10(max(2.0 * abs(tau), 10.0))
    return int(60 + 1.5 * cancel)


def kappa_and_b(tau_grid: np.ndarray, gmax: int,
                extra_dps: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """``kappa[g-1] = K_{USp(2g)}`` and ``bsq[g-1] = b_g(tau)^2``.

    ``kappa`` has shape ``(gmax, len(tau_grid))``; ``bsq`` has shape
    ``(gmax - 1, len(tau_grid))`` and holds ``b_g(tau)^2`` for ``g = 1..gmax-1``
    (``b_g^2 = H_{g+1} H_{g-1} / H_g^2``).
    """

    tau_grid = np.ascontiguousarray(np.asarray(tau_grid, dtype=np.float64))
    CACHE.mkdir(exist_ok=True)
    key = hashlib.sha1(f"kappa{gmax}:{extra_dps}".encode()
                       + tau_grid.tobytes()).hexdigest()[:20]
    path = CACHE / f"kappa_{gmax}_{key}.npz"
    if path.exists():
        z = np.load(path)
        return z["kappa"], z["bsq"]

    h0 = hankels(0.0, gmax + 1, working_dps(gmax + 1, 1.0) + extra_dps)
    kappa = np.empty((gmax, tau_grid.size))
    bsq = np.empty((gmax, tau_grid.size))
    for i, t in enumerate(tau_grid):
        dps = working_dps(gmax + 1, float(t)) + extra_dps
        h = hankels(float(t), gmax + 1, dps)
        saved = mp.dps
        mp.dps = dps
        try:
            for g in range(1, gmax + 1):
                kappa[g - 1, i] = float(mp.log(h[g] / h0[g]))
            for g in range(1, gmax + 1):
                bsq[g - 1, i] = float(h[g + 1] * h[g - 1] / h[g] ** 2)
        finally:
            mp.dps = saved
    np.savez(path, kappa=kappa, bsq=bsq)
    return kappa, bsq


def psi_from_kappa(kappa: np.ndarray, tau_grid: np.ndarray,
                   parts: tuple[int, ...]) -> np.ndarray:
    """``Psi`` of ``prod_i USp(2 parts_i)`` on the grid."""

    return sum(kappa[p - 1] for p in parts) / tau_grid
