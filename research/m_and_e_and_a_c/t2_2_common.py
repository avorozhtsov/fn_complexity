#!/usr/bin/env python3
"""Shared machinery for T2.2: what the exchange rates know about {a_c}.

For a map ``f : A^2 -> A^1`` over ``F_q`` the *signature* is the multiset of
non-zero fiber counts ``N_c = #f^{-1}(c)``.  Writing ``a_c = q - N_c`` for the
Frobenius trace of the fiber, the normalised trace moments are

    m_k = (1/q) * sum_c a_c^k / q^{k/2}

so that ``m_k = O(1)`` under Hasse-Weil.  Exactly ``sum_c N_c = q^2`` forces
``m_1 = 0`` whenever every fiber is non-empty.

The partition function is ``Z_f(beta) = sum_c N_c^beta`` and the exchange rate
is ``C(g -> f) = inf_{beta in [0, inf]} log Z_g(beta) / log Z_f(beta)``.

Two facts drive everything in this note and are re-derived here:

*   ``log Z_f(beta) = (beta+1) log q + C(beta,2) m_2 / q
                      - C(beta,3) m_3 / q^{3/2} + O(q^{-2})``
    for maps all of whose fibers are non-empty (``C(beta,k)`` = binomial).
*   Against the flat reference ``L = (q,...,q)`` the ratio
    ``log Z_f / log Z_L = 1 + C(beta,2) m_2 / ((beta+1) q log q) + ...``
    is minimised at ``beta -> sqrt(2) - 1`` with value
    ``1 - (3 - 2 sqrt 2) m_2 / (2 q log q)``.

Signatures here have length ~q, so the pure-Python solver in ``fn_complexity``
is used only to *certify* headline numbers; bulk work uses the compressed
(value, multiplicity) vectorisation below, which agrees with it to ~1e-14.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate_result  # noqa: E402

SQRT2 = math.sqrt(2.0)
BETA_STAR = SQRT2 - 1.0          # 0.41421356...
KAPPA = 3.0 - 2.0 * SQRT2        # 0.17157288... = -min_beta C(beta,2)*2/(beta+1)


# --------------------------------------------------------------- signatures


@dataclass(frozen=True)
class Sig:
    """A compressed signature: distinct fiber sizes with multiplicities."""

    values: np.ndarray   # int64, strictly decreasing, all >= 1
    mults: np.ndarray    # int64, positive
    name: str = ""

    @property
    def n_fibers(self) -> int:
        return int(self.mults.sum())

    @property
    def total(self) -> int:
        return int((self.values * self.mults).sum())

    @property
    def max_fiber(self) -> int:
        return int(self.values[0])

    def expand(self) -> tuple[int, ...]:
        return tuple(int(v) for v, m in zip(self.values, self.mults) for _ in range(m))

    def log_z(self, beta: np.ndarray) -> np.ndarray:
        """``log Z(beta)`` evaluated on an array of betas (stable)."""
        lv = np.log(self.values.astype(float))
        lm = np.log(self.mults.astype(float))
        e = np.outer(np.asarray(beta, dtype=float), lv) + lm[None, :]
        return _logsumexp(e, axis=1)


def _logsumexp(e: np.ndarray, axis: int) -> np.ndarray:
    top = np.max(e, axis=axis, keepdims=True)
    return (top + np.log(np.sum(np.exp(e - top), axis=axis, keepdims=True))).squeeze(axis)


def sig_from_counts(counts: np.ndarray, name: str = "") -> Sig:
    """Compress a raw ``bincount`` array (zeros dropped: fibers must be non-empty)."""
    nz = counts[counts > 0]
    vals, mult = np.unique(nz, return_counts=True)
    order = np.argsort(-vals)
    return Sig(vals[order].astype(np.int64), mult[order].astype(np.int64), name)


def sig_from_list(pairs: list[tuple[int, int]], name: str = "") -> Sig:
    vals = np.array([p[0] for p in pairs], dtype=np.int64)
    mult = np.array([p[1] for p in pairs], dtype=np.int64)
    order = np.argsort(-vals)
    return Sig(vals[order], mult[order], name)


# ------------------------------------------------------------------- rates


def _ratio(g: Sig, f: Sig, beta: np.ndarray) -> np.ndarray:
    return g.log_z(beta) / f.log_z(beta)


def _horizon(g: Sig, f: Sig) -> float:
    """Beta beyond which every non-maximal Gibbs weight is below 3e-16.

    Same rule as ``fn_complexity._search_horizon``: past this point both ratios
    have saturated to their ``beta = inf`` values, so scanning further only
    manufactures floating-point plateaus.
    """
    gaps: list[float] = []
    for s in (g, f):
        top = math.log(s.max_fiber)
        gaps.extend(top - math.log(v) for v in s.values[1:])
    if not gaps:
        return 64.0
    return min(36.0 / min(gaps), 1.0e7)


def rate(g: Sig, f: Sig, grid: int = 4096) -> tuple[float, float]:
    """Return ``(C(g -> f), argmin beta)``, matching ``fn_complexity``.

    ``g`` is the implementer, ``f`` the implemented signature, exactly as in
    ``exchange_rate(implementer, implemented)``.
    """
    # endpoints
    cand: list[tuple[float, float]] = []
    lz_g0 = math.log(g.n_fibers)
    lz_f0 = math.log(f.n_fibers)
    cand.append((lz_g0 / lz_f0 if lz_f0 > 0 else math.inf, 0.0))
    lmg, lmf = math.log(g.max_fiber), math.log(f.max_fiber)
    cand.append((lmg / lmf if lmf > 0 else math.inf, math.inf))

    beta_max = _horizon(g, f)
    betas = np.concatenate(([0.0], np.geomspace(min(1e-10, beta_max / grid),
                                               beta_max, grid)))
    vals = _ratio(g, f, betas)
    # Strict on the left so a floating-point plateau yields one candidate, not
    # one per grid point.
    interior = np.where((vals[1:-1] < vals[:-2]) & (vals[1:-1] <= vals[2:]))[0] + 1
    for i in interior:
        b, v = _golden(lambda x: float(_ratio(g, f, np.array([x]))[0]),
                       float(betas[i - 1]), float(betas[i + 1]))
        cand.append((v, b))
    v, b = min(cand, key=lambda t: t[0])
    return v, b


def _golden(fn, lo: float, hi: float, tol: float = 1e-13) -> tuple[float, float]:
    inv = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = hi - inv * (hi - lo)
    x2 = lo + inv * (hi - lo)
    y1, y2 = fn(x1), fn(x2)
    for _ in range(200):
        if hi - lo <= tol * (1.0 + abs(x1) + abs(x2)):
            break
        if y1 <= y2:
            hi, x2, y2 = x2, x1, y1
            x1 = hi - inv * (hi - lo)
            y1 = fn(x1)
        else:
            lo, x1, y1 = x1, x2, y2
            x2 = lo + inv * (hi - lo)
            y2 = fn(x2)
    return (x1, y1) if y1 <= y2 else (x2, y2)


def classify(g: Sig, f: Sig, tol: float = 1e-12) -> str:
    """Say whether ``C(g -> f)`` is an endpoint value or a genuine interior tangency.

    A reported finite ``beta`` can still be spurious: once ``beta`` is past the
    saturation horizon both partition functions are affine, the ratio is flat to
    machine precision and the "interior minimum" equals the ``beta = inf`` value.
    Such probes carry no information beyond the endpoint.
    """
    v, b = rate(g, f)
    v0 = math.log(g.n_fibers) / math.log(f.n_fibers)
    vinf = math.log(g.max_fiber) / math.log(f.max_fiber)
    if v >= min(v0, vinf) - tol:
        return "beta=0" if v0 <= vinf else "beta=inf"
    return "interior"


def certified_rate(g: Sig, f: Sig) -> tuple[float, float]:
    """Exact repo solver (slow, O(len) per beta); used to certify headline numbers."""
    r = exchange_rate_result(g.expand(), f.expand())
    return r.rate, r.beta


# ----------------------------------------------------------------- moments


def moments(counts: np.ndarray, q: int, kmax: int = 6) -> dict[str, float]:
    """Normalised trace moments of a full length-q fiber-count vector."""
    a = (q - counts).astype(float)
    out: dict[str, float] = {}
    for k in range(1, kmax + 1):
        out[f"m{k}"] = float((a ** k).sum()) / (q * q ** (k / 2.0))
    out["amax"] = float(a.max())          # most negative trace direction: smallest fiber
    out["amin"] = float(a.min())          # = q - N_max
    out["a_absmax"] = float(np.abs(a).max())
    out["n_image"] = int((counts > 0).sum())
    out["max_fiber"] = int(counts.max())
    return out


# -------------------------------------------------------- reference family


def references(q: int) -> dict[str, Sig]:
    """A small fixed reference family over ``F_q`` (all have q^2 points)."""
    refs: dict[str, Sig] = {}
    refs["L"] = sig_from_list([(q, q)], "L")                                  # flat
    refs["Xsplit"] = sig_from_list([(2 * q - 1, 1), (q - 1, q - 1)], "Xsplit")  # xy
    refs["Xaniso"] = sig_from_list([(q + 1, q - 1), (1, 1)], "Xaniso")          # x^2-dy^2
    refs["Sq"] = sig_from_list([(2 * q, (q - 1) // 2), (q, 1)], "Sq")           # f = x^2
    if q % 3 == 1:
        refs["Cb"] = sig_from_list([(3 * q, (q - 1) // 3), (q, 1)], "Cb")       # f = x^3
    return refs


def flat(n: int, m: int, name: str = "") -> Sig:
    """The flat reference with ``n`` fibers of size ``m``: ``log Z = log n + beta log m``."""
    return Sig(np.array([m], dtype=np.int64), np.array([n], dtype=np.int64),
               name or f"flat({n},{m})")


# ------------------------------------------------------- asymptotic model


def predicted_rate_f_to_L(q: int, m2: float, m3: float = 0.0) -> float:
    """``C(f -> L)`` to order ``q^{-3/2}`` for a map with all fibers non-empty.

    ``C(f -> L) = inf_beta [ (beta+1) log q + C(beta,2) m2/q - C(beta,3) m3/q^{3/2} ]
                          / [ (beta+1) log q ]``.
    """
    b = np.linspace(1e-6, 4.0, 400001)
    num = (b + 1) * math.log(q) + b * (b - 1) / 2 * m2 / q \
        - b * (b - 1) * (b - 2) / 6 * m3 / q ** 1.5
    den = (b + 1) * math.log(q)
    return float(np.min(num / den))
