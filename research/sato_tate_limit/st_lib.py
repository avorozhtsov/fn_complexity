#!/usr/bin/env python3
"""Sato--Tate trace measures, their cumulant generating functions, and ``Psi``.

The object of session brief F is

    Psi_mu(tau) = K_mu(tau)/tau,      K_mu(tau) = log E_mu[e^{tau alpha}],

for ``mu`` the limiting distribution of the normalised Frobenius traces
``alpha_c = (N_c - q)/sqrt(q)`` of a family of curves over ``F_q``.  Two
boundary values are forced:

    Psi_mu(0)   = K'_mu(0) = E_mu[alpha] = 0        (sum_c a_c = 0 identically)
    Psi_mu(inf) = ess sup supp(mu) = alpha_max

and ``Psi_mu`` is non-decreasing on ``(0, inf)``: ``K_mu`` is convex with
``K_mu(0) = 0``, so ``K(tau)/tau`` is the slope of a chord from the origin.

Everything here comes from the **Weyl integration formula**, in the form that
turns the ``g``-dimensional integral into a ``g x g`` determinant.  For a
classical compact group of rank ``N`` whose eigenvalues are ``e^{+-i theta_j}``
together with fixed eigenvalues contributing ``epsilon`` to the trace, Weyl's
formula in the variables ``x_j = cos theta_j`` is

    d mu  proportional to  prod_{i<j} (x_i - x_j)^2  prod_j w(x_j) dx_j

with ``w`` a Jacobi weight depending on the group.  Andreief's identity then
gives, exactly,

    E[e^{tau tr}] = e^{tau epsilon} det( L_tau[x^{i+j}] ) / det( L_0[x^{i+j}] ),
    L_tau[p] = integral p(x) e^{2 tau x} w(x) dx,        0 <= i, j < N.

The three weights that occur are handled in the Chebyshev basis, where every
moment is a Bessel function of argument ``2 tau``:

    Sp(2N):      w = (2/pi) sqrt(1-x^2)        L[T_k] = I_k - (I_{k-2}+I_{k+2})/2
    SO(2N):      w = (1/pi) / sqrt(1-x^2)      L[T_k] = I_k
    SO(2N+1):    w = (1/pi) (1-x)/sqrt(1-x^2)  L[T_k] = I_k - (I_{k-1}+I_{k+1})/2

``validate_library.py`` checks the resulting MGFs against brute-force Weyl
quadrature in ``g`` dimensions and against exactly known moments.

**Composition rules, all exact.**

* *product* (a Jacobian that splits into independent factors): the MGFs
  multiply, so ``K`` and hence ``Psi`` **add**, and ``alpha_max``, the variance
  and the edge exponent all add;
* *multiplicity* (``k`` isogenous copies of one factor, ``alpha = k alpha_G``):
  ``K(tau) = K_G(k tau)``, so ``alpha_max`` scales by ``k`` and the variance by
  ``k^2`` while the edge exponent is unchanged;
* *atom* (a proportion ``p`` of fibres degenerating to ``alpha = 0``):
  ``K(tau) = log(p + (1-p) e^{K_G(tau)})``, which lowers ``Psi`` pointwise.

**Edge data.**  If ``P(alpha_max - alpha < eps) ~ c eps^t`` then

    M(tau) ~ A e^{tau alpha_max} tau^{-t},
    Psi(tau) = alpha_max - (t log tau - log A)/tau + O(tau^{-2}),

with ``t`` and ``log A`` additive over products.  For a rank-``N`` Weyl measure
with weight behaving like ``(1-x)^a`` at the top edge, the density
``prod (u_i-u_j)^2 prod u_i^a`` is homogeneous of degree ``N(N-1) + aN`` in
``u_i = 1-x_i``, and there are ``N`` differentials, so

    t = N^2 + a N,

which is ``g(2g+1)/2 = dim USp(2g)/2`` for ``Sp(2g)`` (``a = 1/2``), ``N^2-N/2``
for ``SO(2N)`` (``a = -1/2``) and ``N^2+N/2`` for ``SO(2N+1)``.

``t`` decides the large-``tau`` order of two measures with the same
``alpha_max``; the variance decides the small-``tau`` order.  Both are recorded.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from pathlib import Path

import numpy as np
from mpmath import besseli, matrix, mp, mpf

mp.dps = 60

HERE = Path(__file__).resolve().parent
CACHE = HERE / ".kcache"


# ---------------------------------------------------------------- Weyl moments


def _cheb_expansion(n: int) -> list[tuple[int, float]]:
    """``x^n = sum_k c_k T_k(x)``; coefficients are exact in float64 for n <= 25."""

    poly = np.zeros(n + 1)
    poly[n] = 1.0
    coeffs = np.polynomial.chebyshev.poly2cheb(poly)
    return [(k, float(c)) for k, c in enumerate(coeffs) if c != 0.0]


_CHEB: dict[int, list[tuple[int, float]]] = {}


def cheb_expansion(n: int) -> list[tuple[int, float]]:
    if n not in _CHEB:
        _CHEB[n] = _cheb_expansion(n)
    return _CHEB[n]


def _T_moment(kind: str, k: int, two_tau):
    """``L_tau[T_k]`` for the three Weyl weights."""

    if kind == "sp":
        return besseli(k, two_tau) - (
            besseli(abs(k - 2), two_tau) + besseli(k + 2, two_tau)
        ) / 2
    if kind == "so_even":
        return besseli(k, two_tau)
    if kind == "so_odd":
        return besseli(k, two_tau) - (
            besseli(abs(k - 1), two_tau) + besseli(k + 1, two_tau)
        ) / 2
    raise ValueError(kind)


def _x_moment(kind: str, n: int, tau):
    two_tau = 2 * tau
    total = mpf(0)
    for k, c in cheb_expansion(n):
        total += mpf(c) * _T_moment(kind, k, two_tau)
    return total


def working_precision(rank: int, tau: float) -> int:
    """Digits needed for the Hankel determinant at ``tau``.

    As ``tau -> inf`` the measure concentrates at the top edge and the moment
    matrix becomes near-singular: the determinant is smaller than the product of
    its entries by a factor ``~ tau^{-rank(rank-1)}``, so that many digits
    cancel.  Budget them, plus a margin.
    """

    cancel = rank * (rank - 1) * math.log10(max(2.0 * abs(tau), 10.0))
    return int(40 + cancel)


def mgf_classical(kind: str, rank: int, epsilon: int, tau: float):
    """``E[e^{tau tr}]`` for a classical group, via the Andreief determinant."""

    if rank == 0:
        return mp.e ** (mpf(epsilon) * mpf(tau))
    mp.dps = max(mp.dps, working_precision(rank, float(tau)))
    t = mpf(tau)
    num = matrix(rank, rank)
    den = matrix(rank, rank)
    for i in range(rank):
        for j in range(rank):
            num[i, j] = _x_moment(kind, i + j, t)
            den[i, j] = _x_moment(kind, i + j, mpf(0))
    return mp.e ** (mpf(epsilon) * t) * mp.det(num) / mp.det(den)


# --------------------------------------------------------------------- groups


@dataclass(frozen=True)
class Group:
    """One irreducible monodromy block."""

    name: str
    kind: str
    rank: int
    epsilon: int
    alpha_max: float
    variance: float
    tail: float
    realisable: bool
    note: str = ""


def edge_exponent(rank: int, a: float) -> float:
    """``t = N^2 + a N`` -- see the module docstring."""

    return rank * rank + a * rank


GROUPS: dict[str, Group] = {}

for _g in (1, 2, 3, 4, 5, 6):
    GROUPS[f"USp{2 * _g}"] = Group(
        name=f"USp{2 * _g}", kind="sp", rank=_g, epsilon=0,
        alpha_max=2.0 * _g, variance=1.0, tail=edge_exponent(_g, 0.5),
        realisable=True, note="big monodromy, generic hyperelliptic pencil",
    )
for _n in (3, 4, 5, 6, 7):
    if _n % 2 == 0:
        GROUPS[f"SO{_n}"] = Group(
            name=f"SO{_n}", kind="so_even", rank=_n // 2, epsilon=0,
            alpha_max=float(_n), variance=1.0, tail=edge_exponent(_n // 2, -0.5),
            realisable=False, note="orthogonal: not the H^1 of a curve family",
        )
    else:
        GROUPS[f"SO{_n}"] = Group(
            name=f"SO{_n}", kind="so_odd", rank=_n // 2, epsilon=1,
            alpha_max=float(_n), variance=1.0, tail=edge_exponent(_n // 2, 0.5),
            realisable=False, note="orthogonal: not the H^1 of a curve family",
        )

GROUPS["SU2"] = Group(
    "SU2", "sp", 1, 0, 2.0, 1.0, 1.5, True, "non-CM elliptic pencil (Sato-Tate)"
)
GROUPS["U1"] = Group(
    "U1", "so_even", 1, 0, 2.0, 2.0, 0.5, True, "torus monodromy (arcsine)"
)


def _group_K_raw(group: Group, tau: np.ndarray) -> np.ndarray:
    saved = mp.dps
    try:
        if group.kind == "sp" and group.rank == 1:
            return np.array([float(mp.log(besseli(1, 2 * mpf(t)) / mpf(t)))
                             for t in tau])
        if group.kind == "so_even" and group.rank == 1:
            return np.array([float(mp.log(besseli(0, 2 * mpf(t)))) for t in tau])
        out = np.empty(tau.size)
        for i, t in enumerate(tau):
            mp.dps = saved
            out[i] = float(mp.log(mgf_classical(group.kind, group.rank,
                                                group.epsilon, t)))
        return out
    finally:
        mp.dps = saved


def group_K(name: str, tau: np.ndarray) -> np.ndarray:
    """Disk-cached ``log E[e^{tau tr}]`` on a grid, computed at 60 digits."""

    group = GROUPS[name]
    tau = np.ascontiguousarray(np.asarray(tau, dtype=np.float64))
    if tau.size < 64:            # single-point probes: never worth a cache file
        return _group_K_raw(group, tau)
    CACHE.mkdir(exist_ok=True)
    key = hashlib.sha1(name.encode() + tau.tobytes()).hexdigest()[:20]
    path = CACHE / f"{name}_{key}.npy"
    if path.exists():
        return np.load(path)
    out = _group_K_raw(group, tau)
    np.save(path, out)
    return out


# ------------------------------------------------------------------- measures


@dataclass(frozen=True)
class Factor:
    """One isogeny block of the Jacobian.

    ``alpha`` contributed is ``k * tr(g)`` for ``g`` Haar on the group, except
    on a set of fibres of density ``atom`` where the block degenerates and
    contributes ``0`` (supersingular / non-connected monodromy component).
    """

    group: str
    multiplicity: int = 1
    atom: float = 0.0

    @property
    def alpha_max(self) -> float:
        return self.multiplicity * GROUPS[self.group].alpha_max

    @property
    def variance(self) -> float:
        return (1.0 - self.atom) * self.multiplicity ** 2 * GROUPS[self.group].variance

    @property
    def tail(self) -> float:
        return GROUPS[self.group].tail

    def K(self, tau: np.ndarray) -> np.ndarray:
        k = group_K(self.group, self.multiplicity * np.asarray(tau, dtype=np.float64))
        if self.atom > 0.0:
            k = np.logaddexp(math.log(self.atom), math.log1p(-self.atom) + k)
        return k

    def label(self) -> str:
        g = self.group
        if self.atom == 0.5 and g == "U1":
            g = "CM"
        elif self.atom > 0.0:
            g = f"{g}[{self.atom:g}]"
        return g if self.multiplicity == 1 else f"{self.multiplicity}.{g}"


@dataclass(frozen=True)
class Measure:
    """A product of factors -- a Jacobian that splits into independent blocks.

    The MGFs of independent factors multiply, so ``K`` and ``Psi`` **add**.
    """

    factors: tuple[Factor, ...]
    name: str = ""

    @property
    def label(self) -> str:
        return self.name or " x ".join(f.label() for f in self.factors)

    @property
    def alpha_max(self) -> float:
        return sum(f.alpha_max for f in self.factors)

    @property
    def variance(self) -> float:
        return sum(f.variance for f in self.factors)

    @property
    def tail(self) -> float:
        return sum(f.tail for f in self.factors)

    @property
    def realisable(self) -> bool:
        return all(GROUPS[f.group].realisable for f in self.factors)

    def K(self, tau: np.ndarray) -> np.ndarray:
        total = np.zeros_like(np.asarray(tau, dtype=np.float64))
        for f in self.factors:
            total = total + f.K(tau)
        return total

    def Psi(self, tau: np.ndarray) -> np.ndarray:
        return self.K(tau) / np.asarray(tau, dtype=np.float64)


def mixture_K(weights: np.ndarray, atom_K: np.ndarray) -> np.ndarray:
    """``K`` of a convex combination of measures, from their ``K`` rows.

    ``E[e^{tau alpha}]`` is linear in the measure, so the MGFs mix; the CGFs
    combine by ``log-sum-exp``.  Components of zero weight are dropped.
    """

    w = np.asarray(weights, dtype=np.float64)
    mask = w > 0.0
    shifted = np.log(w[mask])[:, None] + atom_K[mask]
    m = shifted.max(axis=0)
    return m + np.log(np.exp(shifted - m).sum(axis=0))


# ------------------------------------------------------------ the comparison


def mid_and_extrema(psi_a: np.ndarray, psi_b: np.ndarray,
                    amax_a: float, amax_b: float) -> tuple[float, float, float]:
    """``(mid, sup, inf)`` of ``D = Psi_a - Psi_b`` over ``[0, inf]``.

    The two endpoints are supplied analytically: ``D(0) = 0`` because both
    measures have mean zero, and ``D(inf) = alpha_max(a) - alpha_max(b)``.
    """

    d = np.concatenate([[0.0], psi_a - psi_b, [amax_a - amax_b]])
    hi = float(d.max())
    lo = float(d.min())
    return 0.5 * (hi + lo), hi, lo


def sign_changes(d_interior: np.ndarray, tol: float | None = None) -> int:
    """Number of sign changes of ``D`` on the open interval ``(0, inf)``.

    The default tolerance is relative to the size of ``D``: a pointwise
    domination whose infimum touches zero must not be counted as a crossing.
    """

    if tol is None:
        tol = 1e-9 * max(float(np.abs(d_interior).max()), 1e-300)
    s = np.sign(np.where(np.abs(d_interior) < tol, 0.0, d_interior))
    s = s[s != 0]
    if s.size == 0:
        return 0
    return int(np.count_nonzero(np.diff(s) != 0))


def tau_grid(low: float = 1e-4, high: float = 1e4, points: int = 1600) -> np.ndarray:
    return np.geomspace(low, high, points)
