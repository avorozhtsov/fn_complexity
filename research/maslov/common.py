#!/usr/bin/env python3
"""Shared machinery for session brief M — temperature as a measure.

Conventions
-----------
``s = log beta``.  ``Z_a(beta) = sum_i a_i^beta``, ``F_a = log Z_a``,
``u_a(s) = log F_a(e^s)``.  For an ordered pair ``(a, b)`` put

    f_ab(s) = u_a(s) - u_b(s)

**and this is the sign convention of ``lambda_seed.py``**, i.e. the reverse of
brief D / brief E, where ``A(a,b) = mid(u_b - u_a)``.  Every quantity this
session reports -- critical lambdas, curl fractions, order-agreement fractions,
cycle counts -- is invariant under the flip, which only reverses the
orientation of every cycle and the sign of every potential.  The flip is
applied explicitly in :func:`hodge` so that the reported potential matches the
brief-E one.

The family
----------
    softmax_L(f) = (1/L) log INT e^{Lf} rho
    softmin_L(f) = -(1/L) log INT e^{-Lf} rho
    A_L(a,b)     = (softmax_L(f_ab) + softmin_L(f_ab)) / 2

with ``A_inf = mid = (max+min)/2`` (the framework's comparison) and
``A_0 = INT f_ab rho = Psi(a) - Psi(b)``, ``Psi(a) = INT u_a rho``.

Endpoints
---------
``f_ab`` extends continuously to ``s = -inf`` and ``s = +inf``:

    f_ab(-inf) = log log r_a - log log r_b        (r = number of fibers)
    f_ab(+inf) = log log M_a - log log M_b        (M = largest fiber)

so the *tropical* comparison must include those two values; a soft comparison
must not, because ``rho`` is a probability measure on the finite ``s``-line
(unless one deliberately puts atoms at the ends -- see ``ENDPOINT_PRIOR``).
"""

from __future__ import annotations

import math

import numpy as np

TIE = 1e-10  # the package's tie threshold, per the brief


# --------------------------------------------------------------------- curves


def u_of(sig, s):
    """``u_a(s) = log log Z_a(e^s)``, computed stably.  Same code as the seed."""

    beta = np.exp(np.asarray(s, dtype=float))
    logs = np.log(np.array(sig, dtype=float))
    z = np.outer(beta, logs)
    m = z.max(axis=1)
    return np.log(m + np.log(np.exp(z - m[:, None]).sum(axis=1)))


def u_matrix(sigs, s):
    """``U[i, k] = u_{sigs[i]}(s[k])`` plus the two exact endpoint columns."""

    s = np.asarray(s, dtype=float)
    U = np.empty((len(sigs), s.size))
    for i, a in enumerate(sigs):
        U[i] = u_of(a, s)
    u_minus = np.array([math.log(math.log(len(a))) for a in sigs])
    u_plus = np.array([math.log(math.log(max(a))) for a in sigs])
    return U, u_minus, u_plus


def sigma_of(sig):
    """``sigma = log(R/Lambda) = log(log r / log M)``."""

    return math.log(math.log(len(sig)) / math.log(max(sig)))


def psi_endpoint(sig):
    """``psi = 1/2 log(log r * log M) = 1/2 log phi`` -- brief D(c)/E potential."""

    return 0.5 * math.log(math.log(len(sig)) * math.log(max(sig)))


# ------------------------------------------------------------- soft functionals


def softmax(f, w, lam):
    """``(1/lam) log sum_k w_k e^{lam f_k}``, overflow-free.  ``lam > 0``."""

    f = np.asarray(f, dtype=float)
    m = f.max(axis=-1, keepdims=True)
    return (np.log((w * np.exp(lam * (f - m))).sum(axis=-1)) + lam * m[..., 0]) / lam


def softmin(f, w, lam):
    """``-(1/lam) log sum_k w_k e^{-lam f_k}``, overflow-free."""

    return -softmax(-np.asarray(f, dtype=float), w, lam)


def soft_mid(f, w, lam):
    """``A_lam`` of one function; ``lam=None`` gives the tropical midrange."""

    f = np.asarray(f, dtype=float)
    if lam is None:
        return 0.5 * (f.max(axis=-1) + f.min(axis=-1))
    return 0.5 * (softmax(f, w, lam) + softmin(f, w, lam))


# --------------------------------------------------------------- comparisons


def tropical_matrix(U, u_minus, u_plus):
    """``A_inf[i, j] = mid_s(u_i - u_j)`` with the two endpoints included."""

    n = U.shape[0]
    A = np.empty((n, n))
    for i in range(n):
        D = U[i][None, :] - U  # (n, grid): f_ij
        e0 = u_minus[i] - u_minus
        ei = u_plus[i] - u_plus
        hi = np.maximum(D.max(axis=1), np.maximum(e0, ei))
        lo = np.minimum(D.min(axis=1), np.minimum(e0, ei))
        A[i] = 0.5 * (hi + lo)
    A = 0.5 * (A - A.T)
    np.fill_diagonal(A, 0.0)
    return A


def soft_matrix(U, w, lam):
    """``A_lam[i, j]`` for a prior ``w`` supported on the grid of ``U``.

    ``lam=None`` gives the midrange **over the support of the prior only**
    (no endpoint atoms); that is the correct ``lam -> infinity`` limit of the
    soft family for a compactly supported prior.
    """

    n = U.shape[0]
    A = np.empty((n, n))
    for i in range(n):
        D = U[i][None, :] - U
        if lam is None:
            A[i] = 0.5 * (D.max(axis=1) + D.min(axis=1))
        else:
            A[i] = 0.5 * (softmax(D, w, lam) + softmin(D, w, lam))
    A = 0.5 * (A - A.T)
    np.fill_diagonal(A, 0.0)
    return A


def hodge(A):
    """Least-squares gradient/curl split of an antisymmetric ``A``.

    ``psi_opt`` is reported in the brief-E orientation (``-rowmean`` of the
    brief-E ``A``, which is ``+rowmean`` of ours)."""

    psi = A.mean(axis=1)
    psi = psi - psi.mean()
    grad = psi[:, None] - psi[None, :]
    curl = A - grad
    return {
        "psi": psi,
        "grad": grad,
        "curl": curl,
        "grad_energy": float((grad**2).sum() / (A**2).sum()),
        "curl_energy": float((curl**2).sum() / (A**2).sum()),
        "curl_fraction": float(np.sqrt((curl**2).sum() / (A**2).sum())),
        "curl_sup": float(np.abs(curl).max()),
    }


def three_cycles(A, tie=TIE):
    """Number of directed 3-cycles of the tournament ``sign(A)``."""

    T = (A > tie).astype(np.float64)
    return int(round(float(np.einsum("ij,jk,ki->", T, T, T)) / 3.0))


def order_agreement(scalar, A, tie=TIE):
    """Fraction of ordered pairs on which ``scalar`` reproduces ``sign(A)``.

    Ties in either object count as failures, exactly as brief E counts them.
    ``scalar`` is compared in the same orientation as ``A``: ``A[i,j] > 0``
    should mean ``scalar[i] > scalar[j]``.
    """

    x = np.asarray(scalar, dtype=float)
    d = x[:, None] - x[None, :]
    n = A.shape[0]
    iu = np.triu_indices(n, 1)
    a, b = A[iu], d[iu]
    decided = np.abs(a) > tie
    right = decided & (np.abs(b) > 0) & (np.sign(a) == np.sign(b))
    return float(right.sum()) / float(a.size)


# ------------------------------------------------------------------- priors


def uniform_prior(lo, hi, n):
    s = np.linspace(lo, hi, n)
    w = np.ones(n) / n
    return s, w


def gaussian_prior(mu, sd, n, halfwidth=8.0):
    s = np.linspace(mu - halfwidth * sd, mu + halfwidth * sd, n)
    w = np.exp(-0.5 * ((s - mu) / sd) ** 2)
    w /= w.sum()
    return s, w


def logistic_prior(mu, scale, n, halfwidth=20.0):
    s = np.linspace(mu - halfwidth * scale, mu + halfwidth * scale, n)
    z = (s - mu) / scale
    w = np.exp(-z) / (1.0 + np.exp(-z)) ** 2
    w /= w.sum()
    return s, w


# ------------------------------------------------------------------- pools


CYCLE = [(6, 3, 3), (7, 2, 1), (6, 5, 1)]


def f11_pool():
    """The 296 genus-two pencil signatures over ``F_11`` of brief E.

    Rebuilt here from scratch (third independent implementation) so nothing is
    imported from a cached file.
    """

    return _pencil_pool(11)


def _pencil_pool(q):
    import itertools

    squares = {(x * x) % q for x in range(1, q)}
    chi = np.array([0] + [1 if x in squares else -1 for x in range(1, q)], dtype=np.int64)
    u = np.arange(q)
    K = chi[(u[None, :] + u[:, None]) % q]  # K[c, v] = chi(v + c)
    x = np.arange(q, dtype=np.int64)
    sigs = set()
    for degree in (5, 6):
        powers = np.stack([(x**k) % q for k in range(1, degree + 1)])
        coeffs = np.array(list(itertools.product(range(q), repeat=degree - 1)), dtype=np.int64)
        for start in range(0, coeffs.shape[0], 4096):
            block = coeffs[start : start + 4096]
            m = block.shape[0]
            vals = np.tile(powers[degree - 1], (m, 1))
            for j in range(degree - 1):
                vals = vals + block[:, j : j + 1] * powers[j][None, :]
            vals %= q
            flat = (vals + q * np.arange(m)[:, None]).ravel()
            hist = np.bincount(flat, minlength=q * m).reshape(m, q)
            N = q + hist @ K.T
            N = np.sort(N, axis=1)[:, ::-1]
            for row in N:
                if row.min() > 0:
                    sigs.add(tuple(int(v) for v in row))
    return sorted(sigs)


# ------------------------------------------------------------------- mpmath


def mp_u(sig, s, dps=40):
    import mpmath as mp

    with mp.workdps(dps):
        beta = mp.exp(mp.mpf(s))
        z = mp.fsum([mp.power(mp.mpf(int(a)), beta) for a in sig])
        return mp.log(mp.log(z))


def mp_soft_mid(sig_a, sig_b, s, w, lam, dps=40):
    """40-digit ``A_lam`` for one pair on a given grid."""

    import mpmath as mp

    with mp.workdps(dps):
        f = [mp_u(sig_a, sk, dps) - mp_u(sig_b, sk, dps) for sk in s]
        if lam is None:
            return mp.mpf(max(f)) + mp.mpf(min(f)) >> 1
        L = mp.mpf(lam)
        fmax, fmin = max(f), min(f)
        smax = (mp.log(mp.fsum([wk * mp.exp(L * (fk - fmax)) for wk, fk in zip(w, f)])) + L * fmax) / L
        smin = -(mp.log(mp.fsum([wk * mp.exp(-L * (fk - fmin)) for wk, fk in zip(w, f)])) - L * fmin) / L
        return (smax + smin) / 2
