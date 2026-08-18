"""Session brief I -- signatures held in log form, so multiplicities may be huge.

A signature is stored as two real vectors

    c_i = log(multiplicity of the i-th distinct atom)   (>= 0)
    x_i = log(the i-th distinct atom value)             (>= 0)

so that   F(beta) = log Z(beta) = logsumexp_i ( c_i + beta x_i ),
          R = F(0) = logsumexp(c),   Lambda = max x,
          U(theta) = log F(e^theta),
          d(a,b) = osc over theta in [-inf, inf] of U_b - U_a.

This is exactly common.Sig with mults = exp(c); keeping c in log form lets the
multiplicities be as large as the construction wants (they are e^{O(10^3)} in
the C_4 witness) without any floating-point strain, and rounding exp(c) to an
integer perturbs c by less than e^{-c}, i.e. by nothing.
"""
from __future__ import annotations

import math

import numpy as np


class LogSig:
    __slots__ = ("c", "x")

    def __init__(self, c, x):
        c = np.asarray(c, float)
        x = np.asarray(x, float)
        o = np.argsort(-x)
        self.c, self.x = c[o], x[o]

    @property
    def R(self):
        m = self.c.max()
        return m + math.log(float(np.exp(self.c - m).sum()))

    @property
    def Lam(self):
        return float(self.x.max())

    @property
    def sigma(self):
        return math.log(self.R) - math.log(self.Lam)

    def F(self, beta):
        beta = np.atleast_1d(np.asarray(beta, float))
        v = self.c[None, :] + np.multiply.outer(beta, self.x)
        m = v.max(axis=1, keepdims=True)
        return (m[:, 0] + np.log(np.exp(v - m).sum(axis=1)))

    def U(self, theta):
        return np.log(self.F(np.exp(np.asarray(theta, float))))


def osc(a: LogSig, b: LogSig, lo=-40.0, hi=40.0, step=2e-3, refine=True):
    """(max, min, argmax, argmin) of U_b - U_a; grid + parabolic refinement.

    |(U_b - U_a)'| <= 1 (structure theorem S2), so a grid of step h brackets the
    global extremum within h/2 and the parabolic refinement is O(h^3).
    """
    g = np.arange(lo, hi + 0.5 * step, step)
    phi = b.U(g) - a.U(g)
    e0 = math.log(b.R) - math.log(a.R)
    e1 = math.log(b.Lam) - math.log(a.Lam)

    def f(t):
        return float(b.U(np.array([t]))[0] - a.U(np.array([t]))[0])

    def one(idx, want_max):
        if 0 < idx < len(g) - 1 and refine:
            y0, y1, y2 = phi[idx - 1], phi[idx], phi[idx + 1]
            den = y0 - 2 * y1 + y2
            if (want_max and den < 0) or ((not want_max) and den > 0):
                dl = 0.5 * (y0 - y2) / den
                if abs(dl) <= 1.0:
                    t = float(g[idx] + dl * step)
                    return t, f(t)
        return float(g[idx]), float(phi[idx])

    ti, vi = one(int(np.argmax(phi)), True)
    tj, vj = one(int(np.argmin(phi)), False)
    hi_v, s_hi = (vi, ti) if vi >= max(e0, e1) else (max(e0, e1),
                                                     -math.inf if e0 >= e1
                                                     else math.inf)
    lo_v, s_lo = (vj, tj) if vj <= min(e0, e1) else (min(e0, e1),
                                                     -math.inf if e0 <= e1
                                                     else math.inf)
    return hi_v, lo_v, s_hi, s_lo


def d_of(a: LogSig, b: LogSig, **kw):
    hi, lo, _, _ = osc(a, b, **kw)
    return hi - lo


def dmatrix(sigs, **kw):
    n = len(sigs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = d_of(sigs[i], sigs[j], **kw)
    return D


def to_integers(sigs, target_log=60.0):
    """Exact integer (multiplicity, value) data.

    Multiplicities: m_i = floor(e^{c_i}) -- changes c_i by < e^{-c_i}.
    Values: raise every atom of every signature to the common power p, which is
    the exact reparametrisation beta -> beta/p and leaves d invariant, then
    round; that changes log a_i by at most 1/(2 a_i^p) <= e^{-target_log}/2.
    """
    xmin = min(min(v for v in s.x if v > 0) for s in sigs)
    p = target_log / xmin
    out = []
    for s in sigs:
        rows = []
        for ci, xi in zip(s.c, s.x):
            m = int(math.floor(math.exp(ci))) if ci < 700 else None
            if m is None:                       # keep it exact as 10**k * ...
                m = int(math.floor(math.exp(math.log(10) *
                                            (ci / math.log(10) % 1.0))
                                   * 10 ** int(ci / math.log(10))))
            v = 1 if xi <= 0 else int(round(math.exp(min(p * xi, 5000.0))))
            rows.append((max(m, 1), max(v, 1)))
        out.append(rows)
    return out, p
