"""Session brief I -- the exchange metric as a Hilbert projective metric.

THE REFORMULATION (proved in OBSTRUCTION.md Sec. 2).

    d(a,b) = osc_beta log( F_b(beta) / F_a(beta) )
           = log sup(F_b/F_a) + log sup(F_a/F_b)

is exactly the *Hilbert projective metric* between F_a = log Z_a and
F_b = log Z_b, taken over beta in [0, infinity] (the endpoint beta = infinity
compares the asymptotic slopes Lambda).  Cartesian powers act as F -> kF, which
is precisely the projective rescaling that a Hilbert metric quotients out --
this reproves FINDINGS Sec. 1.6.

The set of achievable F is  { log sum_i m_i e^{beta x_i} : m_i >= 1 integers,
x_i = log a_i >= 0 }.  Its closure under positive rescaling and locally uniform
limits is the *tropical cone*

    C = { Phi : [0,inf) -> (0,inf) convex, nondecreasing, Phi(beta) >= Lambda*beta }

with Lambda = lim Phi'.  The proof is the tropical (zero-temperature) limit
    (1/K) log sum_i m_i^K e^{K beta x_i}  ->  max_i ( log m_i + beta x_i ),
uniform on [0,inf) with error <= log(#lines)/K, so distances converge with the
same O(1/K).  Since  Phi = max_i (c_i + beta x_i)  with c_i = log m_i >= 0, the
condition Phi >= Lambda beta is exactly the "c_i >= 0" that integrality forces,
and it is precisely the lower sandwich  F >= max(R, beta Lambda)  of FINDINGS
Sec. 1.1.

This module implements the tropical cone exactly: Phi is a finite max of lines,
Phi_b/Phi_a is monotone between consecutive breakpoints (ratio of affine
functions), so the oscillation is a MAXIMUM OVER A FINITE SET -- no grid, no
Lipschitz bracket, exact to machine precision.
"""
from __future__ import annotations

import math

import numpy as np


class Trop:
    """Phi(beta) = max_j (c_j + beta x_j), c_j >= 0, x_j >= 0, max x_j > 0."""

    __slots__ = ("c", "x", "bps")

    def __init__(self, c, x):
        c = np.asarray(c, dtype=float)
        x = np.asarray(x, dtype=float)
        if c.min() < -1e-15 or x.min() < -1e-15:
            raise ValueError("c, x must be >= 0")
        c = np.maximum(c, 0.0)
        x = np.maximum(x, 0.0)
        if x.max() <= 0.0:
            raise ValueError("need a positive slope (excludes all-ones)")
        if c.max() <= 0.0:
            raise ValueError("need a positive intercept (excludes r = 1)")
        # keep only the lines on the upper envelope, in increasing slope order
        order = np.lexsort((-c, x))
        c, x = c[order], x[order]
        keep_c, keep_x = [], []
        for ci, xi in zip(c, x):
            if keep_x and xi == keep_x[-1]:
                continue                      # same slope, lower intercept
            while len(keep_x) >= 1:
                if len(keep_x) == 1:
                    if ci >= keep_c[-1]:      # dominates outright
                        keep_c.pop(); keep_x.pop()
                        continue
                    break
                # breakpoint of last two vs breakpoint with the new line
                b_old = (keep_c[-2] - keep_c[-1]) / (keep_x[-1] - keep_x[-2])
                b_new = (keep_c[-1] - ci) / (xi - keep_x[-1])
                if b_new <= b_old:
                    keep_c.pop(); keep_x.pop()
                else:
                    break
            if len(keep_x) == 1 and ci >= keep_c[-1]:
                keep_c.pop(); keep_x.pop()
            keep_c.append(ci); keep_x.append(xi)
        self.c = np.array(keep_c)
        self.x = np.array(keep_x)
        bps = []
        for i in range(len(self.c) - 1):
            bps.append((self.c[i] - self.c[i + 1]) / (self.x[i + 1] - self.x[i]))
        self.bps = np.array(bps) if bps else np.zeros(0)

    # ---- basic quantities ------------------------------------------------
    def val(self, beta):
        beta = np.asarray(beta, dtype=float)
        return (self.c[None, :] + np.multiply.outer(beta, self.x)).max(axis=-1)

    @property
    def Lam(self):
        return float(self.x.max())

    @property
    def R(self):
        return float(self.c.max())

    @property
    def sigma(self):
        return math.log(self.R) - math.log(self.Lam)

    def normalised(self):
        """Projective representative with Phi(0) = 1."""
        return Trop(self.c / self.R, self.x / self.R)


def candidates(a: Trop, b: Trop):
    """Every beta at which log(Phi_b/Phi_a) can have an interior extremum."""
    bs = np.concatenate([a.bps, b.bps])
    bs = bs[bs > 0]
    return np.unique(np.concatenate([[0.0], bs]))


def phi(a: Trop, b: Trop, beta):
    return np.log(b.val(beta)) - np.log(a.val(beta))


def hilbert(a: Trop, b: Trop):
    """d(a,b) = osc over [0, inf] of log(Phi_b/Phi_a).  Exact."""
    bs = candidates(a, b)
    v = phi(a, b, bs)
    e_inf = math.log(b.Lam) - math.log(a.Lam)
    hi = max(float(v.max()), e_inf)
    lo = min(float(v.min()), e_inf)
    return hi - lo


def mid(a: Trop, b: Trop):
    bs = candidates(a, b)
    v = phi(a, b, bs)
    e_inf = math.log(b.Lam) - math.log(a.Lam)
    hi = max(float(v.max()), e_inf)
    lo = min(float(v.min()), e_inf)
    return 0.5 * (hi + lo)


def dmatrix(fs):
    n = len(fs)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = hilbert(fs[i], fs[j])
    return D


def amatrix(fs):
    n = len(fs)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = mid(fs[i], fs[j])
            A[j, i] = -A[i, j]
    return A


# ---------------------------------------------------------------------------
# real signatures realising a tropical point, and the O(1/K) convergence
# ---------------------------------------------------------------------------

def to_signature(t: Trop, K: float, scale: float = 1.0):
    """Integer signature whose (1/K)*F approximates ``scale * t``.

    Line (c, x) becomes  round(exp(K*scale*c)) copies of the atom
    round(exp(K*scale*x)).  Both are exact integers; the error in F is at most
    log(#lines) + rounding, hence O(1/K) in (1/K)F.
    """
    vals, mults = [], []
    for ci, xi in zip(t.c, t.x):
        m = int(round(math.exp(K * scale * ci)))
        v = int(round(math.exp(K * scale * xi)))
        if m < 1:
            m = 1
        if v < 1:
            v = 1
        vals.append(v)
        mults.append(m)
    return vals, mults


# ---------------------------------------------------------------------------
# The FREE parametrisation: centres and breakpoints.
#
# On the j-th linear piece Phi = c_j + x_j beta, so in theta = log beta
#     y = log Phi = log x_j + log(e^{sigma_j} + e^theta),  sigma_j = log(c_j/x_j)
#     y'(theta) = sigmoid(theta - sigma_j).
# Along the upper envelope c decreases and x increases, so sigma_1 > ... > sigma_k
# strictly, and the breakpoints theta_1 < ... < theta_{k-1} are UNCONSTRAINED:
# given sigma_j > sigma_{j+1} and any beta_j > 0,
#     x_{j+1}/x_j = (e^{sigma_j} + beta_j) / (e^{sigma_{j+1}} + beta_j)  >  1 .
# Hence, up to the projective scaling that the Hilbert metric quotients out,
#
#     C  =  { (sigma_1 > ... > sigma_k ;  theta_1 < ... < theta_{k-1}) }
#
# with no further constraints, and  y_a'(theta) = sigmoid(theta - S_a(theta))
# where S_a is the nonincreasing step function equal to sigma_j on the j-th
# piece.  In particular
#
#     sign (U_b - U_a)'(theta)  =  sign ( S_a(theta) - S_b(theta) ) ,
#
# which is the whole combinatorial content of the oscillation problem.
# ---------------------------------------------------------------------------

def from_centers(sigmas, thetas):
    """Build Phi from centres sigma_1 > ... > sigma_k and breakpoints."""
    sig = np.asarray(sigmas, float)
    th = np.asarray(thetas, float)
    k = len(sig)
    assert len(th) == k - 1
    x = np.empty(k)
    x[0] = 1.0
    for j in range(k - 1):
        b = math.exp(th[j])
        x[j + 1] = x[j] * (math.exp(sig[j]) + b) / (math.exp(sig[j + 1]) + b)
    c = x * np.exp(sig)
    return Trop(c, x)


def centers_of(t: Trop):
    """(sigma_j) and (theta_j) of an already-reduced Trop."""
    sig = np.log(t.c) - np.log(t.x)
    th = np.log(t.bps) if len(t.bps) else np.zeros(0)
    return sig, th
