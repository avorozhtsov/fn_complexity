"""Session brief L -- the quantum exchange rate: core numerics.

Two objects live here.

(1) THE SPECTRAL RATE.  For positive operators A, B on a finite-dimensional
    Hilbert space,  Z_A(beta) = Tr A^beta,  F_A = log Z_A, and

        C(A -> B) = inf_beta  F_A(beta) / F_B(beta).

    This is the verbatim extension of the classical exchange rate.  It is a
    function of spec(A) and spec(B) only -- see q_spectral.py.

(2) THE SANDWICHED RATE.  A *quantum signature* is a pair (A, S) of positive
    definite operators on the same space ("resource" and "background").  With
    the sandwiched Renyi quantity

        Qt_b(A||S) = Tr[ (S^{(1-b)/2b} A S^{(1-b)/2b})^b ]
        F_{(A,S)}(b) = log Qt_b(A||S)

    the exchange rate is  Ct((A,S) -> (B,T)) = inf_b F_{(A,S)}/F_{(B,T)}.

    S = I recovers (1) exactly.  [A,S] = 0 recovers the classical signature
    theory with real multiplicities.  [A,S] != 0 does not factor through
    spectra.

NUMERICS.  The literal expression S^{(1-b)/2b} A S^{(1-b)/2b} overflows in
double precision as b -> 0.  Use instead the similarity

    S^{-p} (S^p A S^p) S^p = A S^{2p},        p = (1-b)/(2b),

so the mu_k(b) are the generalised eigenvalues of the pencil (A, S^t) with

    t = (b-1)/b  in  [-1, 1)  for  b in [1/2, oo).

On b >= 1/2 -- exactly the range on which the sandwiched divergence is a
monotone -- the pencil is perfectly conditioned and plain LAPACK suffices.
Below b = 1/2 we fall back to mpmath at adaptive precision (QSigMP).

All matrices here are real symmetric; that is enough to break commutativity and
keeps every eigen-decomposition real.
"""
from __future__ import annotations

import math

import numpy as np

# ---------------------------------------------------------------- basic tools


def logsumexp(v):
    m = float(np.max(v))
    return m + math.log(float(np.sum(np.exp(v - m))))


def sym(M):
    return 0.5 * (M + M.T)


# ------------------------------------------------------- the sandwiched F

class QSig:
    """A quantum signature (A, S):  A, S real symmetric positive definite.

    F(b) = log Tr[(S^{(1-b)/2b} A S^{(1-b)/2b})^b], with the analytic values
    F(0) = log Tr S and F(b)/b -> Lam = log lambda_max(S^{-1/2} A S^{-1/2}).
    """

    __slots__ = ("A", "S", "_w", "_V", "_logs", "_K", "_Lam", "_F0", "_bsafe",
                 "_logc", "_logsc")

    def __init__(self, A, S):
        A = sym(np.asarray(A, dtype=float))
        S = sym(np.asarray(S, dtype=float))
        self.A, self.S = A, S
        w, V = np.linalg.eigh(S)
        if w.min() <= 0:
            raise ValueError("S must be positive definite")
        if np.linalg.eigvalsh(A).min() <= 0:
            raise ValueError("A must be positive definite")
        self._w, self._V, self._logs = w, V, np.log(w)
        Sm = (V * np.power(w, -0.5)) @ V.T
        self._K = sym(Sm @ A @ Sm)                       # S^{-1/2} A S^{-1/2}
        self._Lam = math.log(float(np.linalg.eigvalsh(self._K).max()))
        self._F0 = math.log(float(np.trace(S)))
        # S is used only through S^t; rescale S by its geometric mean so that
        # the exponents stay centred.  Qt_b(A||cS) = c^{1-b} Qt_b(A||S), so the
        # rescaling contributes the exact term (1-b) log c to F.
        self._logc = float(self._logs.mean())
        self._logsc = self._logs - self._logc
        # |t| * spread(log s) must stay well under log(1/eps) ~ 36 for LAPACK.
        spread = float(self._logs.max() - self._logs.min())
        self._bsafe = 1.0 / (1.0 + 18.0 / max(spread, 1e-12))

    # --- named quantities ------------------------------------------------
    @property
    def Lam(self):
        """lim F'(b) = log lambda_max(S^{-1/2} A S^{-1/2}).  ('log max fiber')"""
        return self._Lam

    @property
    def R(self):
        """F(0) = log Tr S.  ('log number of fibers')"""
        return self._F0

    @property
    def beta_safe(self):
        """Smallest b at which the double-precision pencil is trustworthy."""
        return self._bsafe

    def mu(self, beta):
        """Eigenvalues of the pencil (A, (S/c)^t), c the geometric mean of S.

        These are the mu_k of S^{(1-b)/2b} A S^{(1-b)/2b} rescaled by c^{b-1}.
        """
        t = (beta - 1.0) / beta
        B = (self._V * np.exp(t * self._logsc)) @ self._V.T
        L = np.linalg.cholesky(sym(B))
        Li = np.linalg.inv(L)
        return np.linalg.eigvalsh(sym(Li @ self.A @ Li.T))

    def F(self, beta):
        if beta <= 0.0:
            return self._F0
        if beta == math.inf:
            return math.inf
        m = self.mu(beta)
        return ((1.0 - beta) * self._logc
                + logsumexp(beta * np.log(np.maximum(m, 1e-300))))

    def U(self, s):
        """U(s) = log F(e^s)."""
        return math.log(self.F(math.exp(s)))

    # --- vectorised profiles (the search inner loop) ----------------------
    def mu_grid(self, betas):
        """(n, r) array of mu_k(beta) for a whole beta array, batched LAPACK."""
        betas = np.asarray(betas, float)
        t = (betas - 1.0) / betas
        wt = np.exp(t[:, None] * self._logsc[None, :])         # (n, r)
        B = np.einsum("ij,nj,kj->nik", self._V, wt, self._V)
        B = 0.5 * (B + np.transpose(B, (0, 2, 1)))
        L = np.linalg.cholesky(B)
        Li = np.linalg.inv(L)
        M = Li @ self.A @ np.transpose(Li, (0, 2, 1))
        return np.linalg.eigvalsh(0.5 * (M + np.transpose(M, (0, 2, 1))))

    def F_grid(self, betas):
        betas = np.asarray(betas, float)
        lm = betas[:, None] * np.log(np.maximum(self.mu_grid(betas), 1e-300))
        m = lm.max(axis=1)
        return ((1.0 - betas) * self._logc
                + m + np.log(np.exp(lm - m[:, None]).sum(axis=1)))

    def U_grid(self, sgrid):
        return np.log(self.F_grid(np.exp(np.asarray(sgrid, float))))

    # --- classical shadows ------------------------------------------------
    def pinched(self):
        """Decohere A in the eigenbasis of S: the canonical classical shadow."""
        V = self._V
        D = np.diag(np.diag(V.T @ self.A @ V))
        return QSig(V @ D @ V.T, self.S)

    def classical_atoms(self):
        """For a commuting pair: (multiplicities m_i, atoms exp(x_i))."""
        V = self._V
        d = np.diag(V.T @ self.A @ V)
        return self._w.copy(), d / self._w

    def coherence(self):
        """Relative Frobenius weight of the off-diagonal part of A in S's basis."""
        B = self._V.T @ self.A @ self._V
        off = B - np.diag(np.diag(B))
        return float(np.linalg.norm(off) / np.linalg.norm(B))

    def commutator_norm(self):
        C = self.A @ self.S - self.S @ self.A
        return float(np.linalg.norm(C) / (np.linalg.norm(self.A) * np.linalg.norm(self.S)))


def spectral_QSig(A):
    """The spectral object: (A, I).  F(b) = log Tr A^b."""
    return QSig(A, np.eye(np.asarray(A).shape[0]))


# ------------------------------------------------- oscillation / rate on a grid

def s_grid(smin=math.log(0.5), smax=math.log(2000.0), n=2001):
    """Default log-beta grid: the monotone range b in [1/2, 2000].

    The upper end is past the brief's 10^3 horizon; the b = oo endpoint is
    supplied analytically by Lam, so the grid never has to reach it.
    """
    return np.linspace(smin, smax, n)


def refine_extremum(fun, x0, x2, maximize, tol=1e-14):
    """Golden-section refinement on a bracketed extremum of a 1-D function."""
    g = (math.sqrt(5.0) - 1.0) / 2.0
    lo, hi = x0, x2
    a = hi - g * (hi - lo)
    b = lo + g * (hi - lo)
    sgn = -1.0 if maximize else 1.0
    fa, fb = sgn * fun(a), sgn * fun(b)
    for _ in range(200):
        if hi - lo <= tol * (1.0 + abs(a) + abs(b)):
            break
        if fa <= fb:
            hi, b, fb = b, a, fa
            a = hi - g * (hi - lo)
            fa = sgn * fun(a)
        else:
            lo, a, fa = a, b, fb
            b = lo + g * (hi - lo)
            fb = sgn * fun(b)
    return (a, sgn * fa) if fa <= fb else (b, sgn * fb)


def extrema(fun, grid, ends=()):
    """(max, min) of fun over grid u {ends}, interior extrema golden-refined."""
    vals = np.array([fun(s) for s in grid])
    hi = [float(vals.max())] + [e for e in ends if e is not None]
    lo = [float(vals.min())] + [e for e in ends if e is not None]
    # Strict on the left: on a numerically flat stretch a non-strict test fires
    # at every point and the golden refinement then dominates the run time.
    for i in range(1, len(grid) - 1):
        if vals[i] > vals[i - 1] and vals[i] >= vals[i + 1]:
            hi.append(refine_extremum(fun, grid[i - 1], grid[i + 1], True)[1])
        if vals[i] < vals[i - 1] and vals[i] <= vals[i + 1]:
            lo.append(refine_extremum(fun, grid[i - 1], grid[i + 1], False)[1])
    return max(hi), min(lo)


def endpoint_values(qa, qb, include_zero=True):
    """phi(-oo) = log R_b - log R_a,  phi(+oo) = log Lam_b - log Lam_a."""
    out = []
    if include_zero:
        out.append(math.log(qb.R) - math.log(qa.R))
    if qa.Lam > 0 and qb.Lam > 0:
        out.append(math.log(qb.Lam) - math.log(qa.Lam))
    return tuple(out)


def osc_mid(qa, qb, grid=None, include_zero=True):
    """(d, A) = (range, midrange) of phi = U_b - U_a."""
    if grid is None:
        grid = s_grid()
    fun = lambda s: qb.U(s) - qa.U(s)
    P, Q = extrema(fun, grid, endpoint_values(qa, qb, include_zero))
    return P - Q, 0.5 * (P + Q)


def osc_mid_fast(qa, qb, grid, Ua, Ub, include_zero=True, refine=True):
    """osc_mid from precomputed profiles Ua, Ub on `grid` (the search loop)."""
    phi = Ub - Ua
    ends = endpoint_values(qa, qb, include_zero)
    hi = [float(phi.max())] + list(ends)
    lo = [float(phi.min())] + list(ends)
    if refine:
        fun = lambda s: qb.U(s) - qa.U(s)
        for i in range(1, len(grid) - 1):
            if phi[i] > phi[i - 1] and phi[i] >= phi[i + 1]:
                hi.append(refine_extremum(fun, grid[i - 1], grid[i + 1], True)[1])
            if phi[i] < phi[i - 1] and phi[i] <= phi[i + 1]:
                lo.append(refine_extremum(fun, grid[i - 1], grid[i + 1], False)[1])
    P, Q = max(hi), min(lo)
    return P - Q, 0.5 * (P + Q)


def rate(qa, qb, grid=None, include_zero=True):
    """C(a -> b) = inf_b F_a / F_b  (a = implementer, b = implemented)."""
    if grid is None:
        grid = s_grid()
    fun = lambda s: qa.U(s) - qb.U(s)
    _, Q = extrema(fun, grid, endpoint_values(qb, qa, include_zero))
    return math.exp(Q)


def rand_spd(rng, r, lo=1.0, hi=6.0):
    """Random real symmetric positive definite matrix, eigenvalues in [lo,hi]."""
    X = rng.normal(size=(r, r))
    Qm, _ = np.linalg.qr(X)
    w = rng.uniform(lo, hi, size=r)
    return (Qm * w) @ Qm.T


def rand_admissible(rng, r, smax=6.0, kmax=30.0):
    """A random ADMISSIBLE quantum signature:  A >= S >= I.

    S = I + (psd),  A = S^{1/2} K S^{1/2} with K >= I; the two eigenbases are
    independent, so [A, S] != 0 generically.  This is the operator form of
    'multiplicities >= 1, fiber sizes >= 1', the standing hypothesis of
    FINDINGS Sec. 1.1.
    """
    S = rand_spd(rng, r, 1.0, smax)
    K = rand_spd(rng, r, 1.0, kmax)
    w, V = np.linalg.eigh(S)
    Sh = (V * np.sqrt(w)) @ V.T
    return QSig(sym(Sh @ K @ Sh), S)
