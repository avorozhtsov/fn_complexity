"""Core utilities for session brief G (realizability and rigidity).

Conventions (PLAN.md, brief D Part 0):

    Z_a(beta) = sum_i a_i^beta          F_a = log Z_a          u_a = log F_a
    L(a,b)    = -log C(a->b) = sup_beta (u_b - u_a)
    d(a,b)    = L(a,b) + L(b,a) = osc_beta (u_b - u_a)
    A(a,b)    = (L(a,b) - L(b,a))/2 = mid_beta (u_b - u_a)     (a < b iff A>0)

Everything is done in the variable s = log beta, on which the structure of
u is transparent:

    U_a(s) := u_a(e^s) = log(Lam_a) + max(sigma_a, s) + w_a(s)

with  R_a = log r_a,  Lam_a = log max_i a_i,  sigma_a = log(R_a/Lam_a),
and w_a the "bump" of the structure theorem in FINDINGS.md.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np

LOG2 = math.log(2.0)


@dataclass(frozen=True)
class Sig:
    """A signature, stored compressed as (distinct value, multiplicity) pairs.

    ``atoms`` is the expanded multiset when it is small enough to be useful
    (it is only used for printing and for the package cross-check); ``vals``
    and ``mults`` carry the arithmetic, so r may be astronomically large.
    """

    xs: tuple            # log of the distinct atom values, strictly decreasing
    mults: tuple         # positive multiplicities, same length
    exact: tuple = ()    # the exact integer values, when the signature is one

    @staticmethod
    def of(atoms) -> "Sig":
        from collections import Counter
        c = Counter(atoms)
        vals = sorted(c, key=lambda v: math.log(v), reverse=True)
        exact = tuple(vals) if all(isinstance(v, int) for v in vals) else ()
        return Sig(tuple(math.log(v) for v in vals),
                   tuple(float(c[v]) for v in vals), exact)._checked()

    @staticmethod
    def compressed(vals, mults) -> "Sig":
        pairs = sorted(((math.log(v), float(m), v) for v, m in zip(vals, mults)
                        if m > 0), reverse=True)
        exact = tuple(v for _, _, v in pairs)
        if not all(isinstance(v, int) for v in exact):
            exact = ()
        return Sig(tuple(x for x, _, _ in pairs),
                   tuple(m for _, m, _ in pairs), exact)._checked()

    @staticmethod
    def from_logs(xs, mults=None) -> "Sig":
        xs = list(xs)
        mults = [1.0] * len(xs) if mults is None else list(mults)
        agg = {}
        for x, m in zip(xs, mults):
            agg[round(float(x), 14)] = agg.get(round(float(x), 14), 0.0) + float(m)
        ks = sorted(agg, reverse=True)
        return Sig(tuple(ks), tuple(agg[k] for k in ks), ())._checked()

    def _checked(self):
        if not self.xs:
            raise ValueError("empty signature")
        if sum(self.mults) < 2:
            raise ValueError("need at least two atoms (r >= 2)")
        if self.xs[0] <= 0.0:
            raise ValueError("need max atom > 1 (excludes (1,) and all-ones)")
        if self.xs[-1] < 0.0:
            raise ValueError("atoms must be >= 1")
        return self

    @property
    def vals(self):
        return self.exact if self.exact else tuple(math.exp(x) for x in self.xs)

    @property
    def atoms(self):
        out = []
        for v, m in zip(self.vals, self.mults):
            k = int(round(m))
            if k > 4096 or len(out) > 4096:
                return tuple(self.vals)  # too big to expand; values only
            out.extend([v] * k)
        return tuple(out)

    # --- basic invariants -------------------------------------------------
    @property
    def x(self):
        return np.asarray(self.xs, dtype=float)

    @property
    def m(self):
        return np.asarray(self.mults, dtype=float)

    @property
    def r(self):
        return sum(self.mults)

    @property
    def R(self):
        return math.log(sum(self.mults))          # F(0)

    @property
    def Lam(self):
        return self.xs[0]                         # F'(infinity)

    @property
    def tau(self):
        return self.R / self.Lam

    @property
    def sigma(self):
        return math.log(self.R) - math.log(self.Lam)

    @property
    def psi(self):
        """psi = 1/2 log(log r * log max) : the exact-part potential."""
        return 0.5 * (math.log(self.R) + math.log(self.Lam))

    # --- the analytic content --------------------------------------------
    def F(self, beta):
        """log Z(beta) = log sum_i m_i v_i^beta, stable."""
        beta = np.asarray(beta, dtype=float)
        x, m = self.x, self.m
        top = beta * x[0]
        rest = (np.exp(np.multiply.outer(beta, x - x[0])) * m).sum(axis=-1)
        return top + np.log(rest)

    def U(self, s):
        """u(e^s) = log log Z(e^s)."""
        return np.log(self.F(np.exp(np.asarray(s, dtype=float))))

    def w(self, s):
        """the bump: U(s) - log(Lam) - max(sigma, s)."""
        s = np.asarray(s, dtype=float)
        return self.U(s) - math.log(self.Lam) - np.maximum(self.sigma, s)


def envelope(s, sigma):
    """log(1 + e^{-|s-sigma|}) : the proved upper envelope for w."""
    return np.log1p(np.exp(-np.abs(np.asarray(s, float) - sigma)))


# ---------------------------------------------------------------------------
# max / min of phi = U_b - U_a over beta in [0, infinity]
# ---------------------------------------------------------------------------

def endpoints(a: Sig, b: Sig):
    """(phi(-inf), phi(+inf)) = (log(R_b/R_a), log(Lam_b/Lam_a))."""
    return (math.log(b.R) - math.log(a.R), math.log(b.Lam) - math.log(a.Lam))


def make_grid(lo=-30.0, hi=30.0, step=0.002):
    return np.arange(lo, hi + 0.5 * step, step)


_GRID = make_grid()


def _refine(f, i, g, phi, want_max):
    y0, y1, y2 = phi[i - 1], phi[i], phi[i + 1]
    den = y0 - 2.0 * y1 + y2
    if (want_max and den >= 0.0) or ((not want_max) and den <= 0.0):
        return float(g[i]), float(y1)
    delta = 0.5 * (y0 - y2) / den
    if abs(delta) > 1.0:
        return float(g[i]), float(y1)
    h = g[1] - g[0]
    sx = float(g[i] + delta * h)
    return sx, float(f(sx))


def extrema(a: Sig, b: Sig, grid=None):
    """Return (max phi, min phi, argmax s, argmin s) over beta in [0, inf].

    The endpoints beta = 0 and beta = infinity are included exactly; the
    interior is screened on a grid in s = log beta and refined parabolically
    (phi is real-analytic in s, so this is valid).
    """
    g = _GRID if grid is None else grid
    phi = b.U(g) - a.U(g)
    e0, e1 = endpoints(a, b)

    def f(s):
        return float(b.U(np.array([s]))[0] - a.U(np.array([s]))[0])

    hi, lo = max(e0, e1), min(e0, e1)
    s_hi = -math.inf if e0 >= e1 else math.inf
    s_lo = -math.inf if e0 <= e1 else math.inf
    i = int(np.argmax(phi))
    if 0 < i < len(g) - 1:
        sx, yx = _refine(f, i, g, phi, True)
    else:
        sx, yx = float(g[i]), float(phi[i])
    if yx > hi:
        hi, s_hi = yx, sx
    j = int(np.argmin(phi))
    if 0 < j < len(g) - 1:
        sx, yx = _refine(f, j, g, phi, False)
    else:
        sx, yx = float(g[j]), float(phi[j])
    if yx < lo:
        lo, s_lo = yx, sx
    return hi, lo, s_hi, s_lo


def d_and_A(a: Sig, b: Sig, grid=None):
    hi, lo, _, _ = extrema(a, b, grid)
    return hi - lo, 0.5 * (hi + lo)


def parts(a: Sig, b: Sig, grid=None):
    """d, A, the exact part dpsi, the defect D, the overshoots P and Q."""
    hi, lo, s_hi, s_lo = extrema(a, b, grid)
    e0, e1 = endpoints(a, b)
    P = hi - max(e0, e1)
    Q = min(e0, e1) - lo
    return {
        "d": hi - lo,
        "A": 0.5 * (hi + lo),
        "dpsi": b.psi - a.psi,
        "D": 0.5 * (hi + lo) - (b.psi - a.psi),
        "P": P,
        "Q": Q,
        "eps": P + Q,
        "dsigma": b.sigma - a.sigma,
        "s_hi": s_hi,
        "s_lo": s_lo,
    }


# ---------------------------------------------------------------------------
# matrices over a family, from one shared U table
# ---------------------------------------------------------------------------

def matrices(sigs, grid=None):
    g = _GRID if grid is None else grid
    n = len(sigs)
    tab = np.vstack([s.U(g) for s in sigs])
    e_R = np.array([math.log(s.R) for s in sigs])
    e_L = np.array([math.log(s.Lam) for s in sigs])
    D = np.zeros((n, n))
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            phi = tab[j] - tab[i]
            e0, e1 = e_R[j] - e_R[i], e_L[j] - e_L[i]
            hi, lo = max(e0, e1), min(e0, e1)
            k = int(np.argmax(phi))
            if 0 < k < len(g) - 1:
                y0, y1, y2 = phi[k - 1], phi[k], phi[k + 1]
                den = y0 - 2 * y1 + y2
                yx = y1 - 0.125 * (y0 - y2) ** 2 / den if den < 0 else y1
            else:
                yx = phi[k]
            hi = max(hi, yx)
            k = int(np.argmin(phi))
            if 0 < k < len(g) - 1:
                y0, y1, y2 = phi[k - 1], phi[k], phi[k + 1]
                den = y0 - 2 * y1 + y2
                yx = y1 - 0.125 * (y0 - y2) ** 2 / den if den > 0 else y1
            else:
                yx = phi[k]
            lo = min(lo, yx)
            D[i, j] = D[j, i] = hi - lo
            A[i, j] = 0.5 * (hi + lo)
            A[j, i] = -A[i, j]
    return D, A


# ---------------------------------------------------------------------------
# certified extrema.  |phi'(s)| <= 1 because U'(s) = beta F'/F in [0,1] for
# every signature (structure theorem S3), so a grid of step h in s brackets
# the global max within h/2.  That makes a rigorous global search possible.
# ---------------------------------------------------------------------------

def certified_extrema(a: Sig, b: Sig, tol=1e-12, lo=-40.0, hi=40.0, step=0.001):
    """(max, min) of phi = U_b - U_a over beta in [0, inf], with a rigorous
    Lipschitz bracket: every returned value is within ``tol`` of the truth."""
    g = np.arange(lo, hi + 0.5 * step, step)
    phi = b.U(g) - a.U(g)
    e0, e1 = endpoints(a, b)

    def f(s):
        return float(b.U(np.array([s]))[0] - a.U(np.array([s]))[0])

    def golden(l, r, sign):
        gr = (math.sqrt(5.0) - 1.0) / 2.0
        x1, x2 = r - gr * (r - l), l + gr * (r - l)
        y1, y2 = sign * f(x1), sign * f(x2)
        for _ in range(200):
            if r - l < 1e-13:
                break
            if y1 <= y2:
                r, x2, y2 = x2, x1, y1
                x1 = r - gr * (r - l)
                y1 = sign * f(x1)
            else:
                l, x1, y1 = x1, x2, y2
                x2 = l + gr * (r - l)
                y2 = sign * f(x2)
        return (x1, sign * y1) if y1 <= y2 else (x2, sign * y2)

    out = []
    for sign in (-1.0, +1.0):          # sign=-1 -> maximise, +1 -> minimise
        best = min(sign * e0, sign * e1)
        best_s = -math.inf if sign * e0 <= sign * e1 else math.inf
        # every grid cell whose value is within `step` of the best can host a
        # better point (Lipschitz constant 1); refine one point per contiguous
        # run of such cells.
        vals = sign * phi
        cut = min(vals.min(), best) + step
        mask = vals <= cut
        edges = np.flatnonzero(np.diff(mask.astype(np.int8)))
        starts = list(edges[::2] + 1) if not mask[0] else [0] + list(edges[1::2] + 1)
        ends = list(edges[1::2] + 1) if not mask[0] else list(edges[::2] + 1) + [len(mask)]
        for s0, s1 in zip(starts, ends):
            k = s0 + int(np.argmin(vals[s0:s1]))
            l = g[max(k - 1, 0)]
            r = g[min(k + 1, len(g) - 1)]
            s_, y_ = golden(l, r, sign)
            if sign * y_ < best:
                best, best_s = sign * y_, s_
        out.append((sign * best, best_s))
    (mx, s_mx), (mn, s_mn) = out
    return mx, mn, s_mx, s_mn


def mp_extrema(a: Sig, b: Sig, s_hi, s_lo, dps=40):
    """Re-evaluate the two located extrema at ``dps`` digits and return
    (max, min) as mpmath floats.  Endpoints are exact closed forms."""
    from mpmath import mp, mpf, log as mplog, exp as mpexp

    mp.dps = dps

    def mpx(sig):
        """Exact log-atoms at working precision (exact integers when known)."""
        if sig.exact:
            return [mplog(mpf(v)) if v < 2 ** 53 else mplog(mpf(str(v)))
                    for v in sig.exact]
        return [mpf(str(x)) for x in sig.xs]

    xa, xb = mpx(a), mpx(b)

    def mpU(sig, xs, s):
        beta = mpexp(mpf(s))
        z = sum(mpf(m) * mpexp(beta * x) for x, m in zip(xs, sig.mults))
        return mplog(mplog(z))

    def phi(s):
        if s == -math.inf:
            return mplog(mplog(mpf(b.r))) - mplog(mplog(mpf(a.r)))
        if s == math.inf:
            return mplog(xb[0]) - mplog(xa[0])
        return mpU(b, xb, s) - mpU(a, xa, s)

    def refine(s0, want_max):
        if math.isinf(s0):
            return phi(s0)
        h = mpf(10) ** (-8)
        s = mpf(s0)
        for _ in range(60):
            f1 = (phi(s + h) - phi(s - h)) / (2 * h)
            f2 = (phi(s + h) - 2 * phi(s) + phi(s - h)) / (h * h)
            if f2 == 0:
                break
            stepn = f1 / f2
            if abs(stepn) > 1:
                stepn = mpf(1) if stepn > 0 else mpf(-1)
            s = s - stepn
            if abs(stepn) < mpf(10) ** (-(dps // 2)):
                break
        cand = [phi(s), phi(mpf(s0))]
        return max(cand) if want_max else min(cand)

    return refine(s_hi, True), refine(s_lo, False)


def certified_A_d(a: Sig, b: Sig, dps=40):
    """Return (A, d) computed in double and re-verified at ``dps`` digits."""
    mx, mn, s_hi, s_lo = certified_extrema(a, b)
    mmx, mmn = mp_extrema(a, b, s_hi, s_lo, dps)
    A_mp = (mmx + mmn) / 2
    d_mp = mmx - mmn
    return (0.5 * (mx + mn), mx - mn, float(A_mp), float(d_mp),
            abs(float(A_mp) - 0.5 * (mx + mn)))


def hodge(A):
    """HodgeRank split: gradient fraction and residual fraction of ||A||."""
    psi = -A.mean(axis=1)
    G = psi[None, :] - psi[:, None]
    nk = np.linalg.norm(A)
    return np.linalg.norm(G) / nk, np.linalg.norm(A - G) / nk


def three_cycles(A, tol=0.0):
    n = A.shape[0]
    out = []
    for i, j, k in itertools.combinations(range(n), 3):
        if (A[i, j] > tol) == (A[j, k] > tol) == (A[k, i] > tol):
            out.append((i, j, k))
    return out


def tournament(A, tol=1e-10):
    """Score sequence (out-degrees for the relation a < b iff A(a,b)>0)."""
    n = A.shape[0]
    if np.min(np.abs(A[np.triu_indices(n, 1)])) < tol:
        return None
    return tuple(sorted(int(sum(1 for j in range(n) if A[i, j] > 0)) for i in range(n)))
