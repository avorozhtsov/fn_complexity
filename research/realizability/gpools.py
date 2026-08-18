"""Pool construction and a fast A/d evaluator shared by the G scripts."""
from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def integer_pool(n, seed=11, rmax=7, vmax=40, rmin=2):
    """Random integer signatures; (1,) and all-ones are excluded by Sig."""
    rng = random.Random(seed)
    pool = set()
    while len(pool) < n:
        k = rng.randint(rmin, rmax)
        t = tuple(sorted((rng.randint(1, vmax) for _ in range(k)), reverse=True))
        if t[0] > 1:
            pool.add(t)
    return [C.Sig.of(t) for t in sorted(pool)]


class Batch:
    """A pool with a cached U table, giving vectorised A and d matrices."""

    def __init__(self, sigs, lo=-22.0, hi=22.0, step=0.005):
        self.sigs = list(sigs)
        self.g = np.arange(lo, hi + 0.5 * step, step)
        self.step = step
        self.tab = np.vstack([s.U(self.g) for s in self.sigs])
        self.eR = np.array([math.log(s.R) for s in self.sigs])
        self.eL = np.array([math.log(s.Lam) for s in self.sigs])

    def hi_lo(self):
        """(max, min) matrices of phi_ij = U_j - U_i, coarse (error ~step^2)."""
        n = len(self.sigs)
        HI = np.empty((n, n))
        LO = np.empty((n, n))
        for i in range(n):
            diff = self.tab - self.tab[i]
            HI[i] = diff.max(axis=1)
            LO[i] = diff.min(axis=1)
        e0 = self.eR[None, :] - self.eR[:, None]
        e1 = self.eL[None, :] - self.eL[:, None]
        HI = np.maximum(HI, np.maximum(e0, e1))
        LO = np.minimum(LO, np.minimum(e0, e1))
        return HI, LO

    def A_d(self):
        HI, LO = self.hi_lo()
        A = 0.5 * (HI + LO)
        d = HI - LO
        np.fill_diagonal(A, 0.0)
        np.fill_diagonal(d, 0.0)
        return A, d


def refine_A(sigs, idx=None):
    """Exact (certified + 40-digit) A matrix on a small family."""
    ss = sigs if idx is None else [sigs[i] for i in idx]
    n = len(ss)
    A = np.zeros((n, n))
    worst = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            _, _, A_mp, _, err = C.certified_A_d(ss[i], ss[j])
            A[i, j], A[j, i] = A_mp, -A_mp
            worst = max(worst, err)
    return A, worst


def score_seq(A):
    n = A.shape[0]
    return tuple(sorted(int(sum(1 for j in range(n) if A[i, j] > 0)) for i in range(n)))


def canon(A):
    """Canonical form of the tournament: lexicographically minimal adjacency
    bit-string over all vertex permutations.  Distinguishes iso classes."""
    import itertools
    n = A.shape[0]
    M = (A > 0)
    best = None
    for p in itertools.permutations(range(n)):
        bits = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    bits = (bits << 1) | int(M[p[i], p[j]])
        if best is None or bits < best:
            best = bits
    return best


def all_tournaments(n):
    """Canonical forms of all isomorphism classes of tournaments on n nodes."""
    import itertools
    pairs = list(itertools.combinations(range(n), 2))
    seen = {}
    for bits in range(1 << len(pairs)):
        M = np.zeros((n, n), dtype=bool)
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                M[i, j] = True
            else:
                M[j, i] = True
        c = canon(M.astype(float) - 0.5 * (~M))
        if c not in seen:
            seen[c] = M.copy()
    return seen
