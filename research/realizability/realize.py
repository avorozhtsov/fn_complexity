"""Targeted realisation of a prescribed sign pattern for A.

Given a target tournament ``T`` on n vertices (T[i,j] = True meaning we want
A(i,j) > 0), search for n signatures with r atoms each attaining it with the
largest possible margin  min_{i<j} sgn_ij * A(i,j).

Atoms are searched over the reals (legitimate for the structural question,
brief G "Traps").  A witness is converted to integers exactly: raising every
atom of every signature to a common power p is the reparametrisation
beta -> p*beta, under which d and A are *exactly* invariant, so
a_i -> round(a_i^p) with a_i^p >= 10^9 changes log a_i by < 5e-10.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

GRID = np.arange(-18.0, 18.0 + 1e-9, 0.01)   # search grid; witnesses are
# re-certified against mpmath afterwards, so this only has to rank candidates


def _A_matrix(sigs, grid=GRID):
    tab = np.vstack([s.U(grid) for s in sigs])
    eR = np.array([math.log(s.R) for s in sigs])
    eL = np.array([math.log(s.Lam) for s in sigs])
    n = len(sigs)
    A = np.zeros((n, n))
    for i in range(n):
        diff = tab - tab[i]
        hi = np.maximum(diff.max(axis=1), np.maximum(eR - eR[i], eL - eL[i]))
        lo = np.minimum(diff.min(axis=1), np.minimum(eR - eR[i], eL - eL[i]))
        A[i] = 0.5 * (hi + lo)
    np.fill_diagonal(A, 0.0)
    return A


XSNAP = 0.30      # log-atoms below this snap to 0, i.e. to a fiber of size 1


def _sigs_from(x, n, r):
    x = np.asarray(x, float).reshape(n, r)
    out = []
    for row in x:
        xs = np.where(row < XSNAP, 0.0, np.clip(row, 0.0, 12.0))
        if xs.max() <= 0.0:
            xs[int(np.argmax(row))] = XSNAP
        out.append(C.Sig.from_logs(xs))
    return out


def margin(x, n, r, S):
    """-(min_{i<j} S_ij A_ij).  Minimising it maximises the realisation margin.

    NOTE: the all-equal family gives A == 0 and objective 0, which beats every
    partly-wrong configuration; used alone this objective collapses onto that
    degenerate manifold.  It is therefore only used as a *polisher*, started
    from a point that already realises the tournament.
    """
    try:
        A = _A_matrix(_sigs_from(x, n, r))
    except ValueError:
        return 1e3
    iu = np.triu_indices(n, 1)
    return -float(np.min(S[iu] * A[iu]))


def soft(x, n, r, S, kappa=2e-3):
    """-sum tanh(S_ij A_ij / kappa): rewards correct signs, and the degenerate
    A == 0 scores 0, strictly worse than any correct configuration."""
    try:
        A = _A_matrix(_sigs_from(x, n, r))
    except ValueError:
        return 1e3
    iu = np.triu_indices(n, 1)
    return -float(np.sum(np.tanh(S[iu] * A[iu] / kappa)))


def realise(T, r=4, seed=0, maxiter=300, popsize=20, restarts=2, xmax=5.0,
            kappa=2e-3):
    """Search for n signatures with r atoms each realising the tournament T."""
    n = T.shape[0]
    S = np.where(T, 1.0, -1.0)
    S = np.triu(S, 1) - np.triu(S, 1).T
    bounds = [(0.0, xmax)] * (n * r)
    best_x, best_m = None, -np.inf
    for k in range(restarts):
        x, _ = differential_evolution(soft, bounds, args=(n, r, S, kappa),
                                      seed=seed + 977 * k, maxiter=maxiter,
                                      popsize=popsize, F=(0.3, 1.2), CR=0.9)
        x, _ = pattern_search(soft, x, args=(n, r, S, kappa), step=0.2,
                              min_step=1e-6, maxiter=6000, bounds=bounds)
        m = -margin(x, n, r, S)
        if m > 0:                                   # polish the true margin
            for step in (0.1, 0.01, 1e-3):
                x, f = pattern_search(margin, x, args=(n, r, S), step=step,
                                      min_step=1e-8, maxiter=6000, bounds=bounds)
            m = -margin(x, n, r, S)
        if m > best_m:
            best_x, best_m = x, m
    return best_x, best_m


def to_integers(sigs, target=10 ** 10):
    """Round a real witness to an integer one via the exact power symmetry.

    Raising every atom of every signature to a common power p is the
    reparametrisation beta -> p*beta, under which A and d are exactly
    invariant.  p is chosen so the smallest atom above 1 becomes >= ``target``;
    rounding then perturbs each log-atom by at most 1/(2*target).
    """
    xmin = min(min(x for x in s.xs if x > 0) for s in sigs)
    p = math.log(target) / xmin
    out = []
    for s in sigs:
        atoms = []
        for x, m in zip(s.xs, s.mults):
            v = 1 if x <= 0 else int(round(math.exp(min(p * x, 2000.0))))
            atoms.extend([max(v, 1)] * int(round(m)))
        out.append(tuple(sorted(atoms, key=lambda t: math.log(t), reverse=True)))
    return out, p


def certified_matrix(atom_tuples, dps=40):
    """A matrix of an integer family, each entry re-verified at dps digits."""
    ss = [C.Sig.of(t) for t in atom_tuples]
    n = len(ss)
    A = np.zeros((n, n))
    worst = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            _, _, A_mp, _, err = C.certified_A_d(ss[i], ss[j], dps)
            A[i, j], A[j, i] = A_mp, -A_mp
            worst = max(worst, err)
    return A, worst
