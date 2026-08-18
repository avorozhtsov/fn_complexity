"""Dependency-free global/local minimisers (scipy is unusable in this env).

``differential_evolution``  rand/1/bin DE with box bounds.
``pattern_search``          compass search with shrinking step; robust on the
                            non-smooth  min_ij  objectives used here.
"""
from __future__ import annotations

import numpy as np


def differential_evolution(fun, bounds, args=(), seed=0, maxiter=300,
                           popsize=20, F=(0.4, 1.0), CR=0.9, tol=1e-14,
                           init=None, callback=None):
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in bounds], float)
    hi = np.array([b[1] for b in bounds], float)
    dim = len(bounds)
    npop = max(8, min(popsize * dim, 120))
    pop = lo + rng.random((npop, dim)) * (hi - lo)
    if init is not None:
        pop[0] = np.clip(init, lo, hi)
    val = np.array([fun(p, *args) for p in pop])
    best = int(np.argmin(val))
    for _ in range(maxiter):
        f = rng.uniform(*F)
        idx = rng.integers(0, npop, size=(npop, 3))
        a, b, c = pop[idx[:, 0]], pop[idx[:, 1]], pop[idx[:, 2]]
        mutant = np.clip(a + f * (b - c), lo, hi)
        cross = rng.random((npop, dim)) < CR
        jrand = rng.integers(0, dim, size=npop)
        cross[np.arange(npop), jrand] = True
        trial = np.where(cross, mutant, pop)
        tval = np.array([fun(t, *args) for t in trial])
        take = tval < val
        pop[take] = trial[take]
        val[take] = tval[take]
        best = int(np.argmin(val))
        if callback is not None:
            callback(pop[best], val[best])
        if val.max() - val.min() < tol:
            break
    return pop[best], float(val[best])


def pattern_search(fun, x0, args=(), step=0.25, shrink=0.5, min_step=1e-11,
                   maxiter=200000, bounds=None):
    x = np.array(x0, float)
    if bounds is not None:
        lo = np.array([b[0] for b in bounds], float)
        hi = np.array([b[1] for b in bounds], float)
        x = np.clip(x, lo, hi)
    fx = fun(x, *args)
    dim = len(x)
    n = 0
    while step > min_step and n < maxiter:
        improved = False
        for i in range(dim):
            for sgn in (+1.0, -1.0):
                y = x.copy()
                y[i] += sgn * step
                if bounds is not None:
                    y = np.clip(y, lo, hi)
                fy = fun(y, *args)
                n += 1
                if fy < fx:
                    x, fx, improved = y, fy, True
                    break
        if not improved:
            step *= shrink
    return x, float(fx)
