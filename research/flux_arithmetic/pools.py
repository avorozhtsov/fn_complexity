#!/usr/bin/env python3
"""Signature pools for session brief E.

Two kinds of pool live here.

**Arithmetic pools.**  The complete enumeration of genus-two pencils
``y^2 = P(x) + c`` over ``F_q`` with ``P`` monic of degree 5 or 6 and
``P(0) = 0``.  ``P(0)`` may be normalised to zero because shifting it only
permutes the fibers.  For the affine plane curve the fiber count is

    N_c = #{(x,y) : y^2 = P(x) + c} = q + sum_x chi(P(x) + c)
        = q + sum_u n_u chi(u + c),      n_u = #{x : P(x) = u}

so the signature depends on ``P`` only through its value-multiplicity vector
``n``, and the whole enumeration is one histogram plus one ``q x q`` matrix
multiply per batch of polynomials.  This is a second, independent
implementation of the pool of ``research/curve_family_cycles/common.py``;
the two agree (296 signatures at ``q = 11``, 698 at ``q = 13``).

**Matched random controls.**  Signatures that are *not* curve signatures but
share the structural constraints of one: ``q`` positive integer entries summing
to ``q^2`` (equivalently ``q`` traces ``a_c = q - N_c`` summing to zero), with a
prescribed spread.  Three levels of matching are offered, see ``control_pool``.

Run directly to rebuild and cache the pools.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

# --------------------------------------------------------------- arithmetic


def quadratic_character(q: int) -> np.ndarray:
    """``chi`` on ``F_q`` as a length-``q`` array, ``chi(0) = 0``."""

    squares = {(x * x) % q for x in range(1, q)}
    return np.array([0] + [1 if x in squares else -1 for x in range(1, q)], dtype=np.int64)


def _chi_matrix(q: int) -> np.ndarray:
    """``K[c, u] = chi(u + c)``."""

    chi = quadratic_character(q)
    u = np.arange(q)
    return chi[(u[None, :] + u[:, None]) % q]


def _value_histograms(q: int, degree: int, batch: int = 4096) -> np.ndarray:
    """Value-multiplicity vectors ``n_u`` of every monic ``P``, ``P(0) = 0``.

    Returns an array of shape ``(q**(degree-1), q)``.  The row order is the
    lexicographic order of the coefficient tuple ``(a_1, ..., a_{degree-1})``.
    """

    x = np.arange(q, dtype=np.int64)
    powers = np.stack([(x**k) % q for k in range(1, degree + 1)])  # (degree, q)
    count = q ** (degree - 1)
    out = np.empty((count, q), dtype=np.int16)
    coeff_grid = np.array(list(itertools.product(range(q), repeat=degree - 1)), dtype=np.int64)
    for start in range(0, count, batch):
        block = coeff_grid[start : start + batch]  # (m, degree-1)
        m = block.shape[0]
        values = np.tile(powers[degree - 1], (m, 1))  # x**degree
        for j in range(degree - 1):
            values = values + block[:, j : j + 1] * powers[j][None, :]
        values %= q
        flat = (values + q * np.arange(m)[:, None]).ravel()
        out[start : start + m] = np.bincount(flat, minlength=q * m).reshape(m, q)
    return out


def arithmetic_pool(q: int, degrees=(5, 6)):
    """Complete pool of genus-two pencil signatures over ``F_q``.

    Returns ``(signatures, records)`` where ``signatures`` is an ``n x q``
    integer array of sorted (descending) signatures with every entry positive,
    and ``records`` is a list of dicts giving, per signature, the number of
    enumerated polynomials of each degree that realise it and one witness
    coefficient tuple.
    """

    K = _chi_matrix(q)
    index: dict[tuple[int, ...], int] = {}
    counts: list[dict[int, int]] = []
    witness: list[tuple[int, tuple[int, ...]]] = []
    for degree in degrees:
        hist = _value_histograms(q, degree)
        N = q + hist.astype(np.int64) @ K.T  # (count, q): N[p, c]
        N.sort(axis=1)
        N = N[:, ::-1]
        coeff_grid = np.array(
            list(itertools.product(range(q), repeat=degree - 1)), dtype=np.int64
        )
        for row, coeffs in zip(N, coeff_grid):
            key = tuple(int(v) for v in row)
            if key not in index:
                index[key] = len(index)
                counts.append({})
                witness.append((degree, tuple(int(c) for c in coeffs)))
            slot = counts[index[key]]
            slot[degree] = slot.get(degree, 0) + 1
    keys = sorted(index, key=lambda k: index[k])
    signatures = np.array(keys, dtype=np.int64)
    records = [
        {"counts": counts[i], "witness_degree": witness[i][0], "witness_coeffs": witness[i][1]}
        for i in range(len(keys))
    ]
    keep = (signatures > 0).all(axis=1)
    signatures = signatures[keep]
    records = [r for r, k in zip(records, keep) if k]
    return signatures, records


def sampled_pool(q: int, degrees=(5, 6), draws: int = 20000, seed: int = 20260818):
    """Random sample of genus-two pencil signatures over ``F_q``.

    For ``q`` where the complete enumeration is out of reach.  Same
    normalisation as :func:`arithmetic_pool`.
    """

    rng = np.random.default_rng(seed)
    K = _chi_matrix(q)
    x = np.arange(q, dtype=np.int64)
    seen: dict[tuple[int, ...], None] = {}
    for degree in degrees:
        powers = np.stack([(x**k) % q for k in range(1, degree + 1)])
        for start in range(0, draws, 2048):
            m = min(2048, draws - start)
            block = rng.integers(0, q, size=(m, degree - 1))
            values = np.tile(powers[degree - 1], (m, 1))
            for j in range(degree - 1):
                values = values + block[:, j : j + 1] * powers[j][None, :]
            values %= q
            flat = (values + q * np.arange(m)[:, None]).ravel()
            hist = np.bincount(flat, minlength=q * m).reshape(m, q)
            N = q + hist @ K.T
            N.sort(axis=1)
            N = N[:, ::-1]
            for row in N:
                if row.min() > 0:
                    seen.setdefault(tuple(int(v) for v in row), None)
    return np.array(sorted(seen), dtype=np.int64)


# ------------------------------------------------------------- the controls


def _repair(a: np.ndarray, q: int, rng) -> np.ndarray | None:
    """Centre ``a`` on the integer lattice so that ``sum a = 0`` exactly.

    The correction is spread evenly -- subtract the floor of the mean from every
    coordinate, then take one further unit off ``sum a mod q`` coordinates
    chosen uniformly.  A repair that instead hammers the extreme coordinate
    (the obvious implementation) destroys the second moment: the typical excess
    is ``sqrt(q * var)``, so it would shave that many units off the largest
    entries and shrink ``m2`` by tens of percent.
    """

    a = np.asarray(a, dtype=np.int64).copy()
    a -= int(np.floor(a.sum() / q))
    excess = int(a.sum())
    if excess:
        a[rng.choice(q, excess, replace=False)] -= 1
    assert a.sum() == 0
    if (q - a).min() < 1:
        return None
    return a


def control_pool(q: int, n: int, kind: str, reference: np.ndarray | None = None,
                 seed: int = 7):
    """A random control pool of ``n`` signatures with ``q`` positive entries
    summing to ``q**2``.

    Every control satisfies the two structural identities a curve signature
    satisfies -- ``q`` fibers and ``sum_c N_c = q^2``, equivalently ``q`` traces
    summing to zero -- so the ``beta = 0`` endpoint is a tie across the pool
    exactly as it is for the arithmetic pools.  They differ in how much of the
    trace *distribution* they are told to match.

    ``"loose"``
        traces iid uniform on the genus-two Weil box ``[-4 sqrt q, 4 sqrt q]``.
        Deliberately over-spread: ``m2`` about 5 against the pool's 1.33.
    ``"m2matched"``
        traces iid Gaussian, rescaled so that each control signature's ``m2``
        equals an ``m2`` drawn from the arithmetic pool's empirical ``m2``
        distribution.  Matches the spread, not the shape.
    ``"marginal"``
        traces drawn iid from the pooled empirical distribution of all traces
        of ``reference``.  Matches the one-point trace law -- the Sato-Tate
        distribution of the pool -- and destroys only *which* traces co-occur
        inside one pencil.
    ``"sigshuffle"``
        each order statistic drawn independently from that order statistic's
        empirical distribution over ``reference``.
    ``"maxmatched"``
        the tightest control.  The multiset of largest fibers ``max_c N_c`` is
        *identical* to the arithmetic pool's, so the endpoint potential
        ``psi = (1/2) log phi`` -- which E2 shows is the whole gradient part --
        agrees signature by signature; the remaining traces come from the
        pooled marginal law truncated to be compatible.  Any difference in curl
        is then interior structure and nothing else.
    """

    rng = np.random.default_rng(seed)
    if kind == "maxmatched":
        return _max_matched(q, n, reference, rng)
    out: dict[tuple[int, ...], None] = {}
    traces = m2_pool = cols = None
    if kind in ("marginal",):
        assert reference is not None
        traces = (q - reference).ravel()
    if kind in ("m2matched",):
        assert reference is not None
        a = q - reference.astype(float)
        m2_pool = (a**2).sum(axis=1) / q**2
    if kind == "sigshuffle":
        assert reference is not None
        cols = [reference[:, j].copy() for j in range(q)]

    guard = 0
    while len(out) < n and guard < 2000 * n:
        guard += 1
        if kind == "loose":
            hi = 2 * 2 * math.sqrt(q)
            a = np.rint(rng.uniform(-hi, hi, size=q))
        elif kind == "m2matched":
            target = float(rng.choice(m2_pool))
            z = rng.normal(0.0, 1.0, size=q)
            z -= z.mean()
            scale = math.sqrt(q**2 * target / max((z**2).sum(), 1e-12))
            a = np.rint(scale * z)
        elif kind == "marginal":
            a = rng.choice(traces, size=q, replace=True)
        elif kind == "sigshuffle":
            a = q - np.array([rng.choice(c) for c in cols])
        else:
            raise ValueError(kind)
        fixed = _repair(a, q, rng)
        if fixed is None:
            continue
        sig = np.sort(q - fixed)[::-1]
        out.setdefault(tuple(int(v) for v in sig), None)
    return np.array(sorted(out), dtype=np.int64)[:n]


def _max_matched(q: int, n: int, reference: np.ndarray, rng) -> np.ndarray:
    """Controls whose multiset of largest fibers equals the reference's."""

    assert reference is not None
    traces = (q - reference).ravel()
    targets = list(reference.max(axis=1))
    rng.shuffle(targets)
    out: dict[tuple[int, ...], None] = {}
    for target in (targets * (1 + n // max(len(targets), 1) + 1))[: 60 * n]:
        if len(out) >= n:
            break
        a_min = q - int(target)  # the trace of the largest fiber
        for _ in range(400):
            rest = rng.choice(traces[traces >= a_min], size=q - 1, replace=True)
            a = np.concatenate([[a_min], rest])
            a -= int(np.floor(a.sum() / q))
            excess = int(a.sum())
            if excess:
                movable = np.flatnonzero(a[1:] > a_min) + 1
                if len(movable) < excess:
                    continue
                a[rng.choice(movable, excess, replace=False)] -= 1
            if a.sum() != 0 or a.min() != a_min or (q - a).min() < 1:
                continue
            out.setdefault(tuple(int(v) for v in np.sort(q - a)[::-1]), None)
            break
    return np.array(sorted(out), dtype=np.int64)[:n]


if __name__ == "__main__":
    for q in (11, 13):
        sig, rec = arithmetic_pool(q)
        total = sum(sum(r["counts"].values()) for r in rec)
        print(f"q={q}: {len(sig)} positive signatures from {total} enumerated pencils")
        np.save(f"pool_q{q}.npy", sig)
