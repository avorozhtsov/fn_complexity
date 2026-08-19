#!/usr/bin/env python3
"""Superadditivity of the exchange rate under Cartesian products.

    C(a (x) b -> c)  >=  C(a->c) + C(b->c)                       always,
    with equality iff log Z_a/log Z_c and log Z_b/log Z_c share a minimiser.

Reproduces the tables of research/synergy/FINDINGS.md.

    python research/synergy/synergy.py
"""
from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fn_complexity import exchange_rate, exchange_rate_result  # noqa: E402

TOL = 1e-9  # the project's tie threshold is 1e-10; stay an order above it


def tensor(a, b):
    """Signature of the Cartesian product: fiber sizes multiply."""
    return tuple(sorted((x * y for x in a for y in b), reverse=True))


def gap(a, b, c):
    return exchange_rate(tensor(a, b), c) - exchange_rate(a, c) - exchange_rate(b, c)


def classes(q, degree):
    """Signatures of all maps F_q^2 -> F_q of total degree <= `degree`."""
    mons = [(i, j) for i in range(degree + 1) for j in range(degree + 1) if i + j <= degree]
    seen = set()
    for co in itertools.product(range(q), repeat=len(mons)):
        counts = {}
        for x in range(q):
            for y in range(q):
                v = sum(c * pow(x, i, q) * pow(y, j, q) for c, (i, j) in zip(co, mons)) % q
                counts[v] = counts.get(v, 0) + 1
        sig = tuple(sorted(counts.values(), reverse=True))
        if len(sig) > 1:  # drop the constant maps
            seen.add(sig)
    return sorted(seen, reverse=True)


def report(pool, label):
    print(f"\n=== {label} ===")
    print(f"signatures: {pool}")
    strict = []
    worst = float("inf")
    for a, b, c in itertools.product(pool, repeat=3):
        g = gap(a, b, c)
        worst = min(worst, g)
        if g > TOL:
            ra, rb = exchange_rate_result(a, c), exchange_rate_result(b, c)
            strict.append((g, a, b, c, ra, rb))
    strict.sort(reverse=True)
    print(f"triples: {len(pool)**3}   strictly superadditive: {len(strict)}")
    print(f"minimum gap (superadditivity, must be >= 0): {worst:.3e}")
    for g, a, b, c, ra, rb in strict[:2]:
        print(f"  gap={g:.6f}  {a} (x) {b} -> {c}")
        print(f"      C(a->c)={ra.rate:.6f} @ beta={ra.beta}"
              f"   C(b->c)={rb.rate:.6f} @ beta={rb.beta}")


def interior_search(trials=4000, seed=9):
    """Strict cases with BOTH contacts interior -- not an endpoint artefact."""
    rng = random.Random(seed)

    def rs():
        n = rng.randint(3, 7)
        return tuple(sorted((rng.randint(2, 30) for _ in range(n)), reverse=True))

    pool = list({rs() for _ in range(80)})
    finite = lambda t: t not in (0.0, float("inf"))
    both, flat = [], 0
    for _ in range(trials):
        a, b, c = rng.sample(pool, 3)
        ra, rb = exchange_rate_result(a, c), exchange_rate_result(b, c)
        g = exchange_rate(tensor(a, b), c) - ra.rate - rb.rate
        if finite(ra.beta) and finite(rb.beta):
            if g > TOL:
                both.append((g, a, b, c, ra.beta, rb.beta))
            elif abs(ra.beta - rb.beta) > 1.0:
                flat += 1
    both.sort(reverse=True)
    print(f"\n=== random integer signatures, {trials} triples ===")
    print(f"strict with BOTH contacts interior: {len(both)}")
    for g, a, b, c, ba, bb in both[:3]:
        print(f"  gap={g:.6f}  {a} (x) {b} -> {c}   beta_a={ba:.4f}  beta_b={bb:.4f}")
    print(f"contacts far apart but gap below {TOL:g}: {flat}"
          "  <- flat minima and spurious argmins; test the gap, not the argmin")


if __name__ == "__main__":
    report(classes(3, 2), "F_3 quadratic classes")
    report(classes(5, 2), "F_5 quadratic classes")
    report(classes(3, 3), "F_3 cubic classes")
    interior_search()
