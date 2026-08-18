#!/usr/bin/env python3
"""Independent checks, at the level of the original definition ``f = b o g o a``.

Nothing here uses the reduction lemma; the lemma is instead *tested* by check 1,
which compares brute-force affine implementation against ``N_1 = 1`` for all
14 x 14 ordered pairs of classes.
"""

from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from affine_atoms import P, enumerate_classes  # noqa: E402

POINTS2 = tuple(product(range(P), repeat=2))


def brute_force_implements(g_table, f_table) -> bool:
    """True iff ``f = b o g o a`` for some affine ``a : F_3^2 -> F_3^2``, ``b``."""

    for entries in product(range(P), repeat=4):
        for shift in product(range(P), repeat=2):
            a, b_, c, d = entries
            s, t = shift
            composed = []
            for x, y in POINTS2:
                u = (a * x + b_ * y + s) % P
                v = (c * x + d * y + t) % P
                composed.append(g_table[u * P + v])
            for scale in range(P):
                for translate in range(P):
                    if all((scale * value + translate) % P == target
                           for value, target in zip(composed, f_table)):
                        return True
    return False


def check_one_step() -> None:
    classes = enumerate_classes()
    data = json.loads((Path(__file__).resolve().parent / "n1_matrix.json").read_text())
    mismatches = 0
    for g in classes:
        for f in classes:
            direct = brute_force_implements(g.table, f.table)
            reduced = data["n1"][g.key][f.key] in (0, 1)
            if direct != reduced:
                mismatches += 1
                print(f"  MISMATCH {g.key} -> {f.key}: brute {direct}, N_1 {reduced}")
    print(f"check 1: brute-force affine implementation vs (N_1 <= 1) on "
          f"{len(classes) ** 2} ordered pairs: {mismatches} mismatches")


def check_parabolic_implements_rank1() -> None:
    """x^2 = b(g(a(x,y))) with g = x^2+y, a(x,y) = (x,0), b = identity."""

    ok = all(((x * x) % P) == (((x * x) % P + 0) % P) for x, _ in POINTS2)
    values = [((x * x + 0) % P) for x, y in POINTS2]
    target = [(x * x) % P for x, y in POINTS2]
    assert ok and values == target
    print("check 2: x^2 = (x^2+y) o a with a(x,y) = (x,0) — verified on all 9 points")


def check_k2_witness() -> None:
    """f^{x2} <=_aff g^{x3} for g = x^2+y^2, f = x, from the computed witness.

    alpha_1(x) = (x_2, x_1), alpha_2(x) = (x_2, x_1+1), alpha_3(x) = (x_2+1, x_1);
    x_1 = h_1 + 2 h_2 + 1 and x_2 = h_1 + 2 h_3 + 1.
    """

    def g(u, v):
        return (u * u + v * v) % P

    bad = 0
    for x1, y1, x2, y2 in product(range(P), repeat=4):
        h1 = g(x2, x1)
        h2 = g(x2, (x1 + 1) % P)
        h3 = g((x2 + 1) % P, x1)
        if (h1 + 2 * h2 + 1) % P != x1 % P:
            bad += 1
        if (h1 + 2 * h3 + 1) % P != x2 % P:
            bad += 1
    print(f"check 3: x^{{x2}} <=_aff (x^2+y^2)^{{x3}} witness — {bad} failures "
          f"over all 81 points")


def check_k2_witness_split() -> None:
    """f^{x2} <=_aff g^{x3} for g = x^2+y^2, f = xy."""

    def g(u, v):
        return (u * u + v * v) % P

    bad = 0
    for x1, y1, x2, y2 in product(range(P), repeat=4):
        h1 = g((x2 + y2) % P, (x1 + y1) % P)
        h2 = g((x2 + y2) % P, (x1 + 2 * y1) % P)
        h3 = g((x2 + 2 * y2) % P, (x1 + y1) % P)
        if (h1 + 2 * h2) % P != (x1 * y1) % P:
            bad += 1
        if (h1 + 2 * h3) % P != (x2 * y2) % P:
            bad += 1
    print(f"check 4: (xy)^{{x2}} <=_aff (x^2+y^2)^{{x3}} witness — {bad} failures "
          f"over all 81 points")


def check_tensor_synergy() -> None:
    """c = x^2+y is implemented in one shot by (a (x) b) with a = x^2, b = x."""

    bad = 0
    for x, y in POINTS2:
        # alpha(x,y) = (x, 0, y, 0) into F_3^2 x F_3^2, then B = (1,1)
        first = (x * x) % P            # a-component  x^2 evaluated at (x, 0)
        second = y % P                 # b-component  x     evaluated at (y, 0)
        if (first + second) % P != (x * x + y) % P:
            bad += 1
    print(f"check 5: (x^2+y) <=_aff (x^2 (x) x) witness — {bad} failures over 9 points")


if __name__ == "__main__":
    check_one_step()
    check_parabolic_implements_rank1()
    check_k2_witness()
    check_k2_witness_split()
    check_tensor_synergy()
