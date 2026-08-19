#!/usr/bin/env python3
"""Explicit witnesses for the computed values of ``N_2``.

For a pair ``(g, f)`` and a budget ``t`` this finds affine processors
``alpha_1, ..., alpha_t : F_3^4 -> F_3^2`` and coefficients ``lambda, mu`` with

    f(x_1) = sum_j lambda_j  g(alpha_j(x)) + c_1 ,
    f(x_2) = sum_j mu_j      g(alpha_j(x)) + c_2 ,

which is exactly an affine implementation ``f^{x2} <=_aff g^{xt}``.  Each
witness is re-verified from scratch on all 81 points of ``F_3^4``.
"""

from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from affine_atoms import (  # noqa: E402
    P, Subspace, block_points, enumerate_classes, target_tables,
)

VARIABLES = ("x1", "y1", "x2", "y2")


def labelled_atoms(table, k: int):
    """Atom value tables together with the affine map that produced them."""

    points = block_points(k).astype(np.int16)
    lookup = np.array(table, dtype=np.int8)
    seen: dict[bytes, tuple] = {}
    rows, labels = [], []
    for entries in product(range(P), repeat=4 * k):
        matrix = np.array(entries, dtype=np.int16).reshape(2, 2 * k)
        image = (points @ matrix.T) % P
        for shift0 in range(P):
            for shift1 in range(P):
                u = (image[:, 0] + shift0) % P
                v = (image[:, 1] + shift1) % P
                values = lookup[u * P + v]
                key = values.tobytes()
                if key in seen:
                    continue
                seen[key] = (matrix.copy(), (shift0, shift1))
                rows.append(values)
                labels.append((matrix.copy(), (shift0, shift1)))
    return np.array(rows, dtype=np.int8), labels


def form(row, shift) -> str:
    parts = [f"{int(c)}{v}" for c, v in zip(row, VARIABLES) if c]
    if shift:
        parts.append(str(int(shift)))
    return " + ".join(parts) if parts else "0"


def search_witness(atoms, targets, budget, restrict=None):
    width = atoms.shape[1]
    order = np.arange(atoms.shape[0]) if restrict is None else np.asarray(restrict)
    pool = atoms[order]
    start = Subspace(width).extended(np.ones(width, dtype=np.int8))

    def uncovered(basis):
        return basis.extended_many(targets).dimension - basis.dimension

    def search(basis, first, budget, chosen):
        if uncovered(basis) == 0:
            return list(chosen)
        if uncovered(basis) > budget:
            return None
        residual = basis.reduce(pool[first:])
        useful = np.flatnonzero(residual.any(axis=1)) + first
        for j in useful:
            found = search(basis.extended(pool[j]), int(j) + 1, budget - 1,
                           chosen + [int(j)])
            if found is not None:
                return found
        return None

    found = search(start, 0, budget, [])
    return None if found is None else [int(order[j]) for j in found]


def hull_restricted_witness(atoms, targets, budget):
    """Budget-3 search restricted to the only hulls that can carry a solution."""

    width = atoms.shape[1]
    keys = {row.tobytes(): j for j, row in enumerate(atoms)}
    hull = Subspace(width).extended(np.ones(width, dtype=np.int8)).extended_many(targets)
    spaces = [hull]
    seen = set()
    for atom in atoms:
        if hull.contains(atom):
            continue
        bigger = hull.extended(atom)
        key = bigger.rows.tobytes()
        if key not in seen:
            seen.add(key)
            spaces.append(bigger)
    for space in spaces:
        elements = (np.array(list(product(range(P), repeat=space.rows.shape[0])),
                             dtype=np.int8) @ space.rows) % P
        pool = [keys[row.tobytes()] for row in elements if row.tobytes() in keys]
        if not pool:
            continue
        found = search_witness(atoms, targets, budget, restrict=sorted(set(pool)))
        if found is not None:
            return found
    return None


def solve_coefficients(atoms, indices, target):
    """Solve ``target = sum_j c_j atoms[j] + c_0`` over F_3 by brute force."""

    for coefficients in product(range(P), repeat=len(indices)):
        combination = np.zeros(atoms.shape[1], dtype=np.int8)
        for c, j in zip(coefficients, indices):
            combination = (combination + c * atoms[j]) % P
        difference = (target - combination) % P
        if len(set(difference.tolist())) == 1:
            return coefficients, int(difference[0])
    raise AssertionError("no coefficients")


def report(g_key: str, f_key: str, budget: int) -> None:
    classes = {c.key: c for c in enumerate_classes()}
    atoms, labels = labelled_atoms(classes[g_key].table, 2)
    targets = target_tables(classes[f_key].table, 2)
    indices = (hull_restricted_witness(atoms, targets, budget) if budget == 3
               else search_witness(atoms, targets, budget))
    if indices is None:
        print(f"no witness with {budget} atoms for {g_key} -> {f_key}")
        return
    print(f"N_2({g_key} = {classes[g_key].representative} -> "
          f"{f_key} = {classes[f_key].representative}) <= {len(indices)}")
    for slot, j in enumerate(indices, start=1):
        matrix, shift = labels[j]
        print(f"  alpha_{slot}(x) = ( {form(matrix[0], shift[0])} ,"
              f" {form(matrix[1], shift[1])} )")
    for row, name in zip(targets, ("f(x_1)", "f(x_2)")):
        coefficients, constant = solve_coefficients(atoms, indices, row)
        terms = " + ".join(f"{c}*h_{s}" for s, c in enumerate(coefficients, start=1) if c)
        print(f"  {name} = {terms}" + (f" + {constant}" if constant else ""))
        # independent re-verification on all 81 points
        check = np.full(atoms.shape[1], constant, dtype=np.int8)
        for c, j in zip(coefficients, indices):
            check = (check + c * atoms[j]) % P
        assert np.array_equal(check, row % P), "witness failed re-verification"
    print("  re-verified on all 81 points of F_3^4")


if __name__ == "__main__":
    for g_key, f_key, budget in (
        ("anisotropic", "linear", 3),
        ("anisotropic", "split", 3),
        ("rank1", "linear", 4),
        ("rank1", "split", 4),
        ("rank1", "anisotropic", 4),
    ):
        report(g_key, f_key, budget)
        print()
