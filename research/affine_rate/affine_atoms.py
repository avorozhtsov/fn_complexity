#!/usr/bin/env python3
"""Affine-implementation rate ``C_aff`` for polynomial maps ``F_3^2 -> F_3``.

Everything here rests on the reduction lemma proved in ``FINDINGS.md``:

    f^{x k}  <=_aff  g^{x r}
      <=>  there are affine maps  alpha_1, ..., alpha_r : F_3^{2k} -> F_3^2
           such that every  f(x_i),  i = 1..k,  lies in
           span_{F_3}{ g o alpha_1, ..., g o alpha_r }  +  F_3 . 1 .

So the whole search collapses to linear algebra over F_3 inside the space of
functions F_3^{2k} -> F_3.  We call ``g o alpha`` a *g-atom* and define

    N_k(g -> f) = min { |T| : T a set of g-atoms on F_3^{2k},
                        f(x_1), ..., f(x_k)  in  span(T) + F_3.1 } ,

so that  k_{g->f}(r) = max{ k : N_k <= r }  and, N_k being subadditive,

    C_aff(g -> f) = 1 / inf_k ( N_k / k ) = sup_k k / N_k .
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from functools import lru_cache

import numpy as np

P = 3

# ---------------------------------------------------------------------------
# degree-at-most-three functions F_3^2 -> F_3, in the repository's basis
# ---------------------------------------------------------------------------

# x^3 = x and y^3 = y as functions, so these eight monomials are a basis of the
# degree-at-most-three polynomial functions (x^2 y^2 is the missing reduced
# monomial and has degree four).
CUBIC_BASIS: tuple[tuple[int, int], ...] = (
    (0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (2, 1), (1, 2),
)

POINTS2 = tuple(product(range(P), repeat=2))


def values_of_coefficients(coefficients, basis=CUBIC_BASIS) -> tuple[int, ...]:
    """Value table of a polynomial on F_3^2, in the fixed point order."""

    table = []
    for x, y in POINTS2:
        total = 0
        for coefficient, (i, j) in zip(coefficients, basis):
            total += coefficient * pow(x, i) * pow(y, j)
        table.append(total % P)
    return tuple(table)


def signature(table) -> tuple[int, ...]:
    counts = [0] * P
    for value in table:
        counts[value] += 1
    return tuple(sorted((c for c in counts if c), reverse=True))


# ---------------------------------------------------------------------------
# affine group of F_3^2 and the affine equivalence classes
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def invertible_affine_2() -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    result = []
    for entries in product(range(P), repeat=4):
        a, b, c, d = entries
        if (a * d - b * c) % P == 0:
            continue
        for shift in product(range(P), repeat=2):
            result.append((entries, shift))
    return tuple(result)


def compose_table(table, matrix_entries, shift) -> tuple[int, ...]:
    """Value table of ``h o a`` where ``a(x) = matrix x + shift`` on F_3^2."""

    a, b, c, d = matrix_entries
    s, t = shift
    out = []
    for x, y in POINTS2:
        u = (a * x + b * y + s) % P
        v = (c * x + d * y + t) % P
        out.append(table[u * P + v])
    return tuple(out)


def affine_class_of(table) -> tuple[int, ...]:
    """Canonical representative of the affine left-right class of ``table``."""

    best = None
    for matrix_entries, shift in invertible_affine_2():
        moved = compose_table(table, matrix_entries, shift)
        for scale in (1, 2):
            for translate in range(P):
                candidate = tuple((scale * v + translate) % P for v in moved)
                if best is None or candidate < best:
                    best = candidate
    return best


@dataclass(frozen=True)
class MapClass:
    key: str
    representative: str
    table: tuple[int, ...]
    size: int
    signature: tuple[int, ...]
    degree: int


REPRESENTATIVE_NAMES = {
    # coefficient tuples in CUBIC_BASIS -> (key, printed representative)
    (0, 0, 0, 0, 0, 0, 0, 0): ("constant", "0"),
    (0, 1, 0, 0, 0, 0, 0, 0): ("linear", "x"),
    (0, 0, 0, 1, 0, 0, 0, 0): ("rank1", "x^2"),
    (0, 0, 0, 0, 1, 0, 0, 0): ("split", "xy"),
    (0, 0, 0, 1, 0, 1, 0, 0): ("anisotropic", "x^2+y^2"),
    (0, 0, 1, 1, 0, 0, 0, 0): ("parabolic", "x^2+y"),
    (0, 0, 0, 0, 1, 0, 0, 1): ("cubic-711", "xy+xy^2"),
    (0, 0, 1, 0, 0, 1, 1, 0): ("cubic-63", "y+y^2+x^2y"),
    (0, 0, 0, 0, 0, 1, 1, 0): ("cubic-522-a", "y^2+x^2y"),
    (0, 0, 0, 0, 0, 0, 0, 1): ("cubic-522-b", "xy^2"),
    (0, 0, 0, 0, 1, 1, 1, 0): ("cubic-441-a", "xy+y^2+x^2y"),
    (0, 0, 1, 0, 1, 0, 0, 1): ("cubic-441-b", "y+xy+xy^2"),
    (0, 0, 1, 0, 0, 0, 1, 0): ("cubic-333-a", "y+x^2y"),
    (0, 1, 1, 0, 0, 1, 1, 0): ("cubic-333-b", "x+y+y^2+x^2y"),
}


def polynomial_degree(coefficients, basis=CUBIC_BASIS) -> int:
    return max((i + j for coefficient, (i, j) in zip(coefficients, basis) if coefficient),
               default=0)


def enumerate_classes() -> list[MapClass]:
    """All 14 affine classes of degree-at-most-three functions F_3^2 -> F_3."""

    canonical_to_members: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for coefficients in product(range(P), repeat=len(CUBIC_BASIS)):
        table = values_of_coefficients(coefficients)
        canonical = affine_class_of(table)
        canonical_to_members.setdefault(canonical, []).append(coefficients)

    classes = []
    for canonical, members in canonical_to_members.items():
        named = None
        for coefficients in members:
            if coefficients in REPRESENTATIVE_NAMES:
                named = coefficients
                break
        if named is None:
            raise AssertionError(f"unnamed class with canonical form {canonical}")
        key, printed = REPRESENTATIVE_NAMES[named]
        table = values_of_coefficients(named)
        classes.append(MapClass(
            key=key,
            representative=printed,
            table=table,
            size=len(members),
            signature=signature(table),
            degree=polynomial_degree(named),
        ))
    classes.sort(key=lambda c: (c.degree, -max(c.signature), c.key))
    return classes


# ---------------------------------------------------------------------------
# linear algebra over F_3, with pivots carried alongside
# ---------------------------------------------------------------------------

INVERSE = {1: 1, 2: 2}


class Subspace:
    """Row space over F_3 in reduced row echelon form."""

    __slots__ = ("rows", "pivots", "width")

    def __init__(self, width: int):
        self.width = width
        self.rows = np.zeros((0, width), dtype=np.int8)
        self.pivots: list[int] = []

    def copy(self) -> "Subspace":
        other = Subspace(self.width)
        other.rows = self.rows.copy()
        other.pivots = list(self.pivots)
        return other

    @property
    def dimension(self) -> int:
        return len(self.pivots)

    def reduce(self, matrix: np.ndarray) -> np.ndarray:
        """Reduce every row of ``matrix`` modulo the subspace."""

        result = np.array(matrix, dtype=np.int8, copy=True) % P
        if result.ndim == 1:
            result = result.reshape(1, -1)
        for row, column in zip(self.rows, self.pivots):
            coefficient = result[:, column].copy()
            nonzero = coefficient != 0
            if nonzero.any():
                result[nonzero] = (result[nonzero]
                                   - coefficient[nonzero, None] * row[None, :]) % P
        return result

    def contains(self, vector: np.ndarray) -> bool:
        return not self.reduce(vector).any()

    def extended(self, vector: np.ndarray) -> "Subspace":
        """New subspace with ``vector`` adjoined (assumed independent)."""

        residual = self.reduce(vector)[0]
        column = int(np.flatnonzero(residual)[0])
        residual = (residual * INVERSE[int(residual[column])]) % P
        other = Subspace(self.width)
        rows = self.rows.copy()
        if rows.shape[0]:
            coefficient = rows[:, column].copy()
            nonzero = coefficient != 0
            if nonzero.any():
                rows[nonzero] = (rows[nonzero]
                                 - coefficient[nonzero, None] * residual[None, :]) % P
        order = np.argsort(self.pivots + [column])
        other.rows = np.vstack([rows, residual.reshape(1, -1)])[order]
        other.pivots = sorted(self.pivots + [column])
        return other

    def extended_many(self, matrix: np.ndarray) -> "Subspace":
        other = self
        for row in np.atleast_2d(matrix):
            if not other.contains(row):
                other = other.extended(row)
        return other


# ---------------------------------------------------------------------------
# atoms and the quantity N_k
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def block_points(k: int) -> np.ndarray:
    return np.array(tuple(product(range(P), repeat=2 * k)), dtype=np.int8)


def atom_tables(table, k: int) -> np.ndarray:
    """Distinct value tables of ``g o alpha``, alpha affine F_3^{2k} -> F_3^2."""

    points = block_points(k).astype(np.int16)
    lookup = np.array(table, dtype=np.int8)
    seen: set[bytes] = set()
    out: list[np.ndarray] = []
    for entries in product(range(P), repeat=4 * k):
        matrix = np.array(entries, dtype=np.int16).reshape(2, 2 * k)
        image = (points @ matrix.T) % P
        for shift0 in range(P):
            u = (image[:, 0] + shift0) % P
            for shift1 in range(P):
                v = (image[:, 1] + shift1) % P
                values = lookup[u * P + v]
                key = values.tobytes()
                if key not in seen:
                    seen.add(key)
                    out.append(values)
    return np.array(out, dtype=np.int8)


def target_tables(table, k: int) -> np.ndarray:
    """Value tables of ``f(x_1), ..., f(x_k)`` on F_3^{2k}."""

    points = block_points(k)
    lookup = np.array(table, dtype=np.int8)
    rows = [lookup[points[:, 2 * i] * P + points[:, 2 * i + 1]] for i in range(k)]
    return np.array(rows, dtype=np.int8)


def minimal_atom_count(atoms: np.ndarray, targets: np.ndarray,
                       limit: int) -> int | None:
    """Smallest ``t <= limit`` with the targets inside ``span(t atoms) + F_3.1``.

    Returns ``None`` when no such ``t <= limit`` exists.  The search is an
    exhaustive depth-first enumeration of atom subsets in increasing index
    order with two exact prunings, so the returned value is a proved minimum
    and a returned ``None`` proves ``N_k > limit``.
    """

    width = atoms.shape[1]
    start = Subspace(width).extended(np.ones(width, dtype=np.int8))

    def uncovered(basis: Subspace) -> int:
        return basis.extended_many(targets).dimension - basis.dimension

    def search(basis: Subspace, first: int, budget: int) -> bool:
        need = uncovered(basis)
        if need == 0:
            return True
        if need > budget:
            return False
        hull = basis.extended_many(targets)
        if need == budget:
            # span(basis) + span(chosen) must equal span(basis) + span(targets),
            # so every chosen atom lies in the hull; greedy spanning is exact.
            residual = hull.reduce(atoms[first:])
            candidates = np.flatnonzero(~residual.any(axis=1)) + first
            current = basis
            for j in candidates:
                if not current.contains(atoms[j]):
                    current = current.extended(atoms[j])
                    if current.dimension - basis.dimension == need:
                        return True
            return False
        residual = basis.reduce(atoms[first:])
        useful = np.flatnonzero(residual.any(axis=1)) + first
        for j in useful:
            if search(basis.extended(atoms[j]), int(j) + 1, budget - 1):
                return True
        return False

    for t in range(uncovered(start), limit + 1):
        if search(start, 0, t):
            return t
    return None


# ---------------------------------------------------------------------------
# exact N_1 by breadth-first search inside the 3^9 function space
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _digit_table():
    powers = (P ** np.arange(9)).astype(np.int64)
    codes = np.arange(P ** 9, dtype=np.int64)
    digits = (codes[:, None] // powers[None, :]) % P
    return digits.astype(np.int8), powers


def encode9(vector) -> int:
    _, powers = _digit_table()
    return int((np.asarray(vector, dtype=np.int64) * powers).sum())


def minimal_atom_count_k1(atoms, target, limit: int = 9):
    """Exact ``N_1(g -> f)`` by BFS over all 3^9 functions F_3^2 -> F_3.

    ``atoms`` are the value tables of the g-atoms; the search starts from the
    three constant functions and adds one scaled atom per step.  Exhaustive,
    so ``None`` proves ``N_1 > limit`` (and with ``limit = 9``, ``N_1`` infinite).
    """

    digits, powers = _digit_table()
    shifts = []
    seen = set()
    for atom in atoms:
        for scale in (1, 2):
            vector = (np.asarray(atom, dtype=np.int64) * scale) % P
            key = vector.tobytes()
            if key in seen or not vector.any():
                continue
            seen.add(key)
            shifts.append((((digits + vector[None, :].astype(np.int8)) % P
                            ).astype(np.int64) @ powers))

    target_code = encode9(target)
    reached = np.zeros(P ** 9, dtype=bool)
    for constant in range(P):
        reached[encode9(np.full(9, constant, dtype=np.int8))] = True
    if reached[target_code]:
        return 0
    frontier = reached.copy()
    for step in range(1, limit + 1):
        new = np.zeros_like(reached)
        indices = np.flatnonzero(frontier)
        for permutation in shifts:
            new[permutation[indices]] = True
        new &= ~reached
        if not new.any():
            return None
        reached |= new
        if reached[target_code]:
            return step
        frontier = new
    return None
