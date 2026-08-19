#!/usr/bin/env python3
"""The diagonal ansatz over F_3, brute-forced against its proof.

Fix the target's own basis ``e_1,...,e_{2k}`` = ``x_1,y_1,...,x_k,y_k`` and
restrict attention to solutions all of whose atom quadratic parts are diagonal
in it.  A diagonal binary form ``alpha e_s^2 + beta e_t^2`` is

* a **product** (a ``xy``-atom) iff ``alpha beta`` is a non-residue, i.e. over
  ``F_3`` iff it is ``a (e_s^2 - e_t^2)``;
* a **``Q``-form** for ``Q`` anisotropic iff ``alpha beta`` is a residue, i.e.
  over ``F_3`` iff it is ``a (e_s^2 + e_t^2)``;

and every rank-one diagonal form ``a e_s^2`` is both.  So in coefficient space
``F_3^{2k}`` the two atom families are

    products : support <= 1, or support 2 with entries (a, -a)   ("edges")
    Q-forms  : support <= 1, or support 2 with entries (a,  a)

and the targets are

    x_i^2 + y_i^2  ->  e_{2i-1} + e_{2i}     (for the anisotropic target)
    x_i y_i        ->  e_{2i-1} - e_{2i}     (for the hyperbolic target,
                                              since 4 x y = (x+y)^2 - (x-y)^2)

Proved in ISOTROPY.md section 3.4: in the first row the minimum is exactly
``2k`` (every connected component of the edge graph needs a rank-one atom
because the target has nonzero coordinate sum there, so
``r >= (2k - c) + c``).  In the second row the "sum" atoms already have nonzero
coordinate sum, a spanning forest suffices, and the minimum is ``ceil(3k/2)``.
This script checks both by exhaustive subset search.
"""

from __future__ import annotations

import sys
from itertools import combinations
from math import ceil

P = 3


def rank(rows, width):
    rows = [row[:] for row in rows]
    r = 0
    for column in range(width):
        pivot = None
        for i in range(r, len(rows)):
            if rows[i][column] % P:
                pivot = i
                break
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inverse = pow(rows[r][column], P - 2, P)
        rows[r] = [(v * inverse) % P for v in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][column] % P:
                factor = rows[i][column]
                rows[i] = [(rows[i][j] - factor * rows[r][j]) % P
                           for j in range(width)]
        r += 1
    return r


def spans(chosen, targets, width):
    return rank(chosen + targets, width) == rank(chosen, width)


def atoms(width, sign):
    """``sign = -1`` gives the products, ``sign = +1`` gives the Q-forms."""

    out = []
    for s in range(width):
        for a in (1, 2):
            row = [0] * width
            row[s] = a
            out.append(row)
    for s in range(width):
        for t in range(s + 1, width):
            for a in (1, 2):
                row = [0] * width
                row[s] = a
                row[t] = (sign * a) % P
                out.append(row)
    return out


def minimum(width, sign, targets, limit):
    pool = atoms(width, sign)
    for r in range(len(targets), limit + 1):
        for combo in combinations(range(len(pool)), r):
            if spans([pool[j] for j in combo], targets, width):
                return r
    return None


def main(kmax: int = 3) -> None:
    print("diagonal ansatz, resource xy (products), target x^2+y^2:")
    for k in range(1, kmax + 1):
        width = 2 * k
        targets = []
        for i in range(k):
            row = [0] * width
            row[2 * i] = 1
            row[2 * i + 1] = 1
            targets.append(row)
        value = minimum(width, -1, targets, 2 * k)
        print(f"  k={k}: minimum = {value}   predicted 2k = {2 * k}"
              f"   {'ok' if value == 2 * k else 'MISMATCH'}")

    print()
    print("diagonal ansatz, resource x^2+y^2 (Q-forms), target xy:")
    for k in range(1, kmax + 1):
        width = 2 * k
        targets = []
        for i in range(k):
            row = [0] * width
            row[2 * i] = 1
            row[2 * i + 1] = P - 1
            targets.append(row)
        predicted = ceil(3 * k / 2)
        value = minimum(width, 1, targets, 2 * k)
        print(f"  k={k}: minimum = {value}   predicted ceil(3k/2) = {predicted}"
              f"   {'ok' if value == predicted else 'MISMATCH'}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
