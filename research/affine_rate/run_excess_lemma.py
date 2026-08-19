#!/usr/bin/env python3
"""The one-step case of Conjectures E1 and E2 of ISOTROPY.md section 8.

For a set ``S_i`` of rank-at-most-two forms ``q_j`` with supports ``W_j`` whose
span contains a form ``F``, define the excess

    psi = 2 |S_i| - dim  sum_j W_j .

* ``F`` of **rank one**: ``psi >= 1``.  Proved: if the ``W_j`` were independent
  the rank of ``sum c_j q_j`` would be the sum of the ranks, which is at least
  ``|S_i|`` and equal to one only for a single rank-one term -- and then
  ``psi = 2 - 1 = 1`` anyway.
* ``F`` of **rank two and anisotropic**, the ``q_j`` products of linear forms:
  ``psi >= 2``.  The case ``|S_i| = 2`` is the only one not settled by the rank
  count, and it is what this script checks exhaustively over ``F_3``: no two
  products whose supports together span three dimensions have a linear
  combination that is an anisotropic rank-two form.

A ``None`` counterexample is the assertion; anything else refutes E2's one-step
case and would have to be reported.
"""

from __future__ import annotations

from itertools import product

import numpy as np

P = 3
N = 3
VECTORS = [np.array(v, dtype=int) for v in product(range(P), repeat=N)]
HALF = pow(2, P - 2, P)


def rank(rows, width=N):
    A = [[int(x) % P for x in row] for row in rows]
    r = 0
    for column in range(width):
        pivot = None
        for i in range(r, len(A)):
            if A[i][column]:
                pivot = i
                break
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        inverse = pow(A[r][column], P - 2, P)
        A[r] = [(v * inverse) % P for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][column]:
                factor = A[i][column]
                A[i] = [(A[i][j] - factor * A[r][j]) % P for j in range(width)]
        r += 1
    return r


def gram_product(l, m):
    return ((np.outer(l, m) + np.outer(m, l)) * HALF) % P


def anisotropic_rank_two(G):
    if rank(G) != 2:
        return False
    for v in VECTORS:
        if not v.any():
            continue
        if not (G @ v % P).any():        # v lies in the radical
            continue
        if int(v @ G @ v) % P == 0:
            return False
    return True


def main() -> None:
    products = [(l, m, gram_product(l, m))
                for l in VECTORS for m in VECTORS if l.any() and m.any()]
    print(f"products of two nonzero linear forms on F_3^3: {len(products)}")
    counterexample = None
    for l1, m1, G1 in products:
        for l2, m2, G2 in products:
            if rank([l1, m1, l2, m2]) != 3:
                continue
            for c1 in (1, 2):
                for c2 in (1, 2):
                    if anisotropic_rank_two((c1 * G1 + c2 * G2) % P):
                        counterexample = (l1.tolist(), m1.tolist(),
                                          l2.tolist(), m2.tolist(), c1, c2)
                        break
                if counterexample:
                    break
            if counterexample:
                break
        if counterexample:
            break
    print("counterexample to the one-step case of E2:", counterexample)
    print("(None means the exhaustive search confirms it.)")


if __name__ == "__main__":
    main()
