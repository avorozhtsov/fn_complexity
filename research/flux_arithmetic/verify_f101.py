#!/usr/bin/env python3
"""Independent re-verification of the ``F_101`` cycle of brief B's addendum 2.

Signatures recomputed here by point count, flux computed on this thread's grid
engine, rates cross-checked against the package solver.
"""

from __future__ import annotations

import sys

import numpy as np

import curves
import flux

sys.path.insert(0, "../../src")
from fn_complexity import exchange_rate_result  # noqa: E402

P = 101
PENCILS = {
    "f_1": {5: 1, 4: 70, 3: 28, 2: 15, 1: 11},
    "f_2": {5: 1, 4: 42, 3: 32, 2: 74, 1: 96},
    "f_3": {5: 1, 4: 72, 3: 21, 2: 2, 1: 6},
}


def main() -> None:
    names = list(PENCILS)
    sigs = {}
    for name, coeffs in PENCILS.items():
        fibers = curves.pencil_fibers(coeffs, P, 5)
        sig = tuple(sorted((f.n_affine for f in fibers), reverse=True))
        sigs[name] = sig
        mult = sum(1 for v in sig if v == sig[0])
        print(f"{name}: sum {sum(sig)} (= q^2 = {P*P}), max {sig[0]}, "
              f"multiplicity {mult}, min {sig[-1]}, "
              f"smooth fibers {sum(f.smooth for f in fibers)}/{P}")
    S = [sigs[n] for n in names]
    A = flux.flux_matrix(S, flux.beta_grid(P, 2400, 13600))
    print()
    for i in range(3):
        for j in range(i + 1, 3):
            print(f"A({names[i]},{names[j]}) = {A[i,j]:+.9e}   "
                  f"(> 0 means {names[i]} < {names[j]})")
    curl = A[0, 1] + A[1, 2] + A[2, 0]
    total = abs(A[0, 1]) + abs(A[1, 2]) + abs(A[2, 0])
    print(f"curl = {curl:+.9e}   sum|A| = {total:.9e}   r = {abs(curl)/total:.12f}")
    print()
    for i in range(3):
        for j in range(3):
            if i != j:
                print(f"  C({names[i]}->{names[j]}) = "
                      f"{exchange_rate_result(S[i], S[j]).rate:.12f}")


if __name__ == "__main__":
    main()
