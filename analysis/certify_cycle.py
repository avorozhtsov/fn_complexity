#!/usr/bin/env python3
"""Interval-arithmetic check of the three numerical comparisons in Appendix C.

Every claim in Appendix C is elementary except three finite computations.
This script verifies them with mpmath's interval arithmetic, so the printed
verdicts are rigorous rather than floating point.

    python analysis/certify_cycle.py
"""
from mpmath import iv, mp

iv.dps = 50
mp.dps = 50

L2, L3, L5, L6 = (iv.log(iv.mpf(n)) for n in (2, 3, 5, 6))
c = L2 / L3
theta = L5 / L3
alpha = L5 / L6
inv_c = L3 / L2

ok = True


def check(name, holds, detail):
    global ok
    ok = ok and holds
    print(f"  [{'ok' if holds else 'FAIL'}] {name}: {detail}")


print("Appendix C, the three numerical comparisons\n")

# 1. Lemma C.3: A(12) > 1 and B(12) > 1.
s = iv.mpf(12)
A = s ** (inv_c - 1) - s ** (theta - 1)
B = (2 / c) * s ** (inv_c - 2)
check("A(12) > 1", A.a > 1, f"A(12) in [{A.a}, {A.b}]")
check("B(12) > 1", B.a > 1, f"B(12) in [{B.a}, {B.b}]")

# ... and the exponent ordering 1/c - 1 > theta - 1, i.e. (log 3)^2 > log 2 log 5.
d = L3 ** 2 - L2 * L5
check("(log 3)^2 > log 2 log 5", d.a > 0, f"difference in [{d.a}, {d.b}]")

# 2. Proposition C.4, edge f2 < f3: gamma < c, i.e. (log 3)^2 < log 2 log 6.
d = L2 * L6 - L3 ** 2
check("(log 3)^2 < log 2 log 6", d.a > 0, f"difference in [{d.a}, {d.b}]")

# 3. Proposition C.4, edge f3 < f1: the single evaluation at beta = 9/20.
b = iv.mpf(9) / 20
ratio = iv.log(iv.mpf(6) ** b + 1) / iv.log(iv.mpf(5) ** b + iv.mpf(3) ** b)
gap = alpha - ratio
check("R(9/20) < log 5 / log 6", gap.a > 0,
      f"R(9/20) in [{ratio.a}, {ratio.b}], margin >= {gap.a}")

print("\n" + ("all comparisons verified" if ok else "VERIFICATION FAILED"))
raise SystemExit(0 if ok else 1)
