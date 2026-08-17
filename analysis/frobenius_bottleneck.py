#!/usr/bin/env python3
"""The two directions of the exchange rate against a linear map over F_q.

Let ``L = (q, ..., q)`` be the flat signature of a linear map and ``f`` any map
``A^2 -> A^1`` over ``F_q``, with fiber counts ``N_c`` and Frobenius traces
``a_c = q - N_c``.  The two directed rates read completely different things.

FORWARD, and it is an endpoint always::

    C(L -> f) = log q / log(max_c N_c),   attained at beta = infinity.

Proof: ``Z_f(beta) = sum_c N_c^beta <= q (max N)^beta``, so with ``A = log q``
and ``B = log max N``, ``R(beta) = (1+beta)A / log Z_f >= (1+beta)A/(A+beta B)``
which exceeds ``A/B`` for every finite beta when ``B > A``, while
``R(infinity) = A/B``.  So this rate sees the largest fiber and nothing else.

REVERSE, with a universal bottleneck.  Since ``sum_c a_c = 0`` identically
(because ``sum_c N_c = q^2``), the first moment drops out and

    1 - C(f -> L) = (m_2 / (2 q log q)) * min_beta beta(beta-1)/(beta+1) + O(q^-3/2)

with ``m_2 = q^-2 sum_c a_c^2``.  The derivative of ``beta(beta-1)/(beta+1)`` is
``(beta^2 + 2 beta - 1)/(beta+1)^2``, so the minimiser is the positive root of
``beta^2 + 2 beta - 1 = 0``::

    beta* = sqrt(2) - 1,     1 - C(f -> L) = (3 - 2 sqrt 2) m_2 / (2 q log q).

Both constants are independent of the family and of the genus.

FLATNESS.  For ``f = y^2 - P(x)`` the signature is flat -- so both rates are
exactly 1 -- precisely when ``P`` is a permutation polynomial of ``F_q``.  This
subsumes the ``q = 2 mod 3`` phenomenon, which is the case ``P = x^3``, since
``x^3`` permutes ``F_q`` iff ``gcd(3, q-1) = 1``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate, exchange_rate_result  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "frobenius_bottleneck.csv"

PRIMES = [211, 503, 1009]
FAMILIES = {
    "g=1  y^2 = x^3 + x + c": [1, 0, 1, 0],
    "g=2  y^2 = x^5 + x + c": [1, 0, 0, 0, 1, 0],
    "g=3  y^2 = x^7 + x + c": [1, 0, 0, 0, 0, 0, 1, 0],
}
PERMUTATION_TESTS = {"x^3": 3, "x^5": 5, "x^7": 7, "x^2": 2}


def fiber_counts(q: int, coefficients: list[int]) -> np.ndarray:
    """Counts of the fibers of f(x,y) = y^2 - P(x), P given by Horner coefficients."""

    x = np.arange(q)
    y = np.arange(q)
    values = np.zeros(q, dtype=np.int64)
    for coefficient in coefficients:
        values = (values * x + coefficient) % q
    fibers = (y[:, None] ** 2 - values[None, :]) % q
    return np.bincount(fibers.ravel(), minlength=q)


def signature(counts: np.ndarray) -> tuple[int, ...]:
    return tuple(sorted(counts[counts > 0].tolist(), reverse=True))


def golden_minimum(objective, low: float, high: float) -> float:
    for _ in range(400):
        first = low + (high - low) / 3
        second = high - (high - low) / 3
        if objective(first) < objective(second):
            high = second
        else:
            low = first
    return (low + high) / 2


def main() -> int:
    argument = golden_minimum(lambda b: b * (b - 1) / (b + 1), 0.0, 2.0)
    kappa = -argument * (argument - 1) / (argument + 1)
    print("analytic constants: minimise beta(beta-1)/(beta+1)")
    print(f"   beta*     {argument:.12f}   sqrt(2) - 1 = {math.sqrt(2) - 1:.12f}")
    print(f"   -minimum  {kappa:.12f}   3 - 2 sqrt 2 = {3 - 2 * math.sqrt(2):.12f}\n")

    rows = []
    print("forward rate is an endpoint:  C(L -> f) = log q / log(max_c N_c)")
    print(f"   {'q':>5} {'family':>24} {'max N':>7} {'C(L->f)':>14} {'log q/log maxN':>16} {'beta*':>6}")
    for q in PRIMES:
        flat = tuple([q] * q)
        for name, coefficients in FAMILIES.items():
            sig = signature(fiber_counts(q, coefficients))
            result = exchange_rate_result(implemented=sig, implementer=flat)
            predicted = math.log(q) / math.log(max(sig))
            tag = "inf" if math.isinf(result.beta) else f"{result.beta:.2f}"
            print(f"   {q:>5} {name:>24} {max(sig):>7} {result.rate:>14.10f} "
                  f"{predicted:>16.10f} {tag:>6}")
            rows.append([q, name, "C(L->f)", f"{result.rate:.15f}", f"{predicted:.15f}", tag])

    print("\nreverse rate has a universal bottleneck at beta* = sqrt(2) - 1")
    print(f"   {'q':>5} {'family':>24} {'beta* observed':>15} {'(1-C) 2q log q / m2':>21}")
    for q in PRIMES:
        flat = tuple([q] * q)
        for name, coefficients in FAMILIES.items():
            counts = fiber_counts(q, coefficients)
            sig = signature(counts)
            traces = (q - counts).astype(float)
            second_moment = float((traces**2).sum()) / q / q
            result = exchange_rate_result(implemented=flat, implementer=sig)
            scaled = (1 - result.rate) * 2 * q * math.log(q) / second_moment
            print(f"   {q:>5} {name:>24} {result.beta:>15.6f} {scaled:>21.6f}")
            rows.append([q, name, "C(f->L)", f"{result.rate:.15f}",
                         f"{second_moment:.15f}", f"{result.beta:.6f}"])
    print(f"   targets: beta* = {math.sqrt(2) - 1:.6f}, ratio = {3 - 2 * math.sqrt(2):.6f}")

    print("\nflatness: sig(y^2 - P(x)) is flat  <=>  P is a permutation polynomial")
    print(f"   {'q':>5} {'P':>6} {'permutation':>12} {'flat':>6} {'C(f->L)':>13} {'C(L->f)':>13}")
    violations = 0
    for q in (13, 17, 19, 23, 31, 37, 41, 43, 211):
        x = np.arange(q)
        flat = tuple([q] * q)
        for name, power in PERMUTATION_TESTS.items():
            values = (x.astype(object) ** power) % q
            sig = signature(fiber_counts(q, [1] + [0] * power))
            is_permutation = len(set(values.tolist())) == q
            is_flat = len(set(sig)) == 1 and sig[0] == q
            forward = exchange_rate(sig, flat)
            backward = exchange_rate(flat, sig)
            if is_permutation != is_flat:
                violations += 1
            if q in (13, 31, 211):
                print(f"   {q:>5} {name:>6} {str(is_permutation):>12} {str(is_flat):>6} "
                      f"{forward:>13.10f} {backward:>13.10f}")
            rows.append([q, name, "flatness", str(is_permutation), str(is_flat),
                         f"{forward:.15f}"])
    print(f"   disagreements between 'permutation' and 'flat': {violations}")
    print("   the q = 2 mod 3 case is P = x^3, since x^3 permutes F_q iff gcd(3,q-1) = 1")

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["q", "family", "quantity", "value", "reference", "contact"])
        writer.writerows(rows)
    print(f"\nwritten to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
