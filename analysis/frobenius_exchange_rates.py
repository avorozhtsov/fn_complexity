#!/usr/bin/env python3
"""Exchange rates of elliptic-fibration signatures over finite fields.

For a map ``f : A^2 -> A^1`` over ``F_q`` the fiber signature is the list of
point counts of the fibers, so it is a list of Weil numbers.  Two families are
compared here against the linear map, against each other, and against the split
conic:

    E1 : f(x,y) = y^2 - x^3 - x     fibers  y^2 = x^3 + x + c   (large monodromy)
    E0 : f(x,y) = y^2 - x^3         fibers  y^2 = x^3 + c       (sextic twists, CM)

Writing ``N_c`` for the affine point count of the fiber over ``c`` and
``a_c = q - N_c`` for the trace of Frobenius, the script checks:

* ``sum_c a_c = 0`` exactly, because ``sum_c N_c = q^2``;
* the normalised second moment ``q^-2 sum_c a_c^2``, which is ``1`` for the
  large-monodromy family (the semicircle value), ``2`` for the CM family when
  ``q = 1 mod 3``, and ``0`` when ``q = 2 mod 3`` where every fiber is
  supersingular;
* ``Z_f(k) = #(X x_Y ... x_Y X)(F_q)``, the point count of the k-fold fiber
  power, for integer ``k``;
* the Weil-scale expansion ``1 - C(L -> f) ~ 2g / (sqrt(q) log q)``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate_result  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "frobenius_exchange_rates.csv"
PRIMES = [101, 211, 401, 1009, 2003]


def fiber_counts(q: int, linear_coefficient: int) -> np.ndarray:
    """Point counts of the fibers of y^2 - x^3 - a x over F_q."""

    x = np.arange(q)
    y = np.arange(q)
    cubic = (x**3 + linear_coefficient * x) % q
    values = (y[:, None] ** 2 - cubic[None, :]) % q
    return np.bincount(values.ravel(), minlength=q)


def signature(counts: np.ndarray) -> tuple[int, ...]:
    return tuple(sorted(counts.tolist(), reverse=True))


def linear_signature(q: int) -> tuple[int, ...]:
    return tuple([q] * q)


def split_conic_signature(q: int) -> tuple[int, ...]:
    return tuple(sorted([2 * q - 1] + [q - 1] * (q - 1), reverse=True))


def rate(implementer: tuple[int, ...], implemented: tuple[int, ...]):
    return exchange_rate_result(implemented=implemented, implementer=implementer)


def main() -> int:
    rows = []
    print("Frobenius statistics")
    print(
        f"{'q':>5} {'q%3':>4} {'family':>4} {'sum a_c':>9} {'moment2':>9} "
        f"{'min a_c':>9} {'-2 sqrt q':>10} {'distinct':>9}"
    )
    signatures: dict[tuple[int, str], tuple[int, ...]] = {}
    for q in PRIMES:
        for name, coefficient in (("E1", 1), ("E0", 0)):
            counts = fiber_counts(q, coefficient)
            traces = q - counts
            signatures[(q, name)] = signature(counts)
            print(
                f"{q:>5} {q % 3:>4} {name:>4} {int(traces.sum()):>9} "
                f"{float((traces.astype(float) ** 2).sum()) / q / q:>9.4f} "
                f"{int(traces.min()):>9} {-2 * math.sqrt(q):>10.2f} "
                f"{len(set(counts.tolist())):>9}"
            )

    # The partition function at an integer counts points of the fiber power.
    q = PRIMES[0]
    counts = fiber_counts(q, 1)
    for power in (2, 3):
        direct = int((counts.astype(object) ** power).sum())
        print(
            f"\nZ_f({power}) = {direct} = number of {power}-tuples of points of "
            f"A^2(F_{q}) with a common image"
        )

    print("\nExchange rates")
    header = f"{'q':>5} {'pair':>10} {'rate':>13} {'contact':>10} {'(1-C) sqrt(q) log q':>21}"
    print(header)
    for q in PRIMES:
        linear = linear_signature(q)
        conic = split_conic_signature(q)
        pairs = {
            "L->E1": (linear, signatures[(q, "E1")]),
            "E1->L": (signatures[(q, "E1")], linear),
            "L->E0": (linear, signatures[(q, "E0")]),
            "E0->L": (signatures[(q, "E0")], linear),
            "E1->E0": (signatures[(q, "E1")], signatures[(q, "E0")]),
            "E0->E1": (signatures[(q, "E0")], signatures[(q, "E1")]),
            "E1->X": (signatures[(q, "E1")], conic),
            "X->E1": (conic, signatures[(q, "E1")]),
        }
        for label, (implementer, implemented) in pairs.items():
            result = rate(implementer, implemented)
            contact = (
                "0" if result.beta == 0
                else "inf" if math.isinf(result.beta)
                else f"{result.beta:.3f}"
            )
            scaled = (1 - result.rate) * math.sqrt(q) * math.log(q)
            print(f"{q:>5} {label:>10} {result.rate:>13.9f} {contact:>10} {scaled:>21.5f}")
            rows.append([q, q % 3, label, f"{result.rate:.15f}", contact, f"{scaled:.9f}"])
        print()

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["q", "q_mod_3", "pair", "rate", "contact_beta", "one_minus_rate_scaled"])
        writer.writerows(rows)
    print(f"written to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
