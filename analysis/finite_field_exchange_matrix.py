#!/usr/bin/env python3
"""Exchange matrix of the quadratic-map classes over a finite field.

For every prime power ``q`` the degree-at-most-two maps ``F_q^2 -> F_q`` fall
into six classes for odd ``q`` and seven for even ``q``.  Their fiber
signatures are point counts, so each is a Weil number of the corresponding
conic:

    constant                  {q^2}
    linear x                  {q, ..., q}          (q fibers)
    parabolic x^2 + y         {q, ..., q}          (q fibers)  -- same as linear
    pure square x^2   (odd)   {2q, ..., 2q, q}     ((q-1)/2 fibers of size 2q)
    x^2   (even)              {q, ..., q}          -- Frobenius is a bijection
    x^2 + x   (even)          {2q, ..., 2q}        (q/2 fibers)
    split xy                  {2q-1, q-1, ..., q-1}
    anisotropic               {q+1, ..., q+1, 1}

Distinct signatures are therefore ``K`` (constant), ``L`` (linear/parabolic,
and ``x^2`` in even characteristic), ``S`` (``x^2`` for odd ``q``, ``x^2+x``
for even ``q``), ``X`` (split) and ``A`` (anisotropic).  ``K`` is degenerate:
it has a single fiber, so ``log Z_K(0) = 0``.

The script writes the four-by-four matrix of the non-degenerate classes for a
range of ``q``, marks which entries are attained at an endpoint, and checks the
closed forms and asymptotic expansions claimed in
``paper_finite_fields_maps/docs/finite_field_exchange_matrix.md``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate_result  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "finite_field_exchange_matrix.csv"

CLASS_ORDER = ["S", "L", "A", "X"]
PRIME_POWERS = [3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, 31, 49, 64, 101, 121, 256, 509, 1009]


def signatures(q: int) -> dict[str, tuple[int, ...]]:
    """Fiber signatures of the non-degenerate quadratic-map classes over F_q."""

    if q % 2:
        pure_square = tuple(sorted([2 * q] * ((q - 1) // 2) + [q], reverse=True))
    else:
        pure_square = tuple([2 * q] * (q // 2))
    return {
        "S": pure_square,
        "L": tuple([q] * q),
        "X": tuple(sorted([2 * q - 1] + [q - 1] * (q - 1), reverse=True)),
        "A": tuple(sorted([q + 1] * (q - 1) + [1], reverse=True)),
    }


def endpoint_tag(beta: float) -> str:
    if beta == 0.0:
        return "0"
    if math.isinf(beta):
        return "inf"
    return f"{beta:.4f}"


def golden_maximum(objective, low: float, high: float) -> tuple[float, float]:
    for _ in range(400):
        first = low + (high - low) / 3
        second = high - (high - low) / 3
        if objective(first) < objective(second):
            low = first
        else:
            high = second
    argument = (low + high) / 2
    return objective(argument), argument


def closed_forms(q: int) -> dict[tuple[str, str], float]:
    """Entries with an exact endpoint value."""

    log = math.log
    fiber_count_S = (q - 1) // 2 + 1 if q % 2 else q // 2
    return {
        ("L", "S"): log(q) / log(2 * q),
        ("L", "X"): log(q) / log(2 * q - 1),
        ("L", "A"): log(q) / log(q + 1),
        ("S", "L"): log(fiber_count_S) / log(q),
        ("S", "X"): log(fiber_count_S) / log(q),
        ("S", "A"): log(fiber_count_S) / log(q),
        ("A", "S"): log(q + 1) / log(2 * q),
        ("A", "X"): log(q + 1) / log(2 * q - 1),
    }


def complexity_index(signature: tuple[int, ...]) -> float:
    """``log(#fibers) * log(max fiber)``.

    In the endpoint regime -- both directed rates attained at ``beta = 0`` or
    ``beta = infinity`` -- one has ``a < b`` exactly when this number is
    smaller for ``a``, so the comparison is a total preorder and no cycle can
    occur.  The quadratic and cubic classes lie close enough to that regime for
    the index to reproduce every comparison, including the flip at ``q = 3``.
    """

    return math.log(len(signature)) * math.log(max(signature))


CUBIC_F3_SIGNATURES = [(3, 3, 3), (6, 3), (4, 4, 1), (5, 2, 2), (7, 1, 1)]


def report_complexity_index() -> None:
    print("\ncomplexity index phi = log(#fibers) * log(max fiber)")
    for q in (3, 4, 5, 8, 9, 16, 25, 101):
        sig = signatures(q)
        ranked = sorted(CLASS_ORDER, key=lambda key: complexity_index(sig[key]))
        computed = sorted(
            CLASS_ORDER,
            key=lambda name: sum(
                1
                for other in CLASS_ORDER
                if exchange_rate_result(
                    implemented=sig[other], implementer=sig[name]
                ).rate
                < exchange_rate_result(
                    implemented=sig[name], implementer=sig[other]
                ).rate
            ),
            reverse=True,
        )
        flag = "matches" if ranked == computed else "DIFFERS"
        print(f"   q={q:>4}  phi order {' < '.join(ranked)}   {flag}")

    print("\n   cubic maps over F_3, five distinct signatures")
    ranked = sorted(CUBIC_F3_SIGNATURES, key=complexity_index)
    print("   phi order " + " < ".join(str(s) for s in ranked))
    violations = [
        (a, b)
        for index, a in enumerate(CUBIC_F3_SIGNATURES)
        for b in CUBIC_F3_SIGNATURES[index + 1 :]
        if (exchange_rate_result(implemented=b, implementer=a).rate
            < exchange_rate_result(implemented=a, implementer=b).rate)
        != (complexity_index(a) < complexity_index(b))
    ]
    print(f"   phi violations: {violations or 'none'}")


def main() -> int:
    lam, lam_beta = golden_maximum(lambda b: (b + 1 - 2.0**b) / (b + 1), 1e-9, 1.0)
    kappa, kappa_beta = golden_maximum(lambda b: (2 * b - 2.0**b) / (b + 1), 0.5, 3.0)
    print(f"lambda = {lam:.12f} at beta = {lam_beta:.12f}   (X -> L)")
    print(f"kappa  = {kappa:.12f} at beta = {kappa_beta:.12f}   (X -> A)\n")

    rows = []
    worst_closed_form = 0.0
    order_exceptions = []
    for q in PRIME_POWERS:
        sig = signatures(q)
        exact = closed_forms(q)
        rates = {}
        for g in CLASS_ORDER:
            for f in CLASS_ORDER:
                result = exchange_rate_result(implemented=sig[f], implementer=sig[g])
                rates[(g, f)] = result
                rows.append(
                    [
                        q,
                        g,
                        f,
                        f"{result.rate:.15f}",
                        endpoint_tag(result.beta),
                        f"{exact[(g, f)]:.15f}" if (g, f) in exact else "",
                    ]
                )
                if (g, f) in exact:
                    worst_closed_form = max(
                        worst_closed_form, abs(result.rate - exact[(g, f)])
                    )

        order = sorted(
            CLASS_ORDER,
            key=lambda name: sum(
                1 for other in CLASS_ORDER if rates[(name, other)].rate < rates[(other, name)].rate
            ),
            reverse=True,
        )
        if order != ["S", "L", "A", "X"]:
            order_exceptions.append((q, order))

        print(f"q = {q}")
        header = "      " + "".join(f"{name:>14}" for name in CLASS_ORDER)
        print(header)
        for g in CLASS_ORDER:
            cells = "".join(f"{rates[(g, f)].rate:>14.8f}" for f in CLASS_ORDER)
            print(f"  {g} ->{cells}")
        print(
            "        contact: "
            + ", ".join(
                f"{g}->{f}:{endpoint_tag(rates[(g, f)].beta)}"
                for g in CLASS_ORDER
                for f in CLASS_ORDER
                if g != f
            )
        )
        print(f"        comparison order (least first): {' < '.join(order)}\n")

    print(f"max deviation from the closed forms: {worst_closed_form:.3e}")
    if order_exceptions:
        print("comparison order differs from S < L < A < X at:", order_exceptions)
    else:
        print("comparison order is S < L < A < X at every tabulated q")

    report_complexity_index()

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["q", "implementer", "implemented", "rate", "contact_beta", "closed_form"])
        writer.writerows(rows)
    print(f"written to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
