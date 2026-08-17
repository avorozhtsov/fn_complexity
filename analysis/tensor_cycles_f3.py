#!/usr/bin/env python3
"""Cycles in the exchange comparison of homogeneous tensor classes over F_3.

The quadratic and cubic *map* classes over a finite field are totally ordered by
the exchange comparison, and the order agrees with the endpoint index

    phi(a) = log(#fibers) * log(max fiber),

which is a total preorder and therefore admits no cycle.  The index is exact
whenever both directed rates of a pair are attained at an endpoint, so a cycle
can only be built from a pair whose rate is attained at an interior temperature.

This script searches the six homogeneous tensor families over F_3 for such
cycles.  Five of them have none.  Case 3 -- quadratic homogeneous maps
F_3^3 -> F_3^3, the fifty-orbit family -- has seven distinct strict three-cycles,
and every one of them is closed by exactly one phi-violating edge carrying an
interior contact, as predicted.
"""

from __future__ import annotations

from collections import Counter
import csv
import itertools
import math
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate, exchange_rate_result  # noqa: E402
from fn_complexity.homogeneous_tensor_maps import (  # noqa: E402
    TENSOR_CASES,
    _evaluate_form,
    compute_tensor_poset,
    tensor_label,
)

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "tensor_cycles_f3.csv"
CASES = (1, 2, 3, 4, 5)
TOLERANCE = 1e-9


def complexity_index(signature: tuple[int, ...]) -> float:
    return math.log(len(signature)) * math.log(max(signature))


def orbit_signatures(number: int) -> dict[tuple[int, ...], list]:
    """Fiber signature of every orbit representative, grouped by signature."""

    case = TENSOR_CASES[number]
    grouped: dict[tuple[int, ...], list] = {}
    for orbit in compute_tensor_poset(number).orbits:
        counts: Counter = Counter()
        for point in itertools.product(range(3), repeat=case.input_dimension):
            value = tuple(
                _evaluate_form(form, point, case.degree)
                for form in orbit.representative
            )
            counts[value] += 1
        signature = tuple(sorted(counts.values(), reverse=True))
        grouped.setdefault(signature, []).append(orbit)
    return grouped


def non_special(grouped: dict[tuple[int, ...], list]) -> list[tuple[int, ...]]:
    """Drop the constant class and the all-singleton class, as in the first paper."""

    return sorted(s for s in grouped if len(s) >= 2 and s[0] > 1)


def contact(implementer: tuple[int, ...], implemented: tuple[int, ...]) -> str:
    beta = exchange_rate_result(implemented=implemented, implementer=implementer).beta
    if beta == 0.0:
        return "0"
    if math.isinf(beta):
        return "inf"
    return f"{beta:.4f}"


def distinct_cycles(signatures: list[tuple[int, ...]]) -> list[tuple]:
    beats = {
        (a, b): exchange_rate(a, b) < exchange_rate(b, a) - TOLERANCE
        for a, b in itertools.permutations(signatures, 2)
    }
    found: list[tuple] = []
    seen: set = set()
    for triple in itertools.combinations(signatures, 3):
        for ordered in itertools.permutations(triple):
            first, second, third = ordered
            if beats[(first, second)] and beats[(second, third)] and beats[(third, first)]:
                key = min(ordered[index:] + ordered[:index] for index in range(3))
                if key not in seen:
                    seen.add(key)
                    found.append(key)
                break
    return found


def main() -> int:
    rows = []
    for number in CASES:
        grouped = orbit_signatures(number)
        signatures = non_special(grouped)
        cycles = distinct_cycles(signatures)
        violations = [
            (a, b)
            for a, b in itertools.combinations(signatures, 2)
            if abs(complexity_index(a) - complexity_index(b)) > 1e-12
            and abs(exchange_rate(a, b) - exchange_rate(b, a)) > TOLERANCE
            and (exchange_rate(a, b) < exchange_rate(b, a))
            != (complexity_index(a) < complexity_index(b))
        ]
        print(f"case {number}: {TENSOR_CASES[number].title}")
        print(
            f"   {len(compute_tensor_poset(number).orbits)} orbits, "
            f"{len(signatures)} non-special signatures, "
            f"{len(cycles)} distinct strict 3-cycles, "
            f"{len(violations)} phi-violating pairs"
        )

        for cycle in cycles:
            for index in range(3):
                source, target = cycle[index], cycle[(index + 1) % 3]
                forward, backward = exchange_rate(source, target), exchange_rate(target, source)
                rows.append(
                    [
                        number,
                        " -> ".join("{" + ",".join(map(str, s)) + "}" for s in cycle),
                        "{" + ",".join(map(str, source)) + "}",
                        "{" + ",".join(map(str, target)) + "}",
                        f"{forward:.15f}",
                        f"{backward:.15f}",
                        contact(source, target),
                        contact(target, source),
                        f"{abs(forward - backward):.15f}",
                    ]
                )

        if cycles:
            widest = max(
                cycles,
                key=lambda c: min(
                    abs(exchange_rate(c[i], c[(i + 1) % 3]) - exchange_rate(c[(i + 1) % 3], c[i]))
                    for i in range(3)
                ),
            )
            print("   widest cycle:")
            for signature in widest:
                orbit = grouped[signature][0]
                print(
                    f"     {signature}  phi={complexity_index(signature):.6f}  "
                    f"{orbit.key}: {tensor_label(orbit, TENSOR_CASES[number])}"
                )
            for index in range(3):
                source, target = widest[index], widest[(index + 1) % 3]
                forward, backward = exchange_rate(source, target), exchange_rate(target, source)
                flag = (
                    "phi-violating"
                    if (complexity_index(source) < complexity_index(target)) is False
                    else "phi-consistent"
                )
                print(
                    f"     {'{' + ','.join(map(str, source)) + '}'} < "
                    f"{'{' + ','.join(map(str, target)) + '}'}  "
                    f"C={forward:.9f} [{contact(source, target)}] / "
                    f"{backward:.9f} [{contact(target, source)}]  "
                    f"margin={abs(forward - backward):.3e}  {flag}"
                )
        print()

    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "case",
                "cycle",
                "implementer",
                "implemented",
                "rate",
                "reverse_rate",
                "contact_beta",
                "reverse_contact_beta",
                "margin",
            ]
        )
        writer.writerows(rows)
    print(f"written to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
