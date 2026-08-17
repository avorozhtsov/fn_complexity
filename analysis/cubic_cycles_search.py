#!/usr/bin/env python3
"""Cycle search in two families too large to classify: cubic maps over F_8 and
cubic homogeneous maps F_3^3 -> F_3^3.

Neither family has a computed orbit list --- the cubic ``P^2`` case has at least
1,632,040 classes --- but the exchange comparison only sees fiber signatures,
which are far coarser than orbits.  This script therefore samples the map spaces
directly, keeps one witness map per distinct signature, and searches the strict
comparison graph for cycles.

Sampling mixes uniform coefficient vectors with sparse ones, because uniform
sampling alone almost never reaches the structured maps that sit in small
orbits.

Finding a cycle is a proof: every signature analysed is realised by the witness
map printed alongside it.  Finding none is evidence, not a theorem.
"""

from __future__ import annotations

from collections import Counter
import csv
import itertools
import math
from pathlib import Path
import random
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate, exchange_rate_result  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "cubic_cycles_search.csv"
SAMPLES = 150_000
SIGNATURE_CAP = 140
TOLERANCE = 1e-9
SEED = 17

# F_8 = F_2[a]/(a^3 + a + 1), as in the cubic F_8 note.
GF8_POLY = 0b1011
BIVARIATE_EXPONENTS = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0), (2, 1), (1, 2), (0, 3)]
MONOMIAL_NAMES = ["1", "x", "y", "x^2", "xy", "y^2", "x^3", "x^2y", "xy^2", "y^3"]


def gf8_multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 0b1000:
            left ^= GF8_POLY
    return result


GF8 = np.array([[gf8_multiply(a, b) for b in range(8)] for a in range(8)], dtype=np.int64)


def coefficient_vector(rng: random.Random, length: int, modulus: int) -> list[int]:
    """Uniform half the time, sparse the other half."""

    if rng.random() < 0.5:
        return [rng.randrange(modulus) for _ in range(length)]
    support = rng.sample(range(length), rng.randint(1, 4))
    vector = [0] * length
    for position in support:
        vector[position] = rng.randrange(1, modulus)
    return vector


def gf8_cubic_witnesses(samples: int) -> dict[tuple[int, ...], str]:
    rng = random.Random(SEED)
    points = [(x, y) for x in range(8) for y in range(8)]

    def power(value: int, exponent: int) -> int:
        result = 1
        for _ in range(exponent):
            result = GF8[result, value]
        return result

    table = np.array(
        [[GF8[power(x, i), power(y, j)] for (i, j) in BIVARIATE_EXPONENTS] for (x, y) in points],
        dtype=np.int64,
    )
    witnesses: dict[tuple[int, ...], str] = {}
    for _ in range(samples):
        coefficients = coefficient_vector(rng, 10, 8)
        values = np.zeros(len(points), dtype=np.int64)
        for index, coefficient in enumerate(coefficients):
            if coefficient:
                values ^= GF8[coefficient, table[:, index]]
        counts = np.bincount(values, minlength=8)
        signature = tuple(sorted((int(v) for v in counts if v), reverse=True))
        witnesses.setdefault(
            signature,
            " + ".join(
                f"{c}·{MONOMIAL_NAMES[i]}" for i, c in enumerate(coefficients) if c
            )
            or "0",
        )
    return witnesses


def f3_cubic_p2_witnesses(samples: int) -> dict[tuple[int, ...], str]:
    rng = random.Random(SEED + 1)
    points = list(itertools.product(range(3), repeat=3))
    exponents = [e for e in itertools.product(range(4), repeat=3) if sum(e) == 3]
    names = ["".join(v * e for v, e in zip("xyz", exp)) for exp in exponents]
    table = np.array(
        [[(p[0] ** e[0] * p[1] ** e[1] * p[2] ** e[2]) % 3 for e in exponents] for p in points],
        dtype=np.int64,
    )
    witnesses: dict[tuple[int, ...], str] = {}
    for _ in range(samples):
        rows = [coefficient_vector(rng, 10, 3) for _ in range(3)]
        coefficients = np.array(rows, dtype=np.int64)
        values = (table @ coefficients.T) % 3
        keys = values[:, 0] * 9 + values[:, 1] * 3 + values[:, 2]
        counts = np.bincount(keys, minlength=27)
        signature = tuple(sorted((int(v) for v in counts if v), reverse=True))
        label = "⟨" + "; ".join(
            " + ".join(f"{c}{names[i]}" if c > 1 else names[i]
                       for i, c in enumerate(row) if c) or "0"
            for row in rows
        ) + "⟩"
        witnesses.setdefault(signature, label)
    return witnesses


def complexity_index(signature: tuple[int, ...]) -> float:
    return math.log(len(signature)) * math.log(max(signature))


def contact(implementer: tuple[int, ...], implemented: tuple[int, ...]) -> str:
    beta = exchange_rate_result(implemented=implemented, implementer=implementer).beta
    if beta == 0.0:
        return "0"
    if math.isinf(beta):
        return "inf"
    return f"{beta:.4f}"


def analyse(name: str, witnesses: dict[tuple[int, ...], str], rows: list) -> None:
    candidates = [s for s in witnesses if len(s) >= 2 and s[0] > 1]
    signatures = sorted(candidates[:SIGNATURE_CAP])
    rate = {
        (a, b): exchange_rate(a, b) for a, b in itertools.permutations(signatures, 2)
    }
    violations = [
        (a, b)
        for a, b in itertools.combinations(signatures, 2)
        if abs(rate[(a, b)] - rate[(b, a)]) > TOLERANCE
        and abs(complexity_index(a) - complexity_index(b)) > 1e-12
        and (rate[(a, b)] < rate[(b, a)]) != (complexity_index(a) < complexity_index(b))
    ]
    cycles = []
    for triple in itertools.combinations(signatures, 3):
        for ordered in itertools.permutations(triple):
            first, second, third = ordered
            if (
                rate[(first, second)] < rate[(second, first)] - TOLERANCE
                and rate[(second, third)] < rate[(third, second)] - TOLERANCE
                and rate[(third, first)] < rate[(first, third)] - TOLERANCE
            ):
                cycles.append(ordered)
                break
    print(f"\n{name}")
    print(
        f"  {len(witnesses)} distinct signatures sampled, {len(signatures)} non-special "
        f"analysed, {len(violations)} phi-violating pairs, {len(cycles)} strict 3-cycles"
    )
    if not cycles:
        print("  no cycle found (evidence, not a theorem)")
        return

    def margin(cycle):
        return min(
            abs(rate[(cycle[i], cycle[(i + 1) % 3])] - rate[(cycle[(i + 1) % 3], cycle[i])])
            for i in range(3)
        )

    widest = max(cycles, key=margin)
    print(f"  widest cycle, minimum margin {margin(widest):.3e}:")
    for signature in widest:
        print(f"    {signature}")
        print(f"      phi={complexity_index(signature):.6f}   witness {witnesses[signature]}")
    for index in range(3):
        source, target = widest[index], widest[(index + 1) % 3]
        forward, backward = rate[(source, target)], rate[(target, source)]
        consistent = (complexity_index(source) < complexity_index(target))
        print(
            f"    edge {index + 1}: C={forward:.9f} [{contact(source, target)}] / "
            f"{backward:.9f} [{contact(target, source)}]  margin={abs(forward - backward):.3e}"
            f"  {'phi-consistent' if consistent else 'phi-VIOLATING'}"
        )
        rows.append(
            [
                name,
                "{" + ",".join(map(str, source)) + "}",
                "{" + ",".join(map(str, target)) + "}",
                f"{forward:.15f}",
                f"{backward:.15f}",
                contact(source, target),
                contact(target, source),
                f"{abs(forward - backward):.15f}",
                "consistent" if consistent else "violating",
                witnesses[source],
            ]
        )
    single = sum(
        1
        for cycle in cycles
        if sum(
            1
            for i in range(3)
            if (rate[(cycle[i], cycle[(i + 1) % 3])] < rate[(cycle[(i + 1) % 3], cycle[i])])
            != (complexity_index(cycle[i]) < complexity_index(cycle[(i + 1) % 3]))
        )
        == 1
    )
    print(f"  cycles closed by exactly one phi-violating edge: {single} of {len(cycles)}")


def main() -> int:
    rows: list = []
    analyse("cubic maps over F_8", gf8_cubic_witnesses(SAMPLES), rows)
    analyse("cubic homogeneous maps F_3^3 -> F_3^3", f3_cubic_p2_witnesses(SAMPLES), rows)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "family",
                "implementer",
                "implemented",
                "rate",
                "reverse_rate",
                "contact_beta",
                "reverse_contact_beta",
                "margin",
                "phi",
                "witness",
            ]
        )
        writer.writerows(rows)
    print(f"\nwritten to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
