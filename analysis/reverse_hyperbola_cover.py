#!/usr/bin/env python3
"""Analyze C((3,1) | (2,2)) through n < 100000.

For this pair, (2,2)^k has 2^k fibers of size 2^k.  The fibers of
(3,1)^n have sizes 3^j with multiplicity binomial(n, j), so feasibility is

    sum(binomial(n, j), j >= ceil(k log(2)/log(3))) >= 2^k.

The logarithmic tail evaluation is certified with an exact integer fallback
whenever the floating-point comparison is close.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from fractions import Fraction
import math
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import (  # noqa: E402
    continued_fraction_convergents,
    exchange_rate_result,
    k_max,
)

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
LN2 = math.log(2)
LOG2_OVER_LOG3 = LN2 / math.log(3)


def threshold(k: int) -> int:
    """Return the least j with 3**j >= 2**k."""

    estimate = k * LOG2_OVER_LOG3
    nearest = round(estimate)
    if abs(estimate - nearest) > 1e-10:
        return math.ceil(estimate)
    # This path is rare.  Exact integer comparison avoids a bad ceil near an
    # exceptionally close rational approximation to log(2)/log(3).
    candidate = math.floor(estimate)
    return candidate if pow(3, candidate) >= pow(2, k) else candidate + 1


def log_binomial_tail(n: int, first_j: int) -> float:
    """Return log(sum(C(n,j), j >= first_j)) with a stable relative sum."""

    if first_j > n:
        return -math.inf
    if first_j <= 0:
        return n * LN2
    log_first = (
        math.lgamma(n + 1)
        - math.lgamma(first_j + 1)
        - math.lgamma(n - first_j + 1)
    )
    relative_sum = 1.0
    relative_term = 1.0
    for j in range(first_j, n):
        relative_term *= (n - j) / (j + 1)
        relative_sum += relative_term
        if relative_term < relative_sum * 1e-16:
            break
    return log_first + math.log(relative_sum)


def feasible(n: int, k: int) -> bool:
    """Test whether (2,2)^k is implemented by (3,1)^n."""

    first_j = threshold(k)
    margin = log_binomial_tail(n, first_j) - k * LN2
    if abs(margin) > 1e-8:
        return margin > 0
    exact_tail = sum(math.comb(n, j) for j in range(first_j, n + 1))
    return exact_tail >= 1 << k


def finite_values(n_max: int) -> list[tuple[int, int]]:
    """Compute all (n, k_max(n)) sequentially."""

    values: list[tuple[int, int]] = []
    k = 0
    for n in range(1, n_max + 1):
        while feasible(n, k + 1):
            k += 1
        while not feasible(n, k):
            k -= 1
        values.append((n, k))
    return values


def formula(e: Fraction, d: Fraction) -> str:
    """Format the standard c_n=e-d/n convention."""

    if d < 0:
        return f"c_n = {e} + {-d}/n"
    return f"c_n = {e} - {d}/n"


def write_best_convergents(
    values: list[tuple[int, int]], convergents: tuple[Fraction, ...], path: Path
) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "order",
                "e",
                "d",
                "formula",
                "point_count",
                "first_n",
                "last_n",
                "n_values",
            ]
        )
        for order, e in enumerate(convergents, 1):
            groups: dict[Fraction, list[int]] = defaultdict(list)
            for n, k in values:
                groups[e * n - k].append(n)
            d, ns = max(
                groups.items(),
                key=lambda item: (
                    len(item[1]),
                    -abs(item[0]),
                    -item[0],
                ),
            )
            writer.writerow(
                [
                    order,
                    str(e),
                    str(d),
                    formula(e, d),
                    len(ns),
                    ns[0],
                    ns[-1],
                    " ".join(map(str, ns)),
                ]
            )


def write_e1_partition(values: list[tuple[int, int]], path: Path) -> None:
    groups: dict[int, list[int]] = defaultdict(list)
    for n, k in values:
        groups[n - k].append(n)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "m",
                "standard_formula",
                "plus_formula_d",
                "point_count",
                "first_n",
                "last_n",
                "n_values",
            ]
        )
        for m in sorted(groups):
            ns = groups[m]
            writer.writerow(
                [
                    m,
                    f"c_n = 1 - {m}/n",
                    -m,
                    len(ns),
                    ns[0],
                    ns[-1],
                    " ".join(map(str, ns)),
                ]
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n-max",
        type=int,
        default=99999,
        help="last included n (99999 means n < 100000)",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.n_max < 1:
        raise SystemExit("--n-max must be positive")

    values = finite_values(args.n_max)
    for n, computed in values[: min(250, args.n_max)]:
        exact = k_max((3, 1), (2, 2), n)
        if computed != exact:
            raise RuntimeError(
                f"binomial-tail computation disagrees at n={n}: "
                f"{computed} != {exact}"
            )

    rate = exchange_rate_result((3, 1), (2, 2)).rate
    convergents = continued_fraction_convergents(
        rate, maximum_denominator=1_000_000
    )[:8]
    suffix = f"n-1-{args.n_max}"
    best_path = (
        OUTPUT_DIRECTORY
        / f"best_convergent_hyperbolas_3-1_over_2-2_{suffix}.csv"
    )
    partition_path = (
        OUTPUT_DIRECTORY / f"e1_partition_3-1_over_2-2_{suffix}.csv"
    )
    write_best_convergents(values, convergents, best_path)
    write_e1_partition(values, partition_path)

    print(f"C((3,1) | (2,2)) = {rate:.15f}")
    print(f"convergents: {', '.join(map(str, convergents))}")
    print(f"k_max({args.n_max}) = {values[-1][1]}")
    print(best_path)
    print(partition_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
