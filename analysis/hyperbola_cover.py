#!/usr/bin/env python3
"""Generate the complete convergent-hyperbola catalog and cover analysis."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import (  # noqa: E402
    COVER_COLORS,
    continued_fraction_convergents,
    enumerate_hyperbolas,
    first_uncovered_greedy_cover,
    k_max,
    minimum_curve_cover,
    point_curve_map,
)

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CATALOG_PATH = OUTPUT_DIRECTORY / "hyperbolas_2-2_over_3-1.csv"
POINT_MAP_PATH = OUTPUT_DIRECTORY / "point_hyperbola_map_2-2_over_3-1.csv"
COVER_PATH = OUTPUT_DIRECTORY / "minimum_cover_2-2_over_3-1.csv"
GREEDY_PATH = OUTPUT_DIRECTORY / "first_uncovered_greedy_2-2_over_3-1.csv"


def fraction_text(value) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def formula(curve) -> str:
    return f"c_n = {fraction_text(curve.slope)} - ({fraction_text(curve.offset)})/n"


def main() -> int:
    values = [(n, k_max((2, 2), (3, 1), n)) for n in range(1, 101)]
    rate = math.log(2) / math.log(3)
    convergents = continued_fraction_convergents(rate, maximum_denominator=100)
    curves = enumerate_hyperbolas(values, convergents)
    mapping = point_curve_map(values, curves)
    minimum = minimum_curve_cover(values, curves)
    greedy = first_uncovered_greedy_cover(values, curves)
    minimum_index = {curve.identifier: index for index, curve in enumerate(minimum.curves)}

    with CATALOG_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "curve_id",
                "e_numerator",
                "e_denominator",
                "d_numerator",
                "d_denominator",
                "formula",
                "point_count",
                "n_values",
                "in_minimum_cover",
            ]
        )
        for curve in curves:
            writer.writerow(
                [
                    curve.identifier,
                    curve.slope.numerator,
                    curve.slope.denominator,
                    curve.offset.numerator,
                    curve.offset.denominator,
                    formula(curve),
                    len(curve.points),
                    " ".join(map(str, curve.points)),
                    curve.identifier in minimum_index,
                ]
            )

    with POINT_MAP_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "n",
                "k_max",
                "c_n",
                "candidate_curve_count",
                "candidate_curve_ids",
                "minimum_cover_curve_id",
                "color",
            ]
        )
        for n, k in values:
            minimum_curve = next(curve for curve in minimum.curves if n in curve.points)
            color = COVER_COLORS[minimum_index[minimum_curve.identifier] % len(COVER_COLORS)]
            writer.writerow(
                [
                    n,
                    k,
                    f"{k / n:.15f}",
                    len(mapping[n]),
                    ";".join(curve.identifier for curve in mapping[n]),
                    minimum_curve.identifier,
                    color,
                ]
            )

    with COVER_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["order", "color", "curve_id", "formula", "point_count", "n_values"])
        for index, curve in enumerate(minimum.curves):
            writer.writerow(
                [
                    index + 1,
                    COVER_COLORS[index % len(COVER_COLORS)],
                    curve.identifier,
                    formula(curve),
                    len(curve.points),
                    " ".join(map(str, curve.points)),
                ]
            )

    remaining = set(range(1, 101))
    with GREEDY_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "order",
                "curve_id",
                "formula",
                "new_point_count",
                "new_n_values",
                "all_curve_n_values",
            ]
        )
        for index, curve in enumerate(greedy):
            newly_covered = sorted(remaining.intersection(curve.points))
            remaining.difference_update(newly_covered)
            writer.writerow(
                [
                    index + 1,
                    curve.identifier,
                    formula(curve),
                    len(newly_covered),
                    " ".join(map(str, newly_covered)),
                    " ".join(map(str, curve.points)),
                ]
            )

    print(f"convergents: {', '.join(map(str, convergents))}")
    print(f"candidate curves: {len(curves)}")
    print(f"first-uncovered greedy curves: {len(greedy)}")
    print(
        f"minimum cover: {len(minimum.curves)} curves "
        f"({minimum.search_nodes} branch-and-bound nodes)"
    )
    for path in (CATALOG_PATH, POINT_MAP_PATH, COVER_PATH, GREEDY_PATH):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

