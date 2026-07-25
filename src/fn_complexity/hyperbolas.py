"""Exact hyperbola structure in finite exchange-rate sequences."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Iterable, Sequence

COVER_COLORS = (
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#3B82A0",
    "#6A3D9A",
    "#A6761D",
    "#1B9E77",
    "#E7298A",
    "#7570B3",
    "#4D7C0F",
    "#B7791F",
    "#A50F15",
    "#084594",
    "#54278F",
)


@dataclass(frozen=True)
class Hyperbola:
    """An exact curve ``c_n = slope - offset / n`` and its matching n values."""

    slope: Fraction
    offset: Fraction
    points: tuple[int, ...]

    def value(self, n: int) -> Fraction:
        return self.slope - self.offset / n

    @property
    def formula(self) -> str:
        return f"c_n = {self.slope} - ({self.offset})/n"

    @property
    def identifier(self) -> str:
        def fraction_slug(value: Fraction) -> str:
            sign = "m" if value < 0 else "p"
            magnitude = abs(value)
            return f"{sign}{magnitude.numerator}_{magnitude.denominator}"

        return f"e_{fraction_slug(self.slope)}__d_{fraction_slug(self.offset)}"


@dataclass(frozen=True)
class CoverResult:
    """A minimum-cardinality cover and search diagnostics."""

    curves: tuple[Hyperbola, ...]
    search_nodes: int


def continued_fraction_convergents(
    value: float, *, maximum_denominator: int
) -> tuple[Fraction, ...]:
    """Return positive continued-fraction convergents within a denominator bound."""

    if not math.isfinite(value) or value <= 0:
        raise ValueError("value must be finite and positive")
    if maximum_denominator < 1:
        raise ValueError("maximum_denominator must be positive")

    x = value
    p_previous_previous, p_previous = 0, 1
    q_previous_previous, q_previous = 1, 0
    result: list[Fraction] = []
    for _ in range(64):
        coefficient = math.floor(x)
        numerator = coefficient * p_previous + p_previous_previous
        denominator = coefficient * q_previous + q_previous_previous
        if denominator > maximum_denominator:
            break
        convergent = Fraction(numerator, denominator)
        if convergent > 0 and convergent not in result:
            result.append(convergent)
        remainder = x - coefficient
        if abs(remainder) < 1e-15:
            break
        p_previous_previous, p_previous = p_previous, numerator
        q_previous_previous, q_previous = q_previous, denominator
        x = 1.0 / remainder
    return tuple(result)


def enumerate_hyperbolas(
    values: Sequence[tuple[int, int]], slopes: Iterable[Fraction]
) -> tuple[Hyperbola, ...]:
    """Enumerate every exact curve induced by the supplied rational slopes."""

    if not values:
        raise ValueError("values cannot be empty")
    curves: list[Hyperbola] = []
    for slope in slopes:
        groups: dict[Fraction, list[int]] = defaultdict(list)
        for n, k in values:
            if n <= 0:
                raise ValueError("n must be positive")
            groups[slope * n - k].append(n)
        curves.extend(
            Hyperbola(slope, offset, tuple(ns))
            for offset, ns in sorted(groups.items())
        )
    return tuple(curves)


def point_curve_map(
    values: Sequence[tuple[int, int]], curves: Sequence[Hyperbola]
) -> dict[int, tuple[Hyperbola, ...]]:
    """Map every n to all enumerated curves containing its point."""

    mapping: dict[int, list[Hyperbola]] = {n: [] for n, _ in values}
    for curve in curves:
        for n in curve.points:
            mapping[n].append(curve)
    return {n: tuple(matches) for n, matches in mapping.items()}


def maximum_coverage_curve_per_slope(
    curves: Sequence[Hyperbola],
) -> tuple[Hyperbola, ...]:
    """Select one longest curve per slope, preferring the smallest absolute d."""

    by_slope: dict[Fraction, list[Hyperbola]] = defaultdict(list)
    for curve in curves:
        by_slope[curve.slope].append(curve)
    return tuple(
        max(
            slope_curves,
            key=lambda curve: (
                len(curve.points),
                -abs(curve.offset),
                -curve.offset,
            ),
        )
        for slope_curves in by_slope.values()
    )


def first_uncovered_greedy_cover(
    values: Sequence[tuple[int, int]], curves: Sequence[Hyperbola]
) -> tuple[Hyperbola, ...]:
    """Coverage-first version of the proposed first-uncovered-point algorithm."""

    remaining = {n for n, _ in values}
    selected: list[Hyperbola] = []
    while remaining:
        first = min(remaining)
        choices = [curve for curve in curves if first in curve.points]
        if not choices:
            raise ValueError(f"no curve covers n={first}")
        chosen = max(
            choices,
            key=lambda curve: (
                len(remaining.intersection(curve.points)),
                len(curve.points),
                -abs(curve.offset),
                -curve.slope.denominator,
            ),
        )
        selected.append(chosen)
        remaining.difference_update(chosen.points)
    return tuple(selected)


def minimum_curve_cover(
    values: Sequence[tuple[int, int]], curves: Sequence[Hyperbola]
) -> CoverResult:
    """Find and certify a minimum-cardinality cover by branch and bound."""

    point_order = tuple(n for n, _ in values)
    point_index = {n: index for index, n in enumerate(point_order)}
    all_points_mask = (1 << len(point_order)) - 1
    masks = [
        sum(1 << point_index[n] for n in curve.points)
        for curve in curves
    ]
    covering_curves: list[list[int]] = [[] for _ in point_order]
    for curve_index, curve in enumerate(curves):
        for n in curve.points:
            covering_curves[point_index[n]].append(curve_index)

    # Every fixed slope partitions the data.  The best such partition is a
    # strong deterministic initial upper bound.
    by_slope: dict[Fraction, list[int]] = defaultdict(list)
    for index, curve in enumerate(curves):
        by_slope[curve.slope].append(index)
    best = min(
        by_slope.values(),
        key=lambda indices: (len(indices), curves[indices[0]].slope.denominator),
    )
    best = list(best)
    memo: dict[int, int] = {}
    search_nodes = 0

    def search(covered: int, selected: list[int]) -> None:
        nonlocal best, search_nodes
        search_nodes += 1
        if covered == all_points_mask:
            if len(selected) < len(best):
                best = selected.copy()
            return
        if len(selected) >= len(best) - 1:
            return

        remaining = all_points_mask ^ covered
        maximum_gain = max((mask & remaining).bit_count() for mask in masks)
        optimistic_more = (remaining.bit_count() + maximum_gain - 1) // maximum_gain
        if len(selected) + optimistic_more >= len(best):
            return
        previous_depth = memo.get(covered)
        if previous_depth is not None and previous_depth <= len(selected):
            return
        memo[covered] = len(selected)

        uncovered_indices = [
            index
            for index in range(len(point_order))
            if remaining & (1 << index)
        ]
        pivot = min(
            uncovered_indices,
            key=lambda index: (
                len(covering_curves[index]),
                -max(
                    (masks[curve_index] & remaining).bit_count()
                    for curve_index in covering_curves[index]
                ),
            ),
        )
        options = sorted(
            covering_curves[pivot],
            key=lambda index: (masks[index] & remaining).bit_count(),
            reverse=True,
        )
        for curve_index in options:
            search(covered | masks[curve_index], selected + [curve_index])

    search(0, [])
    selected_curves = sorted(
        (curves[index] for index in best),
        key=lambda curve: (
            min(curve.points),
            curve.slope.denominator,
            curve.offset,
        ),
    )
    return CoverResult(tuple(selected_curves), search_nodes)
