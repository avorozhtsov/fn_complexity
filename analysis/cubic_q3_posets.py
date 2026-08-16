#!/usr/bin/env python3
"""Exhaustively verify the two cubic-map posets over F_3."""

from __future__ import annotations

from itertools import product
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity.cubic_field_maps import (  # noqa: E402
    AFFINE_INPUT_CLASSES,
    AFFINE_INPUT_COVERS,
    CUBIC_Q3_BASIS,
    QUADRATIC_INPUT_CLASSES,
    QUADRATIC_INPUT_COVERS,
    cubic_q3_map_count,
)


P = 3
POINTS = tuple(product(range(P), repeat=2))
POINT_INDEX = {point: index for index, point in enumerate(POINTS)}


def polynomial_values(
    coefficients: tuple[int, ...],
    exponents: tuple[tuple[int, int], ...],
) -> tuple[int, ...]:
    return tuple(
        sum(
            coefficient * pow(x, a, P) * pow(y, b, P)
            for coefficient, (a, b) in zip(coefficients, exponents)
        )
        % P
        for x, y in POINTS
    )


def scalar_polynomials(
    exponents: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        polynomial_values(coefficients, exponents)
        for coefficients in product(range(P), repeat=len(exponents))
    )


def input_point_maps(
    scalar_functions: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            POINT_INDEX[(first[index], second[index])]
            for index in range(len(POINTS))
        )
        for first in scalar_functions
        for second in scalar_functions
    )


def orbit_partition(
    functions: tuple[tuple[int, ...], ...],
    invertible_inputs: tuple[tuple[int, ...], ...],
) -> tuple[list[set[tuple[int, ...]]], dict[tuple[int, ...], int]]:
    unseen = set(functions)
    orbits: list[set[tuple[int, ...]]] = []
    orbit_index: dict[tuple[int, ...], int] = {}
    while unseen:
        representative = min(unseen)
        orbit: set[tuple[int, ...]] = set()
        for input_map in invertible_inputs:
            pulled_back = tuple(representative[index] for index in input_map)
            for scale in (1, 2):
                for shift in range(P):
                    orbit.add(
                        tuple((scale * value + shift) % P for value in pulled_back)
                    )
        index = len(orbits)
        for function in orbit:
            orbit_index[function] = index
        unseen.difference_update(orbit)
        orbits.append(orbit)
    return orbits, orbit_index


def transitive_closure(adjacency: list[set[int]]) -> list[set[int]]:
    closure: list[set[int]] = []
    for start in range(len(adjacency)):
        seen = {start}
        stack = [start]
        while stack:
            source = stack.pop()
            for target in adjacency[source]:
                if target not in seen:
                    seen.add(target)
                    stack.append(target)
        closure.append(seen)
    return closure


def strongly_connected_components(closure: list[set[int]]) -> list[set[int]]:
    unused = set(range(len(closure)))
    components: list[set[int]] = []
    while unused:
        source = min(unused)
        component = {
            target
            for target in unused
            if target in closure[source] and source in closure[target]
        }
        components.append(component)
        unused.difference_update(component)
    return components


def hasse_covers(reachability: list[set[int]]) -> set[tuple[int, int]]:
    strict = {
        (source, target)
        for source, targets in enumerate(reachability)
        for target in targets
        if source != target
    }
    return {
        (source, target)
        for source, target in strict
        if not any(
            (source, middle) in strict and (middle, target) in strict
            for middle in range(len(reachability))
            if middle not in (source, target)
        )
    }


def main() -> int:
    cubic_functions = scalar_polynomials(CUBIC_Q3_BASIS)
    assert len(cubic_functions) == cubic_q3_map_count()
    assert len(set(cubic_functions)) == cubic_q3_map_count()

    affine_basis = ((0, 0), (1, 0), (0, 1))
    affine_scalars = scalar_polynomials(affine_basis)
    affine_inputs = input_point_maps(affine_scalars)
    invertible_affine_inputs = tuple(
        input_map for input_map in affine_inputs if len(set(input_map)) == 9
    )
    assert len(affine_inputs) == 729
    assert len(invertible_affine_inputs) == 432

    orbits, orbit_index = orbit_partition(
        cubic_functions, invertible_affine_inputs
    )
    assert len(orbits) == 14

    key_to_orbit: dict[str, int] = {}
    orbit_to_key: dict[int, str] = {}
    for item in AFFINE_INPUT_CLASSES:
        representative = polynomial_values(item.coefficients, CUBIC_Q3_BASIS)
        index = orbit_index[representative]
        assert index not in orbit_to_key
        assert len(orbits[index]) == item.size
        key_to_orbit[item.key] = index
        orbit_to_key[index] = item.key
    assert len(orbit_to_key) == len(orbits)

    affine_reachability: list[set[int]] = []
    for orbit in orbits:
        representative = next(iter(orbit))
        affine_reachability.append(
            {
                orbit_index[tuple(representative[index] for index in input_map)]
                for input_map in affine_inputs
            }
        )
    affine_closure = transitive_closure(affine_reachability)
    affine_covers = {
        (orbit_to_key[source], orbit_to_key[target])
        for source, target in hasse_covers(affine_closure)
    }
    assert affine_covers == set(AFFINE_INPUT_COVERS)

    quadratic_basis = (
        (0, 0),
        (1, 0),
        (0, 1),
        (2, 0),
        (1, 1),
        (0, 2),
    )
    quadratic_scalars = scalar_polynomials(quadratic_basis)
    assert len(quadratic_scalars) == 729

    quadratic_reachability = [{key_to_orbit["constant"]} for _ in orbits]
    representatives = [next(iter(orbit)) for orbit in orbits]
    for source, representative in enumerate(representatives):
        targets = quadratic_reachability[source]
        for first in quadratic_scalars:
            for second in quadratic_scalars:
                pulled_back = tuple(
                    representative[POINT_INDEX[(first[index], second[index])]]
                    for index in range(9)
                )
                target = orbit_index.get(pulled_back)
                if target is not None:
                    targets.add(target)

    quadratic_closure = transitive_closure(quadratic_reachability)
    components = strongly_connected_components(quadratic_closure)
    orbit_to_component = {
        orbit: component_index
        for component_index, component in enumerate(components)
        for orbit in component
    }

    component_names: dict[int, str] = {}
    for component_index, component in enumerate(components):
        keys = {orbit_to_key[orbit] for orbit in component}
        if keys == {"constant"}:
            component_names[component_index] = "constant"
        elif keys == {"rank1", "cubic-63"}:
            component_names[component_index] = "two-valued"
        else:
            assert keys == set(key_to_orbit) - {
                "constant",
                "rank1",
                "cubic-63",
            }
            component_names[component_index] = "surjective"

    component_sizes = {
        component_names[index]: sum(len(orbits[orbit]) for orbit in component)
        for index, component in enumerate(components)
    }
    assert component_sizes == {
        item.key: item.size for item in QUADRATIC_INPUT_CLASSES
    }
    for component_index, component in enumerate(components):
        image_sizes = {
            len(set(function))
            for orbit in component
            for function in orbits[orbit]
        }
        expected_image_sizes = {
            "constant": {1},
            "two-valued": {2},
            "surjective": {3},
        }
        assert image_sizes == expected_image_sizes[component_names[component_index]]

    component_reachability = [set() for _ in components]
    for source, targets in enumerate(quadratic_closure):
        component_source = orbit_to_component[source]
        component_reachability[component_source].update(
            orbit_to_component[target] for target in targets
        )
    quadratic_covers = {
        (component_names[source], component_names[target])
        for source, target in hasse_covers(component_reachability)
    }
    assert quadratic_covers == set(QUADRATIC_INPUT_COVERS)

    print(f"cubic functions: {len(cubic_functions):,}")
    print(
        f"affine inputs: {len(affine_inputs):,} "
        f"({len(invertible_affine_inputs):,} invertible)"
    )
    print(f"affine classes: {len(orbits)}; Hasse covers: {len(affine_covers)}")
    print(f"quadratic inputs: {len(quadratic_scalars) ** 2:,}")
    print(f"generated-preorder classes: {len(components)}; sizes: {component_sizes}")
    print("verification: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
