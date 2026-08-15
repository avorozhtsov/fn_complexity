#!/usr/bin/env python3
"""Generate the 97-signature cached exchange matrix and cycle analysis."""

from __future__ import annotations

import csv
from itertools import combinations_with_replacement
import math
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import ExchangeRateCache  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "exchange_matrix_extended.csv"
MARKDOWN_PATH = OUTPUT_DIRECTORY / "exchange_matrix_extended.md"
SPECIFICATION = ((1, 6), (2, 6), (3, 5), (4, 4))


def label(signature: tuple[int, ...]) -> str:
    return "{" + ",".join(map(str, signature)) + "}"


def signatures() -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []
    for length, maximum in SPECIFICATION:
        result.extend(
            tuple(reversed(values))
            for values in combinations_with_replacement(
                range(1, maximum + 1), length
            )
        )
    return result


def formatted_rate(value: float, digits: int) -> str:
    return "∞" if math.isinf(value) else f"{value:.{digits}f}"


def comparison_graph(
    all_signatures: list[tuple[int, ...]],
    rates: dict[tuple[tuple[int, ...], tuple[int, ...]], float],
    *,
    tolerance: float = 1e-10,
) -> tuple[
    dict[tuple[int, ...], set[tuple[int, ...]]],
    int,
]:
    """Return edges b -> a iff C(a -> b) > C(b -> a), plus tie count."""

    outgoing = {signature: set() for signature in all_signatures}
    tie_count = 0
    for index, a in enumerate(all_signatures):
        for b in all_signatures[index + 1 :]:
            a_over_b = rates[a, b]
            b_over_a = rates[b, a]
            if a_over_b > b_over_a + tolerance:
                outgoing[b].add(a)
            elif b_over_a > a_over_b + tolerance:
                outgoing[a].add(b)
            else:
                tie_count += 1
    return outgoing, tie_count


def strongly_connected_components(
    all_signatures: list[tuple[int, ...]],
    outgoing: dict[tuple[int, ...], set[tuple[int, ...]]],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Return deterministic strongly connected components (Tarjan)."""

    index = 0
    indices: dict[tuple[int, ...], int] = {}
    lowlinks: dict[tuple[int, ...], int] = {}
    stack: list[tuple[int, ...]] = []
    on_stack: set[tuple[int, ...]] = set()
    components: list[tuple[tuple[int, ...], ...]] = []

    def visit(node: tuple[int, ...]) -> None:
        nonlocal index
        indices[node] = index
        lowlinks[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)

        for successor in sorted(
            outgoing[node], key=lambda signature: (len(signature), signature)
        ):
            if successor not in indices:
                visit(successor)
                lowlinks[node] = min(lowlinks[node], lowlinks[successor])
            elif successor in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[successor])

        if lowlinks[node] == indices[node]:
            component: list[tuple[int, ...]] = []
            while True:
                member = stack.pop()
                on_stack.remove(member)
                component.append(member)
                if member == node:
                    break
            components.append(
                tuple(sorted(component, key=lambda value: (len(value), value)))
            )

    for signature in sorted(
        all_signatures, key=lambda value: (len(value), value)
    ):
        if signature not in indices:
            visit(signature)
    return tuple(
        sorted(
            components,
            key=lambda component: (
                len(component[0]),
                component[0],
                len(component),
            ),
        )
    )


def condensation_layers(
    components: tuple[tuple[tuple[int, ...], ...], ...],
    outgoing: dict[tuple[int, ...], set[tuple[int, ...]]],
) -> tuple[tuple[int, ...], ...]:
    """Topologically layer the DAG obtained by contracting each SCC."""

    component_of = {
        signature: component_index
        for component_index, component in enumerate(components)
        for signature in component
    }
    component_outgoing = {index: set() for index in range(len(components))}
    indegree = {index: 0 for index in range(len(components))}
    for source, successors in outgoing.items():
        source_component = component_of[source]
        for successor in successors:
            target_component = component_of[successor]
            if (
                source_component != target_component
                and target_component not in component_outgoing[source_component]
            ):
                component_outgoing[source_component].add(target_component)
                indegree[target_component] += 1

    remaining = set(range(len(components)))
    layers: list[tuple[int, ...]] = []
    while remaining:
        layer = tuple(
            sorted(
                (
                    component_index
                    for component_index in remaining
                    if indegree[component_index] == 0
                ),
                key=lambda component_index: (
                    len(components[component_index][0]),
                    components[component_index][0],
                ),
            )
        )
        if not layer:
            raise AssertionError("the SCC condensation graph must be acyclic")
        layers.append(layer)
        for component_index in layer:
            remaining.remove(component_index)
            for successor in component_outgoing[component_index]:
                indegree[successor] -= 1
    return tuple(layers)


def directed_triangle(
    component: tuple[tuple[int, ...], ...],
    outgoing: dict[tuple[int, ...], set[tuple[int, ...]]],
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]] | None:
    """Return a representative directed 3-cycle from an SCC, if present."""

    members = set(component)
    for first in component:
        for second in sorted(
            outgoing[first] & members,
            key=lambda signature: (len(signature), signature),
        ):
            for third in sorted(
                outgoing[second] & members,
                key=lambda signature: (len(signature), signature),
            ):
                if first in outgoing[third]:
                    return first, second, third
    return None


def main() -> int:
    all_signatures = signatures()
    cache = ExchangeRateCache()
    seeded = cache.seed_from_matrix_csv(CSV_PATH)
    rates = {
        (g, f): cache.get(g, f)
        for g in all_signatures
        for f in all_signatures
    }
    cache.save()
    outgoing, tie_count = comparison_graph(all_signatures, rates)
    components = strongly_connected_components(all_signatures, outgoing)
    component_layers = condensation_layers(components, outgoing)
    cyclic_components = tuple(
        component for component in components if len(component) > 1
    )
    ordered = [
        signature
        for layer in component_layers
        for component_index in layer
        for signature in components[component_index]
    ]
    rank = {
        signature: index
        for index, signature in enumerate(ordered, 1)
    }
    layer_number = {
        signature: index
        for index, layer in enumerate(component_layers, 1)
        for component_index in layer
        for signature in components[component_index]
    }
    component_number = {
        signature: index
        for index, component in enumerate(components, 1)
        for signature in component
    }
    component_size = {
        signature: len(component)
        for component in components
        for signature in component
    }
    edge_count = sum(len(successors) for successors in outgoing.values())

    with CSV_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "implementer g / implemented f",
                "topological_rank",
                "topological_layer",
                "scc",
                "scc_size",
            ]
            + [label(signature) for signature in ordered]
        )
        for g in ordered:
            writer.writerow(
                [
                    label(g),
                    rank[g],
                    layer_number[g],
                    component_number[g],
                    component_size[g],
                ]
                + [
                    "inf" if math.isinf(rates[g, f])
                    else f"{rates[g, f]:.15f}"
                    for f in ordered
                ]
            )

    with MARKDOWN_PATH.open("w", encoding="utf-8") as stream:
        stream.write("# Extended exchange-rate matrix\n\n")
        stream.write(
            "Rows are implementers `g`, columns are implemented signatures "
            "`f`, and each cell is `C(g -> f)`.\n\n"
        )
        stream.write(
            f"There are **{len(ordered)} signatures**, **{edge_count} strict "
            f"comparison edges**, **{tie_count} tied pairs**, and "
            f"**{len(component_layers)} condensation-DAG layers**. "
            + (
                f"The strict graph contains **{len(cyclic_components)} "
                "cyclic strongly connected component(s)**.\n\n"
                if cyclic_components
                else "The strict graph is acyclic.\n\n"
            )
        )
        if cyclic_components:
            stream.write("## Cyclic strongly connected components\n\n")
            for index, component in enumerate(cyclic_components, 1):
                stream.write(
                    f"{index}. "
                    + ", ".join(
                        f"`{label(signature)}`" for signature in component
                    )
                    + "\n"
                )
                triangle = directed_triangle(component, outgoing)
                if triangle is not None:
                    closed_triangle = triangle + (triangle[0],)
                    stream.write(
                        "   Representative strict 3-cycle: "
                        + " → ".join(
                            f"`{label(signature)}`"
                            for signature in closed_triangle
                        )
                        + ".\n"
                    )
            stream.write("\n")

        stream.write("## Condensation-DAG layers, low to high\n\n")
        for index, layer in enumerate(component_layers, 1):
            layer_parts = []
            for component_index in layer:
                component = components[component_index]
                rendered = ", ".join(
                    f"`{label(signature)}`" for signature in component
                )
                layer_parts.append(
                    f"({rendered})" if len(component) > 1 else rendered
                )
            stream.write(
                f"{index}. "
                + "; ".join(layer_parts)
                + "\n"
            )

        stream.write("\n## Matrix\n\n")
        headers = ["g \\ f"] + [label(signature) for signature in ordered]
        stream.write("| " + " | ".join(headers) + " |\n")
        stream.write("| " + " | ".join(["---"] * len(headers)) + " |\n")
        for g in ordered:
            cells = [label(g)] + [
                formatted_rate(rates[g, f], 6) for f in ordered
            ]
            stream.write("| " + " | ".join(cells) + " |\n")

    print(f"signatures: {len(ordered)}")
    print(f"strict comparison edges: {edge_count}")
    print(f"tied pairs: {tie_count}")
    print(f"condensation-DAG layers: {len(component_layers)}")
    print(f"cyclic SCCs: {len(cyclic_components)}")
    for component in cyclic_components:
        print("  cycle SCC:", ", ".join(label(value) for value in component))
    print(
        f"cache: {cache.hits} hits, {cache.misses} misses, "
        f"{seeded} values seeded, {len(cache)} total"
    )
    print(CSV_PATH)
    print(MARKDOWN_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
