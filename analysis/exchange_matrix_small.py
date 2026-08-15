#!/usr/bin/env python3
"""Generate the small exchange-rate matrix requested in the investigation."""

from __future__ import annotations

import csv
from itertools import combinations_with_replacement
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate_result  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIRECTORY / "exchange_matrix_small.csv"
MARKDOWN_PATH = OUTPUT_DIRECTORY / "exchange_matrix_small.md"


def label(signature: tuple[int, ...]) -> str:
    return "{" + ",".join(map(str, signature)) + "}"


def signatures() -> list[tuple[int, ...]]:
    length_three = [
        tuple(reversed(values))
        for values in combinations_with_replacement(range(1, 4), 3)
    ]
    # Positive n_i < 2 leaves only n_i=1.
    return length_three + [(1, 1, 1, 1)]


def main() -> int:
    all_signatures = signatures()
    rates = {
        (g, f): exchange_rate_result(
            g, f, grid_size=2048, tolerance=1e-13
        ).rate
        for g in all_signatures
        for f in all_signatures
    }

    # Orient an edge b -> a when a is more complex than b:
    # C(a -> b) > C(b -> a).  Kahn's algorithm then lists low complexity first.
    tolerance = 1e-12
    outgoing = {signature: set() for signature in all_signatures}
    indegree = {signature: 0 for signature in all_signatures}
    for index, a in enumerate(all_signatures):
        for b in all_signatures[index + 1 :]:
            if rates[a, b] > rates[b, a] + tolerance:
                outgoing[b].add(a)
                indegree[a] += 1
            elif rates[b, a] > rates[a, b] + tolerance:
                outgoing[a].add(b)
                indegree[b] += 1

    remaining = set(all_signatures)
    layers: list[tuple[tuple[int, ...], ...]] = []
    while remaining:
        layer = tuple(
            sorted(
                (signature for signature in remaining if indegree[signature] == 0),
                key=lambda signature: (len(signature), signature),
            )
        )
        if not layer:
            raise RuntimeError(
                "the pairwise complexity relation contains a directed cycle"
            )
        layers.append(layer)
        for signature in layer:
            remaining.remove(signature)
            for successor in outgoing[signature]:
                indegree[successor] -= 1
    ordered = [signature for layer in layers for signature in layer]
    rank = {
        signature: index
        for index, signature in enumerate(ordered, 1)
    }
    layer_number = {
        signature: index
        for index, layer in enumerate(layers, 1)
        for signature in layer
    }

    with CSV_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "implementer g / implemented f",
                "topological_rank",
                "topological_layer",
            ]
            + [label(signature) for signature in ordered]
        )
        for g in ordered:
            writer.writerow(
                [
                    label(g),
                    rank[g],
                    layer_number[g],
                ]
                + [f"{rates[g, f]:.15f}" for f in ordered]
            )

    with MARKDOWN_PATH.open("w", encoding="utf-8") as stream:
        stream.write("# Small exchange-rate matrix\n\n")
        stream.write(
            "Rows are implementers `g`, columns are implemented signatures "
            "`f`, and each cell is `C(g -> f)`.\n\n"
        )
        stream.write(
            "The order is a topological sort of the relation `b -> a` when "
            "`C(a -> b) > C(b -> a)`. Lower-complexity signatures are shown "
            "first.\n\n"
        )
        headers = ["g \\ f"] + [label(signature) for signature in ordered]
        stream.write("| " + " | ".join(headers) + " |\n")
        stream.write("| " + " | ".join(["---"] * len(headers)) + " |\n")
        for g in ordered:
            cells = [label(g)] + [
                f"{rates[g, f]:.6f}" for f in ordered
            ]
            stream.write("| " + " | ".join(cells) + " |\n")

        stream.write("\n## Topological order\n\n")
        for index, layer in enumerate(layers, 1):
            stream.write(
                f"{index}. "
                + ", ".join(f"`{label(signature)}`" for signature in layer)
                + "\n"
            )
        edge_count = sum(len(successors) for successors in outgoing.values())
        stream.write(
            f"\nThe graph has `{edge_count}` strict comparison edges. "
            f"There are `{len(layers)}` topological layers.\n"
        )

    print(f"signatures: {len(ordered)}")
    print(f"strict comparison edges: {sum(len(x) for x in outgoing.values())}")
    print(f"topological layers: {len(layers)}")
    print("order:", ", ".join(label(signature) for signature in ordered))
    print(CSV_PATH)
    print(MARKDOWN_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
