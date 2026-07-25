#!/usr/bin/env python3
"""Verify stabilization of Appendix B's first 99 signatures.

This is the large screening calculation behind
``analysis/appendix_b_signatures.py``.  It requires NumPy and SciPy and uses a
vectorized power-sum grid before contracting the strict-comparison graph.
The final displayed rates are computed separately by the project's
high-accuracy persistent cache.
"""

from __future__ import annotations

import argparse
from itertools import combinations_with_replacement
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "analysis"))

from appendix_b_signatures import STABLE_FIRST_99  # noqa: E402


def candidate_signatures(budget: int) -> list[tuple[int, ...]]:
    """Return the non-special signatures with ``a_1 + 2 len(a) <= budget``."""

    result: list[tuple[int, ...]] = []
    for length in range(2, (budget - 1) // 2 + 1):
        maximum = budget - 2 * length
        result.extend(
            tuple(reversed(values))
            for values in combinations_with_replacement(
                range(1, maximum + 1), length
            )
            if values[-1] > 1
        )
    return sorted(result, key=lambda signature: (len(signature), signature))


def approximate_rate_matrix(
    signatures: list[tuple[int, ...]],
    *,
    budget: int,
    grid_size: int,
) -> np.ndarray:
    """Evaluate the power-sum infimum on a shared geometric beta grid."""

    beta = np.concatenate(
        (
            np.array([0.0]),
            np.geomspace(1.0e-10, 36.0 * budget, grid_size - 1),
        )
    )
    log_partitions = np.empty((len(signatures), grid_size))
    for index, signature in enumerate(signatures):
        weighted_logs = beta[:, None] * np.log(
            np.asarray(signature, dtype=float)
        )
        largest = weighted_logs.max(axis=1)
        log_partitions[index] = largest + np.log(
            np.exp(weighted_logs - largest[:, None]).sum(axis=1)
        )

    rates = np.full((len(signatures), len(signatures)), np.inf)
    for beta_index in range(grid_size):
        np.minimum(
            rates,
            (
                log_partitions[:, beta_index, None]
                / log_partitions[None, :, beta_index]
            ),
            out=rates,
        )

    log_maxima = np.log(
        np.asarray([signature[0] for signature in signatures], dtype=float)
    )
    np.minimum(
        rates,
        log_maxima[:, None] / log_maxima[None, :],
        out=rates,
    )
    return rates


def condensation_order(
    signatures: list[tuple[int, ...]],
    rates: np.ndarray,
    *,
    tolerance: float,
) -> list[tuple[int, ...]]:
    """Return the deterministic SCC-condensation order."""

    # Edge i -> j means signature j wins the directed comparison with i.
    adjacency = csr_matrix(rates.T > rates + tolerance)
    component_count, old_labels = connected_components(
        adjacency,
        directed=True,
        connection="strong",
    )

    old_components: list[list[tuple[int, ...]]] = [
        [] for _ in range(component_count)
    ]
    for signature, label in zip(signatures, old_labels):
        old_components[label].append(signature)
    old_component_tuples = [tuple(component) for component in old_components]
    old_indices = sorted(
        range(component_count),
        key=lambda index: (
            len(old_component_tuples[index][0]),
            old_component_tuples[index][0],
            len(old_component_tuples[index]),
        ),
    )
    components = [old_component_tuples[index] for index in old_indices]
    new_index = np.empty(component_count, dtype=int)
    for new, old in enumerate(old_indices):
        new_index[old] = new
    labels = new_index[old_labels]

    membership = csr_matrix(
        (
            np.ones(len(signatures), dtype=np.int32),
            (np.arange(len(signatures)), labels),
        ),
        shape=(len(signatures), component_count),
    )
    component_adjacency = (
        membership.T @ adjacency.astype(np.int32) @ membership
    ).astype(bool).tocsr()
    component_adjacency.setdiag(False)
    component_adjacency.eliminate_zeros()

    outgoing = [
        set(
            component_adjacency.indices[
                component_adjacency.indptr[index]
                : component_adjacency.indptr[index + 1]
            ]
        )
        for index in range(component_count)
    ]
    indegree = (
        np.asarray(component_adjacency.sum(axis=0))
        .ravel()
        .astype(int)
        .tolist()
    )
    remaining = set(range(component_count))
    ordered: list[tuple[int, ...]] = []
    while remaining:
        layer = sorted(
            (
                index
                for index in remaining
                if indegree[index] == 0
            ),
            key=lambda index: (
                len(components[index][0]),
                components[index][0],
            ),
        )
        if not layer:
            raise AssertionError("the SCC condensation graph must be acyclic")
        for index in layer:
            ordered.extend(components[index])
        for index in layer:
            remaining.remove(index)
            for successor in outgoing[index]:
                indegree[successor] -= 1
    return ordered


def first_signatures(
    budget: int,
    *,
    count: int,
    grid_size: int,
    tolerance: float,
) -> list[tuple[int, ...]]:
    signatures = candidate_signatures(budget)
    rates = approximate_rate_matrix(
        signatures,
        budget=budget,
        grid_size=grid_size,
    )
    return condensation_order(
        signatures,
        rates,
        tolerance=tolerance,
    )[:count]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--budgets", nargs="+", type=int, default=[18, 19])
    parser.add_argument("--count", type=int, default=99)
    parser.add_argument("--grid-size", type=int, default=256)
    parser.add_argument("--tolerance", type=float, default=2.0e-6)
    args = parser.parse_args()

    orders = []
    for budget in args.budgets:
        order = first_signatures(
            budget,
            count=args.count,
            grid_size=args.grid_size,
            tolerance=args.tolerance,
        )
        orders.append(order)
        print(
            f"B={budget}: {len(candidate_signatures(budget))} candidates; "
            f"tail={order[-5:]}"
        )

    if any(order != orders[0] for order in orders[1:]):
        raise AssertionError("the requested budget shells have not stabilized")
    if args.count == len(STABLE_FIRST_99):
        if tuple(orders[0]) != STABLE_FIRST_99:
            raise AssertionError(
                "screened order differs from the stored Appendix B cutoff"
            )
    print(f"stable first {args.count} across budgets {args.budgets}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
