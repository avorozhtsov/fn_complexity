"""Finite-shell searches for comparison clusters of signatures.

The non-special signatures are exhausted by

    S_B = {a : len(a) >= 2, a[0] > 1, a[0] + 2 len(a) <= B}.

For signatures ``a`` and ``b``, the comparison graph contains the arrow

    a -> b  iff  C(a -> b) >= C(b -> a).

This module uses the same vectorized partition-function screening calculation as the
Appendix B stabilization analysis.  NumPy is imported lazily so the rest of
the library remains dependency-free.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement
from typing import Callable

from .core import Signature, normalize_signature


@dataclass(frozen=True)
class ClusterSearchResult:
    """Result of searching the component of one signature in finite shells."""

    target: Signature
    shell_budget: int
    candidate_count: int
    cluster_size: int
    members: tuple[Signature, ...]
    discovered_at: tuple[int, ...]
    requested_maximum: int | None

    @property
    def truncated(self) -> bool:
        """Whether the component had more members than were returned."""

        return self.cluster_size > len(self.members)

    @property
    def reached_requested_maximum(self) -> bool:
        """Whether the search found the requested number of members."""

        return (
            self.requested_maximum is not None
            and len(self.members) == self.requested_maximum
        )


def shell_budget(signature: Signature) -> int:
    """Return the least ``B`` for which a non-special signature is in ``S_B``."""

    return signature[0] + 2 * len(signature)


def candidate_signatures(budget: int) -> list[Signature]:
    """Return all non-special signatures in the finite shell ``S_B``."""

    result: list[Signature] = []
    for length in range(2, max(2, (budget - 1) // 2 + 1)):
        maximum = budget - 2 * length
        if maximum < 2:
            continue
        result.extend(
            tuple(reversed(values))
            for values in combinations_with_replacement(
                range(1, maximum + 1), length
            )
            if values[-1] > 1
        )
    return sorted(result, key=lambda signature: (len(signature), signature))


def approximate_rate_matrix(
    signatures: list[Signature],
    *,
    budget: int,
    grid_size: int,
):
    """Evaluate all directed exchange rates on a shared inverse-T grid."""

    try:
        import numpy as np
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "cluster search requires NumPy; install fn-complexity[analysis]"
        ) from error

    if grid_size < 16:
        raise ValueError("grid_size must be at least 16")
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


def component_for_target(
    signatures: list[Signature],
    rates,
    target: Signature,
    *,
    tolerance: float,
    strict: bool = False,
) -> tuple[Signature, ...]:
    """Return the target's SCC in the strict or non-strict comparison graph."""

    try:
        import numpy as np
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "cluster search requires NumPy; install fn-complexity[analysis]"
        ) from error

    if tolerance < 0:
        raise ValueError("tolerance must be non-negative")
    try:
        target_index = signatures.index(target)
    except ValueError as error:
        raise ValueError(f"target {target} is not in the candidate shell") from error

    # rates[i, j] = C(i -> j).  Arrows point toward the less-complex
    # signature.  In non-strict mode near-equalities are
    # conservatively represented by arrows in both directions.  In strict
    # mode an inequality must exceed the numerical tolerance.
    adjacency = (
        rates > rates.T + tolerance
        if strict
        else rates >= rates.T - tolerance
    )

    def reachable(matrix) -> "np.ndarray":
        seen = np.zeros(len(signatures), dtype=bool)
        frontier = np.zeros(len(signatures), dtype=bool)
        frontier[target_index] = True
        while frontier.any():
            seen |= frontier
            successors = matrix[frontier].any(axis=0)
            frontier = successors & ~seen
        return seen

    strongly_connected = reachable(adjacency) & reachable(adjacency.T)
    return tuple(
        signature for signature, included in zip(signatures, strongly_connected)
        if included
    )


def search_cluster(
    target: Signature,
    *,
    n_max: int | None = None,
    max_b: int = 18,
    grid_size: int = 256,
    tolerance: float = 2.0e-6,
    strict: bool = False,
    progress: Callable[[int, int], None] | None = None,
) -> ClusterSearchResult:
    """Search ``Cl(target)`` in nested finite shells.

    If ``strict`` is false, arrows use ``C(a -> b) >= C(b -> a)``; otherwise they
    use strict inequality.  If ``n_max`` is supplied, shells are examined in
    increasing order and the search stops once that many distinct component
    members have been found.  Otherwise only ``S_max_b`` is evaluated and its
    complete component is returned.
    """

    normalized_target = normalize_signature(target)
    if len(normalized_target) < 2 or normalized_target[0] == 1:
        raise ValueError(
            "cluster shells contain only non-special signatures: "
            "at least two fibers and a largest fiber greater than one"
        )
    if n_max is not None and n_max < 1:
        raise ValueError("n_max must be positive")
    if max_b < 1:
        raise ValueError("max_b must be positive")

    first_budget = shell_budget(normalized_target)
    if first_budget > max_b:
        raise ValueError(
            f"target first belongs to S_{first_budget}, beyond --max-b {max_b}"
        )

    budgets = (
        range(first_budget, max_b + 1)
        if n_max is not None
        else (max_b,)
    )
    discovered: list[Signature] = []
    discovered_set: set[Signature] = set()
    discovered_at: list[int] = []
    last_component: tuple[Signature, ...] = ()
    last_candidates: list[Signature] = []
    last_budget = first_budget

    for budget in budgets:
        candidates = candidate_signatures(budget)
        if progress is not None:
            progress(budget, len(candidates))
        rates = approximate_rate_matrix(
            candidates,
            # A common grid makes the induced graphs genuinely nested even
            # when the search evaluates several successive shells.
            budget=max_b,
            grid_size=grid_size,
        )
        component = component_for_target(
            candidates,
            rates,
            normalized_target,
            tolerance=tolerance,
            strict=strict,
        )
        new_members = sorted(
            (member for member in component if member not in discovered_set),
            key=lambda signature: (
                signature != normalized_target,
                shell_budget(signature),
                len(signature),
                signature,
            ),
        )
        for member in new_members:
            if n_max is not None and len(discovered) >= n_max:
                break
            discovered.append(member)
            discovered_set.add(member)
            discovered_at.append(budget)

        last_component = component
        last_candidates = candidates
        last_budget = budget
        if n_max is not None and len(discovered) >= n_max:
            break

    if n_max is None:
        discovered = sorted(
            last_component,
            key=lambda signature: (
                signature != normalized_target,
                shell_budget(signature),
                len(signature),
                signature,
            ),
        )
        discovered_at = [last_budget] * len(discovered)

    return ClusterSearchResult(
        target=normalized_target,
        shell_budget=last_budget,
        candidate_count=len(last_candidates),
        cluster_size=len(last_component),
        members=tuple(discovered),
        discovered_at=tuple(discovered_at),
        requested_maximum=n_max,
    )
