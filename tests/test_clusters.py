import pytest
import numpy as np

from fn_complexity.clusters import (
    candidate_signatures,
    component_for_target,
    search_cluster,
    shell_budget,
)


def test_candidate_shell_exhaustion_and_budget():
    assert candidate_signatures(5) == []
    assert candidate_signatures(6) == [(2, 1), (2, 2)]
    assert shell_budget((3, 1, 1)) == 9


def test_cluster_search_stops_at_requested_count():
    result = search_cluster(
        (3, 1, 1),
        n_max=3,
        max_b=10,
        grid_size=64,
    )

    assert result.target == (3, 1, 1)
    assert len(result.members) == 3
    assert result.reached_requested_maximum
    assert result.shell_budget <= 10
    assert result.target in result.members


def test_cluster_search_rejects_target_outside_limit():
    with pytest.raises(ValueError, match="beyond --max-b"):
        search_cluster((3, 1, 1), max_b=8)


def test_known_cluster_size_in_shell_12():
    result = search_cluster((3, 1, 1), max_b=12)

    assert result.candidate_count == 129
    assert result.cluster_size == 112
    assert len(result.members) == 112
    assert not result.truncated


def test_strict_mode_does_not_turn_equal_rates_into_arrows():
    signatures = [(2, 1), (2, 2)]
    equal_rates = np.ones((2, 2))

    assert component_for_target(
        signatures,
        equal_rates,
        (2, 1),
        tolerance=1.0e-9,
    ) == tuple(signatures)
    assert component_for_target(
        signatures,
        equal_rates,
        (2, 1),
        tolerance=1.0e-9,
        strict=True,
    ) == ((2, 1),)
