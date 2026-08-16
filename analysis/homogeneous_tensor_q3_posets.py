#!/usr/bin/env python3
"""Exhaustively verify the tractable homogeneous-tensor posets over F_3."""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity.homogeneous_tensor_maps import (  # noqa: E402
    TENSOR_CASES,
    compute_tensor_poset,
    minimum_orbit_counts_for_case6,
    rank_counts_for_case6,
)


EXPECTED = {
    1: (5, 5),
    2: (7, 8),
    3: (50, 210),
    4: (26, 66),
    5: (19, 32),
}


def main() -> int:
    print("case\ttensors\tclasses\tcovers\tcheck")
    for case_number, (expected_classes, expected_covers) in EXPECTED.items():
        case = TENSOR_CASES[case_number]
        poset = compute_tensor_poset(case_number)
        assert len(poset.orbits) == expected_classes
        assert len(poset.covers) == expected_covers
        assert sum(orbit.size for orbit in poset.orbits) == case.tensor_count
        print(
            f"{case_number}\t{case.tensor_count}\t{len(poset.orbits)}\t"
            f"{len(poset.covers)}\tok"
        )

    counts = rank_counts_for_case6()
    assert sum(counts.values()) == TENSOR_CASES[6].tensor_count
    minima = minimum_orbit_counts_for_case6()
    rank1_exact = len(compute_tensor_poset(4).orbits) - 1
    lower_bound = minima[3] + minima[2] + rank1_exact + 1
    assert lower_bound == 1_632_040
    print(
        f"6\t{TENSOR_CASES[6].tensor_count}\t>={lower_bound}\tN/A\t"
        "rank-stratum lower bound verified"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
