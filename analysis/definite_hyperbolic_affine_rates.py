#!/usr/bin/env python3
"""Affine exchange rates between the definite and the hyperbolic binary form.

The pair is g_+ = x^2 + y^2 and g_- = x^2 - y^2: definite versus indefinite
over R, anisotropic versus split over F_q with q = 3 mod 4.  Three computations
back the corresponding section of paper_finite_fields_maps/main.tex:

1. an exact three-source/two-target certificate with linear processors, so
   C_aff(g_+ -> g_-) >= 2/3 over every field;
2. an exhaustive search over F_3 for k target blocks out of k+1 source copies,
   k = 2, 3, in both directions;
3. a numerical feasibility scan over R of the same span condition.

Items 1 and 2 are exact.  Item 3 is a search: a nonzero residual is evidence
that no conversion exists at that block size, not a proof.  The one proved
upper bound, C_aff(g_- -> g_+) = 1/2, comes from the inertia monotone in the
paper; the scan only illustrates it and is never promoted to a proof.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = PROJECT_ROOT / "analysis"

# The certificate: three copies of x^2 + y^2 implement two copies of x^2 - y^2.
# Input processor: copy i reads the two coordinates CERTIFICATE_INPUTS[i] of
# z = (x1, y1, x2, y2).  Output processor: the rows of CERTIFICATE_OUTPUT.
CERTIFICATE_INPUTS = ((1, 3), (0, 2), (0, 3))
CERTIFICATE_OUTPUT = ((-1, 0, 1), (0, 1, -1))

DEFINITE = ((1, 0), (0, 1))
HYPERBOLIC = ((1, 0), (0, -1))

REAL_SCAN_CASES = (
    ("definite", "hyperbolic", 1, 1),
    ("definite", "hyperbolic", 1, 2),
    ("definite", "hyperbolic", 2, 2),
    ("definite", "hyperbolic", 2, 3),
    ("definite", "hyperbolic", 3, 4),
    ("definite", "hyperbolic", 4, 5),
    ("definite", "hyperbolic", 5, 7),
    ("definite", "hyperbolic", 4, 6),
    ("hyperbolic", "definite", 1, 1),
    ("hyperbolic", "definite", 1, 2),
    ("hyperbolic", "definite", 2, 3),
    ("hyperbolic", "definite", 2, 4),
    ("hyperbolic", "definite", 3, 5),
    ("hyperbolic", "definite", 3, 6),
)

FORMS = {"definite": DEFINITE, "hyperbolic": HYPERBOLIC}


@dataclass(frozen=True)
class FiniteSearch:
    """One exhaustive F_q search for k target blocks out of k+1 source copies."""

    source: str
    target: str
    k: int
    r: int
    pullbacks: int
    target_span_is_pullback_free: bool
    best_difference_rank: int
    found: bool


@dataclass(frozen=True)
class RealScan:
    """One numerical feasibility probe of the span condition over R."""

    source: str
    target: str
    k: int
    r: int
    residual: float
    feasible: bool


def block_matrix(form, k: int, j: int) -> np.ndarray:
    """The matrix of the j-th target block of f^{x k} as a form on K^{2k}."""

    matrix = np.zeros((2 * k, 2 * k), dtype=np.int64)
    matrix[2 * j:2 * j + 2, 2 * j:2 * j + 2] = np.array(form, dtype=np.int64)
    return matrix


def verify_certificate(bound: int = 4, moduli: tuple[int, ...] = (3, 5, 7, 11)) -> bool:
    """Check the three-to-two identity over Z on a box, then over small fields."""

    def evaluate(z, modulus=None):
        outputs = [
            z[a] * z[a] + z[b] * z[b] for a, b in CERTIFICATE_INPUTS
        ]
        got = [
            sum(row[i] * outputs[i] for i in range(3)) for row in CERTIFICATE_OUTPUT
        ]
        want = [z[0] * z[0] - z[1] * z[1], z[2] * z[2] - z[3] * z[3]]
        if modulus is not None:
            got = [value % modulus for value in got]
            want = [value % modulus for value in want]
        return got == want

    if not all(evaluate(z) for z in product(range(-bound, bound + 1), repeat=4)):
        return False
    return all(
        evaluate(z, modulus)
        for modulus in moduli
        for z in product(range(modulus), repeat=4)
    )


def pullback_set(form, n: int, p: int) -> set[bytes]:
    """Every M = A^T G A with A in Mat_{2 x n}(F_p), keyed by its bytes."""

    digits = np.array(list(product(range(p), repeat=2 * n)), dtype=np.int64)
    matrices = digits.reshape(-1, 2, n)
    G = np.array(form, dtype=np.int64)
    pullbacks = np.einsum("iax,ab,iby->ixy", matrices, G, matrices) % p
    return {matrix.astype(np.int8).tobytes() for matrix in pullbacks}


def rank_mod_p(rows, p: int) -> int:
    """Rank over F_p of a list of vectors."""

    matrix = np.array(rows, dtype=np.int64) % p
    if matrix.size == 0:
        return 0
    pivot_row = 0
    row_count, column_count = matrix.shape
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row, column] % p),
            None,
        )
        if pivot is None:
            continue
        matrix[[pivot_row, pivot]] = matrix[[pivot, pivot_row]]
        inverse = pow(int(matrix[pivot_row, column]), p - 2, p)
        matrix[pivot_row] = (matrix[pivot_row] * inverse) % p
        for row in range(row_count):
            if row != pivot_row and matrix[row, column] % p:
                matrix[row] = (matrix[row] - matrix[row, column] * matrix[pivot_row]) % p
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def finite_field_search(source: str, target: str, k: int, p: int = 3) -> FiniteSearch:
    """Exhaustive test of f^{x k} <= g^{x (k+1)} over F_p, homogeneous part.

    With r = k+1 the span W of the pullbacks has dimension at most k+1, so it is
    V + <P> for V the span of the k target blocks and P any pullback outside V.
    The pullback set is closed under nonzero scaling, hence every generator can
    be normalized to P + v with v in V, and the search reduces to: is there a
    pullback P whose set {v in V : P + v is a pullback} spans V?  The degenerate
    alternative, W = V, is excluded by checking that V contains no nonzero
    pullback; that check is reported rather than assumed.
    """

    n = 2 * k
    r = k + 1
    keys = pullback_set(FORMS[source], n, p)
    blocks = [block_matrix(FORMS[target], k, j) % p for j in range(k)]
    span = [
        sum(
            (c * block for c, block in zip(coefficients, blocks)),
            np.zeros((n, n), dtype=np.int64),
        )
        % p
        for coefficients in product(range(p), repeat=k)
    ]
    pullback_free = not any(
        element.any() and element.astype(np.int8).tobytes() in keys
        for element in span
    )
    best_rank = 0
    found = False
    for key in keys:
        candidate = np.frombuffer(key, dtype=np.int8).reshape(n, n).astype(np.int64)
        if not candidate.any():
            continue
        shifts = [
            element
            for element in span
            if ((candidate + element) % p).astype(np.int8).tobytes() in keys
        ]
        if len(shifts) <= k:
            continue
        rank = rank_mod_p([element.ravel() for element in shifts], p)
        best_rank = max(best_rank, rank)
        if rank == k:
            found = True
            break
    return FiniteSearch(
        source=source,
        target=target,
        k=k,
        r=r,
        pullbacks=len(keys),
        target_span_is_pullback_free=pullback_free,
        best_difference_rank=best_rank,
        found=found,
    )


def real_probe(source: str, target: str, k: int, r: int, restarts: int, seed: int) -> RealScan:
    """Minimize the distance from the target blocks to the span of r pullbacks."""

    from scipy.optimize import minimize

    G = np.array(FORMS[source], dtype=float)
    targets = np.array(
        [block_matrix(FORMS[target], k, j).astype(float).ravel() for j in range(k)]
    )

    def loss(parameters: np.ndarray) -> float:
        A = parameters.reshape(r, 2, 2 * k)
        pullbacks = np.einsum("iax,ab,iby->ixy", A, G, A).reshape(r, -1)
        norms = np.linalg.norm(pullbacks, axis=1, keepdims=True)
        pullbacks = pullbacks / np.maximum(norms, 1e-9)
        coefficients, *_ = np.linalg.lstsq(pullbacks.T, targets.T, rcond=None)
        return float(((pullbacks.T @ coefficients - targets.T) ** 2).sum())

    generator = np.random.default_rng(seed)
    best = np.inf
    for _ in range(restarts):
        result = minimize(
            loss,
            generator.normal(size=r * 4 * k),
            method="L-BFGS-B",
            options=dict(maxiter=5000, ftol=1e-16, gtol=1e-12),
        )
        best = min(best, float(result.fun))
        if best < 1e-16:
            break
    return RealScan(source, target, k, r, best, best < 1e-14)


def write_finite_searches(searches, output: Path) -> None:
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "source",
                "target",
                "k",
                "r",
                "rate",
                "pullbacks",
                "target_span_pullback_free",
                "best_difference_rank",
                "conversion_found",
            )
        )
        for search in searches:
            writer.writerow(
                (
                    search.source,
                    search.target,
                    search.k,
                    search.r,
                    f"{search.k / search.r:.6f}",
                    search.pullbacks,
                    str(search.target_span_is_pullback_free).lower(),
                    search.best_difference_rank,
                    str(search.found).lower(),
                )
            )


def write_real_scan(scans, output: Path) -> None:
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("source", "target", "k", "r", "rate", "residual", "feasible"))
        for scan in scans:
            writer.writerow(
                (
                    scan.source,
                    scan.target,
                    scan.k,
                    scan.r,
                    f"{scan.k / scan.r:.6f}",
                    f"{scan.residual:.3e}",
                    str(scan.feasible).lower(),
                )
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-k", type=int, default=3, help="largest F_3 block count")
    parser.add_argument("--restarts", type=int, default=60)
    parser.add_argument("--skip-real", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ANALYSIS_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    verified = verify_certificate()
    print(f"three-to-two certificate verified: {str(verified).lower()}")
    if not verified:
        return 1

    searches = [
        finite_field_search(source, target, k)
        for k in range(2, args.max_k + 1)
        for source, target in (("definite", "hyperbolic"), ("hyperbolic", "definite"))
    ]
    finite_path = args.output_dir / "definite_hyperbolic_f3_search.csv"
    write_finite_searches(searches, finite_path)
    for search in searches:
        print(
            f"F_3 {search.source} -> {search.target}: k={search.k} r={search.r} "
            f"conversion found = {str(search.found).lower()} "
            f"({search.pullbacks} pullbacks, best difference rank "
            f"{search.best_difference_rank} of {search.k})"
        )
    print(finite_path)

    if args.skip_real:
        return 0

    scans = [
        real_probe(source, target, k, r, args.restarts, 97 * k + r)
        for source, target, k, r in REAL_SCAN_CASES
    ]
    real_path = args.output_dir / "definite_hyperbolic_real_scan.csv"
    write_real_scan(scans, real_path)
    for scan in scans:
        print(
            f"R  {scan.source} -> {scan.target}: k={scan.k} r={scan.r} "
            f"rate={scan.k / scan.r:.4f} residual={scan.residual:.2e} "
            f"{'feasible' if scan.feasible else 'not found'}"
        )
    print(real_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
