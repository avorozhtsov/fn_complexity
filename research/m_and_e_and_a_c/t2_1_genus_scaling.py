#!/usr/bin/env python3
"""T2.1 -- the genus scaling law ``1 - C(L -> f) ~ 2g / (sqrt(q) log q)``.

For ``f(x, y) = y^2 - h(x)`` with ``deg h = 2g + 1`` over a prime field ``F_q``
the fibers are the affine hyperelliptic curves ``y^2 = h(x) + c`` of genus ``g``.
Because the degree is odd the smooth projective model has exactly one rational
point at infinity, so the affine count is ``N_c = q - a_c`` with ``a_c`` the
trace of Frobenius, and ``|a_c| <= 2 g sqrt(q)`` (Weil).

Two exact facts do all the work (see the note for proofs):

* ``sum_c a_c = 0`` because ``sum_c N_c = q^2``;
* ``C(L -> f) = log q / log(max_c N_c)`` exactly, the infimum always being at
  ``beta = infinity``, for *any* ``f : A^2 -> A^1``.

So with ``m = max_c(-a_c) = max_c N_c - q`` we have the closed form

    (1 - C) sqrt(q) log q = sqrt(q) log q * log(1 + m/q) / log(q + m)
                          = mu - mu^2 (1/2 + 1/log q) / sqrt(q) + O(1/q),

where ``mu = m / sqrt(q)``.  Testing the "2g law" therefore means testing
``mu -> 2g``, i.e. the *extreme value* statistics of the trace distribution.

Traces are computed exactly (no ``q^2`` loop) from

    N_c = sum_x (1 + chi(h(x) + c)) = q + (m_h * chi)[c],

a circular cross-correlation of the value-multiplicity vector of ``h`` with the
Legendre symbol, evaluated by FFT.  The result is rounded to integers; the FFT
error is ~1e-9 even at q ~ 10^6.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate_result  # noqa: E402

OUTPUT_DIRECTORY = Path(__file__).resolve().parent


# --------------------------------------------------------------------------
# families
# --------------------------------------------------------------------------
# Coefficients are listed from the leading term down; the constant term is
# irrelevant because the family sweeps h(x) + c over all c in F_q.
FAMILIES: dict[str, tuple[int, tuple[int, ...]]] = {
    # genus 0 controls
    "G0a": (0, (1, 0)),                                  # h = x,  f = y^2 - x: parabolas
    "G0b": (0, (1, 1, 0)),                               # h = x^2 + x: split conics
    "G0c": (0, (-1, 0, 0)),                              # h = -x^2, f = x^2 + y^2
    # genus 1
    "E1": (1, (1, 0, 1, 0)),                             # h = x^3 + x (the paper's family)
    "E1b": (1, (1, 1, 0, 0)),                            # h = x^3 + x^2
    "E1c": (1, (1, 0, 3, 0)),                            # h = x^3 + 3x
    "E1d": (1, (1, 2, 1, 0)),                            # h = x^3 + 2x^2 + x
    # genus 2
    "H2": (2, (1, 0, 0, 1, 0, 0)),                       # h = x^5 + x^2
    "H2b": (2, (1, 1, 0, 0, 1, 0)),                      # h = x^5 + x^4 + x
    "H2c": (2, (1, 0, 1, 1, 0, 0)),                      # h = x^5 + x^3 + x^2
    "H2d": (2, (1, 3, 1, 2, 1, 0)),                      # h = x^5 + 3x^4 + x^3 + 2x^2 + x
    # genus 3
    "H3": (3, (1, 0, 0, 0, 0, 1, 0, 0)),                 # h = x^7 + x^2
    "H3b": (3, (1, 0, 1, 0, 1, 0, 1, 0)),                # h = x^7 + x^5 + x^3 + x
    "H3c": (3, (1, 1, 0, 2, 0, 0, 1, 0)),                # h = x^7 + x^6 + 2x^4 + x
    "H3d": (3, (1, 2, 1, 0, 3, 1, 0, 0)),                # h = x^7 + 2x^6 + x^5 + 3x^3 + x^2
    # genus 4, for the exponent fit only
    "H4": (4, (1, 0, 0, 0, 0, 0, 0, 1, 0, 0)),           # h = x^9 + x^2
    "H4b": (4, (1, 1, 0, 2, 0, 0, 1, 0, 1, 0)),          # h = x^9 + x^8 + 2x^6 + x^3 + x
}

SPLIT_CONIC = "XY"  # f = x*y, handled separately (not of the form y^2 - h(x))


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_at_least(n: int) -> int:
    while not is_prime(n):
        n += 1
    return n


def legendre_table(q: int) -> np.ndarray:
    """chi[v] in {-1, 0, 1}: the quadratic character of F_q, chi[0] = 0."""

    chi = -np.ones(q, dtype=np.int64)
    squares = (np.arange(1, q, dtype=np.int64) ** 2) % q
    chi[squares] = 1
    chi[0] = 0
    return chi


def horner_values(coefficients: tuple[int, ...], q: int) -> np.ndarray:
    """h(x) mod q for x = 0, ..., q-1, without overflow."""

    x = np.arange(q, dtype=np.int64)
    value = np.zeros(q, dtype=np.int64)
    for coefficient in coefficients:
        value = (value * x + coefficient) % q
    return value


def traces_hyperelliptic(coefficients: tuple[int, ...], q: int) -> np.ndarray:
    """a_c = q - #{(x,y) : y^2 = h(x) + c} for all c in F_q."""

    multiplicity = np.bincount(horner_values(coefficients, q), minlength=q).astype(float)
    chi = legendre_table(q).astype(float)
    correlation = np.fft.irfft(
        np.conj(np.fft.rfft(multiplicity)) * np.fft.rfft(chi), n=q
    )
    return -np.rint(correlation).astype(np.int64)


def traces_hyperelliptic_direct(coefficients: tuple[int, ...], q: int) -> np.ndarray:
    """Brute-force reference implementation, O(q^2), for validation only."""

    h = horner_values(coefficients, q)
    y = np.arange(q, dtype=np.int64)
    values = (y[:, None] ** 2 - h[None, :]) % q
    return q - np.bincount(values.ravel(), minlength=q)


def traces_split_conic(q: int) -> np.ndarray:
    """a_c for f(x, y) = x y: N_c = q - 1 for c != 0 and N_0 = 2q - 1."""

    traces = np.ones(q, dtype=np.int64)
    traces[0] = -(q - 1)
    return traces


def rate_closed_form(q: int, max_count: int) -> float:
    return math.log(q) / math.log(max_count)


def scaled_deviation(q: int, max_count: int) -> float:
    return (1.0 - rate_closed_form(q, max_count)) * math.sqrt(q) * math.log(q)


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------
def validate(verbose: bool = True) -> None:
    """FFT traces vs brute force, and the closed-form rate vs the repo solver."""

    problems = []
    for q in (101, 211, 401):
        for name, (genus, coefficients) in FAMILIES.items():
            fast = traces_hyperelliptic(coefficients, q)
            slow = traces_hyperelliptic_direct(coefficients, q)
            if not np.array_equal(fast, slow):
                problems.append(f"FFT != brute force for {name} at q={q}")
            if fast.sum() != 0:
                problems.append(f"sum a_c != 0 for {name} at q={q}")
            bound = 2 * genus * math.sqrt(q)
            smooth_violations = int((np.abs(fast) > bound).sum())
            if smooth_violations and verbose:
                print(
                    f"  note: {name} q={q} has {smooth_violations} fibers with "
                    f"|a_c| > 2g sqrt(q) (singular fibers of the pencil)"
                )
    # closed form for C(L -> f) against the numerical infimum solver
    for q in (101, 211, 401, 1009):
        for name in ("E1", "H2", "H3", "H4"):
            genus, coefficients = FAMILIES[name]
            counts = q - traces_hyperelliptic(coefficients, q)
            solver = exchange_rate_result(
                implementer=tuple([q] * q), implemented=tuple(counts.tolist())
            )
            closed = rate_closed_form(q, int(counts.max()))
            if abs(solver.rate - closed) > 1e-12 or not math.isinf(solver.beta):
                problems.append(
                    f"solver {solver.rate!r} beta={solver.beta} vs closed {closed!r} "
                    f"for {name} at q={q}"
                )
    if verbose:
        if problems:
            for problem in problems:
                print("  FAIL:", problem)
        else:
            print("  all validation checks passed")


# --------------------------------------------------------------------------
# main sweep
# --------------------------------------------------------------------------
TARGETS = [100, 200, 400, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000,
           256000, 512000, 1000000]


def main() -> int:
    primes = [prime_at_least(t) for t in TARGETS]

    print("=" * 78)
    print("validation")
    print("=" * 78)
    validate()

    rows: list[list] = []
    print()
    for name, (genus, coefficients) in FAMILIES.items():
        print("=" * 78)
        print(f"family {name}: h = {poly_string(coefficients)}, genus {genus}, "
              f"Weil bound 2g = {2 * genus}")
        print("=" * 78)
        print(f"{'q':>8} {'q%4':>4} {'m=max(-a)':>10} {'mu=m/sqrt q':>12} "
              f"{'max|a|/sqrt q':>14} {'(1-C)sqrt q log q':>18} "
              f"{'2g-mu':>9} {'m2':>7} {'m4':>7} {'#|a|>2g sq':>11}")
        for q in primes:
            traces = traces_hyperelliptic(coefficients, q)
            record = summarize(name, genus, q, traces)
            rows.append(record)
            print(format_row(record))
        print()

    # split conic, genus 0 but with a reducible fiber
    print("=" * 78)
    print("family XY: f = x*y (genus 0 fibers, but the fiber over 0 is reducible)")
    print("=" * 78)
    print(f"{'q':>8} {'q%4':>4} {'m=max(-a)':>10} {'mu=m/sqrt q':>12} "
          f"{'max|a|/sqrt q':>14} {'(1-C)sqrt q log q':>18} "
          f"{'2g-mu':>9} {'m2':>7} {'m4':>7} {'#|a|>2g sq':>11}")
    for q in primes:
        traces = traces_split_conic(q)
        record = summarize(SPLIT_CONIC, 0, q, traces)
        rows.append(record)
        print(format_row(record))
    print()

    corrections(primes)
    endpoint_agreement(rows)

    path = OUTPUT_DIRECTORY / "t2_1_genus_scaling.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["family", "genus", "q", "m", "mu", "mu_abs", "rate", "scaled_deviation",
             "deficit", "moment2", "moment4", "n_over_weil"]
        )
        writer.writerows(rows)
    print(f"written to {path.name}")
    return 0


def corrections(primes: list[int]) -> None:
    """The deterministic gap between the rate quantity and mu = m/sqrt(q).

    (1 - C) sqrt(q) log q = mu - mu^2 (1/2 + 1/log q) / sqrt(q) + O(1/q),
    so the rate quantity always *undershoots* mu, by an amount that is
    itself of order 1/sqrt(q) -- the same order as the genus-1 deficit.
    """

    print("=" * 78)
    print("deterministic gap  mu - (1-C) sqrt(q) log q  vs its leading term")
    print("=" * 78)
    print(f"{'family':>6} {'q':>8} {'mu':>8} {'rate quantity':>14} {'gap':>9} "
          f"{'mu^2(1/2+1/log q)/sqrt q':>25}")
    for name in ("E1", "H2", "H3"):
        genus, coefficients = FAMILIES[name]
        for q in primes:
            if q not in (101, 2003, 32003, 1000003):
                continue
            traces = traces_hyperelliptic(coefficients, q)
            m = int((-traces).max())
            root = math.sqrt(q)
            mu = m / root
            scaled = scaled_deviation(q, q + m)
            leading = mu**2 * (0.5 + 1.0 / math.log(q)) / root
            print(f"{name:>6} {q:>8} {mu:>8.4f} {scaled:>14.4f} {mu - scaled:>9.4f} "
                  f"{leading:>25.4f}")
    print()


def endpoint_agreement(rows: list[list]) -> None:
    """Does the rate quantity track max_c(-a_c) or max_c |a_c|?

    The rate only sees the *largest fiber*, so it tracks max_c(-a_c).  The two
    coincide exactly when the extreme trace happens to be negative.
    """

    print("=" * 78)
    print("max_c(-a_c) versus max_c |a_c| (only the first drives the rate)")
    print("=" * 78)
    print(f"{'genus':>6} {'families x primes':>18} {'agree':>7} {'share':>7} "
          f"{'mean |a| / (-a) - 1':>20}")
    for genus in (1, 2, 3, 4):
        selected = [row for row in rows if row[1] == genus and row[0] != SPLIT_CONIC]
        agree = sum(1 for row in selected if abs(row[4] - row[5]) < 1e-12)
        excess = float(np.mean([row[5] / row[4] - 1.0 for row in selected]))
        print(f"{genus:>6} {len(selected):>18} {agree:>7} "
              f"{agree / len(selected):>7.2f} {excess:>20.4f}")
    print()


def poly_string(coefficients: tuple[int, ...]) -> str:
    degree = len(coefficients) - 1
    parts = []
    for index, coefficient in enumerate(coefficients):
        power = degree - index
        if coefficient == 0:
            continue
        if power == 0:
            parts.append(str(coefficient))
        else:
            head = "" if coefficient == 1 else str(coefficient)
            parts.append(f"{head}x^{power}" if power > 1 else f"{head}x")
    return " + ".join(parts) or "0"


def summarize(name: str, genus: int, q: int, traces: np.ndarray) -> list:
    assert int(traces.sum()) == 0, f"sum a_c != 0 for {name} at q={q}"
    m = int((-traces).max())
    root = math.sqrt(q)
    mu = m / root
    mu_abs = float(np.abs(traces).max()) / root
    counts = q - traces
    rate = rate_closed_form(q, int(counts.max()))
    scaled = (1.0 - rate) * root * math.log(q)
    moment2 = float((traces.astype(float) ** 2).sum()) / q / q
    moment4 = float((traces.astype(float) ** 4).sum()) / q / q**2
    over = int((np.abs(traces) > 2 * genus * root).sum())
    return [name, genus, q, m, mu, mu_abs, rate, scaled, 2 * genus - mu,
            moment2, moment4, over]


def format_row(record: list) -> str:
    (_, genus, q, m, mu, mu_abs, _, scaled, deficit, moment2, moment4, over) = record
    return (f"{q:>8} {q % 4:>4} {m:>10} {mu:>12.4f} {mu_abs:>14.4f} {scaled:>18.4f} "
            f"{deficit:>9.4f} {moment2:>7.3f} {moment4:>7.3f} {over:>11}")


if __name__ == "__main__":
    raise SystemExit(main())
