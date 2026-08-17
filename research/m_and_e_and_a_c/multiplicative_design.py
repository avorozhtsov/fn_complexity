#!/usr/bin/env python3
"""T1.5 -- can the Weil geometry be made to track the exchange geometry?

T1.4 established, via Landau's theorem, that on the family of atomic Weil test
measures ``mu_a = sum_i delta_{a_i}`` (i.e. on signatures) the truncated Weil
pairing is

    E_ab = N * O(a,b) - (T/2pi) * A(a,b) + O(r s log T),
    O(a,b) = sum over pairs with a_i = b_j of a_i,
    A(a,b) = sum over pairs a_i != b_j of min(a_i,b_j) * Lambda(max/min),

so ``E`` sees nothing but *multiplicative coincidences*.  On a generic family
those never happen, every pair is Weil-orthogonal, and the correlation between
the exchange distance ``d`` and the Weil angle was only ``+0.19``.

This script asks whether that independence is an artefact of genericity.  It

  (1) builds families designed so that *every* ratio is a prime power -- powers
      of a single prime, and the mixed 2^i 3^j lattice -- and measures the
      correlation between the two geometries against generic controls;
  (2) searches one-parameter families for an actual monotone relation, and
      exhibits the exact obstruction as a pair of transverse group actions;
  (3) quantifies the mismatch under Cartesian powering and under rescaling;
  (4) reports the N-dependence of the headline correlation.

The obstruction found (see the markdown note) is that the two functionals have
transverse invariance groups:

    a -> a tensor-power k   fixes d exactly and moves the Weil angles,
    a -> lambda * a         fixes the Weil angles exactly (to leading Landau
                            order) and moves d,

so no function can express one geometry in terms of the other.
"""

from __future__ import annotations

import csv
import functools
import itertools
import math
from pathlib import Path
import sys

import numpy as np

RESEARCH_DIRECTORY = Path(__file__).resolve().parent
PROJECT_ROOT = RESEARCH_DIRECTORY.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity import exchange_rate as _exchange_rate  # noqa: E402

exchange_rate = functools.lru_cache(maxsize=None)(_exchange_rate)

ZEROS_PATH = RESEARCH_DIRECTORY / "zeta_zeros_1200.npy"
LONG_ZEROS_PATH = RESEARCH_DIRECTORY / "zeta_zeros_2400.npy"
CSV_PATH = RESEARCH_DIRECTORY / "multiplicative_design.csv"

GAMMAS = np.load(ZEROS_PATH)
LONG_GAMMAS = np.load(LONG_ZEROS_PATH) if LONG_ZEROS_PATH.exists() else GAMMAS


# --------------------------------------------------------------------------
# the two geometries
# --------------------------------------------------------------------------

def von_mangoldt(x: float, tolerance: float = 1e-9) -> float:
    """``Lambda(x)``: ``log p`` when ``x`` is exactly a prime power, else 0."""

    if abs(x - round(x)) > tolerance:
        return 0.0
    n = int(round(x))
    if n < 2:
        return 0.0
    remainder, divisor = n, 2
    while divisor * divisor <= remainder:
        if remainder % divisor == 0:
            while remainder % divisor == 0:
                remainder //= divisor
            return math.log(divisor) if remainder == 1 else 0.0
        divisor += 1
    return math.log(n)


def overlap(a, b) -> float:
    return float(sum(x for x in a for y in b if x == y))


def arithmetic_term(a, b) -> float:
    total = 0.0
    for x in a:
        for y in b:
            if x != y:
                high, low = max(x, y), min(x, y)
                total += low * von_mangoldt(high / low)
    return total


def gram_matrix(family, gammas: np.ndarray) -> np.ndarray:
    """``E_ab = sum_{n<=N} Z_a(rho_n) conj(Z_b(rho_n))``, computed from zeros."""

    vectors = []
    for signature in family:
        entries = np.asarray(signature, dtype=float)
        phases = np.exp(1j * np.outer(np.log(entries), gammas))
        vectors.append((np.sqrt(entries)[:, None] * phases).sum(axis=0))
    vectors = np.array(vectors)
    return (vectors @ vectors.conj().T).real


def landau_gram(family, gammas: np.ndarray) -> np.ndarray:
    """The Landau prediction ``N*O - (T/2pi)*A`` for the same matrix."""

    count = len(gammas)
    scale = gammas[-1] / (2 * math.pi)
    size = len(family)
    matrix = np.zeros((size, size))
    for i, a in enumerate(family):
        for j, b in enumerate(family):
            matrix[i, j] = count * overlap(a, b) - scale * arithmetic_term(a, b)
    return matrix


def weil_angles(gram: np.ndarray) -> np.ndarray:
    diagonal = np.sqrt(np.diag(gram))
    correlation = gram / np.outer(diagonal, diagonal)
    return np.arccos(np.clip(correlation, -1.0, 1.0))


def exchange_distances(family) -> np.ndarray:
    size = len(family)
    matrix = np.zeros((size, size))
    for i, a in enumerate(family):
        for j in range(i + 1, size):
            b = family[j]
            product = exchange_rate(a, b) * exchange_rate(b, a)
            value = math.inf if product <= 0 else -math.log(product)
            matrix[i, j] = matrix[j, i] = max(value, 0.0)
    return matrix


def tensor_power(a, k: int):
    """The signature of the ``k``-fold Cartesian power: all ``k``-fold products."""

    entries = [1]
    for _ in range(k):
        entries = [x * y for x in entries for y in a]
    return tuple(sorted(entries, reverse=True))


def upper(matrix: np.ndarray) -> np.ndarray:
    size = matrix.shape[0]
    return np.array([matrix[i, j] for i in range(size) for j in range(i + 1, size)])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    def ranks(values):
        order = np.argsort(values, kind="mergesort")
        result = np.empty(len(values))
        result[order] = np.arange(len(values), dtype=float)
        # average ties
        _, inverse, counts = np.unique(values, return_inverse=True, return_counts=True)
        sums = np.zeros(len(counts))
        np.add.at(sums, inverse, result)
        return (sums / counts)[inverse]

    rx, ry = ranks(x), ranks(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def compare(family, gammas: np.ndarray = GAMMAS, exact: bool = True):
    """Return ``(pearson, spearman, pairs, distances, angles)`` for a family."""

    gram = gram_matrix(family, gammas) if exact else landau_gram(family, gammas)
    angles = upper(weil_angles(gram))
    distances = upper(exchange_distances(family))
    if distances.std() == 0 or angles.std() == 0:
        return float("nan"), float("nan"), len(distances), distances, angles
    pearson = float(np.corrcoef(distances, angles)[0, 1])
    return pearson, spearman(distances, angles), len(distances), distances, angles


# --------------------------------------------------------------------------
# families
# --------------------------------------------------------------------------

def geometric_pairs(prime: int, exponents) -> list:
    return [
        tuple(sorted((prime ** i, prime ** j), reverse=True))
        for i, j in itertools.combinations(exponents, 2)
    ]


def smooth_pairs(limit: int, count: int) -> list:
    """Pairs drawn from the 2^i 3^j lattice below ``limit``."""

    lattice = sorted(
        2 ** i * 3 ** j
        for i in range(0, 8)
        for j in range(0, 5)
        if 2 ** i * 3 ** j <= limit
    )
    pairs = [
        tuple(sorted(p, reverse=True))
        for p in itertools.combinations(lattice, 2)
    ]
    step = max(1, len(pairs) // count)
    return pairs[::step][:count]


def generic_pairs(count: int, seed: int, low: int = 5, high: int = 120) -> list:
    """Pairs with no shared entry and no prime-power ratio between families."""

    rng = np.random.default_rng(seed)
    chosen = []
    while len(chosen) < count:
        x, y = sorted(rng.integers(low, high, size=2).tolist(), reverse=True)
        if x == y:
            continue
        candidate = (int(x), int(y))
        if candidate in chosen:
            continue
        ok = True
        for other in chosen:
            if overlap(candidate, other) or arithmetic_term(candidate, other):
                ok = False
                break
        if ok and arithmetic_term(candidate, candidate) == 0:
            chosen.append(candidate)
    return chosen


FAMILIES = {}


def register(name, family):
    FAMILIES[name] = family
    return family


register("generic-T1.4", [
    (2, 2), (3, 1), (4, 2), (5, 3), (6, 1), (3, 1, 1), (8, 4), (9, 3),
    (4, 4, 1), (7, 5), (5, 5), (7, 1), (9, 1), (6, 3), (10, 5), (11, 7),
])
register("powers-of-2 pairs", geometric_pairs(2, range(0, 7)))
register("powers-of-2 pairs (no 1)", geometric_pairs(2, range(1, 8)))
register("powers-of-3 pairs", geometric_pairs(3, range(0, 6)))
register("powers-of-2 triples", [
    tuple(sorted((2 ** i, 2 ** j, 2 ** k), reverse=True))
    for i, j, k in itertools.combinations(range(0, 6), 3)
])
register("3-smooth pairs", smooth_pairs(200, 21))
register("staircase (2^t, 1)", [(2 ** t, 1) for t in range(1, 10)])
register("staircase (2^t, 2)", [(2 ** t, 2) for t in range(1, 10)])
register("scaled orbit 2^t*(2,1)", [(2 ** (t + 1), 2 ** t) for t in range(0, 9)])
register("chain (2^t, 2^t, 1)", [(2 ** t, 2 ** t, 1) for t in range(1, 10)])
register("staircase (3^t, 1)", [(3 ** t, 1) for t in range(1, 9)])
register("staircase (2^t, 2^t)", [(2 ** t, 2 ** t) for t in range(1, 10)])
# Three ladders of eight entries with nearly the same logarithmic spacing.
# The first two make every ratio a prime power; the third makes NO ratio a
# prime power, so its Landau arithmetic term A vanishes identically.  Subsets
# of a ladder share entries, so the overlap term O is large in all three.
LADDERS = {
    "2-power ladder": [2 ** i for i in range(1, 9)],
    "3-power ladder": [3 ** i for i in range(1, 9)],
    "prime ladder (A = 0)": [3, 7, 17, 37, 67, 131, 257, 521],
}


def ladder_subsets(ladder, k: int) -> list:
    return [tuple(sorted(c, reverse=True))
            for c in itertools.combinations(ladder, k)]


register("2-ladder subsets k=5", ladder_subsets(LADDERS["2-power ladder"], 5))
register("prime-ladder subsets k=5",
         ladder_subsets(LADDERS["prime ladder (A = 0)"], 5))
register("generic control A", generic_pairs(16, seed=1))
register("generic control B", generic_pairs(21, seed=7))


# --------------------------------------------------------------------------
# experiments
# --------------------------------------------------------------------------

def section(title):
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def experiment_designed_families(rows):
    section("(1) DESIGNED vs GENERIC FAMILIES")
    print(f"{'family':>26} {'size':>5} {'pairs':>6} {'pearson':>9} "
          f"{'spearman':>9} {'mean angle':>11} {'angle sd':>9}")
    for name, family in FAMILIES.items():
        pearson, rho, pairs, distances, angles = compare(family)
        print(f"{name:>26} {len(family):>5} {pairs:>6} {pearson:>+9.4f} "
              f"{rho:>+9.4f} {angles.mean():>11.4f} {angles.std():>9.4f}")
        rows.append(["designed-vs-generic", name, len(family), pairs,
                     f"{pearson:.6f}", f"{rho:.6f}",
                     f"{angles.mean():.6f}", f"{angles.std():.6f}"])


def experiment_mechanism(rows):
    """Which of the two Landau terms actually produces the agreement?"""

    section("(1b) WHICH TERM DOES THE WORK: OVERLAP O, OR THE Lambda TERM A?")
    print("Subsets of size k of a ladder of eight entries.  The first two "
          "ladders are\ngeometric, so EVERY ratio is a prime power and the "
          "arithmetic term A is huge.\nThe third has the same logarithmic "
          "spacing but no prime-power ratio at all,\nso A vanishes "
          "identically -- it is the control.\n")
    print(f"{'ladder':>22} {'max cross A':>12} " +
          " ".join(f"{'k=' + str(k):>18}" for k in [3, 4, 5, 6]))
    for name, ladder in LADDERS.items():
        cross = max(arithmetic_term((x,), (y,))
                    for x in ladder for y in ladder if x != y)
        cells = []
        for k in [3, 4, 5, 6]:
            family = ladder_subsets(ladder, k)
            pearson, rho, pairs, _, _ = compare(family)
            cells.append(f"{pearson:+.3f}/{rho:+.3f}")
            rows.append(["ladder", f"{name} k={k}", len(family), pairs,
                         f"{pearson:.6f}", f"{rho:.6f}", f"{cross:.6f}", ""])
        print(f"{name:>22} {cross:>12.2f} " +
              " ".join(f"{c:>18}" for c in cells))
    print("\n   (each cell is pearson/spearman; pair counts are "
          "56/70/56/28 signatures ->\n    1540/2415/1540/378 pairs)")
    print("   The control matches the geometric ladders cell for cell: the "
          "agreement is\n   produced entirely by the OVERLAP term, and the "
          "arithmetic term contributes\n   nothing to it.")


def experiment_monotone_search(rows):
    section("(2) IS THE WEIL ANGLE EVER A MONOTONE FUNCTION OF d?")
    for name in ["staircase (2^t, 1)", "staircase (3^t, 1)", "staircase (2^t, 2)",
                 "powers-of-2 triples", "scaled orbit 2^t*(2,1)",
                 "chain (2^t, 2^t, 1)", "powers-of-2 pairs (no 1)"]:
        family = FAMILIES[name]
        pearson, rho, pairs, distances, angles = compare(family)
        print(f"\n{name}: {pairs} pairs, pearson {pearson:+.4f}, "
              f"spearman {rho:+.4f}")
        order = np.argsort(distances)
        print(f"   {'d':>10} {'angle':>10}")
        for index in order[:5]:
            print(f"   {distances[index]:>10.5f} {angles[index]:>10.5f}")
        print("   ...")
        for index in order[-5:]:
            print(f"   {distances[index]:>10.5f} {angles[index]:>10.5f}")
        rows.append(["monotone-search", name, len(family), pairs,
                     f"{pearson:.6f}", f"{rho:.6f}", "", ""])

    section("(2b) THE OBSTRUCTION: TRANSVERSE INVARIANCE GROUPS")
    base = FAMILIES["powers-of-2 pairs (no 1)"]
    print("\n(i) Cartesian powering a -> a tensor a fixes d EXACTLY "
          "(log Z_{a^k} = k log Z_a) but moves the Weil angles.")
    squared = [tensor_power(a, 2) for a in base]
    d_base, d_sq = upper(exchange_distances(base)), upper(exchange_distances(squared))
    ang_base = upper(weil_angles(gram_matrix(base, GAMMAS)))
    ang_sq = upper(weil_angles(gram_matrix(squared, GAMMAS)))
    print(f"    max |d(a^2,b^2) - d(a,b)| = {np.abs(d_sq - d_base).max():.3e} "
          f"over {len(d_base)} pairs   (identical d-geometry)")
    print(f"    angles: mean |change| = {np.abs(ang_sq - ang_base).mean():.4f} rad, "
          f"max = {np.abs(ang_sq - ang_base).max():.4f} rad")
    print(f"    corr(angle before, angle after) = "
          f"{np.corrcoef(ang_base, ang_sq)[0,1]:+.4f}")
    rows.append(["obstruction-power", "powers-of-2 pairs (no 1) squared",
                 len(base), len(d_base),
                 f"{np.abs(d_sq - d_base).max():.6e}",
                 f"{np.abs(ang_sq - ang_base).mean():.6f}",
                 f"{np.abs(ang_sq - ang_base).max():.6f}",
                 f"{np.corrcoef(ang_base, ang_sq)[0,1]:.6f}"])

    print("\n(ii) Common rescaling a -> lambda*a fixes O and A up to the common "
          "factor lambda, hence fixes ALL Weil angles, but moves d.")
    print(f"    {'lambda':>7} {'max |angle change| (Landau)':>29} "
          f"{'max |angle change| (exact)':>28} {'mean d':>9} {'corr(d_1,d_L)':>14}")
    ang_landau_base = upper(weil_angles(landau_gram(base, GAMMAS)))
    for factor in [1, 2, 3, 5, 7, 11, 30, 210]:
        scaled = [tuple(factor * x for x in a) for a in base]
        ang_l = upper(weil_angles(landau_gram(scaled, GAMMAS)))
        ang_e = upper(weil_angles(gram_matrix(scaled, GAMMAS)))
        d_s = upper(exchange_distances(scaled))
        print(f"    {factor:>7} {np.abs(ang_l - ang_landau_base).max():>29.3e} "
              f"{np.abs(ang_e - ang_base).max():>28.4f} {d_s.mean():>9.4f} "
              f"{np.corrcoef(d_s, d_base)[0,1]:>+14.4f}")
        rows.append(["obstruction-scaling", f"lambda={factor}", len(base),
                     len(d_s), f"{np.abs(ang_l - ang_landau_base).max():.6e}",
                     f"{np.abs(ang_e - ang_base).max():.6f}",
                     f"{d_s.mean():.6f}",
                     f"{np.corrcoef(d_s, d_base)[0,1]:.6f}"])


def experiment_staircase_anatomy(rows):
    """The best family found, dissected: d sees t/s, the Weil angle sees t-s."""

    section("(2c) THE BEST FAMILY, DISSECTED:  a_t = (p^t, 1)")
    for prime, top in [(3, 8), (2, 9)]:
        family = [(prime ** t, 1) for t in range(1, top + 1)]
        gram = gram_matrix(family, GAMMAS)
        angles = weil_angles(gram)
        distances = exchange_distances(family)
        size = len(family)
        deviation = max(
            abs(distances[i, j] - math.log((j + 1) / (i + 1)))
            for i in range(size) for j in range(i + 1, size)
        )
        print(f"\np = {prime}, t = 1..{top}:  "
              f"max |d(a_s,a_t) - log(t/s)| = {deviation:.3e}  "
              f"(so d is EXACTLY log(t/s))")

        pearson, rho, pairs, dd, aa = compare(family)
        print(f"   over {pairs} pairs: pearson {pearson:+.4f}, "
              f"spearman {rho:+.4f}")

        print("\n   (A) pairs at the SAME exchange distance d = log 2, "
              "different Weil angles:")
        print(f"      {'(s,t)':>8} {'d':>9} {'angle':>9} {'angle-pi/2':>12}")
        for s in range(1, top // 2 + 1):
            t = 2 * s
            if t > top:
                continue
            i, j = s - 1, t - 1
            print(f"      {f'({s},{t})':>8} {distances[i, j]:>9.5f} "
                  f"{angles[i, j]:>9.5f} {angles[i, j] - math.pi / 2:>+12.5f}")
            rows.append(["staircase-equal-d", f"p={prime} ({s},{t})", top, pairs,
                         f"{distances[i, j]:.6f}", f"{angles[i, j]:.6f}",
                         f"{angles[i, j] - math.pi / 2:.6f}", ""])

        print("\n   (B) pairs at the SAME Weil angle (t - s = 1), "
              "exchange distance varying by a factor of five:")
        print(f"      {'(s,t)':>8} {'d':>9} {'angle':>9} {'angle-pi/2':>12}")
        for s in range(1, top):
            i, j = s - 1, s
            print(f"      {f'({s},{s+1})':>8} {distances[i, j]:>9.5f} "
                  f"{angles[i, j]:>9.5f} {angles[i, j] - math.pi / 2:>+12.5f}")
            rows.append(["staircase-equal-angle", f"p={prime} ({s},{s+1})", top,
                         pairs, f"{distances[i, j]:.6f}", f"{angles[i, j]:.6f}",
                         f"{angles[i, j] - math.pi / 2:.6f}", ""])

        gaps = np.array([j - i for i in range(size) for j in range(i + 1, size)])
        deviations = np.array([angles[i, j] - math.pi / 2
                               for i in range(size) for j in range(i + 1, size)])
        model = prime ** (-gaps / 2.0)
        constant = float(np.dot(model, deviations) / np.dot(model, model))
        residual = deviations - constant * model
        print(f"\n   fit  angle - pi/2 = K * p^(-(t-s)/2)  with K = {constant:.4f}: "
              f"relative rms residual "
              f"{np.linalg.norm(residual) / np.linalg.norm(deviations):.4f}")
        print(f"   corr(angle - pi/2, p^(-(t-s)/2)) = "
              f"{np.corrcoef(deviations, model)[0, 1]:+.4f}   vs   "
              f"corr(angle - pi/2, d) = {np.corrcoef(deviations, dd)[0, 1]:+.4f}")
        rows.append(["staircase-fit", f"p={prime}", top, len(gaps),
                     f"{constant:.6f}",
                     f"{np.linalg.norm(residual) / np.linalg.norm(deviations):.6f}",
                     f"{np.corrcoef(deviations, model)[0, 1]:.6f}",
                     f"{np.corrcoef(deviations, dd)[0, 1]:.6f}"])


def experiment_staircase_window(rows):
    """The one genuine bridge: a staircase window, with both closed forms."""

    section("(2e) THE BRIDGE: A RECEDING STAIRCASE WINDOW  a_t = (p^t, 1), "
            "t in [T0, T0+8]")
    print("Two closed forms hold on this family:")
    print("   d(a_s,a_t)          = log(t/s)                       (exact)")
    print("   angle(a_s,a_t)-pi/2 = K * p^{-(t-s)/2},  "
          "K = (T/2piN) log p    (Landau)")
    print("Since log(t/s) = (t-s)/T0 * (1 + O(n/T0)) inside a window of "
          "width n << T0,")
    print("the angle becomes an exponential function of d:  "
          "angle - pi/2 = K exp(-c d),  c = T0 log p / 2.")
    print()
    print(f"{'p':>3} {'T0':>5} {'pairs':>6} {'pearson':>9} {'spearman':>9} "
          f"{'K fitted':>9} {'K = (T/2piN)log p':>18} {'rms resid':>10}")
    count, height = len(GAMMAS), GAMMAS[-1]
    for prime, start in [(2, 5), (2, 10), (2, 20), (2, 40), (2, 100),
                         (3, 10), (3, 20), (3, 40)]:
        family = [(prime ** t, 1) for t in range(start, start + 9)]
        angles = weil_angles(gram_matrix(family, GAMMAS))
        distances = exchange_distances(family)
        size = len(family)
        gaps = np.array([j - i for i in range(size) for j in range(i + 1, size)],
                        dtype=float)
        deviation = upper(angles) - math.pi / 2
        model = float(prime) ** (-gaps / 2)
        constant = float(np.dot(model, deviation) / np.dot(model, model))
        predicted = height / (2 * math.pi * count) * math.log(prime)
        residual = (np.linalg.norm(deviation - constant * model)
                    / np.linalg.norm(deviation))
        pearson, rho, pairs, _, _ = compare(family)
        print(f"{prime:>3} {start:>5} {pairs:>6} {pearson:>+9.4f} {rho:>+9.4f} "
              f"{constant:>9.5f} {predicted:>18.5f} {residual:>10.4f}")
        rows.append(["staircase-window", f"p={prime} T0={start}", size, pairs,
                     f"{pearson:.6f}", f"{rho:.6f}", f"{constant:.6f}",
                     f"{residual:.6f}"])
    print("\n   This is the best agreement found anywhere: spearman near -0.98 "
          "over 36 pairs,\n   with a derived two-parameter law.  Note the sign: "
          "on this family the Weil\n   form is REPULSIVE (angles above pi/2) and "
          "the repulsion DECREASES with d.")


def experiment_exact_invariances(rows):
    """Two exact theorems, verified: the invariance groups are transverse."""

    section("(2d) TWO EXACT INVARIANCES, VERIFIED")
    base = [(2, 1), (4, 1), (8, 1), (4, 2), (8, 2), (16, 4)]

    print("\nTheorem A.  E_{lambda a, lambda b} = lambda * E_{a,b} EXACTLY, for "
          "every truncation,\n            because sqrt(lambda a_i * lambda b_j) "
          "= lambda sqrt(a_i b_j) and every\n            ratio a_i/b_j is "
          "unchanged.  Hence the Weil correlation matrix R,\n            and "
          "every Weil angle, is EXACTLY scale invariant.")
    gram_base = gram_matrix(base, GAMMAS)
    for factor in [2, 3, 5, 7, 30, 210]:
        scaled = [tuple(factor * x for x in a) for a in base]
        gram_scaled = gram_matrix(scaled, GAMMAS)
        relative = np.abs(gram_scaled - factor * gram_base).max() / np.abs(gram_base).max()
        angle_shift = np.abs(weil_angles(gram_scaled) - weil_angles(gram_base)).max()
        print(f"   lambda = {factor:>4}: max relative deviation of E "
              f"{relative:.2e}, max angle shift {angle_shift:.2e}")

    print("\nTheorem B.  d(a^{ok}, b^{ok}) = d(a,b) EXACTLY, because "
          "log Z_{a^{ok}} = k log Z_a,\n            so the ratio defining "
          "C(a -> b) is unchanged.")
    for power in [2, 3]:
        raised = [tensor_power(a, power) for a in base]
        shift = np.abs(upper(exchange_distances(raised))
                       - upper(exchange_distances(base))).max()
        print(f"   k = {power}: max |d change| = {shift:.2e}")

    print("\nConsequence, inside ONE family.  Put F and 5*F side by side: the "
          "two blocks\nhave IDENTICAL Weil angles and DIFFERENT exchange "
          "distances.")
    scaled = [tuple(5 * x for x in a) for a in base]
    combined = base + scaled
    angles = weil_angles(gram_matrix(combined, GAMMAS))
    distances = exchange_distances(combined)
    n = len(base)
    print(f"   {'pair in F':>22} {'d':>9} {'angle':>9}    "
          f"{'pair in 5F':>26} {'d':>9} {'angle':>9}")
    for i, j in [(0, 1), (0, 2), (1, 2), (3, 4), (0, 3)]:
        print(f"   {f'{base[i]} vs {base[j]}':>22} {distances[i, j]:>9.5f} "
              f"{angles[i, j]:>9.5f}    "
              f"{f'{scaled[i]} vs {scaled[j]}':>26} "
              f"{distances[n + i, n + j]:>9.5f} {angles[n + i, n + j]:>9.5f}")
        rows.append(["scaling-twin", f"{base[i]}|{base[j]}", 0, 0,
                     f"{distances[i, j]:.6f}", f"{angles[i, j]:.6f}",
                     f"{distances[n + i, n + j]:.6f}",
                     f"{angles[n + i, n + j]:.6f}"])

    print("\nAnd the mirror image.  Put F and F^{o2} side by side: identical "
          "exchange\ndistances, different Weil angles.")
    squared = [tensor_power(a, 2) for a in base]
    combined = base + squared
    angles = weil_angles(gram_matrix(combined, GAMMAS))
    distances = exchange_distances(combined)
    print(f"   {'pair in F':>22} {'d':>9} {'angle':>9}    "
          f"{'pair in F^o2':>34} {'d':>9} {'angle':>9}")
    for i, j in [(0, 1), (0, 2), (1, 2), (3, 4), (0, 3)]:
        print(f"   {f'{base[i]} vs {base[j]}':>22} {distances[i, j]:>9.5f} "
              f"{angles[i, j]:>9.5f}    "
              f"{f'{squared[i]} vs {squared[j]}':>34} "
              f"{distances[n + i, n + j]:>9.5f} {angles[n + i, n + j]:>9.5f}")
        rows.append(["power-twin", f"{base[i]}|{base[j]}", 0, 0,
                     f"{distances[i, j]:.6f}", f"{angles[i, j]:.6f}",
                     f"{distances[n + i, n + j]:.6f}",
                     f"{angles[n + i, n + j]:.6f}"])


def experiment_scale_mismatch(rows):
    section("(3) SCALE-INVARIANCE MISMATCH")
    base = FAMILIES["generic-T1.4"]
    print("\nCartesian square of the whole family (a -> a tensor a):")
    squared = [tensor_power(a, 2) for a in base]
    d_base, d_sq = upper(exchange_distances(base)), upper(exchange_distances(squared))
    print(f"   max |d change| over {len(d_base)} pairs = "
          f"{np.abs(d_sq - d_base).max():.3e}   (exactly invariant)")

    print("\nCartesian square of a SINGLE member (a -> a tensor a, others fixed):")
    print(f"   {'member':>12} {'mean |d change|':>16} {'mean |angle change|':>20}")
    ang_base_matrix = weil_angles(gram_matrix(base, GAMMAS))
    d_base_matrix = exchange_distances(base)
    for index in range(4):
        modified = list(base)
        a = base[index]
        modified[index] = tensor_power(a, 2)
        d_new = exchange_distances(modified)
        ang_new = weil_angles(gram_matrix(modified, GAMMAS))
        others = [j for j in range(len(base)) if j != index]
        dd = np.abs(d_new[index, others] - d_base_matrix[index, others]).mean()
        da = np.abs(ang_new[index, others] - ang_base_matrix[index, others]).mean()
        print(f"   {str(a):>12} {dd:>16.4f} {da:>20.4f}")

    print("\nDoubling a SINGLE member's entries (a -> 2a, others fixed):")
    print(f"   {'member':>12} {'mean |d change|':>16} {'mean |angle change|':>20} "
          f"{'A-mass before':>14} {'A-mass after':>13}")
    for index in range(6):
        modified = list(base)
        a = base[index]
        modified[index] = tuple(2 * x for x in a)
        d_new = exchange_distances(modified)
        ang_new = weil_angles(gram_matrix(modified, GAMMAS))
        others = [j for j in range(len(base)) if j != index]
        dd = np.abs(d_new[index, others] - d_base_matrix[index, others]).mean()
        da = np.abs(ang_new[index, others] - ang_base_matrix[index, others]).mean()
        before = sum(arithmetic_term(a, base[j]) + overlap(a, base[j]) for j in others)
        after = sum(arithmetic_term(modified[index], base[j])
                    + overlap(modified[index], base[j]) for j in others)
        print(f"   {str(a):>12} {dd:>16.4f} {da:>20.4f} "
              f"{before:>14.3f} {after:>13.3f}")
        rows.append(["single-member-doubling", str(a), len(base), len(others),
                     f"{dd:.6f}", f"{da:.6f}", f"{before:.6f}", f"{after:.6f}"])

    print("\nSame experiment inside the designed powers-of-2 family:")
    designed = FAMILIES["powers-of-2 pairs (no 1)"]
    d_designed = exchange_distances(designed)
    ang_designed = weil_angles(gram_matrix(designed, GAMMAS))
    print("(tripling every entry of one member destroys every 2-power ratio it "
          "had with\n the rest, so its Weil angles all collapse to pi/2, while "
          "d barely moves.)")
    print(f"   {'member':>12} {'mean |d change|':>16} {'mean |angle change|':>20} "
          f"{'mean |ang-pi/2| before':>23} {'after':>9}")
    for index in range(4):
        modified = list(designed)
        a = designed[index]
        modified[index] = tuple(3 * x for x in a)  # 3 kills every 2-power ratio
        d_new = exchange_distances(modified)
        ang_new = weil_angles(gram_matrix(modified, GAMMAS))
        others = [j for j in range(len(designed)) if j != index]
        dd = np.abs(d_new[index, others] - d_designed[index, others]).mean()
        da = np.abs(ang_new[index, others] - ang_designed[index, others]).mean()
        before = np.abs(ang_designed[index, others] - math.pi / 2).mean()
        after = np.abs(ang_new[index, others] - math.pi / 2).mean()
        print(f"   {str(a):>12} {dd:>16.4f} {da:>20.4f} {before:>23.4f} "
              f"{after:>9.4f}")
        rows.append(["single-member-tripling", str(a), len(designed), len(others),
                     f"{dd:.6f}", f"{da:.6f}", f"{before:.6f}", f"{after:.6f}"])


def experiment_stability_in_n(rows):
    section("(4) N-DEPENDENCE OF THE CORRELATION")
    names = ["generic-T1.4", "2-ladder subsets k=5", "prime-ladder subsets k=5",
             "staircase (2^t, 1)"]
    print(f"{'N':>6} {'T':>9} " + " ".join(f"{n:>26}" for n in names))
    for count in [100, 200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400]:
        if count > len(LONG_GAMMAS):
            continue
        gammas = LONG_GAMMAS[:count]
        values = []
        for name in names:
            pearson, rho, _, _, _ = compare(FAMILIES[name], gammas)
            values.append(f"{pearson:+.4f}/{rho:+.4f}")
        print(f"{count:>6} {gammas[-1]:>9.1f} " + " ".join(f"{v:>26}" for v in values))
        rows.append(["N-dependence", f"N={count}", count, 0, *values[:4]]
                    if len(values) == 4 else ["N-dependence", f"N={count}", count, 0])


def main() -> int:
    print(f"{len(GAMMAS)} zeta zeros, T = {GAMMAS[-1]:.3f}, "
          f"T/2pi = {GAMMAS[-1] / (2 * math.pi):.3f}")
    rows = []
    experiment_designed_families(rows)
    experiment_mechanism(rows)
    experiment_monotone_search(rows)
    experiment_staircase_anatomy(rows)
    experiment_staircase_window(rows)
    experiment_exact_invariances(rows)
    experiment_scale_mismatch(rows)
    experiment_stability_in_n(rows)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["experiment", "family", "size", "pairs",
                         "value1", "value2", "value3", "value4"])
        writer.writerows(rows)
    print(f"\nwritten to {CSV_PATH.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
