#!/usr/bin/env python3
"""Cycle search among families of curves.

Three parts.

**Part 1, exhaustive.**  Over ``F_11`` and ``F_13`` every genus-two pencil
``y^2 = P(x) + c`` with ``P`` monic of degree 5 or 6 is enumerated (the constant
term of ``P`` may be set to zero, since shifting it only permutes the fibers).
The exchange comparison is computed on every distinct fiber signature and every
strict three-cycle is found.  This is a complete search, not a sample.

**Part 2, sampled.**  At larger ``q`` a pool mixing hyperelliptic, superelliptic,
quadratic-twist, additive and dense families is sampled, and the search is run
inside each class of equal largest fiber.

**Part 3.**  The observed cycle margins against the predicted scale
``1/(sqrt(q) log q)``.

Every cycle reported here is re-verified pairwise against the package's
``exchange_rate_result``; the widest one is proved by interval arithmetic in
``certify.py``.

    python research/curve_family_cycles/search.py
"""

from __future__ import annotations

import collections
import csv
import itertools
import json
import math
from pathlib import Path
import random
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from common import Engine, beta_grid, build_pool  # noqa: E402
from fn_complexity import exchange_rate_result  # noqa: E402

GRID_POINTS = 30_000
TOLERANCE = 1e-9
EXHAUSTIVE_Q = (11, 13)
SAMPLED_Q = (31, 101, 211)


# ------------------------------------------------------- genus-two pencils


def quadratic_character_table(q: int) -> np.ndarray:
    squares = {(i * i) % q for i in range(1, q)}
    chi = np.array([0 if u == 0 else (1 if u in squares else -1) for u in range(q)])
    return np.array([[chi[(u + c) % q] for c in range(q)] for u in range(q)], dtype=np.int64)


def all_pencils(q: int, degree: int) -> tuple[np.ndarray, np.ndarray]:
    """Signatures of ``y^2 = P(x) + c`` for every monic ``P`` with ``P(0) = 0``.

    ``#f^{-1}(c) = sum_x (1 + chi(P(x) + c)) = q + sum_u n_u chi(u + c)`` with
    ``n_u = #{x : P(x) = u}``, so one matrix product covers the whole degree.
    """

    exponents = np.arange(1, degree + 1)
    powers = np.array([[pow(x, int(e), q) for e in exponents] for x in range(q)], dtype=np.int64)
    free = np.array(list(itertools.product(range(q), repeat=degree - 1)), dtype=np.int64)
    coefficients = np.concatenate([free, np.ones((free.shape[0], 1), dtype=np.int64)], axis=1)
    values = (coefficients @ powers.T) % q
    counts = np.zeros(values.shape, dtype=np.int64)
    np.add.at(
        counts,
        (np.repeat(np.arange(values.shape[0]), q), values.ravel()),
        1,
    )
    fibers = q + counts @ quadratic_character_table(q)
    return np.sort(fibers, axis=1)[:, ::-1], coefficients


def _poly_mod(a: list[int], q: int) -> list[int]:
    while a and a[-1] % q == 0:
        a.pop()
    return [v % q for v in a]


def _poly_gcd(a: list[int], b: list[int], q: int) -> list[int]:
    a, b = _poly_mod(list(a), q), _poly_mod(list(b), q)
    while b:
        inverse = pow(b[-1], q - 2, q)
        shift = len(a) - len(b)
        while shift >= 0 and a:
            factor = (a[-1] * inverse) % q
            for index, coefficient in enumerate(b):
                a[index + shift] = (a[index + shift] - factor * coefficient) % q
            a = _poly_mod(a, q)
            shift = len(a) - len(b)
        a, b = b, a
    return a


def smooth_fibers(coefficients: np.ndarray, q: int) -> int:
    """Number of ``c`` for which ``P(x) + c`` is squarefree over the closure.

    ``P + c`` is squarefree exactly when it is coprime to ``P'``, and ``P'`` does
    not depend on ``c``.
    """

    p = [0] + [int(v) for v in coefficients]
    derivative = [(index * value) % q for index, value in enumerate(p)][1:]
    total = 0
    for c in range(q):
        shifted = list(p)
        shifted[0] = (shifted[0] + c) % q
        if len(_poly_gcd(shifted, list(derivative), q)) <= 1:
            total += 1
    return total


def exhaustive_pool(q: int) -> tuple[dict[tuple[int, ...], dict], list]:
    """Distinct signatures with a witness each, plus the raw enumeration.

    The smoothness of a pencil costs ``q`` polynomial gcds, far too much to pay
    on every one of the ``q^5`` polynomials, so the first witness found is kept
    here and upgraded on demand by ``best_witness``.
    """

    pool: dict[tuple[int, ...], dict] = {}
    tables = []
    for degree in (5, 6):
        fibers, coefficients = all_pencils(q, degree)
        tables.append((degree, fibers, coefficients))
        for row in range(fibers.shape[0]):
            if fibers[row, -1] == 0:
                continue
            signature = tuple(int(v) for v in fibers[row])
            if signature not in pool:
                pool[signature] = {
                    "degree": degree,
                    "coefficients": [int(v) for v in coefficients[row]],
                    "smooth": None,
                }
    return pool, tables


def best_witness(signature: tuple[int, ...], tables: list, q: int, cap: int = 400) -> dict:
    """Return the realising pencil with the most smooth fibers, capped search."""

    target = np.array(signature, dtype=np.int64)
    best = None
    for degree, fibers, coefficients in tables:
        matches = np.flatnonzero((fibers == target).all(axis=1))
        for row in matches[:cap]:
            smooth = smooth_fibers(coefficients[row], q)
            if best is None or smooth > best["smooth"]:
                best = {
                    "degree": degree,
                    "coefficients": [int(v) for v in coefficients[row]],
                    "smooth": smooth,
                }
            if smooth == q:
                return best
    return best


def polynomial_text(coefficients: list[int]) -> str:
    terms = []
    for index, value in enumerate(coefficients, start=1):
        if not value:
            continue
        power = "x" if index == 1 else f"x^{index}"
        terms.append(power if value == 1 else f"{value}{power}")
    return " + ".join(reversed(terms)) or "0"


# ------------------------------------------------------------ comparisons


def comparison(signatures: list[tuple[int, ...]], q: int) -> tuple[np.ndarray, np.ndarray]:
    engine = Engine(signatures, beta_grid(q, points=GRID_POINTS))
    return engine.rate_matrix(chunk=64)


def three_cycles(rates: np.ndarray) -> list[tuple[float, int, int, int]]:
    difference = rates - rates.T
    strict = difference < -TOLERANCE
    size = rates.shape[0]
    successors = [np.flatnonzero(strict[index]) for index in range(size)]
    found = []
    for first in range(size):
        for second in successors[first]:
            if second <= first:
                continue
            for third in successors[second]:
                if third <= first:
                    continue
                if strict[third, first]:
                    margin = min(
                        -difference[first, second],
                        -difference[second, third],
                        -difference[third, first],
                    )
                    found.append((float(margin), first, second, third))
    found.sort(reverse=True)
    return found


def cross_check(signatures: list[tuple[int, ...]], rates: np.ndarray, count: int = 300) -> float:
    """Largest disagreement between the grid engine and the package's solver."""

    rng = np.random.default_rng(7)
    worst = 0.0
    size = len(signatures)
    for _ in range(count):
        i, j = (int(v) for v in rng.integers(0, size, size=2))
        if i == j:
            continue
        exact = exchange_rate_result(
            implementer=signatures[i], implemented=signatures[j], grid_size=16384
        ).rate
        worst = max(worst, abs(exact - rates[i, j]))
    return worst


def family_kind(entry) -> str:
    """The family a pool member came from, with the degree parameters dropped."""

    return entry.family.split(" deg")[0].split(" (")[0]


def edge_status(a: tuple[int, ...], b: tuple[int, ...]) -> str:
    """How the computed edge ``a`` before ``b`` stands against ``phi``."""

    gap = math.log(len(a)) * math.log(max(a)) - math.log(len(b)) * math.log(max(b))
    if abs(gap) <= 1e-12:
        return "blind"
    return "consistent" if gap < 0 else "violating"


def contact(result, implementer: tuple[int, ...], implemented: tuple[int, ...]) -> float:
    """The contact temperature, with a flat approach to ``beta = infinity`` fixed.

    On these signatures the ratio can approach its ``beta = infinity`` limit from
    above so slowly that it is constant to double precision from ``beta`` of a
    few hundred on, and the golden-section refinement then reports a spurious
    interior minimiser.  Whenever the rate agrees with the endpoint value to
    machine precision the contact is at infinity.
    """

    endpoint = math.log(max(implementer)) / math.log(max(implemented))
    if abs(result.rate - endpoint) < 1e-13:
        return math.inf
    return result.beta


def verify_cycle(cycle: list[tuple[int, ...]]) -> list[dict]:
    """Recompute every rate of a cycle with the package, independently."""

    out = []
    for index in range(3):
        a, b = cycle[index], cycle[(index + 1) % 3]
        forward = exchange_rate_result(implementer=a, implemented=b, grid_size=16384)
        backward = exchange_rate_result(implementer=b, implemented=a, grid_size=16384)
        out.append(
            {
                "forward": forward.rate,
                "forward_beta": contact(forward, a, b),
                "backward": backward.rate,
                "backward_beta": contact(backward, b, a),
                "margin": backward.rate - forward.rate,
                "phi": edge_status(a, b),
            }
        )
    return out


# ------------------------------------------------------------------ parts


def part_one(rows: list, headline: dict) -> None:
    print("\n" + "=" * 78)
    print("Part 1 --- exhaustive over genus-two pencils y^2 = P(x) + c")
    print("=" * 78)
    for q in EXHAUSTIVE_Q:
        pool, tables = exhaustive_pool(q)
        signatures = sorted(pool)
        rates, betas = comparison(signatures, q)
        difference = rates - rates.T
        strict = difference < -TOLERANCE
        phi = np.array([math.log(len(s)) * math.log(max(s)) for s in signatures])
        gap = phi[:, None] - phi[None, :]
        blind = strict & (np.abs(gap) <= 1e-12)
        violating = strict & (gap > 1e-12)
        cycles = three_cycles(rates)
        interior = np.isfinite(betas) & (betas > 0)
        print(
            f"\nq = {q}: {len(signatures)} distinct signatures, "
            f"{int(strict.sum())} strict ordered pairs, "
            f"{int((np.abs(difference) <= TOLERANCE).sum() - len(signatures)) // 2} tied pairs"
        )
        print(
            f"    phi-blind strict pairs {int(blind.sum())}, "
            f"phi-violating ordered pairs {int(violating.sum())}, "
            f"rates with interior contact {100 * interior.mean():.1f}%"
        )
        print(
            f"    grid engine against the package solver on 300 random pairs: "
            f"max deviation {cross_check(signatures, rates):.2e}"
        )
        # Inside a phi-class both endpoints give the ratio 1 exactly, so a rate
        # there is either 1 (endpoint) or strictly below it (interior).  Check
        # that the first case never occurs: every comparison inside a class is
        # then decided in the interior, which is the claim of the note.
        tops = np.array([max(s) for s in signatures])
        together = (tops[:, None] == tops[None, :]) & ~np.eye(len(signatures), dtype=bool)
        print(
            f"    pairs inside a phi-class {int(together.sum())}: "
            f"rates equal to 1 {int(((np.abs(rates - 1.0) < 1e-12) & together).sum())}, "
            f"contacts not interior {int((together & ~interior).sum())}"
        )
        print(f"    strict 3-cycles: {len(cycles)}")
        if not cycles:
            continue
        pattern = collections.Counter()
        for _, i, j, k in cycles:
            triple = (signatures[i], signatures[j], signatures[k])
            pattern[
                tuple(
                    sorted(
                        edge_status(triple[t], triple[(t + 1) % 3]) for t in range(3)
                    )
                )
            ] += 1
        for key, count in pattern.most_common():
            print(f"      edges {key}: {count}")

        involved = {index for entry in cycles for index in entry[1:]}
        for index in involved:
            pool[signatures[index]] = best_witness(signatures[index], tables, q)
        smooth_cycles = [
            entry
            for entry in cycles
            if all(pool[signatures[t]]["smooth"] == q for t in entry[1:])
        ]
        print(
            f"      cycles whose three pencils have every fiber smooth: "
            f"{len(smooth_cycles)}"
        )
        for label, chosen in (("widest", cycles[:1]), ("widest all-smooth", smooth_cycles[:1])):
            if not chosen:
                continue
            margin, i, j, k = chosen[0]
            triple = [signatures[i], signatures[j], signatures[k]]
            checked = verify_cycle(triple)
            print(f"\n    {label} cycle, minimum margin {margin:.6e}:")
            for position, index in enumerate((i, j, k)):
                record = pool[signatures[index]]
                name = "ABC"[position]
                print(
                    f"      {name}:  y^2 = {polynomial_text(record['coefficients'])} + c "
                    f"(deg {record['degree']}, {record['smooth']}/{q} fibers smooth)"
                )
                print(
                    f"           sigma = {signatures[index]}   "
                    f"phi = {phi[index]:.6f}   "
                    f"m2 = {sum((q - n) ** 2 for n in signatures[index]) / q ** 2:.6f}"
                )
            for position, entry in enumerate(checked):
                source, target = "ABC"[position], "ABC"[(position + 1) % 3]
                print(
                    f"      {source} < {target}: "
                    f"C({source}->{target}) = {entry['forward']:.12f} "
                    f"[beta = {entry['forward_beta']:.6g}]   "
                    f"C({target}->{source}) = {entry['backward']:.12f} "
                    f"[beta = {entry['backward_beta']:.6g}]   "
                    f"margin {entry['margin']:.4e}   phi-{entry['phi']}"
                )
                rows.append(
                    [
                        f"F_{q} genus-2 pencils ({label})",
                        source,
                        target,
                        "{" + ",".join(map(str, triple[position])) + "}",
                        "{" + ",".join(map(str, triple[(position + 1) % 3])) + "}",
                        f"{entry['forward']:.15f}",
                        f"{entry['forward_beta']:.6g}",
                        f"{entry['backward']:.15f}",
                        f"{entry['backward_beta']:.6g}",
                        f"{entry['margin']:.15f}",
                        entry["phi"],
                    ]
                )
            if label == "widest all-smooth" and q == 11:
                # the cycle certify.py proves, recorded so the two cannot drift
                headline.update(
                    q=q,
                    signatures=[list(s) for s in triple],
                    pencils=[
                        f"y^2 = {polynomial_text(pool[s]['coefficients'])} + c" for s in triple
                    ],
                    degrees=[pool[s]["degree"] for s in triple],
                    smooth_fibers=[pool[s]["smooth"] for s in triple],
                    margins=[entry["margin"] for entry in checked],
                    contacts=[
                        [
                            "inf" if math.isinf(entry["forward_beta"]) else entry["forward_beta"],
                            "inf" if math.isinf(entry["backward_beta"]) else entry["backward_beta"],
                        ]
                        for entry in checked
                    ],
                    phi_status=[entry["phi"] for entry in checked],
                )


def part_two(rows: list, scaling: list) -> None:
    print("\n" + "=" * 78)
    print("Part 2 --- sampled pools of mixed curve families at larger q")
    print("=" * 78)
    for q in SAMPLED_Q:
        pool = build_pool(q)
        classes = collections.defaultdict(list)
        for entry in pool:
            classes[entry.max_fiber].append(entry)
        ranked = sorted(classes.items(), key=lambda item: -len(item[1]))[:4]
        print(f"\nq = {q}: {len(pool)} distinct signatures over {len(classes)} phi-classes")
        best = (0.0,)
        mixed = (0.0,)
        shuffler = random.Random(3)
        for max_fiber, members in ranked:
            # The pool is generated family by family, so taking a prefix would
            # fill every class with hyperelliptic pencils alone.
            members = list(members)
            shuffler.shuffle(members)
            members = members[:200]
            if len(members) < 20:
                continue
            signatures = [entry.signature for entry in members]
            rates, _ = comparison(signatures, q)
            cycles = three_cycles(rates)
            widest = cycles[0][0] if cycles else 0.0
            print(
                f"    max fiber {max_fiber}: {len(signatures)} signatures, "
                f"{len(cycles)} strict 3-cycles, widest margin {widest:.3e}"
            )
            if cycles and widest > best[0]:
                triple = [signatures[t] for t in cycles[0][1:]]
                best = (widest, max_fiber, triple, [members[t] for t in cycles[0][1:]])
            for margin, *indices in cycles:
                kinds = {family_kind(members[t]) for t in indices}
                if len(kinds) == 3 and margin > mixed[0]:
                    mixed = (
                        margin,
                        max_fiber,
                        [signatures[t] for t in indices],
                        [members[t] for t in indices],
                    )
        if len(best) > 1:
            scaling.append((q, best[0]))
            checked = verify_cycle(best[2])
            print(f"    widest cycle at q = {q}, margin {best[0]:.6e}, max fiber {best[1]}:")
            for entry, member in zip(checked, best[3]):
                print(f"      {member.family}: {member.witness[:64]}")
        if len(mixed) > 1:
            print(
                f"    widest cycle spanning three different family types, "
                f"margin {mixed[0]:.6e}, max fiber {mixed[1]}:"
            )
            for member in mixed[3]:
                print(f"      {member.family}: {member.witness[:64]}")
            for position, entry in enumerate(verify_cycle(mixed[2])):
                rows.append(
                    [
                        f"F_{q} sampled pool (three family types)",
                        "ABC"[position],
                        "ABC"[(position + 1) % 3],
                        "{" + ",".join(map(str, mixed[2][position])) + "}",
                        "{" + ",".join(map(str, mixed[2][(position + 1) % 3])) + "}",
                        f"{entry['forward']:.15f}",
                        f"{entry['forward_beta']:.6g}",
                        f"{entry['backward']:.15f}",
                        f"{entry['backward_beta']:.6g}",
                        f"{entry['margin']:.15f}",
                        entry["phi"],
                    ]
                )
            for position, entry in enumerate(checked):
                rows.append(
                    [
                        f"F_{q} sampled pool",
                        "ABC"[position],
                        "ABC"[(position + 1) % 3],
                        "{" + ",".join(map(str, best[2][position])) + "}",
                        "{" + ",".join(map(str, best[2][(position + 1) % 3])) + "}",
                        f"{entry['forward']:.15f}",
                        f"{entry['forward_beta']:.6g}",
                        f"{entry['backward']:.15f}",
                        f"{entry['backward_beta']:.6g}",
                        f"{entry['margin']:.15f}",
                        entry["phi"],
                    ]
                )


def part_three(scaling: list) -> None:
    print("\n" + "=" * 78)
    print("Part 3 --- margins against the predicted scale 1/(sqrt(q) log q)")
    print("=" * 78)
    print(f"\n    {'q':>5} {'widest margin':>16} {'1/(sqrt q log q)':>18} {'ratio':>10}")
    for q, margin in scaling:
        predicted = 1.0 / (math.sqrt(q) * math.log(q))
        print(f"    {q:>5} {margin:>16.6e} {predicted:>18.6e} {margin / predicted:>10.4f}")


def main() -> int:
    rows: list = []
    headline: dict = {}
    scaling: list = []
    part_one(rows, headline)
    part_two(rows, scaling)
    part_three(scaling)

    with (HERE / "cycles.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "family",
                "source",
                "target",
                "sigma_source",
                "sigma_target",
                "rate_forward",
                "beta_forward",
                "rate_backward",
                "beta_backward",
                "margin",
                "phi_status",
            ]
        )
        writer.writerows(rows)
    (HERE / "headline_cycle.json").write_text(json.dumps(headline, indent=2), encoding="utf-8")
    print("\nwritten: cycles.csv, headline_cycle.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
