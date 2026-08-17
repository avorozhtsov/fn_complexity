#!/usr/bin/env python3
"""T2.3 part 3: a systematic search for exact entries of the exchange matrix.

An entry ``C(g -> f)`` is *exact* if it equals a recognisable closed-form number
(a small rational, or ``log a / log b`` for small integers) rather than merely
converging to one.  The search runs over a pool of structured maps
``f : A^2 -> A^1`` and over the primes ``q <= 500``, and asks three questions.

1.  **Degenerate signatures.**  Which pool members have a signature with only
    one or two distinct fiber sizes, and under what condition on ``q``?
2.  **Collisions.**  Which pairs of pool members have *equal* signatures --
    equivalently ``C(f -> g) = C(g -> f) = 1`` exactly -- and for which ``q``?
    Every collision set is then matched against a congruence ``q = r mod m``.
3.  **Other exact values.**  A full ordered scan of the pairwise rates, testing
    each against small rationals ``p/r`` and against ``log a / log b`` for
    ``a, b <= 24``, at tolerance ``1e-12``.

The headline criterion the search keeps rediscovering is proved in the note:
for ``f = y^2 - P(x)`` the signature is flat exactly when ``P`` is a permutation
polynomial of ``F_q``.  That is checked here directly over the whole pool.

Run:  python research/m_and_e_and_a_c/t2_3_detector_search.py [q_max] [q_pair_max]
"""

from __future__ import annotations

import csv
from fractions import Fraction
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import rate, sig_from_counts  # noqa: E402
import ffmaps as F  # noqa: E402
import t2_3_common as T  # noqa: E402

CSV_COLLIDE = HERE / "t2_3_collisions.csv"
CSV_EXACT = HERE / "t2_3_exact_values.csv"
CSV_DEGEN = HERE / "t2_3_degenerate.csv"


# ------------------------------------------------------------------- the pool


def build_pool(q: int) -> list[tuple[str, np.ndarray, tuple | None]]:
    """``(label, fiber counts, branch polynomial or None)``."""
    out: list[tuple[str, np.ndarray, tuple | None]] = []
    out.append(("L = x", F.flat_map(q), None))
    out.append(("const = 0", T.constant_map(q), None))
    out.append(("Xsplit = xy", F.split_conic(q), None))
    out.append(("Xaniso = x^2-d y^2", F.aniso_conic(q), None))
    for d in range(2, 13):
        P = tuple([0] * d + [1])
        out.append((f"y^2 - x^{d}", F.hyperelliptic(q, list(P)), P))
    for n in range(2, 9):
        P = tuple(c % q for c in T.dickson(n))
        out.append((f"y^2 - D_{n}(x)", F.hyperelliptic(q, list(P)), P))
    for d, a in ((3, 1), (3, 2), (5, 1), (5, 2), (7, 1)):
        P = tuple([0, a] + [0] * (d - 2) + [1])
        out.append((f"y^2 - (x^{d}+{a}x)", F.hyperelliptic(q, list(P)), P))
    for r in (3, 4, 5):
        for d in (2, 3, 4, 5, 6):
            out.append((f"y^{r} - x^{d}", F.superelliptic(q, r, [0] * d + [1]), None))
    for a in (2, 3, 4, 5):
        for b in (2, 3, 4, 5):
            if a <= b:
                out.append((f"x^{a} + y^{b}",
                            F.additive(q, [0] * a + [1], [0] * b + [1]), None))
    out.append(("x^2 (pushforward)", F.bilinear_family(q, _mono(2)), None))
    out.append(("x^3 (pushforward)", F.bilinear_family(q, _mono(3)), None))
    out.append(("x y^2 - x^2 (quartic twists)", T.quartic_twist_map(q), None))
    for label, P in (("x^3+x", [0, 1, 0, 1]), ("x^3+1", [1, 0, 0, 1]),
                     ("x^3+x+1", [1, 1, 0, 1])):
        out.append((f"({label}) y^2 (quadratic twists)",
                    T.quadratic_twist_map(q, P), None))
    return out


def _mono(dx: int) -> np.ndarray:
    A = np.zeros((dx + 1, 1), dtype=np.int64)
    A[dx, 0] = 1
    return A


# ---------------------------------------------------------------- exact values


def _targets() -> list[tuple[float, str]]:
    out: dict[float, str] = {}
    for r in range(1, 13):
        for p in range(0, 3 * r + 1):
            v = Fraction(p, r)
            if 0.0 <= float(v) <= 3.0:
                out.setdefault(float(v), f"{v.numerator}/{v.denominator}")
    for a in range(2, 25):
        for b in range(2, 25):
            v = math.log(a) / math.log(b)
            if 0.0 < v <= 3.0:
                out.setdefault(v, f"log {a}/log {b}")
    return sorted(out.items())


TARGETS = _targets()
TOL = 1e-12


def _nearest(value: float) -> tuple[float, str]:
    lo, hi = 0, len(TARGETS)
    while lo < hi:
        mid = (lo + hi) // 2
        if TARGETS[mid][0] < value:
            lo = mid + 1
        else:
            hi = mid
    best = min(TARGETS[max(lo - 2, 0):lo + 2], key=lambda t: abs(t[0] - value),
               default=(float("inf"), ""))
    return best


def recognise_loose(value: float, tol: float = 1e-5) -> str | None:
    v, name = _nearest(value)
    return name if abs(v - value) <= tol else None


def recognise(value: float) -> str | None:
    lo, hi = 0, len(TARGETS)
    while lo < hi:
        mid = (lo + hi) // 2
        if TARGETS[mid][0] < value - TOL:
            lo = mid + 1
        else:
            hi = mid
    for v, name in TARGETS[lo:lo + 3]:
        if abs(v - value) <= TOL:
            return name
    return None


# --------------------------------------------------------------------- passes


def classify(primes: list[int]):
    """Per prime, bucket the pool by signature; track class membership."""
    labels: list[str] | None = None
    degen: dict[str, dict[int, int]] = {}
    member: dict[str, dict[int, str]] = {}          # map -> q -> class name
    other_groups: dict[tuple[str, ...], dict[int, bool]] = {}
    perm_check = {"agree": 0, "disagree": 0}
    for q in primes:
        pool = build_pool(q)
        if labels is None:
            labels = [p[0] for p in pool]
        named = {
            T.sig_key(*_sig(F.flat_map(q))): "L",
            T.sig_key(*_sig(F.split_conic(q))): "Xsplit",
            T.sig_key(*_sig(F.aniso_conic(q))): "Xaniso",
            T.sig_key(*_sig(T.constant_map(q))): "const",
        }
        buckets: dict[tuple, list[str]] = {}
        for label, counts, P in pool:
            s_ = sig_from_counts(counts)
            key = T.sig_key(s_.values, s_.mults)
            degen.setdefault(label, {})[q] = len(s_.values)
            member.setdefault(label, {})[q] = named.get(key, "-")
            buckets.setdefault(key, []).append(label)
            if P is not None:
                flat = named.get(key) == "L"
                if flat == T.is_permutation_polynomial(list(P), q):
                    perm_check["agree"] += 1
                else:
                    perm_check["disagree"] += 1
        for key, members in buckets.items():
            if len(members) >= 2 and key not in named:
                other_groups.setdefault(tuple(sorted(members)), {})[q] = True
    for grp in other_groups:
        for q in primes:
            other_groups[grp].setdefault(q, False)
    return labels, degen, member, other_groups, perm_check


def _sig(counts):
    s = sig_from_counts(counts)
    return s.values, s.mults


def _condition(primes, flag) -> str:
    n = sum(flag[q] for q in primes)
    if n == 0:
        return "never"
    if n == len(primes):
        return "always"
    cls = T.residue_classes(primes, flag)
    return f"q = {cls[1]} mod {cls[0]}" if cls else f"not a congruence (n={n})"


def report_degenerate(primes, degen, rows):
    print("1. signatures with one or two distinct fiber sizes")
    print(f"   {'map':30s} {'#q flat':>8} {'condition (flat)':32s} "
          f"{'#q 2-valued':>12} {'condition (2-valued)':32s}")
    for label, per_q in degen.items():
        flat = {q: per_q[q] == 1 for q in primes}
        two = {q: per_q[q] == 2 for q in primes}
        nf, nt = sum(flat.values()), sum(two.values())
        if nf == 0 and nt == 0:
            continue
        sf, st = _condition(primes, flat), _condition(primes, two)
        print(f"   {label:30s} {nf:>8} {sf:32s} {nt:>12} {st:32s}")
        rows.append([label, nf, sf, nt, st])


def report_collisions(primes, member, other_groups, rows):
    print("\n2. exact signature coincidences with the named references")
    print("   (each is C = 1 in both directions between the map and the reference)")
    print(f"   {'map':30s} {'= L (flat)':34s} {'= Xsplit':34s} {'= Xaniso':30s}")
    for label, per_q in member.items():
        conds = []
        for name in ("L", "Xsplit", "Xaniso"):
            flag = {q: per_q[q] == name for q in primes}
            conds.append(_condition(primes, flag))
        if all(c == "never" for c in conds):
            continue
        print(f"   {label:30s} {conds[0]:34s} {conds[1]:34s} {conds[2]:30s}")
        rows.append([label, "L", conds[0]])
        rows.append([label, "Xsplit", conds[1]])
        rows.append([label, "Xaniso", conds[2]])

    print("\n   coincidence classes NOT equal to a named reference "
          "(genuinely new shared signatures)")
    shown = 0
    for grp, flag in sorted(other_groups.items(), key=lambda kv: -sum(kv[1].values())):
        cond = _condition(primes, flag)
        if cond == "never":
            continue
        shown += 1
        if shown > 25:
            continue
        print(f"      {cond:34s} {' == '.join(grp)}")
        rows.append(["|".join(grp), "shared", cond])
    print(f"   {shown} such classes (first 25 shown)")


def report_exact_rates(primes, rows):
    print("\n3. full ordered scan of pairwise rates against closed-form targets")
    found: dict[tuple[str, str, str], list[int]] = {}
    n_pairs = 0
    for q in primes:
        pool = build_pool(q)
        sigs = [(lbl, sig_from_counts(c)) for lbl, c, _ in pool]
        for i, (la, sa) in enumerate(sigs):
            for j, (lb, sb) in enumerate(sigs):
                if i == j:
                    continue
                n_pairs += 1
                # screen on a coarse grid, then re-solve at full resolution
                v, _ = rate(sa, sb, grid=192)
                if recognise_loose(v) is None:
                    continue
                v, _ = rate(sa, sb)
                name = recognise(v)
                if name is not None:
                    found.setdefault((la, lb, name), []).append(q)
    print(f"   scanned {n_pairs} ordered rates over q in "
          f"[{primes[0]}, {primes[-1]}]  (tolerance {TOL:g})")
    non_one = {k: v for k, v in found.items() if k[2] != "1/1"}
    print(f"   hits with value 1 (signature collisions or one-sided): "
          f"{len(found) - len(non_one)} distinct (ordered pair, value) combinations")
    single = sum(1 for qs in non_one.values() if len(qs) == 1)
    print(f"   hits with a value other than 1: {len(non_one)}, of which {single} occur at "
          f"a single small q (accidents of log a/log b with a, b <= 24 when q is tiny)")
    print(f"   persistent hits (at least two primes):")
    print(f"   {'implementer':32s} {'implemented':32s} {'value':14s} {'#q':>4} {'condition':30s}")
    seen_values = set()
    for (la, lb, name), qs in sorted(non_one.items()):
        flag = {q: q in set(qs) for q in primes}
        cls = T.residue_classes(primes, flag) if len(qs) < len(primes) else None
        cond = ("always" if len(qs) == len(primes) else
                (f"q = {cls[1]} mod {cls[0]}" if cls else "not a congruence"))
        rows.append([la, lb, name, len(qs), cond])
        if len(qs) < 2:
            continue
        seen_values.add(name)
        print(f"   {la:32s} {lb:32s} {name:14s} {len(qs):>4} {cond:30s}")
    print(f"   distinct persistent values found: {sorted(seen_values)}")


def main(argv: list[str]) -> int:
    q_max = int(argv[0]) if argv else 500
    q_pair_max = int(argv[1]) if len(argv) > 1 else 200
    primes = T.primes_upto(q_max, lo=11)
    labels, degen, member, other_groups, perm = classify(primes)
    drows: list[list] = []
    report_degenerate(primes, degen, drows)
    crows: list[list] = []
    report_collisions(primes, member, other_groups, crows)
    print(f"\n   permutation-polynomial criterion 'sig(y^2 - P) = L  <=>  P permutes F_q':"
          f"  {perm['agree']} agreements, {perm['disagree']} violations")
    erows: list[list] = []
    report_exact_rates([q for q in primes if q <= q_pair_max], erows)
    for path, header, rws in ((CSV_DEGEN, ["map", "n_q_flat", "flat_condition",
                                           "n_q_two_valued", "two_valued_condition"], drows),
                              (CSV_COLLIDE, ["map", "reference", "condition"], crows),
                              (CSV_EXACT, ["implementer", "implemented", "value",
                                           "n_q", "condition"], erows)):
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rws)
    print(f"\nwritten: {CSV_DEGEN.name}, {CSV_COLLIDE.name}, {CSV_EXACT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
