#!/usr/bin/env python3
"""Do the F_3 tensor three-cycles of ``analysis/tensor_cycles_f3.py`` survive C_aff?

The cycles are cycles of *signatures*.  ``C_aff`` is not a function of the
signature, so the question only makes sense after choosing orbits.  For each
cycle we list every triple of orbits realising the three signatures and read off
the one-step linear implementation order between them.

The point of the exercise is the following consequence of the reduction lemma:
for maps with equal output dimension ``m`` one has ``k(r) <= r`` and therefore

    f <= g   (one step)   ==>   C_aff(g -> f) = 1 ,

while ``C_aff(f -> g) < 1`` whenever ``g`` is not reachable from ``f``.  So along
every comparable pair the ``C_aff`` comparison is forced and agrees with the
degeneration order, which is a partial order and admits no cycle.
"""

from __future__ import annotations

import csv
import itertools
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fn_complexity.homogeneous_tensor_maps import (  # noqa: E402
    TENSOR_CASES, _evaluate_form, compute_tensor_poset,
)

CASE = 3
CYCLE_CSV = PROJECT_ROOT / "analysis" / "tensor_cycles_f3.csv"


def orbit_signature(case, representative) -> tuple[int, ...]:
    counts: Counter = Counter()
    for point in itertools.product(range(3), repeat=case.input_dimension):
        value = tuple(_evaluate_form(form, point, case.degree)
                      for form in representative)
        counts[value] += 1
    return tuple(sorted(counts.values(), reverse=True))


def reachability(poset) -> dict[str, set[str]]:
    """``below[a]`` = orbits implementable from ``a`` (transitive closure)."""

    keys = [orbit.key for orbit in poset.orbits]
    direct = {key: set() for key in keys}
    for source, target in poset.covers:
        direct[source].add(target)
    below = {key: set() for key in keys}

    def visit(key: str) -> set[str]:
        if below[key]:
            return below[key]
        result = set()
        for child in direct[key]:
            result.add(child)
            result |= visit(child)
        below[key] = result
        return result

    for key in keys:
        visit(key)
    return below


def parse_signature(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.strip("{}").split(","))


def main() -> None:
    case = TENSOR_CASES[CASE]
    poset = compute_tensor_poset(CASE)
    below = reachability(poset)
    signatures = {orbit.key: orbit_signature(case, orbit.representative)
                  for orbit in poset.orbits}
    by_signature: dict[tuple[int, ...], list[str]] = {}
    for key, sig in signatures.items():
        by_signature.setdefault(sig, []).append(key)

    cycles = []
    with CYCLE_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if int(row["case"]) != CASE:
                continue
            if row["cycle"] not in [c[0] for c in cycles]:
                nodes = tuple(parse_signature(part)
                              for part in row["cycle"].split(" -> "))
                cycles.append((row["cycle"], nodes))

    print(f"case {CASE}: {len(poset.orbits)} orbits, "
          f"{len(by_signature)} distinct signatures")
    print(f"{len(cycles)} distinct signature three-cycles\n")

    total, refuted, undetermined = 0, 0, 0
    for label, nodes in cycles:
        print(f"cycle {label}")
        print("   C_sig orientation:  "
              + " < ".join("{" + ",".join(map(str, s)) + "}" for s in nodes)
              + " < (first)")
        for sig in nodes:
            print(f"   {sig}  realised by {len(by_signature[sig])} orbit(s): "
                  f"{sorted(by_signature[sig])}")
        for triple in itertools.product(*(sorted(by_signature[s]) for s in nodes)):
            if len(set(triple)) < 3:
                continue
            total += 1
            # the C_sig cycle needs triple[i] < triple[i+1] for every i
            broken = []
            facts = []
            for index in range(3):
                a, b = triple[index], triple[(index + 1) % 3]
                # needed: a < b, i.e. C_aff(a->b) < C_aff(b->a)
                if b in below[a]:
                    # a implements b, so C_aff(a->b) = 1 >= C_aff(b->a): NOT a < b
                    facts.append(f"{b} <= {a}")
                    broken.append(f"{a} < {b}")
                elif a in below[b]:
                    facts.append(f"{a} <= {b}")
            if broken:
                refuted += 1
                verdict = ("REFUTED: C_aff cannot have "
                           + ", ".join(broken))
            else:
                undetermined += 1
                verdict = "undetermined"
            detail = ("  [" + ",  ".join(facts) + "]") if facts else ""
            print(f"   triple {triple}: {verdict}{detail}")
        print()

    print(f"orbit triples examined: {total}")
    print(f"  C_sig cycle refuted by the affine order: {refuted}")
    print(f"  undetermined at this level: {undetermined}")


if __name__ == "__main__":
    main()
