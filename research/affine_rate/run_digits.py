#!/usr/bin/env python3
"""Forty-digit values of the headline ``C_sig`` numbers.

Every class in the pool has fiber sizes summing to 9, so ``Z(1) = 9`` and the
ratio ``log Z_g / log Z_f`` equals 1 at ``beta = 1`` for every pair: the infimum
is therefore at most 1 and is either at ``beta = 0`` (ratio ``log #fibers_g /
log #fibers_f``), at ``beta = infinity`` (ratio ``log max_g / log max_f``), or at
an interior contact.  Interior contacts are located by ternary search at 80
working digits and cross-checked against a 4001-point grid.
"""

from __future__ import annotations

from mpmath import mp, mpf, log, exp, nstr

mp.dps = 80

PAIRS = (
    ("(3,3,3)", (3, 3, 3), "(6,3)", (6, 3)),
    ("(6,3)", (6, 3), "(3,3,3)", (3, 3, 3)),
    ("(6,3)", (6, 3), "(5,2,2)", (5, 2, 2)),
    ("(6,3)", (6, 3), "(4,4,1)", (4, 4, 1)),
    ("(5,2,2)", (5, 2, 2), "(6,3)", (6, 3)),
    ("(4,4,1)", (4, 4, 1), "(6,3)", (6, 3)),
    ("(4,4,1)", (4, 4, 1), "(3,3,3)", (3, 3, 3)),
    ("(5,2,2)", (5, 2, 2), "(3,3,3)", (3, 3, 3)),
    ("(7,1,1)", (7, 1, 1), "(6,3)", (6, 3)),
    ("(7,1,1)", (7, 1, 1), "(3,3,3)", (3, 3, 3)),
)


def ratio(source, target, beta):
    numerator = log(sum(exp(beta * log(mpf(n))) for n in source))
    denominator = log(sum(exp(beta * log(mpf(n))) for n in target))
    return numerator / denominator


def infimum(source, target):
    zero = log(mpf(len(source))) / log(mpf(len(target)))
    infinity = log(mpf(max(source))) / log(mpf(max(target)))
    grid = [mpf(10) ** (mpf(e) / 500) for e in range(-2000, 2001)]
    values = [(ratio(source, target, b), b) for b in grid]
    values.sort(key=lambda pair: pair[0])
    best, beta = values[0]
    position = sorted(grid).index(beta)
    low = sorted(grid)[max(position - 1, 0)]
    high = sorted(grid)[min(position + 1, len(grid) - 1)]
    for _ in range(600):
        c = low + (high - low) / 3
        d = high - (high - low) / 3
        if ratio(source, target, c) <= ratio(source, target, d):
            high = d
        else:
            low = c
    interior = ratio(source, target, (low + high) / 2)
    contact = (low + high) / 2
    candidates = [(zero, mpf(0), "beta = 0"), (infinity, mp.inf, "beta = inf"),
                  (interior, contact, "interior")]
    candidates.sort(key=lambda triple: triple[0])
    return candidates[0]


def main() -> None:
    for source_name, source, target_name, target in PAIRS:
        value, contact, where = infimum(source, target)
        location = where if where != "interior" else f"beta = {nstr(contact, 12)}"
        print(f"C_sig({source_name} -> {target_name}) = {nstr(value, 40)}")
        print(f"    contact: {location}")
    print()
    print(f"log 3 / log 6 = {nstr(log(mpf(3)) / log(mpf(6)), 40)}")
    print(f"log 2 / log 3 = {nstr(log(mpf(2)) / log(mpf(3)), 40)}")
    print(f"log 4 / log 6 = {nstr(log(mpf(4)) / log(mpf(6)), 40)}")


if __name__ == "__main__":
    main()
