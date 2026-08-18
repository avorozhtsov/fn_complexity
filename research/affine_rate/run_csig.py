#!/usr/bin/env python3
"""Signature rate ``C_sig`` on the 14 F_3 map classes, at 40 digits.

``C_sig(g -> f) = inf_beta log Z_g(beta) / log Z_f(beta)``.  Every class here
has fiber sizes summing to 9, so ``Z(1) = 9`` for all of them and the ratio is
exactly 1 at ``beta = 1``: hence ``C_sig <= 1`` throughout this pool, and the
infimum is attained at an endpoint whenever the ratio is monotone on each side
of 1.  The values are computed with ``mpmath`` at 60 working digits on a
logarithmic grid plus the two endpoints, then refined by golden section.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from mpmath import mp, mpf, log, exp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from affine_atoms import enumerate_classes  # noqa: E402

mp.dps = 60
OUTPUT = Path(__file__).resolve().parent / "csig_matrix.json"


def log_partition(signature, beta):
    if beta == mp.inf:
        raise ValueError
    return log(sum(exp(beta * log(mpf(n))) for n in signature))


def ratio(source, target, beta):
    numerator = log_partition(source, beta)
    denominator = log_partition(target, beta)
    return numerator / denominator


def endpoint_zero(source, target):
    return log(mpf(len(source))) / log(mpf(len(target)))


def endpoint_infinity(source, target):
    return log(mpf(max(source))) / log(mpf(max(target)))


def csig(source, target):
    """Infimum over ``beta in [0, inf]`` of ``log Z_source / log Z_target``."""

    if len(target) == 1:            # constant target: produced for free
        return mp.inf
    if len(source) == 1:            # constant resource: log Z_source(0) = 0
        return mpf(0)
    candidates = [endpoint_zero(source, target), endpoint_infinity(source, target)]
    grid = [mpf(10) ** (mpf(e) / 24) for e in range(-24 * 4, 24 * 4 + 1)]
    values = [(ratio(source, target, b), b) for b in grid if b != 1]
    values.append((mpf(1), mpf(1)))
    values.sort()
    best_value, best_beta = values[0]
    candidates.append(best_value)
    # golden-section refinement around the best interior grid point
    index = [b for _, b in sorted(values, key=lambda t: t[1])]
    position = index.index(best_beta)
    low = index[max(position - 1, 0)]
    high = index[min(position + 1, len(index) - 1)]
    if low < high:
        phi = (mpf(5) ** mpf("0.5") - 1) / 2
        for _ in range(300):
            c = high - phi * (high - low)
            d = low + phi * (high - low)
            if ratio(source, target, c) < ratio(source, target, d):
                high = d
            else:
                low = c
        candidates.append(ratio(source, target, (low + high) / 2))
    return min(candidates)


def main() -> None:
    classes = [c for c in enumerate_classes() if c.key != "constant"]
    keys = [c.key for c in classes]
    signatures = {c.key: c.signature for c in classes}

    table = {}
    for g in keys:
        row = {}
        for f in keys:
            value = csig(signatures[g], signatures[f])
            row[f] = "inf" if value == mp.inf else mp.nstr(value, 40)
        table[g] = row

    OUTPUT.write_text(json.dumps({
        "signatures": {k: list(v) for k, v in signatures.items()},
        "csig": table,
    }, indent=1))

    print("C_sig(g -> f), rows g:")
    print(f"{'':<14}" + "".join(f"{k[:11]:>12}" for k in keys))
    for g in keys:
        cells = []
        for f in keys:
            v = table[g][f]
            cells.append("inf" if v == "inf" else f"{float(v):.6f}")
        print(f"{g:<14}" + "".join(f"{c:>12}" for c in cells))


if __name__ == "__main__":
    main()
