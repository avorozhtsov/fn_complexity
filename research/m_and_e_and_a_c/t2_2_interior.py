#!/usr/bin/env python3
"""T2.2 part (3): what the interior of the rate curve adds to the endpoints.

Two things are checked.

A.  *Tangency reading.*  If ``C(f -> g)`` is attained at an interior ``beta0``
    then ``log Z_f(beta0) = C * log Z_g(beta0)`` exactly, so one rate is one
    exact evaluation of the Mellin-type transform ``beta -> log Z_f(beta)``.
    Since ``log Z_f`` is convex, the rates against *all* flat references
    ``flat(n, m)`` reconstruct ``log Z_f`` on ``[0, inf)`` as an upper envelope
    of lines, hence determine the whole multiset ``{N_c}``.  A finite reference
    family gives only finitely many tangent lines.

B.  *Where the tangency sits.*  For any signature whose fibers are all ``q(1 +
    O(q^{-1/2}))`` the tangency against ``L = (q, ..., q)`` is pinned at
    ``beta -> sqrt(2) - 1``, independently of f.  So the interior effectively
    contributes a single functional, and higher moments enter only through the
    ``q^{-1/2}``-suppressed corrections at that one point.

Run:  python research/m_and_e_and_a_c/t2_2_interior.py [q]
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import (  # noqa: E402
    BETA_STAR, KAPPA, flat, moments, rate, references, sig_from_counts,
)
import ffmaps as F  # noqa: E402


def tangency_reading(q: int) -> None:
    refs = references(q)
    counts = F.hyperelliptic(q, [0, 1, 0, 1])
    f = sig_from_counts(counts)
    print("A. one rate = one exact evaluation of log Z_f")
    print(f"   {'probe':14s} {'C':>18s} {'beta0':>12s} {'log Z_f(beta0) from C':>24s} "
          f"{'direct':>22s} {'err':>10s}")
    for name, g in refs.items():
        for role in ("impd", "impl"):
            c, b = rate(f, g) if role == "impd" else rate(g, f)
            if math.isinf(b) or b == 0.0:
                continue
            lzg = float(g.log_z(np.array([b]))[0])
            recovered = c * lzg if role == "impd" else lzg / c
            direct = float(f.log_z(np.array([b]))[0])
            tag = f"C(f->{name})" if role == "impd" else f"C({name}->f)"
            print(f"   {tag:14s} {c:18.15f} {b:12.6f} {recovered:24.15f} "
                  f"{direct:22.15f} {abs(recovered - direct):10.1e}")


def flat_envelope(q: int) -> None:
    """Reconstruct log Z_f from rates against flat references only."""
    counts = F.hyperelliptic(q, [0, 1, 0, 1])
    f = sig_from_counts(counts)
    print("\n   reconstruction of log Z_f from flat references flat(n, m):")
    print(f"   {'n':>7s} {'m':>7s} {'C(f->flat)':>18s} {'beta0':>10s} "
          f"{'reconstructed':>20s} {'true':>20s}")
    for n, m in ((q, q), (q, 2 * q), (q, 5 * q), (q * q, q), (3, q), (q ** 3, q)):
        g = flat(n, m)
        c, b = rate(f, g)
        if math.isinf(b):
            print(f"   {n:7d} {m:7d} {c:18.12f} {'inf':>10s} "
                  f"{'(endpoint: max fiber)':>20s}")
            continue
        if b == 0.0:
            print(f"   {n:7d} {m:7d} {c:18.12f} {'0':>10s} "
                  f"{'(endpoint: #fibers)':>20s}")
            continue
        lzg = math.log(n) + b * math.log(m)
        print(f"   {n:7d} {m:7d} {c:18.12f} {b:10.5f} {c * lzg:20.15f} "
              f"{float(f.log_z(np.array([b]))[0]):20.15f}")


def tangency_pinned(q: int, rng: np.random.Generator) -> None:
    L = references(q)["L"]
    print("\nB. the interior tangency against L is pinned at sqrt(2)-1")
    print(f"   {'family':14s} {'m2':>12s} {'beta0':>12s} {'beta0-(sqrt2-1)':>18s} "
          f"{'1-C(f->L)':>14s} {'kappa m2/(2q log q)':>21s}")
    fams = [("y^2=x^3+x+c", F.hyperelliptic(q, [0, 1, 0, 1])),
            ("y^2=x^3+c", F.hyperelliptic(q, [0, 0, 0, 1])),
            ("y^2=x^5+x+c", F.hyperelliptic(q, [0, 1, 0, 0, 0, 1])),
            ("y^2=x^7+x+c", F.hyperelliptic(q, [0, 1, 0, 0, 0, 0, 0, 1])),
            ("y^3=x^4+c", F.superelliptic(q, 3, [0, 0, 0, 0, 1])),
            ("y^2=x^9+x+c", F.hyperelliptic(q, [0, 1] + [0] * 7 + [1]))]
    for _ in range(4):
        c = [int(v) for v in rng.integers(0, q, size=11)] + [1]
        fams.append((f"y^2=P11(x)+c", F.hyperelliptic(q, c)))
    for name, counts in fams:
        if (counts > 0).sum() != q:
            name += " (*)"
        s = sig_from_counts(counts)
        m = moments(counts, q)
        c, b = rate(s, L)
        pred = KAPPA * m["m2"] / (2 * q * math.log(q))
        bs = "inf" if math.isinf(b) else f"{b:.9f}"
        db = "-" if math.isinf(b) else f"{b - BETA_STAR:+.9f}"
        print(f"   {name:14s} {m['m2']:12.6f} {bs:>12s} {db:>18s} "
              f"{1 - c:14.6e} {pred:21.6e}")
    print("   (*) marks maps whose image is not all of F_q")


def endpoint_blindness(q: int) -> None:
    """The rate curve reads the largest fiber exactly and the smallest not at all."""
    print("\nC. sensitivity of log Z_f(beta) to a single fiber: d/dN_c ~ N_c^(beta-1)")
    for b in (BETA_STAR, 1.0, 2.0, 8.0):
        big, small = q + int(2 * math.sqrt(q)), q - int(2 * math.sqrt(q))
        print(f"   beta = {b:6.3f}:  weight(largest fiber)/weight(smallest fiber) = "
              f"{(big / small) ** (b - 1):.6f}")
    print("   beta < 1 weights small fibers *more* per unit change, but no beta < 0 is")
    print("   available, so min_c N_c (= q - max_c a_c) is never isolated by an endpoint;")
    print("   only max_c N_c (= q - min_c a_c) is, at beta = infinity.")


def main(q: int) -> None:
    rng = np.random.default_rng(7)
    print(f"q = {q}\n")
    tangency_reading(q)
    flat_envelope(q)
    tangency_pinned(q, rng)
    endpoint_blindness(q)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 211)
