#!/usr/bin/env python3
"""T2.2 part (3), sharpest form: pairs on which the *whole* rate vector nearly
collides while the trace distribution genuinely differs.

Searches hyperelliptic maps ``f_P(x, y) = y^2 - P(x)`` with all fibers non-empty,
groups them by the two endpoint readouts ``(|image|, N_max)``, and reports the
pair minimising the sup-norm distance between the full 8-probe rate vectors
among pairs whose ``(m2, m3)`` differ.

The winner has identical image size, identical largest fiber and identical m2,
but a *different smallest fiber*: the rate curve reads ``max_c N_c`` exactly
(beta = infinity) and has no endpoint at all for ``min_c N_c``, so the positive
extreme trace only leaks in through the ``q^{-1/2}``-suppressed higher moments.

Run:  python research/m_and_e_and_a_c/t2_2_vector.py [q] [n_per_degree]
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import (  # noqa: E402
    BETA_STAR, certified_rate, moments, rate, references, sig_from_counts,
)
import ffmaps as F  # noqa: E402


def binom(b: float, k: int) -> float:
    v = 1.0
    for i in range(k):
        v *= (b - i)
    return v / math.factorial(k)


def coefficient(q: int, k: int) -> float:
    """``c_k`` in ``1 - C(f->L) = sum_{k>=2} c_k m_k + O(q^-2)``."""
    return ((-1) ** (k + 1)) * binom(BETA_STAR, k) / (
        (BETA_STAR + 1) * math.log(q) * q ** (k / 2))


def main(q: int, n_each: int) -> None:
    rng = np.random.default_rng(555)
    refs = references(q)
    probes = [("C(f->L)", refs["L"], 0), ("C(L->f)", refs["L"], 1),
              ("C(f->Xsplit)", refs["Xsplit"], 0), ("C(Xsplit->f)", refs["Xsplit"], 1),
              ("C(f->Xaniso)", refs["Xaniso"], 0), ("C(Xaniso->f)", refs["Xaniso"], 1),
              ("C(f->Sq)", refs["Sq"], 0), ("C(Sq->f)", refs["Sq"], 1)]
    pool = []
    for deg in (3, 5, 7, 9):
        for _ in range(n_each):
            c = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
            counts = F.hyperelliptic(q, c)
            if (counts > 0).sum() != q:
                continue
            s = sig_from_counts(counts)
            pool.append({"P": tuple(c), "m": moments(counts, q), "s": s,
                         "r": [rate(s, g)[0] if d == 0 else rate(g, s)[0]
                               for _, g, d in probes]})
    print(f"q = {q}, pool {len(pool)} full-image hyperelliptic maps")

    grp: dict[tuple[int, int], list[dict]] = {}
    for r in pool:
        grp.setdefault((r["m"]["n_image"], r["m"]["max_fiber"]), []).append(r)
    best = None
    n_pairs = 0
    for g in grp.values():
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                a, b = g[i], g[j]
                if (abs(a["m"]["m2"] - b["m"]["m2"]) < 1e-12
                        and abs(a["m"]["m3"] - b["m"]["m3"]) < 1e-9):
                    continue
                n_pairs += 1
                d = max(abs(x - y) for x, y in zip(a["r"], b["r"]))
                if best is None or d < best[0]:
                    best = (d, a, b)
    print(f"pairs sharing (|image|, N_max) with different (m2, m3): {n_pairs}")
    d, a, b = best
    print(f"\ntightest full-vector collision: sup over the 8 probes |dC| = {d:.3e}")
    for tag, r in (("f", a), ("g", b)):
        m = r["m"]
        print(f"  {tag}: P = {r['P']}")
        print(f"     m2 = {m['m2']:.12f}  m3 = {m['m3']:+.9f}  m4 = {m['m4']:.6f}  "
              f"N_max = {m['max_fiber']}  N_min = {q - int(m['amax'])}  "
              f"|image| = {m['n_image']}")
    dm = {k: b["m"][k] - a["m"][k] for k in ("m2", "m3", "m4", "m5")}
    print(f"  dm2 = {dm['m2']:+.4e}  dm3 = {dm['m3']:+.4e}  dm4 = {dm['m4']:+.4e}"
          f"  d(min fiber) = {int(a['m']['amax'] - b['m']['amax'])}")
    for (name, _, _), x, y in zip(probes, a["r"], b["r"]):
        print(f"  {name:14s} {x:.15f} {y:.15f}  |d| = {abs(x - y):.3e}")
    pred = sum(coefficient(q, k) * dm[f"m{k}"] for k in (2, 3, 4, 5))
    print(f"\n  predicted dC(f->L) from the moment ladder = {abs(pred):.4e}")
    print(f"  observed                                  = "
          f"{abs(a['r'][0] - b['r'][0]):.4e}")
    print("  repo-solver certification of C(f->L):")
    for tag, r in (("f", a), ("g", b)):
        v, bt = certified_rate(r["s"], refs["L"])
        print(f"    exact C({tag}->L) = {v:.15f}   beta = {bt:.9f}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 211,
         int(sys.argv[2]) if len(sys.argv) > 2 else 400)
