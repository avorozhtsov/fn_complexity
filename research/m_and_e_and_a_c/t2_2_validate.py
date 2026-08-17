#!/usr/bin/env python3
"""Certify the vectorised rate solver in t2_2_common against fn_complexity."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import certified_rate, rate, references, sig_from_counts  # noqa: E402
import ffmaps as F  # noqa: E402


def main() -> None:
    rng = np.random.default_rng(11)
    worst = 0.0
    worst_beta = 0.0
    n = 0
    t_fast = t_slow = 0.0
    for q in (101, 211):
        refs = references(q)
        pool = []
        for _ in range(15):
            c = [int(v) for v in rng.integers(0, q, size=5)] + [1]
            pool.append(F.hyperelliptic(q, c))
        for _ in range(15):
            A = rng.integers(0, q, size=(4, 4))
            pool.append(F.bilinear_family(q, A))
        for _ in range(10):
            pc = [int(v) for v in rng.integers(0, q, size=4)] + [1]
            qc = [int(v) for v in rng.integers(0, q, size=3)] + [1]
            pool.append(F.additive(q, pc, qc))
        for counts in pool:
            f = sig_from_counts(counts)
            for name, g in refs.items():
                for a, b in ((g, f), (f, g)):
                    t = time.perf_counter()
                    v1, b1 = rate(a, b)
                    t_fast += time.perf_counter() - t
                    t = time.perf_counter()
                    v2, b2 = certified_rate(a, b)
                    t_slow += time.perf_counter() - t
                    n += 1
                    worst = max(worst, abs(v1 - v2))
                    if not (math.isinf(b1) or math.isinf(b2)):
                        worst_beta = max(worst_beta, abs(b1 - b2))
                    elif math.isinf(b1) != math.isinf(b2):
                        print(f"BETA ENDPOINT MISMATCH {name} {b1} {b2}")
    print(f"{n} rate comparisons")
    print(f"max |C_fast - C_exact| = {worst:.3e}")
    print(f"max |beta_fast - beta_exact| (finite betas) = {worst_beta:.3e}")
    print(f"time: fast {t_fast:.2f}s, repo solver {t_slow:.2f}s "
          f"(speedup {t_slow / t_fast:.1f}x)")


if __name__ == "__main__":
    main()
