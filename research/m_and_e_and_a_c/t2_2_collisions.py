#!/usr/bin/env python3
"""T2.2 support: signature collisions among random maps, and a direct fit of

    1 - C(f -> L) = kappa m2 / (2 q log q) + gamma m3 / (q^{3/2} log q) + O(q^-2),
    kappa = 3 - 2 sqrt 2,  gamma = C(beta*, 3) / (beta* + 1),  beta* = sqrt 2 - 1.

The ratio of the two coefficients is ``2 gamma / (kappa sqrt q) = 0.528558 / sqrt q``:
that single number is how much less the rate weighs the third trace moment than
the second.

Run:  python research/m_and_e_and_a_c/t2_2_collisions.py [q ...]
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import BETA_STAR, KAPPA, moments, rate, references, sig_from_counts  # noqa: E402
import ffmaps as F  # noqa: E402

GAMMA = (BETA_STAR * (BETA_STAR - 1) * (BETA_STAR - 2) / 6) / (BETA_STAR + 1)


def key_of(counts: np.ndarray) -> tuple:
    s = sig_from_counts(counts)
    return (tuple(int(v) for v in s.values), tuple(int(m) for m in s.mults))


def collisions(q: int, rng: np.random.Generator, n: int = 400) -> None:
    print(f"\nq = {q}: distinct signatures among {n} random maps per family")
    fams = {
        "y^2 = P_3(x) + c  (elliptic)": lambda: F.hyperelliptic(
            q, [int(v) for v in rng.integers(0, q, size=3)] + [1]),
        "y^2 = P_5(x) + c  (genus 2)": lambda: F.hyperelliptic(
            q, [int(v) for v in rng.integers(0, q, size=5)] + [1]),
        "y^2 = P_7(x) + c  (genus 3)": lambda: F.hyperelliptic(
            q, [int(v) for v in rng.integers(0, q, size=7)] + [1]),
        "dense bidegree (2,2)": lambda: F.bilinear_family(q, rng.integers(0, q, size=(3, 3))),
        "dense bidegree (3,3)": lambda: F.bilinear_family(q, rng.integers(0, q, size=(4, 4))),
        "dense bidegree (5,5)": lambda: F.bilinear_family(q, rng.integers(0, q, size=(6, 6))),
        "bilinear a+bx+cy+dxy": lambda: F.bilinear_family(q, rng.integers(0, q, size=(2, 2))),
        "P_4(x) + Q_4(y)": lambda: F.additive(
            q, [int(v) for v in rng.integers(0, q, size=4)] + [1],
            [int(v) for v in rng.integers(0, q, size=4)] + [1]),
    }
    for name, gen in fams.items():
        keys: dict[tuple, int] = {}
        for _ in range(n):
            k = key_of(gen())
            keys[k] = keys.get(k, 0) + 1
        top = max(keys.values())
        print(f"  {name:30s} {len(keys):5d} distinct   (largest class {top})")
    d = math.gcd(4, q - 1)
    print("  [elliptic] y^2 = x^3 + a x + t and y^2 = x^3 + a u^-4 x + t u^-6 have the")
    print("   same fiber multiset (x -> u^2 x, y -> u^3 y, t -> t u^6 permutes the base),")
    print(f"   so the signature only depends on a mod (F_q^*)^4: gcd(4, q-1) = {d} classes.")


def model_fit(q: int, rng: np.random.Generator, n: int = 600) -> None:
    L = references(q)["L"]
    rows = []
    for deg in (3, 5, 7, 9, 11):
        for _ in range(n // 5):
            c = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
            counts = F.hyperelliptic(q, c)
            if (counts > 0).sum() != q:
                continue
            m = moments(counts, q)
            if m["a_absmax"] > 6 * math.sqrt(q) or m["m2"] < 1e-9:
                continue
            v, _ = rate(sig_from_counts(counts), L)
            rows.append((m["m2"], m["m3"], 1.0 - v))
    m2 = np.array([r[0] for r in rows])
    m3 = np.array([r[1] for r in rows])
    y = np.array([r[2] for r in rows])
    coef, *_ = np.linalg.lstsq(np.vstack([m2, m3]).T, y, rcond=None)
    pred_k = KAPPA / (2 * q * math.log(q))
    pred_g = GAMMA / (q ** 1.5 * math.log(q))
    print(f"\nq = {q}: fit of 1 - C(f->L) on (m2, m3), {len(rows)} Weil-regime maps")
    print(f"  m2 coefficient: fitted {coef[0]:.6e}  predicted {pred_k:.6e}  "
          f"ratio {coef[0] / pred_k:.6f}")
    print(f"  m3 coefficient: fitted {coef[1]:.6e}  predicted {pred_g:.6e}  "
          f"ratio {coef[1] / pred_g:.6f}")
    print(f"  rms residual, m2 only     : {np.sqrt(((y - m2 * pred_k) ** 2).mean()):.3e}")
    print(f"  rms residual, m2 and m3   : "
          f"{np.sqrt(((y - m2 * pred_k - m3 * pred_g) ** 2).mean()):.3e}")
    print(f"  spread of 1 - C(f->L)     : {y.min():.3e} .. {y.max():.3e}")
    print(f"  m3/m2 sensitivity ratio   : {pred_g / pred_k:.6f} "
          f"(= 0.528558 / sqrt q = {0.528558 / math.sqrt(q):.6f})")


def main(qs: list[int]) -> None:
    rng = np.random.default_rng(31337)
    print(f"kappa = {KAPPA:.9f}   gamma = {GAMMA:.9f}   "
          f"2 gamma / kappa = {2 * GAMMA / KAPPA:.9f}")
    for q in qs:
        collisions(q, rng)
        model_fit(q, rng)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [101, 211, 503])
