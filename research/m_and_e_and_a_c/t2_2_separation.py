#!/usr/bin/env python3
"""T2.2 separation: matched pairs of maps and what the rates do or do not see.

For ``f_P(x, y) = y^2 - P(x)`` over ``F_q`` (q odd prime) the fiber over c is the
affine curve ``y^2 = P(x) + c`` and ``a_c = -sum_x chi(P(x) + c)``.  Orthogonality
of the quadratic character gives the exact identity

    sum_c a_c^2 = q * K_P - q^2,      K_P = #{(x, x') : P(x) = P(x')},

so ``m2 = K_P/q - 1`` is an *exact integer-valued* invariant of P.  That makes it
easy to build pairs of maps with identical m2 (identical K_P) and different m3,
and pairs with identical image size and identical largest fiber but different m2.

Run:  python research/m_and_e_and_a_c/t2_2_separation.py [q]
"""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import (  # noqa: E402
    BETA_STAR, KAPPA, Sig, certified_rate, moments, rate, references, sig_from_counts,
)
import ffmaps as F  # noqa: E402


def k_invariant(coeffs, q: int) -> int:
    x = np.arange(q, dtype=np.int64)
    v = F.poly_eval(coeffs, x, q)
    return int((np.bincount(v, minlength=q).astype(np.int64) ** 2).sum())


def make_pool(q: int, deg: int, count: int, rng: np.random.Generator) -> list[dict]:
    pool = []
    for _ in range(count):
        c = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
        counts = F.hyperelliptic(q, c)
        mom = moments(counts, q)
        pool.append({"P": tuple(c), "counts": counts, "mom": mom,
                     "K": k_invariant(c, q), "sig": sig_from_counts(counts)})
    return pool


def probe_row(sig: Sig, refs: dict[str, Sig]) -> dict[str, tuple[float, float]]:
    out = {}
    out["C(f->L)"] = rate(sig, refs["L"])
    out["C(L->f)"] = rate(refs["L"], sig)
    out["C(f->Xsplit)"] = rate(sig, refs["Xsplit"])
    out["C(Xsplit->f)"] = rate(refs["Xsplit"], sig)
    out["C(f->Xaniso)"] = rate(sig, refs["Xaniso"])
    out["C(Xaniso->f)"] = rate(refs["Xaniso"], sig)
    return out


def show_pair(title: str, a: dict, b: dict, q: int, refs: dict[str, Sig]) -> None:
    print(f"\n--- {title} ---")
    for tag, r in (("f", a), ("g", b)):
        m = r["mom"]
        print(f"  {tag}: P={r['P']}  K={r['K']}  m2={m['m2']:.12f}  m3={m['m3']:+.9f}  "
              f"m4={m['m4']:.6f}  N_max={m['max_fiber']}  |image|={m['n_image']}")
    ra, rb = probe_row(a["sig"], refs), probe_row(b["sig"], refs)
    print(f"  {'probe':16s} {'f':>20s} {'g':>20s} {'|diff|':>12s}  {'beta_f':>10s}")
    for k in ra:
        va, vb = ra[k][0], rb[k][0]
        bt = ra[k][1]
        print(f"  {k:16s} {va:20.15f} {vb:20.15f} {abs(va - vb):12.3e}  "
              f"{('inf' if math.isinf(bt) else f'{bt:.6f}'):>10s}")
    return ra, rb


def main(q: int) -> None:
    rng = np.random.default_rng(4242)
    refs = references(q)
    print(f"q = {q}   kappa = {KAPPA:.9f}   beta* = {BETA_STAR:.9f}")

    pool = make_pool(q, 5, 900, rng) + make_pool(q, 7, 900, rng) + make_pool(q, 3, 600, rng)
    print(f"pool: {len(pool)} hyperelliptic maps y^2 = P(x) + c, deg P in {{3,5,7}}")
    sigs = {}
    for r in pool:
        key = (tuple(int(v) for v in r["sig"].values), tuple(int(m) for m in r["sig"].mults))
        sigs.setdefault(key, []).append(r)
    print(f"distinct signatures: {len(sigs)}")

    # ---- (2) same m2 exactly, different m3
    by_k: dict[int, list[dict]] = {}
    for r in pool:
        by_k.setdefault(r["K"], []).append(r)
    best = None
    for K, group in by_k.items():
        if len(group) < 2:
            continue
        group = sorted(group, key=lambda r: r["mom"]["m3"])
        cand = (group[0], group[-1])
        gap = cand[1]["mom"]["m3"] - cand[0]["mom"]["m3"]
        if best is None or gap > best[0]:
            best = (gap, cand)
    if best:
        show_pair("SAME m2 (identical K), maximal m3 gap", best[1][0], best[1][1], q, refs)

    # ---- (3) same |image| AND same N_max, different m2
    by_ext: dict[tuple[int, int], list[dict]] = {}
    for r in pool:
        by_ext.setdefault((r["mom"]["n_image"], r["mom"]["max_fiber"]), []).append(r)
    best = None
    for key, group in by_ext.items():
        if len(group) < 2:
            continue
        group = sorted(group, key=lambda r: r["mom"]["m2"])
        gap = group[-1]["mom"]["m2"] - group[0]["mom"]["m2"]
        if best is None or gap > best[0]:
            best = (gap, (group[0], group[-1]), key)
    if best:
        show_pair(f"SAME image size and SAME largest fiber {best[2]}, different m2",
                  best[1][0], best[1][1], q, refs)

    # ---- (3') same |image|, same N_max, same m2, different m3  (sharpest)
    by_ext3: dict[tuple[int, int, int], list[dict]] = {}
    for r in pool:
        by_ext3.setdefault(
            (r["mom"]["n_image"], r["mom"]["max_fiber"], r["K"]), []).append(r)
    best = None
    for key, group in by_ext3.items():
        if len(group) < 2:
            continue
        group = sorted(group, key=lambda r: r["mom"]["m3"])
        gap = group[-1]["mom"]["m3"] - group[0]["mom"]["m3"]
        if best is None or gap > best[0]:
            best = (gap, (group[0], group[-1]), key)
    if best:
        ra, rb = show_pair(
            "SHARPEST: same image size, same largest fiber, same m2, different m3",
            best[1][0], best[1][1], q, refs)
        print("  certifying with the repo solver:")
        for tag, r in (("f", best[1][0]), ("g", best[1][1])):
            v, b = certified_rate(r["sig"], refs["L"])
            print(f"    exact C({tag}->L) = {v:.15f}  beta = {b:.9f}")

    # ---- converse: different m2 but equal rates?
    print("\n--- converse: distinct signatures with (near-)equal C(f->L) ---")
    print("  the single rate C(f->L) is one real functional of the pair (m2, m3):")
    print(f"     1 - C(f->L) = kappa m2/(2 q log q) + gamma m3/(q^(3/2) log q) + O(q^-2)")
    kf = KAPPA / (2 * q * math.log(q))
    # gamma = C(beta*,3)/(beta*+1); the m3 term enters log Z with a minus sign
    gamma = (BETA_STAR * (BETA_STAR - 1) * (BETA_STAR - 2) / 6) / (BETA_STAR + 1)
    gf = gamma / (q ** 1.5 * math.log(q))
    print(f"     m2 coefficient = {kf:.6e},  m3 coefficient = {gf:.6e},  "
          f"ratio = {gf / kf:.6f}")
    print(f"     => rates are blind to the direction  delta m2 = {-gf / kf:+.6f} delta m3")

    reps = [group[0] for group in sigs.values()]
    vals = sorted(((rate(r["sig"], refs["L"])[0], r) for r in reps),
                  key=lambda t: t[0])
    ties = [(abs(vals[i + 1][0] - vals[i][0]), vals[i][1], vals[i + 1][1])
            for i in range(len(vals) - 1)]
    close = sorted((t for t in ties if abs(t[1]["mom"]["m2"] - t[2]["mom"]["m2"]) > 1e-9),
                   key=lambda t: t[0])
    print(f"  distinct signatures: {len(reps)}; "
          f"adjacent pairs with |Delta C(f->L)| < 1e-10: "
          f"{sum(1 for t in ties if t[0] < 1e-10)}")
    print(f"  {'|dC(f->L)|':>12s} {'m2 (f)':>14s} {'m2 (g)':>14s} {'d m2':>11s} "
          f"{'d m3':>11s} {'d m2 + r d m3':>14s}")
    for d, a, b in close[:5]:
        dm2 = b["mom"]["m2"] - a["mom"]["m2"]
        dm3 = b["mom"]["m3"] - a["mom"]["m3"]
        print(f"  {d:12.3e} {a['mom']['m2']:14.9f} {b['mom']['m2']:14.9f} "
              f"{dm2:+11.3e} {dm3:+11.3e} {dm2 + (gf / kf) * dm3:+14.3e}")
    if close:
        d, a, b = close[0]
        print("\n  certifying the tightest near-collision with the repo solver:")
        for tag, r in (("f", a), ("g", b)):
            v, bt = certified_rate(r["sig"], refs["L"])
            m = r["mom"]
            print(f"    {tag}: P={r['P']}  m2={m['m2']:.12f}  m3={m['m3']:+.9f}  "
                  f"N_max={m['max_fiber']}  exact C({tag}->L) = {v:.15f}")
        # full rate vector on the near-collision pair
        ra, rb = probe_row(a["sig"], refs), probe_row(b["sig"], refs)
        print(f"    {'probe':16s} {'|diff|':>12s}")
        for k in ra:
            print(f"    {k:16s} {abs(ra[k][0] - rb[k][0]):12.3e}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 211)
