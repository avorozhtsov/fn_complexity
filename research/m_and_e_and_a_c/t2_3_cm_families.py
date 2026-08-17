#!/usr/bin/env python3
"""T2.3 part 1: CM families and the exact congruence detectors they carry.

The established example is ``f(x,y) = y^2 - x^3`` (sextic twists, CM by
``Z[zeta_3]``), whose signature collapses to the flat one exactly when
``q = 2 mod 3``.  This script asks how general that is.

Contents
--------
1.  The exact second-moment identity ``sum_c a_c^2 = q (nu(P) - q)`` for
    ``f = y^2 - P(x)``, and its corollary
    ``sum_c a_c^2 = (gcd(d, q-1) - 1) q (q - 1)`` for ``P = x^d``.
2.  The signature of ``y^2 - x^d`` depends only on ``e = gcd(d, q-1)``,
    is flat exactly when ``e = 1`` and is the split conic exactly when ``e = 2``.
3.  The quartic-twist family ``f = x y^2 - x^2`` (CM by ``Z[i]``): the exact
    coincidence at ``q = 3 mod 4`` is with the SPLIT CONIC, not with ``L``.
4.  The genus-2 family ``y^2 - x^5`` (CM by ``Z[zeta_5]``) and ``y^2 - x^7``.
5.  Quadratic-twist families ``f = P(x) y^2`` of a fixed curve ``E``: the
    signature is exactly three-valued and detects ``a_E = 0`` (supersingularity
    of E at q), including for non-CM E where that is not a congruence.

Run:  python research/m_and_e_and_a_c/t2_3_cm_families.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import KAPPA, rate, references, sig_from_counts  # noqa: E402
import ffmaps as F  # noqa: E402
import t2_3_common as T  # noqa: E402

CSV_MAIN = HERE / "t2_3_cm_families.csv"
CSV_TWIST = HERE / "t2_3_quadratic_twists.csv"


def tag_of(sig, refs) -> str:
    k = T.sig_key(sig.values, sig.mults)
    for name in ("L", "Xsplit", "Xaniso"):
        if name in refs and k == T.sig_key(refs[name].values, refs[name].mults):
            return name
    return f"{len(sig.values)}-valued"


# --------------------------------------------------------------- section 1


def check_identities(primes: list[int]) -> None:
    print("1. exact identities")
    rng = np.random.default_rng(20260817)
    worst_sum, bad = 0.0, 0
    for q in primes:
        for deg in (2, 3, 4, 5, 6, 7, 9):
            for _ in range(3):
                P = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
                counts = F.hyperelliptic(q, P)
                a = (q - counts).astype(np.int64)
                worst_sum = max(worst_sum, abs(int(a.sum())))
                if int((a * a).sum()) != q * (T.nu(P, q) - q):
                    bad += 1
    print(f"   sum_c a_c = 0            : worst |sum| over the sample = {worst_sum}")
    print(f"   sum_c a_c^2 = q(nu(P)-q) : mismatches = {bad}")

    bad = 0
    for q in primes:
        for d in range(1, 15):
            P = [0] * d + [1]
            a = (q - F.hyperelliptic(q, P)).astype(np.int64)
            if int((a * a).sum()) != (math.gcd(d, q - 1) - 1) * q * (q - 1):
                bad += 1
    print(f"   P = x^d: sum_c a_c^2 = (gcd(d,q-1)-1) q(q-1) : mismatches = {bad}")


# --------------------------------------------------------------- section 2


def gcd_classification(primes: list[int]) -> None:
    print("\n2. signature of y^2 - x^d depends only on e = gcd(d, q-1)")
    collisions, bad = 0, 0
    for q in primes:
        refs = references(q)
        by_e: dict[int, set] = {}
        for d in range(1, 15):
            s = sig_from_counts(F.hyperelliptic(q, [0] * d + [1]))
            by_e.setdefault(math.gcd(d, q - 1), set()).add(T.sig_key(s.values, s.mults))
        for e, keys in by_e.items():
            if len(keys) != 1:
                bad += 1
            collisions += 1
        # e = 1 <-> flat, e = 2 <-> split conic
        for e, keys in by_e.items():
            k = next(iter(keys))
            if e == 1 and k != T.sig_key(refs["L"].values, refs["L"].mults):
                bad += 1
            if e == 2 and k != T.sig_key(refs["Xsplit"].values, refs["Xsplit"].mults):
                bad += 1
    print(f"   {collisions} (q, e) classes over {len(primes)} primes, "
          f"violations of 'one signature per e' + 'e=1 flat, e=2 split conic': {bad}")


# --------------------------------------------------------------- section 3+4


FAMILIES = [
    ("y^2-x^3  (sextic, CM Z[zeta_3])", lambda q: T.sextic_twist_map(q), 3),
    ("xy^2-x^2 (quartic, CM Z[i])", lambda q: T.quartic_twist_map(q), 4),
    ("y^2-x^4  (CM Z[i], quartic model)", lambda q: F.hyperelliptic(q, [0] * 4 + [1]), 4),
    ("y^2-x^5  (genus 2, CM Z[zeta_5])", lambda q: F.hyperelliptic(q, [0] * 5 + [1]), 5),
    ("y^2-x^7  (genus 3, CM Z[zeta_7])", lambda q: F.hyperelliptic(q, [0] * 7 + [1]), 7),
]


def family_table(primes: list[int], rows: list[list]) -> None:
    print("\n3-4. CM families against the flat map L and the split conic X")
    for label, build, modulus in FAMILIES:
        print(f"\n   {label}    (congruence tested: q mod {modulus})")
        print(f"      {'q':>5} {'q%m':>4} {'sig':>12} {'m2':>9} {'C(f->L)':>16} "
              f"{'C(L->f)':>16} {'C(f->X)':>16} {'C(X->f)':>16}")
        for q in primes:
            counts = build(q)
            sig = sig_from_counts(counts)
            refs = references(q)
            mom = T.trace_moments(counts, q)
            c_fL, _ = rate(sig, refs["L"])
            c_Lf, _ = rate(refs["L"], sig)
            c_fX, _ = rate(sig, refs["Xsplit"])
            c_Xf, _ = rate(refs["Xsplit"], sig)
            tg = tag_of(sig, refs)
            print(f"      {q:>5} {q % modulus:>4} {tg:>12} {mom['m2']:>9.4f} "
                  f"{c_fL:>16.12f} {c_Lf:>16.12f} {c_fX:>16.12f} {c_Xf:>16.12f}")
            rows.append([label, q, q % modulus, tg, f"{mom['m2']:.9f}",
                         f"{c_fL:.15f}", f"{c_Lf:.15f}", f"{c_fX:.15f}", f"{c_Xf:.15f}"])


# ----------------------------------------------------------------- section 5


TWIST_CURVES = [
    ("E: y^2=x^3+x   (CM Z[i])", [0, 1, 0, 1]),
    ("E: y^2=x^3+1   (CM Z[zeta_3])", [1, 0, 0, 1]),
    ("E: y^2=x^3+x+1 (non-CM)", [1, 1, 0, 1]),
]


def quadratic_twist_table(primes: list[int], rows: list[list]) -> None:
    print("\n5. quadratic-twist families f = P(x) y^2: exactly three-valued")
    for label, P in TWIST_CURVES:
        ss = []
        bad = 0
        flags: dict[int, bool] = {}
        for q in primes:
            counts = T.quadratic_twist_map(q, P)
            sig = sig_from_counts(counts)
            x = np.arange(q, dtype=np.int64)
            z = int((F.poly_eval(P, x, q) == 0).sum())
            a_E = q - int(F.hyperelliptic(q, P)[0])
            predicted = sorted({v for v in (q + z * (q - 1), q - z - a_E, q - z + a_E)
                                if v > 0}, reverse=True)
            if [int(v) for v in sig.values] != predicted:
                bad += 1
            flags[q] = (a_E == 0)
            if a_E == 0:
                refs = references(q)
                ss.append((q, z, tag_of(sig, refs)))
                rows.append([label, q, z, tag_of(sig, refs)])
        cls = T.residue_classes(primes, flags)
        print(f"\n   {label}")
        print(f"      exact three-value formula violations over {len(primes)} primes: {bad}")
        print(f"      supersingular q (a_E = 0) up to {primes[-1]}: "
              f"{[s[0] for s in ss]}")
        shapes = sorted({(s[1], s[2]) for s in ss})
        print(f"      signature shapes there (z, tag): {shapes}")
        print(f"      congruence explaining the set: "
              f"{'q = ' + str(cls[1]) + ' mod ' + str(cls[0]) if cls else 'NONE up to mod 60'}")


def main() -> int:
    primes = [q for q in T.primes_upto(240, lo=11)]
    big = [101, 211, 401, 421, 461, 491, 601, 661]
    check_identities([101, 211, 401])
    gcd_classification(primes)
    rows: list[list] = []
    family_table(big, rows)
    with CSV_MAIN.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["family", "q", "q_mod_m", "signature", "m2",
                    "C_f_to_L", "C_L_to_f", "C_f_to_Xsplit", "C_Xsplit_to_f"])
        w.writerows(rows)
    trows: list[list] = []
    quadratic_twist_table(T.primes_upto(500, lo=5), trows)
    with CSV_TWIST.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["curve", "supersingular_q", "z", "signature"])
        w.writerows(trows)
    print(f"\nwritten: {CSV_MAIN.name}, {CSV_TWIST.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
