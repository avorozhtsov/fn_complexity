#!/usr/bin/env python3
"""T2.3 part 2: is the Katz-Sarnak symmetry type visible in the exchange rates?

Set-up.  For ``f : A^2 -> A^1`` over ``F_q`` with full image, ``N_c = q - a_c``.
Two exact facts pin down what a rate can possibly see:

    (i)  sum_c a_c   = 0                    (because sum_c N_c = q^2)
    (ii) sum_c a_c^2 = Z_f(2) - q^3         (definition of the fiber square)

(i) is the structural obstruction of this note: the *first* moment of the
normalised traces is identically zero for every fibration of the affine plane,
and the first moment is exactly the statistic that separates orthogonal from
symplectic symmetry in Katz-Sarnak.  So symmetry type cannot be read off.

(ii) says the leading correction that the exchange rate does see is the point
count of the fiber square, i.e. by Lang-Weil the number of ``F_q``-rational
irreducible components of ``X x_Y X``.  For ``f = y^2 - P(x)`` this is
``nu(P)/q``, and ``nu(P)/q`` converges to the *rank* of the monodromy group of
``P`` (its number of orbits on ordered pairs).  Rank 2 (2-transitive monodromy)
gives the semicircle value ``m2 = 1`` for every genus.

The script therefore reports, at fixed q:

A. the exact vanishing of the first moment across every family, including ones
   that are genuinely orthogonal (quadratic twists) and unitary-flavoured;
B. ``m2`` versus the monodromy rank, for cyclic (x^d), dihedral (Dickson) and
   full-symmetric (random P) branch maps;
C. the genus / monodromy split: ``C(L -> f)`` reads ``2g`` through the extreme
   trace, ``C(f -> L)`` reads ``m2`` through the interior tangency, and the two
   readings order the families differently;
D. matched pairs: same genus + same rank, different symmetry mechanism -> the
   rate vectors agree to the order where the moments agree;
E. the one thing that does separate a twist (orthogonal) family from a
   big-monodromy (symplectic) one: not the value of m2 at a single q, but its
   distribution over q -- deterministic 1 versus Sato-Tate on [0, 4].

Run:  python research/m_and_e_and_a_c/t2_3_symmetry.py
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

CSV_FIXED_Q = HERE / "t2_3_symmetry_fixed_q.csv"
CSV_ST = HERE / "t2_3_twist_sato_tate.csv"


def genus_of_hyperelliptic(deg: int) -> int:
    return (deg - 1) // 2


# ------------------------------------------------------------------ A and B


def first_moment_and_rank(q: int) -> None:
    print(f"A. the first moment is identically zero  (q = {q})")
    rng = np.random.default_rng(11)
    items = []
    for deg in (3, 5, 7):
        P = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
        items.append((f"y^2 = P_{deg}(x) + c (generic, symplectic)", F.hyperelliptic(q, P)))
    items.append(("y^2 = x^3 + c (sextic twists, CM)", T.sextic_twist_map(q)))
    items.append(("x y^2 - x^2 (quartic twists, CM)", T.quartic_twist_map(q)))
    items.append(("P(x) y^2, P = x^3+x+1 (quadratic twists, orthogonal)",
                  T.quadratic_twist_map(q, [1, 1, 0, 1])))
    items.append(("x^3 + y^3 (Fermat pencil)", F.additive(q, [0, 0, 0, 1], [0, 0, 0, 1])))
    items.append(("x y (split conic)", F.split_conic(q)))
    print(f"   {'family':56s} {'sum a_c':>9} {'m2 (exact)':>12} {'Z_f(2)/q^3 - 1':>16}")
    for label, counts in items:
        mom = T.trace_moments(counts, q)
        z2 = float((counts.astype(np.float64) ** 2).sum())
        print(f"   {label:56s} {mom['sum_a']:>9.0f} {mom['m2']:>12.6f} "
              f"{z2 / q ** 3 - 1:>16.6f}")

    print(f"\nB. m2 = nu(P)/q - 1 and the monodromy rank of P  (q = {q})")
    print(f"   {'P':28s} {'monodromy':22s} {'rank':>5} {'nu(P)/q':>10} {'m2':>10}")
    rows = []
    for d in (2, 3, 4, 5, 6):
        e = math.gcd(d, q - 1)
        rows.append((f"x^{d}", f"cyclic Z/{e} (e=gcd(d,q-1))", e, [0] * d + [1]))
    for n in (3, 4, 5, 6):
        g1, g2 = math.gcd(n, q - 1), math.gcd(n, q + 1)
        rank = 1 + (n // 2 if n % 2 == 0 else (n - 1) // 2)
        rows.append((f"D_{n} (Dickson)", f"dihedral, gcds ({g1},{g2})", rank, T.dickson(n)))
    for deg in (4, 5, 6):
        P = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
        rows.append((f"random deg {deg}", f"S_{deg} (generic)", 2, P))
    for label, mono, rank, P in rows:
        nuP = T.nu([c % q for c in P], q)
        print(f"   {label:28s} {mono:22s} {rank:>5} {nuP / q:>10.4f} {nuP / q - 1:>10.4f}")


# ---------------------------------------------------------------------- C, D


def fixed_q_table(q: int, rows: list[list]) -> None:
    print(f"\nC. genus versus monodromy at fixed q = {q}")
    refs = references(q)
    rng = np.random.default_rng(2024)
    fam: list[tuple[str, int, np.ndarray]] = []
    for deg in (3, 5, 7, 9):
        P = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
        fam.append((f"generic g={genus_of_hyperelliptic(deg)}  y^2=P_{deg}(x)+c",
                    genus_of_hyperelliptic(deg), F.hyperelliptic(q, P)))
    fam.append(("CM      g=1  y^2=x^3+c", 1, T.sextic_twist_map(q)))
    fam.append(("CM      g=1  y^2=x^4+c", 1, F.hyperelliptic(q, [0] * 4 + [1])))
    fam.append(("CM      g=2  y^2=x^5+c", 2, F.hyperelliptic(q, [0] * 5 + [1])))
    fam.append(("CM      g=2  y^2=x^6+c", 2, F.hyperelliptic(q, [0] * 6 + [1])))
    fam.append(("CM      g=3  y^2=x^7+c", 3, F.hyperelliptic(q, [0] * 7 + [1])))
    fam.append(("Dickson g=1  y^2=D_4(x)+c", 1,
                F.hyperelliptic(q, [c % q for c in T.dickson(4)])))
    fam.append(("Dickson g=2  y^2=D_5(x)+c", 2,
                F.hyperelliptic(q, [c % q for c in T.dickson(5)])))
    fam.append(("Dickson g=2  y^2=D_6(x)+c", 2,
                F.hyperelliptic(q, [c % q for c in T.dickson(6)])))
    fam.append(("twist   g=1  P(x)y^2, P=x^3+x+1", 1,
                T.quadratic_twist_map(q, [1, 1, 0, 1])))

    print(f"   {'family':32s} {'g':>2} {'m2':>8} {'maxfib-q':>9} {'|a|max/rq':>10} "
          f"{'C(f->L)':>16} {'C(L->f)':>16} {'m2 from C(f->L)':>16} {'2g from C(L->f)':>16}")
    table = []
    for label, g, counts in fam:
        mom = T.trace_moments(counts, q)
        sig = sig_from_counts(counts)
        c_fL, b_fL = rate(sig, refs["L"])
        c_Lf, _ = rate(refs["L"], sig)
        table.append((label, g, mom, c_fL, c_Lf, b_fL))
        print(f"   {label:32s} {g:>2} {mom['m2']:>8.4f} "
              f"{mom['max_fiber'] - q:>9d} {mom['absmax'] / math.sqrt(q):>10.3f} "
              f"{c_fL:>16.12f} {c_Lf:>16.12f} "
              f"{T.m2_from_rate(c_fL, q, KAPPA):>16.4f} {T.weil_scale(c_Lf, q):>16.4f}")
        rows.append([q, label, g, f"{mom['m2']:.9f}", f"{mom['absmax']:.0f}",
                     f"{c_fL:.15f}", f"{c_Lf:.15f}",
                     f"{T.weil_scale(c_Lf, q):.9f}",
                     f"{T.m2_from_rate(c_fL, q, KAPPA):.9f}"])

    order_fL = [t[0] for t in sorted(table, key=lambda t: t[3])]
    order_Lf = [t[0] for t in sorted(table, key=lambda t: t[4])]
    print("\n   rank by C(f->L) (ascending = biggest m2 first) vs by C(L->f) "
          "(ascending = biggest fiber first)")
    for i, (a, b) in enumerate(zip(order_fL, order_Lf), 1):
        mark = "  " if a == b else " *"
        print(f"     {i:2d}.{mark} {a:34s} | {b}")
    disagree = sum(1 for a, b in zip(order_fL, order_Lf) if a != b)
    print(f"   positions where the two readings disagree: {disagree}/{len(order_fL)}")


def matched_pairs(q: int) -> None:
    print(f"\nD. matched pairs at q = {q}")
    refs = references(q)
    probes = [("C(f->L)", "L", "d"), ("C(L->f)", "L", "i"),
              ("C(f->Xsplit)", "Xsplit", "d"), ("C(Xsplit->f)", "Xsplit", "i"),
              ("C(f->Xaniso)", "Xaniso", "d"), ("C(Xaniso->f)", "Xaniso", "i")]

    def vec(counts):
        s = sig_from_counts(counts)
        return [rate(s, refs[r])[0] if role == "d" else rate(refs[r], s)[0]
                for _, r, role in probes]

    # (a) same genus, same nu (hence identical m2), different polynomial
    rng = np.random.default_rng(99)
    buckets: dict[int, list] = {}
    for _ in range(4000):
        P = [int(v) for v in rng.integers(0, q, size=3)] + [1]
        buckets.setdefault(T.nu(P, q), []).append(P)
        cand = buckets[T.nu(P, q)]
        if len(cand) >= 2:
            m3a = T.trace_moments(F.hyperelliptic(q, cand[0]), q)["m3"]
            m3b = T.trace_moments(F.hyperelliptic(q, cand[-1]), q)["m3"]
            if abs(m3a - m3b) > 0.05:
                pair_a, pair_b = cand[0], cand[-1]
                break
    else:
        pair_a = pair_b = None

    pairs = []
    if pair_a is not None:
        pairs.append((f"g=1 y^2=P(x)+c, P={pair_a}", F.hyperelliptic(q, pair_a),
                      f"g=1 y^2=P(x)+c, P={pair_b}", F.hyperelliptic(q, pair_b)))
    pairs.append(("CM quartic model y^2=x^4+c", F.hyperelliptic(q, [0] * 4 + [1]),
                  "quartic twist family x y^2-x^2", T.quartic_twist_map(q)))
    pairs.append(("generic g=1 (symplectic, m2=1)",
                  F.hyperelliptic(q, [int(v) for v in rng.integers(0, q, size=3)] + [1]),
                  "CM g=1 (finite monodromy, m2=2)", T.sextic_twist_map(q)))
    for la, ca, lb, cb in pairs:
        ma, mb = T.trace_moments(ca, q), T.trace_moments(cb, q)
        va, vb = vec(ca), vec(cb)
        print(f"\n   {la}\n   {lb}")
        print(f"      m2 {ma['m2']:.6f} / {mb['m2']:.6f}   "
              f"m3 {ma['m3']:.6f} / {mb['m3']:.6f}   "
              f"m4 {ma['m4']:.6f} / {mb['m4']:.6f}   "
              f"max fiber {ma['max_fiber']} / {mb['max_fiber']}")
        for (name, _, _), x, y in zip(probes, va, vb):
            print(f"      {name:14s} {x:.12f}  {y:.12f}   diff = {abs(x - y):.3e}")


# ------------------------------------------------------------------------ E


def sato_tate(rows: list[list]) -> None:
    print("\nE. distribution of m2 over q: symplectic (deterministic) vs twist (random)")
    P = [1, 1, 0, 1]
    primes = T.primes_upto(2000, lo=11)
    rng = np.random.default_rng(5)
    twist_m2, gen_m2 = [], []
    for q in primes:
        x = np.arange(q, dtype=np.int64)
        z = int((F.poly_eval(P, x, q) == 0).sum())
        if z != 0:
            continue
        a_E = q - int(F.hyperelliptic(q, P)[0])
        # exact: nu = q for an irreducible cubic branch locus is not used here;
        # the twist family's m2 follows from the exact three-value signature.
        m2 = (q - 1) * (z * z + a_E * a_E) / q ** 2 + (z * (q - 1) / q) ** 2
        twist_m2.append((q, a_E, m2))
        rows.append([q, a_E, f"{m2:.9f}"])
    for q in [101, 211, 401, 601, 809, 1009]:
        Pg = [int(v) for v in rng.integers(0, q, size=3)] + [1]
        gen_m2.append(T.nu(Pg, q) / q - 1)
    tm = np.array([t[2] for t in twist_m2])
    print(f"   quadratic-twist family of y^2=x^3+x+1, over the {len(tm)} primes q < 2000")
    print(f"   with irreducible branch cubic (z = 0):")
    print(f"      mean m2 = {tm.mean():.4f}   (Sato-Tate prediction E[4cos^2] = 1)")
    print(f"      range   = [{tm.min():.4f}, {tm.max():.4f}]   (Weil bound: [0, 4))")
    print(f"      std     = {tm.std():.4f}   (Sato-Tate prediction sqrt(E[16cos^4]-1) "
          f"= {math.sqrt(2.0 - 1.0):.4f})")
    edges = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    hist, _ = np.histogram(tm, bins=edges)
    st = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        th = np.linspace(0, math.pi, 20001)
        mask = (4 * np.cos(th) ** 2 >= lo) & (4 * np.cos(th) ** 2 < hi)
        st.append(float(np.trapezoid((2 / math.pi) * np.sin(th) ** 2 * mask, th)) * len(tm))
    print(f"      {'bin':>12} {'observed':>9} {'Sato-Tate':>10}")
    for (lo, hi), h, s in zip(zip(edges[:-1], edges[1:]), hist, st):
        print(f"      [{lo:.1f},{hi:.1f})".rjust(12) + f" {h:>9d} {s:>10.1f}")
    print(f"   generic symplectic families instead give m2 = {np.mean(gen_m2):.4f} "
          f"+- {np.std(gen_m2):.4f} at every q (deterministic, no spread)")


def main() -> int:
    q = 601
    first_moment_and_rank(q)
    rows: list[list] = []
    fixed_q_table(421, rows)
    fixed_q_table(q, rows)
    matched_pairs(q)
    with CSV_FIXED_Q.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["q", "family", "genus", "m2", "max_abs_a", "C_f_to_L", "C_L_to_f",
                    "weil_scale", "m2_recovered"])
        w.writerows(rows)
    strows: list[list] = []
    sato_tate(strows)
    with CSV_ST.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["q", "a_E", "m2_of_twist_family"])
        w.writerows(strows)
    print(f"\nwritten: {CSV_FIXED_Q.name}, {CSV_ST.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
