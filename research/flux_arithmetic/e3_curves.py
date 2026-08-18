#!/usr/bin/env python3
"""E3 -- does the certified cycle descend to the curves?

Four sections.

1. The arithmetic of the eleven fibers of each of the three certified ``F_11``
   pencils: point counts over ``F_11`` and ``F_121``, L-polynomials, splitting
   of the Jacobian, ``p``-rank, and the pencil-level invariants (value
   multiplicities of ``P``, ``nu(P)``, ramification of the branch map).

2. The fibers of the map ``pencil -> signature``, counted in orbits of the
   affine group that acts by isomorphisms of the fibration.

3. Whether two pencils with the *same* signature can have *different* fiberwise
   isogeny data.  This is brief C's crux: if they can, the signature -- and
   therefore the flux, which is a function of the signature -- collapses a
   distinction the curves make.

4. The quantified version of brief C's claim.

Writes ``e3_certified_fibers.csv``, ``e3_signature_fibers.csv``,
``e3_collisions.csv``.
"""

from __future__ import annotations

import collections
import csv
import itertools
import math

import numpy as np

import curves
import flux
import orbits
import pools

COARSE, FINE = 2400, 13600

TRIPLE = {
    "A": (5, {5: 1, 4: 3, 3: 4, 2: 1, 1: 1}),
    "B": (6, {6: 1, 5: 9, 4: 7, 3: 2, 1: 1}),
    "C": (6, {6: 1, 5: 10, 4: 8, 3: 8, 2: 1, 1: 2}),
}


def weil_ok(f: curves.Fiber, p: int) -> bool:
    """Honda-Tate/Ruck admissibility of ``(s1, e2)`` for a genus-two curve."""

    if abs(f.s1) > 4 * math.sqrt(p):
        return False
    return (2 * abs(f.s1) * math.sqrt(p) - 2 * p) <= f.e2 <= f.s1**2 / 4 + 2 * p


def section1(p: int = 11) -> list[dict]:
    rows = []
    print("=" * 78)
    print("E3.1  the three certified pencils, fiber by fiber")
    summary = {}
    for name, (degree, coeffs) in TRIPLE.items():
        fibers = curves.pencil_fibers(coeffs, p, degree)
        sig = tuple(sorted((f.n_affine for f in fibers), reverse=True))
        nu = int((curves.value_multiplicities(coeffs, p) ** 2).sum())
        a = p - np.array([f.n_affine for f in fibers], dtype=float)
        m2 = float((a**2).sum() / p**2)
        ram = collections.Counter(
            curves.factorisation_type(coeffs, p, degree, u) for u in range(p)
        )
        Lmultiset = collections.Counter((f.s1, f.e2) for f in fibers)
        summary[name] = {
            "degree": degree,
            "signature": sig,
            "nu": nu,
            "m2": m2,
            "m2_from_nu": nu / p - 1,
            "splits": sum(f.splits(p) for f in fibers),
            "prank": collections.Counter(f.p_rank(p) for f in fibers),
            "all_smooth": all(f.smooth for f in fibers),
            "weil_ok": all(weil_ok(f, p) for f in fibers),
            "L": Lmultiset,
            "ram": ram,
            "distinct_L": len(Lmultiset),
            "crit": curves.derivative_type(coeffs, p, degree),
            "rational_crit_values": curves.rational_critical_values(coeffs, p, degree),
            "cycle_types": sorted({tuple(sorted(d for d, m in t for _ in range(m)))
                                   for t in ram}),
        }
        for f in fibers:
            rows.append(
                {
                    "pencil": name,
                    "degree": degree,
                    "c": f.c,
                    "N_affine": f.n_affine,
                    "points_F11": f.points_fp,
                    "points_F121": f.points_fp2,
                    "s1": f.s1,
                    "e2": f.e2,
                    "splits": int(f.splits(p)),
                    "p_rank": f.p_rank(p),
                    "smooth": int(f.smooth),
                    "weil_admissible": int(weil_ok(f, p)),
                }
            )
        s = summary[name]
        print(f"\n  pencil {name}  (degree {degree})   signature {sig}")
        print(f"    all fibers smooth: {s['all_smooth']};  every (s1,e2) Weil-admissible: {s['weil_ok']}")
        print(f"    nu(P) = {nu};  m2 = {m2:.6f};  nu(P)/q - 1 = {s['m2_from_nu']:.6f}")
        print(f"    Jacobians isogenous to a product of elliptic curves: {s['splits']}/{p}")
        print(f"    p-rank distribution: {dict(sorted(s['prank'].items()))}")
        print(f"    distinct isogeny classes among the {p} fibers: {s['distinct_L']}")
        print(f"    critical points: P'(x) factors with degrees {s['crit']}; "
              f"{s['rational_crit_values']} critical values lie in F_{p}")
        print(f"    branch map P(x) - u, Frobenius cycle types over F_{p}:")
        for t, k in sorted(ram.items()):
            cycle = tuple(sorted(d for d, m in t for _ in range(m)))
            print(f"        {k:2d} x  cycle type {cycle}")

    print("\n  do the three pencils share any fiber isogeny class?")
    for x, y in itertools.combinations(TRIPLE, 2):
        shared = set(summary[x]["L"]) & set(summary[y]["L"])
        print(f"    {x} & {y}: {len(shared)} shared classes  {sorted(shared)}")
    print("\n  classical invariants that separate the three:")
    for key in ("degree", "nu", "splits", "distinct_L"):
        print(f"    {key:12s} " + "  ".join(f"{k}={summary[k][key]}" for k in TRIPLE))
    print("    p-rank      " + "  ".join(f"{k}={dict(sorted(summary[k]['prank'].items()))}" for k in TRIPLE))
    return rows


def section2(q: int) -> list[dict]:
    print("=" * 78)
    print(f"E3.2  the fibers of  pencil -> signature  at q = {q}")
    K = pools._chi_matrix(q)
    per_signature: dict[tuple[int, ...], dict] = {}
    orbit_reps: dict[tuple[int, ...], list[tuple[int, int]]] = {}
    for degree in (5, 6):
        hist = pools._value_histograms(q, degree)
        N = q + hist.astype(np.int64) @ K.T
        N.sort(axis=1)
        N = N[:, ::-1]
        can = orbits.canonical(q, degree)
        seen = set()
        for i in range(len(can)):
            key = (degree, int(can[i]))
            if key in seen:
                continue
            seen.add(key)
            sig = tuple(int(v) for v in N[i])
            if min(sig) <= 0:
                continue
            per_signature.setdefault(sig, {"orbits": 0, "by_degree": collections.Counter()})
            per_signature[sig]["orbits"] += 1
            per_signature[sig]["by_degree"][degree] += 1
            orbit_reps.setdefault(sig, []).append(key)
        print(f"    degree {degree}: |G0| = {orbits.group_order(q, degree)}, "
              f"{len(np.unique(can))} orbits from {len(can)} polynomials")

    sizes = collections.Counter(v["orbits"] for v in per_signature.values())
    total = sum(v["orbits"] for v in per_signature.values())
    print(f"    {total} pencils (affine-group orbits) -> {len(per_signature)} signatures")
    print(f"    fiber-size distribution of pencil -> signature:")
    for k in sorted(sizes):
        print(f"        {sizes[k]:4d} signatures are hit by {k:3d} pencil(s)")
    print(f"    mean fiber {total/len(per_signature):.2f}, "
          f"max fiber {max(sizes)}, injective on {sizes.get(1,0)}/{len(per_signature)} signatures")
    rows = [
        {
            "q": q,
            "signature": " ".join(str(v) for v in sig),
            "orbits": rec["orbits"],
            "orbits_deg5": rec["by_degree"].get(5, 0),
            "orbits_deg6": rec["by_degree"].get(6, 0),
        }
        for sig, rec in sorted(per_signature.items())
    ]
    return rows, orbit_reps


def section3(q: int, orbit_reps: dict) -> list[dict]:
    print("=" * 78)
    print(f"E3.3  do co-signature pencils differ arithmetically?  q = {q}")
    grids = {d: orbits.coefficient_grid(q, d) for d in (5, 6)}

    def lmultiset(degree: int, index: int):
        coeffs = {degree: 1}
        for j, value in enumerate(grids[degree][index], start=1):
            if value:
                coeffs[j] = int(value)
        fibers = curves.pencil_fibers(coeffs, q, degree)
        return tuple(sorted((f.s1, f.e2) for f in fibers)), coeffs

    rows = []
    same_L = diff_L = 0
    diff_within_degree = 0
    examples = []
    for sig, reps in sorted(orbit_reps.items()):
        if len(reps) < 2:
            continue
        data = {}
        for degree, index in reps:
            L, coeffs = lmultiset(degree, index)
            data[(degree, index)] = (L, coeffs)
        classes = collections.Counter(v[0] for v in data.values())
        if len(classes) == 1:
            same_L += 1
        else:
            diff_L += 1
            per_degree = collections.defaultdict(set)
            for (degree, _), (L, _) in data.items():
                per_degree[degree].add(L)
            if any(len(v) > 1 for v in per_degree.values()):
                diff_within_degree += 1
            if len(examples) < 3:
                examples.append((sig, data))
        rows.append(
            {
                "q": q,
                "signature": " ".join(str(v) for v in sig),
                "pencils": len(reps),
                "distinct_isogeny_multisets": len(classes),
            }
        )
    print(f"    {same_L + diff_L} signatures are realised by 2 or more pencils")
    print(f"        {same_L} of them by pencils all sharing one fiberwise isogeny multiset")
    print(f"        {diff_L} of them by pencils with *different* isogeny multisets")
    print(f"        ({diff_within_degree} of those already differ within a single degree)")
    if examples:
        sig, data = examples[0]
        print(f"\n    example, signature {sig}:")
        for (degree, index), (L, coeffs) in sorted(data.items()):
            terms = " + ".join(
                f"{v}x^{k}" for k, v in sorted(coeffs.items(), reverse=True) if v
            )
            print(f"        deg {degree}  P = {terms}")
            print(f"          (s1,e2) multiset {L}")
    return rows


def section4(q: int) -> dict:
    print("=" * 78)
    print(f"E3.4  the quantified brief-C statement at q = {q}")
    S, _ = pools.arithmetic_pool(q)
    A = flux.flux_matrix(S, flux.beta_grid(q, COARSE, FINE))
    h = flux.hodge(A)
    cyc, tri = flux.cycle_count(A)
    n = len(S)
    off = ~np.eye(n, dtype=bool)
    strict = np.abs(A) > flux.TIE
    print(f"    n = {n} signatures, {n*(n-1)//2} unordered pairs, "
          f"{int((off & ~strict).sum())//2} ties")
    print(f"    energy: gradient {h['grad_energy']*100:.3f}% , curl {h['curl_energy']*100:.3f}%")
    print(f"    the best scalar psi_opt gets the direction right on "
          f"{h['order_agreement']*100:.3f}% of ordered pairs")
    print(f"    strict 3-cycles: {cyc} of {tri} triangles ({cyc/tri*100:.5f}%)")
    return {
        "q": q,
        "n": n,
        "grad_energy": h["grad_energy"],
        "curl_energy": h["curl_energy"],
        "order_agreement": h["order_agreement"],
        "cycles": cyc,
        "triangles": tri,
    }


def main() -> None:
    rows1 = section1()
    with open("e3_certified_fibers.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows1[0]))
        w.writeheader()
        w.writerows(rows1)

    all2, all3, all4 = [], [], []
    for q in (11, 13):
        rows2, reps = section2(q)
        if q == 11:
            print("\n    the three certified signatures, and how many pencils share each:")
            for name, (degree, coeffs) in TRIPLE.items():
                fibers = curves.pencil_fibers(coeffs, q, degree)
                sig = tuple(sorted((f.n_affine for f in fibers), reverse=True))
                print(f"      {name}: {len(reps[sig])} pencil(s) realise {sig}")
        all2.extend(rows2)
        all3.extend(section3(q, reps))
        all4.append(section4(q))
    with open("e3_signature_fibers.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all2[0]))
        w.writeheader()
        w.writerows(all2)
    with open("e3_collisions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all3[0]))
        w.writeheader()
        w.writerows(all3)
    with open("e3_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all4[0]))
        w.writeheader()
        w.writerows(all4)


if __name__ == "__main__":
    main()
