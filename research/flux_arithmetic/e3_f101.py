#!/usr/bin/env python3
"""E3.5 -- the ``F_101`` witness of brief B's addendum 2, where two of the three
pencils agree on every classical statistic the signature exposes.

``f_1`` and ``f_2`` share the largest fiber, ``m_2``, ``nu(P)``, hence ``phi``
and the addendum's ``phi~``.  Their comparison is decided by the interior alone.
This script asks what *does* separate them: the full trace vectors, the
fiberwise L-polynomials and isogeny classes, the splitting and ``p``-rank of the
Jacobians, the ramification and monodromy of the branch map, and whether they
are the same pencil up to an affine coordinate change.

Writes ``e3_f101.csv``.
"""

from __future__ import annotations

import collections
import csv
import math

import numpy as np

import curves

P = 101
PENCILS = {
    "f_1": {5: 1, 4: 70, 3: 28, 2: 15, 1: 11, 0: 31},
    "f_2": {5: 1, 4: 42, 3: 32, 2: 74, 1: 96, 0: 60},
    "f_3": {5: 1, 4: 72, 3: 21, 2: 2, 1: 6, 0: 57},
}


def normalise(coeffs):
    """Drop the constant term: it only relabels the fibers."""

    out = dict(coeffs)
    out.pop(0, None)
    return out


def affine_orbit(coeffs, degree=5, q=P):
    """Least encoding of the ``G0``-orbit of a monic ``P`` with ``P(0) = 0``."""

    c = [coeffs.get(j, 0) % q for j in range(1, degree)]
    squares = {(x * x) % q for x in range(1, q)}
    best = None
    for a in range(1, q):
        if pow(a, degree, q) not in squares:
            continue
        inv = pow(a, q - 2, q)
        for b in range(q):
            new = []
            for j in range(1, degree):
                acc = 0
                for k in range(j, degree + 1):
                    ck = c[k - 1] if k < degree else 1
                    acc += ck * math.comb(k, j) * pow(b, k - j, q)
                new.append(acc % q * pow(a, j, q) % q * pow(inv, degree, q) % q)
            key = tuple(new)
            best = key if best is None or key < best else best
    return best


def main() -> None:
    rows = []
    data = {}
    for name, coeffs in PENCILS.items():
        c = normalise(coeffs)
        fibers = curves.pencil_fibers(c, P, 5)
        sig = tuple(sorted((f.n_affine for f in fibers), reverse=True))
        traces = tuple(sorted(f.s1 for f in fibers))
        a = P - np.array([f.n_affine for f in fibers], dtype=float)
        nu = int((curves.value_multiplicities(c, P) ** 2).sum())
        L = collections.Counter((f.s1, f.e2) for f in fibers)
        ram = collections.Counter(
            tuple(sorted(d for d, m in curves.factorisation_type(c, P, 5, u)
                         for _ in range(m)))
            for u in range(P)
        )
        data[name] = {
            "signature": sig,
            "traces": traces,
            "max_fiber": sig[0],
            "m2": float((a**2).sum() / P**2),
            "nu": nu,
            "splits": sum(f.splits(P) for f in fibers),
            "prank": collections.Counter(f.p_rank(P) for f in fibers),
            "distinct_L": len(L),
            "L": L,
            "smooth": sum(f.smooth for f in fibers),
            "weil_ok": all(
                abs(f.s1) <= 4 * math.sqrt(P)
                and 2 * abs(f.s1) * math.sqrt(P) - 2 * P <= f.e2 <= f.s1**2 / 4 + 2 * P
                for f in fibers
            ),
            "crit": curves.derivative_type(c, P, 5),
            "ram": ram,
            "orbit": affine_orbit(c),
        }
    for name, d in data.items():
        print(f"\n{name}: max fiber {d['max_fiber']}, m2 = {d['m2']:.6f}, "
              f"nu(P) = {d['nu']}")
        print(f"   smooth fibers {d['smooth']}/{P};  "
              f"all (s1,e2) Weil-admissible: {d['weil_ok']}")
        print(f"   split Jacobians {d['splits']}/{P};  "
              f"p-rank {dict(sorted(d['prank'].items()))}")
        print(f"   distinct isogeny classes among the {P} fibers: {d['distinct_L']}")
        print(f"   P'(x) factors with degrees {d['crit']}")
        print(f"   branch-map Frobenius cycle types: {dict(sorted(d['ram'].items()))}")
        rows.append({
            "pencil": name, "max_fiber": d["max_fiber"], "m2": d["m2"],
            "nu": d["nu"], "splits": d["splits"],
            "p_rank_2": d["prank"].get(2, 0), "p_rank_1": d["prank"].get(1, 0),
            "p_rank_0": d["prank"].get(0, 0), "distinct_L": d["distinct_L"],
            "smooth": d["smooth"],
            "signature": " ".join(str(v) for v in d["signature"]),
        })

    print("\n" + "=" * 70)
    print("f_1 against f_2 -- what the signature does NOT expose")
    d1, d2 = data["f_1"], data["f_2"]
    print(f"   same signature?            {d1['signature'] == d2['signature']}")
    print(f"   same multiset of traces?   {d1['traces'] == d2['traces']}")
    print(f"   same G0-orbit (same pencil up to x -> ax+b)?  "
          f"{d1['orbit'] == d2['orbit']}")
    shared = set(d1["L"]) & set(d2["L"])
    print(f"   shared fiber isogeny classes: {len(shared)} of "
          f"{d1['distinct_L']} and {d2['distinct_L']}")
    print(f"   same L multiset?           {d1['L'] == d2['L']}")
    print(f"   split Jacobians:           {d1['splits']} vs {d2['splits']}")
    print(f"   p-rank:                    {dict(sorted(d1['prank'].items()))} vs "
          f"{dict(sorted(d2['prank'].items()))}")
    print(f"   branch ramification:       {d1['crit']} vs {d2['crit']}")
    print(f"   cycle-type multiset equal? {d1['ram'] == d2['ram']}")
    s1, s2 = d1["signature"], d2["signature"]
    for i, (x, y) in enumerate(zip(s1, s2)):
        if x != y:
            print(f"   the two signatures first differ at order statistic {i}: "
                  f"{x} vs {y}")
            break
    v1 = P - np.array(s1, dtype=float)
    v2 = P - np.array(s2, dtype=float)
    for k in (3, 4, 5, 6):
        print(f"   m{k}: {(v1**k).sum()/P**k:+.6f} vs {(v2**k).sum()/P**k:+.6f}")
    mult1 = sum(1 for v in s1 if v == s1[0])
    mult2 = sum(1 for v in s2 if v == s2[0])
    print(f"   multiplicity of the largest fiber: {mult1} vs {mult2}")
    print(f"   smallest fiber: {s1[-1]} vs {s2[-1]}")

    with open("e3_f101.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
