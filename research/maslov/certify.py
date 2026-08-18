#!/usr/bin/env python3
"""40-digit certification of the session's headline numbers.

Three blocks.

1.  The known 3-cycle's ``lambda_c`` under uniform priors with right endpoint
    ``S_+ = 6, 12, 20``: at ``lambda`` slightly below / above the reported
    ``lambda_c`` the sign pattern is certified in 40-digit arithmetic.

2.  The ``F_11`` cycle ``(48, 64, 51)``: its tropical ``A`` and its ``A_lambda``
    on either side of the reported prior-free ``lambda_c = 1158.98``.

3.  The two flat witnesses of ``m3_phase.py`` are certified inside that script;
    the check is repeated here against ``m3_witness.json``.

Writes ``certify.json``.
"""

from __future__ import annotations

import json
import math

import mpmath as mp
import numpy as np

import common as C

DPS = 40


def mp_A(sigs_or_flat, lam, lo, hi, ng, flat=False):
    """40-digit ``A_lambda`` around a triple; returns the three edge values."""

    with mp.workdps(DPS):
        S = [mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * i / (ng - 1) for i in range(ng)]
        U = []
        for a in sigs_or_flat:
            if flat:
                R, L = mp.mpf(a[0]), mp.mpf(a[1])
                U.append([mp.log(R + mp.e**sk * L) for sk in S])
            else:
                U.append([mp.log(mp.log(mp.fsum([mp.power(mp.mpf(int(x)), mp.e**sk)
                                                 for x in a]))) for sk in S])
        wk = mp.mpf(1) / ng
        L = mp.mpf(lam)
        vals = []
        for i in range(3):
            f = [U[i][t] - U[(i + 1) % 3][t] for t in range(ng)]
            fmax, fmin = max(f), min(f)
            smax = (mp.log(mp.fsum([wk * mp.e ** (L * (z - fmax)) for z in f])) + L * fmax) / L
            smin = -(mp.log(mp.fsum([wk * mp.e ** (-L * (z - fmin)) for z in f])) - L * fmin) / L
            vals.append((smax + smin) / 2)
        return vals


def double_A(sigs, lam, lo, hi, ng):
    s, w = C.uniform_prior(lo, hi, ng)
    U = [C.u_of(a, s) for a in sigs]
    return [C.soft_mid(U[i] - U[(i + 1) % 3], w, lam) for i in range(3)]


def main():
    rec = {}
    print("=" * 78)
    print("1.  known 3-cycle: lambda_c depends on the uniform support only through S_+")
    print("=" * 78)
    NG = 2001  # a coarse but honestly declared discrete prior, identical in both arithmetics
    rec["known"] = []
    for lo, hi, lc in ((-20, 6, 184.49), (-8, 6, 184.49), (-12, 12, 286.71), (-8, 20, 349.92)):
        for tag, lam in (("below", lc * 0.9), ("above", lc * 1.1)):
            v = mp_A(C.CYCLE, lam, lo, hi, NG)
            d = double_A(C.CYCLE, lam, lo, hi, NG)
            cyc = all(z > 0 for z in v)
            err = max(abs(float(v[i]) - d[i]) for i in range(3))
            print(f"  supp [{lo},{hi}]  lambda = {lam:9.3f} ({tag:5s})  cycle: {str(cyc):5s}"
                  f"   |double-mp40| <= {err:.2e}")
            print("      " + "  ".join(mp.nstr(z, 22) for z in v))
            rec["known"].append(dict(lo=lo, hi=hi, lam=lam, tag=tag, cycle=cyc,
                                     A=[mp.nstr(z, 30) for z in v], err=err))

    print()
    print("=" * 78)
    print("2.  F_11 cycle (48, 64, 51): prior-free lambda_c ~ 1158.98")
    print("=" * 78)
    pool = C.f11_pool()
    tri = [pool[48], pool[64], pool[51]]
    print("  " + " ; ".join(str(a) for a in tri))
    rec["f11"] = []
    for lo, hi in ((-8, 8), (-24, 24)):
        for tag, lam in (("below", 1000.0), ("above", 1400.0)):
            v = mp_A(tri, lam, lo, hi, NG)
            d = double_A(tri, lam, lo, hi, NG)
            cyc = all(z > 0 for z in v)
            err = max(abs(float(v[i]) - d[i]) for i in range(3))
            print(f"  supp [{lo},{hi}]  lambda = {lam:8.1f} ({tag:5s})  cycle: {str(cyc):5s}"
                  f"   |double-mp40| <= {err:.2e}")
            print("      " + "  ".join(mp.nstr(z, 22) for z in v))
            rec["f11"].append(dict(lo=lo, hi=hi, lam=lam, tag=tag, cycle=cyc,
                                   A=[mp.nstr(z, 30) for z in v], err=err))

    print()
    print("=" * 78)
    print("3.  the flat witnesses of m3_phase.py")
    print("=" * 78)
    wit = json.load(open("m3_witness.json"))
    for tag, d in wit["witnesses"].items():
        print(f"  {tag}: powers {d['powers']}, band {d['cycle_band']}, "
              f"|double-mp40| <= {d['double_vs_mpmath']:.2e}")
        print("      " + "  ".join(d["mpmath40"]))
    rec["flat"] = {k: {"powers": v["powers"], "band": v["cycle_band"],
                       "mpmath40": v["mpmath40"], "err": v["double_vs_mpmath"]}
                   for k, v in wit["witnesses"].items()}

    json.dump(rec, open("certify.json", "w"), indent=1)
    print("\nwrote certify.json")


if __name__ == "__main__":
    main()
