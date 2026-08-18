"""Session brief I -- from the certified cone witness to integer signatures.

i_certify.py certifies, at 40 digits, a point of the cone C whose Hilbert
metric is EXACTLY s.C_4 with s = 0.2019801983...  This script turns that point
into honest integer signatures and measures how the distortion falls.

The translation is the tropical construction: the line (c, x) of Phi becomes
    floor(e^{lam c}) copies of the atom round(e^{lam x p}),
lam a projective scale (Cartesian powers) and p the common exponent that makes
every atom an integer (a_i -> a_i^p is the exact reparametrisation
beta -> beta/p, under which d is invariant).  Writing B = lam max_j c_j, the
smoothing of the tropical corners costs O(log k / B) in U = log F, so

    distortion  =  1 + C / B  + O(B^-2),      C ~ 4.5,

and B = log(number of fibres) up to O(1).  The infimum of the C_4 distortion
over signature 4-tuples is therefore exactly 1: FINDINGS Sec. 4.3's 1.2557 is a
stalled search, not an obstruction.

    python research/realizability/i_witness.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import i_certify as CERT  # noqa: E402
import i_logsig as L  # noqa: E402

HERE = Path(__file__).resolve().parent
C4 = np.array([[0., 1., 2., 1.],
               [1., 0., 1., 2.],
               [2., 1., 0., 1.],
               [1., 2., 1., 0.]])
IU = np.triu_indices(4, 1)


def cone_lines():
    """(c, x) arrays of the four cone elements, from the certified data."""
    theta, S, s, _ = CERT.build(CERT.SEED["z"], CERT.SEED["perm"], dps=50)
    cs, xs = [], []
    for a in range(4):
        c, x = CERT.lines(theta, S[a])
        cs.append(np.array([float(v) for v in c]))
        xs.append(np.array([float(v) for v in x]))
    return np.array(cs), np.array(xs), float(s)


def ladder(cs, xs, budgets):
    print(f"  {'B = max log m':>14} {'log10(fibres)':>14} "
          f"{'distortion':>18} {'scale':>10} {'B(dist-1)':>10}")
    rows = []
    for B in budgets:
        lam = B / cs.max()
        sigs = [L.LogSig(lam * cs[a], lam * xs[a]) for a in range(4)]
        D = L.dmatrix(sigs, step=5e-4, lo=-25.0, hi=25.0)
        r = D[IU] / C4[IU]
        dist = float(r.max() / r.min())
        lr = max(s.R for s in sigs) / math.log(10.0)
        print(f"  {B:14.0f} {lr:14.1f} {dist:18.12f} {r.min():10.6f} "
              f"{B*(dist-1):10.4f}", flush=True)
        rows.append({"B": B, "log10_fibres": lr, "distortion": dist,
                     "scale": float(r.min())})
    return rows


def integerise(cs, xs, B, digits=25):
    """Exact integer (multiplicity, atom) data at multiplicity budget B."""
    lam = B / cs.max()
    c = lam * cs
    x = lam * xs
    p = digits * math.log(10.0) / x[x > 0].min()
    from mpmath import mp, mpf, exp as mexp, floor as mfloor
    mp.dps = 80
    fam = []
    for a in range(4):
        rows = []
        for ci, xi in zip(c[a], x[a]):
            m = int(mfloor(mexp(mpf(float(ci)))))
            v = int(mfloor(mexp(mpf(float(p * xi))) + mpf(1) / 2))
            rows.append((max(m, 1), max(v, 1)))
        fam.append(rows)
    return fam, p


def mp_d(fam, dps=40, grid=2001):
    """d matrix of an exact integer family, evaluated at ``dps`` digits."""
    from mpmath import mp, mpf, log, exp
    mp.dps = dps
    data = []
    for rows in fam:
        cs = [log(mpf(m)) for m, _ in rows]
        xs = [log(mpf(v)) for _, v in rows]
        data.append((cs, xs))

    def F(idx, theta):
        cs, xs = data[idx]
        beta = exp(theta)
        vs = [cs[i] + beta * xs[i] for i in range(len(cs))]
        M = max(vs)
        return M + log(sum(exp(v - M) for v in vs))

    lo, hi = mpf(-30), mpf(30)
    ts = [lo + (hi - lo) * i / (grid - 1) for i in range(grid)]
    U = [[log(F(a, t)) for t in ts] for a in range(4)]
    D = [[mpf(0)] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i + 1, 4):
            v = [U[j][k] - U[i][k] for k in range(grid)]
            e0 = (log(log(sum(mpf(m) for m, _ in fam[j])))
                  - log(log(sum(mpf(m) for m, _ in fam[i]))))
            e1 = (log(log(mpf(max(v for _, v in fam[j]))))
                  - log(log(mpf(max(v for _, v in fam[i])))))
            D[i][j] = D[j][i] = max(max(v), e0, e1) - min(min(v), e0, e1)
    return D


def main():
    cs, xs, s = cone_lines()
    print("=== the certified cone witness, as (intercept, slope) data ===")
    print(f"  scale s = {s:.15f}")
    for a in range(4):
        print(f"   Phi_{a}: c = " + ", ".join(f"{v:.4f}" for v in cs[a])
              + " | x = " + ", ".join(f"{v:.4f}" for v in xs[a]))

    print("\n=== distortion of the induced signature families ===")
    rows = ladder(cs, xs, [25, 50, 100, 200, 400, 800, 1600, 3200, 6400,
                           12800, 25600])

    print("\n=== a small integer family, verified at 40 digits ===")
    fam, p = integerise(cs, xs, 25.0, digits=6)
    D = mp_d(fam, dps=40, grid=1201)
    r = [float(D[i][j]) / C4[i][j] for i, j in zip(*IU)]
    print(f"  common exponent p = {p:.3f}")
    for a, rows_a in enumerate(fam):
        print(f"   sig {a}: " + ", ".join(
            f"[{len(str(v))}-digit atom]^({m})" for m, v in rows_a))
    print("  d (40 digits):")
    for i in range(4):
        print("   " + "  ".join(f"{float(D[i][j]):12.9f}" for j in range(4)))
    print(f"  distortion = {max(r)/min(r):.12f}")

    (HERE / "i_witness.json").write_text(json.dumps(
        {"cone_scale": s, "ladder": rows,
         "integer_family_B25": [[[m, v] for m, v in f] for f in fam],
         "integer_family_distortion": max(r) / min(r)}, indent=1))
    print("\n  written: i_witness.json")


if __name__ == "__main__":
    main()
