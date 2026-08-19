"""G1 -- a certified 3-cycle among *two-fiber* signatures.

The dense sweep of `g1_atoms.py` finds thousands of them; this script extracts
one, maximises its margin, converts it to integers by the exact power symmetry
(a -> a^p is beta -> p*beta, under which A is invariant) and certifies the
three signs against 40-digit mpmath.

Two fibers is the minimum: a one-fiber signature (a,) has u_a(s) = s +
log log a, so u_b - u_a is the constant log(log b / log a), the flow is exact
and every tournament on one-fiber signatures is a total order.

    python research/realizability/g1_r2_witness.py
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from g1_atoms import A_matrix, two_fiber_grid  # noqa: E402
from optimizers import pattern_search  # noqa: E402

HERE = Path(__file__).resolve().parent
GRID = C.make_grid(-14.0, 14.0, 0.004)


def find_cycles(A, limit=200):
    S = A > 0
    n = A.shape[0]
    out = []
    for i in range(n):
        o = np.flatnonzero(S[i])
        inn = np.flatnonzero(S[:, i])
        if len(o) == 0 or len(inn) == 0:
            continue
        sub = S[np.ix_(o, inn)]
        for a, b in np.argwhere(sub):
            t = (i, int(o[a]), int(inn[b]))
            if t[0] < t[1] and t[0] < t[2]:
                out.append(t)
                if len(out) >= limit:
                    return out
    return out


def margin_of(v):
    """v = [x1,d1,x2,d2,x3,d3] -> min of the three A's around the triangle."""
    try:
        ss = [C.Sig.from_logs([max(v[2 * i], 1e-6),
                               max(v[2 * i] - abs(v[2 * i + 1]), 0.0)])
              for i in range(3)]
    except ValueError:
        return -1e3
    As = [C.parts(a, b, GRID)["A"] for a, b in zip(ss, ss[1:] + ss[:1])]
    return min(As)


def neg_margin(v):
    return -margin_of(v)


def main():
    sigs, params = two_fiber_grid(60, 60, 0.25, 4.0)
    A = A_matrix(sigs, C.make_grid(-11.0, 11.0, 0.01))
    cyc = find_cycles(A)
    print(f"two-fiber sweep: {len(cyc)} three-cycles listed (of the full count)")

    # best starting triangle by margin
    best = None
    for (i, j, k) in cyc:
        v = np.array([params[i][0], params[i][1],
                      params[j][0], params[j][1],
                      params[k][0], params[k][1]])
        m = margin_of(v)
        if best is None or m > best[0]:
            best = (m, v)
    print(f"best grid margin: {best[0]:.6e}")

    v, f = pattern_search(neg_margin, best[1], step=0.05, min_step=1e-10,
                          maxiter=40000,
                          bounds=[(1e-3, 12.0), (0.0, 12.0)] * 3)
    m = -f
    print(f"after refinement:  margin = {m:.9e}")
    xs = [(v[2 * i], max(v[2 * i] - abs(v[2 * i + 1]), 0.0)) for i in range(3)]
    for t in xs:
        print(f"   atoms = ({math.exp(t[0]):.9f}, {math.exp(t[1]):.9f})")

    # integerise via the exact power symmetry
    nz = [u for t in xs for u in t if u > 0]
    p = math.log(10 ** 12) / min(nz)
    ints = []
    for t in xs:
        ints.append(tuple(1 if u <= 0 else int(round(math.exp(min(p * u, 2000.0))))
                          for u in t))
    print(f"\npower p = {p:.6f}; integer signatures:")
    for t in ints:
        print(f"   {t}")

    fam = [C.Sig.of(t) for t in ints]
    Ac = np.zeros((3, 3))
    err = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            _, _, A_mp, _, e = C.certified_A_d(fam[i], fam[j], dps=40)
            Ac[i, j], Ac[j, i] = A_mp, -A_mp
            err = max(err, e)
    tri = [Ac[0, 1], Ac[1, 2], Ac[2, 0]]
    is_cycle = all(t > 0 for t in tri) or all(t < 0 for t in tri)
    print(f"\ncertified at 40 digits (max |double - mpmath| = {err:.2e}):")
    print(f"   A(1,2) = {tri[0]:+.12f}")
    print(f"   A(2,3) = {tri[1]:+.12f}")
    print(f"   A(3,1) = {tri[2]:+.12f}")
    print(f"   directed 3-cycle = {is_cycle};  margin = "
          f"{min(abs(t) for t in tri):.6e};  |curl| = {abs(sum(tri)):.9f}"
          f" = sum|A| = {sum(abs(t) for t in tri):.9f}")
    sg = [s.sigma for s in fam]
    ps = [s.psi for s in fam]
    print(f"   sigma = {['%.6f' % t for t in sg]}   spread {max(sg)-min(sg):.6f}")
    print(f"   psi   = {['%.6f' % t for t in ps]}   spread {max(ps)-min(ps):.6f}"
          f"   (bound for a 3-cycle: log 2 = {C.LOG2:.6f})")
    print(f"   mean |A| = {sum(abs(t) for t in tri)/3:.9f}   "
          f"(bound (log2)/2 = {C.LOG2/2:.9f})")

    (HERE / "g1_r2_witness.json").write_text(json.dumps({
        "signatures": [list(map(int, t)) for t in ints],
        "A": [tri[0], tri[1], tri[2]],
        "margin": min(abs(t) for t in tri),
        "mp_error": err, "is_cycle": bool(is_cycle),
        "sigma": sg, "psi": ps,
    }, indent=1))
    with (HERE / "g1_r2_witness.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["signature", "sigma", "psi"])
        for t, s, pv in zip(ints, sg, ps):
            wr.writerow([str(t), f"{s:.9f}", f"{pv:.9f}"])
        wr.writerow([])
        wr.writerow(["edge", "A_certified"])
        for lbl, t in zip(("1->2", "2->3", "3->1"), tri):
            wr.writerow([lbl, f"{t:.12f}"])
    print(f"\nwrote {HERE/'g1_r2_witness.json'}")


if __name__ == "__main__":
    main()
