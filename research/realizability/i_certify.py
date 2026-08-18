"""Session brief I -- the C_4 witness, certified at 40 digits.

The witness is a point of the cone

    C = { Phi convex nondecreasing on [0,inf), Phi > 0, Phi(beta) >= Lam beta },

which is exactly the projective closure of { F_a = log Z_a } (i_cone.py,
i_validate_cone.py).  Each Phi_a is a max of five lines c + x beta, given by
four common breakpoints theta_1 < ... < theta_4 and, on each of the five cells,
a value of the scale function S_a = log((Phi - beta Phi')/Phi').  The data are
reconstructed here IN MPMATH AT 40 DIGITS from the seven numbers
(L_1, L_2, L_3, rho_1, rho_2, rho_3, s) by solving

    delta_{L_k}( S - theta_k ) = u_{a,k},   delta_L(v) = sp(L-v) - sp(-v)

with a 200-step bisection at 60 digits, and then

    d(a,b) = osc over {breakpoints} u {0, infinity} of log(Phi_b/Phi_a)

is evaluated exactly (Phi_b/Phi_a is monotone between consecutive breakpoints,
being a ratio of affine functions).  The certificate reports

  * max_ij | d_ij - s * C4_ij |                      (should be ~1e-35)
  * the margin of every strict inequality used       (should be >= 1e-6)

    python research/realizability/i_certify.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from mpmath import mp, mpf, log, exp, log1p

HERE = Path(__file__).resolve().parent

# the feasible seven-tuple found by i_pattern.py for the theta-ordering
# (beta+, gamma+, gamma-, beta-)  =  A B D C,  refined by compass search
SEED = {
    "perm": [0, 1, 3, 2],
    "z": [1.345131336145283, 0.3294966044895359, 0.6637440044070804,
          0.21189697313675845, 0.6668277893188684, 1.9278507551483397,
          -1.5995856145379987],
}

T0 = [[0, 0, 0, 0],
      [1, 0, 0, 1],
      [2, 1, 0, 1],
      [1, 1, 0, 0]]

C4 = [[0, 1, 2, 1],
      [1, 0, 1, 2],
      [2, 1, 0, 1],
      [1, 2, 1, 0]]


def sp(t):
    return log1p(exp(t)) if t < 0 else t + log1p(exp(-t))


def delta(v, L):
    return sp(L - v) - sp(-v)


def delta_inv(u, L):
    lo, hi = mpf(-1), mpf(1)
    while delta(lo, L) < u:
        lo -= 2 + abs(lo)
    while delta(hi, L) > u:
        hi += 2 + abs(hi)
    for _ in range(400):
        mid = (lo + hi) / 2
        if delta(mid, L) > u:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def build(z, perm, dps=60):
    mp.dps = dps
    Tm = [[T0[a][p] for p in perm] for a in range(4)]
    L = [exp(mpf(repr(z[i]))) for i in range(3)]
    rho = [mpf(repr(z[3 + i])) for i in range(3)]
    s = exp(mpf(repr(z[6])))
    theta = [mpf(0), L[0], L[0] + L[1], L[0] + L[1] + L[2]]
    S = [[None] * 5 for _ in range(4)]
    margins = []
    for k in range(3):
        for a in range(4):
            u = rho[k] + s * (Tm[a][k + 1] - Tm[a][k])
            margins.append(("u in (0,L)", float(min(u, L[k] - u))))
            S[a][k + 1] = theta[k] + delta_inv(u, L[k])
    top = max(S[a][1] for a in range(4)) + mpf(1) / 2
    bot = min(S[a][3] for a in range(4)) - mpf(1) / 2
    for a in range(4):
        S[a][0] = top
        S[a][4] = bot
    return theta, S, s, margins


def lines(theta, Sa):
    """(c_j, x_j) of Phi = max_j (c_j + x_j beta) from centres and breakpoints."""
    k = len(Sa)
    x = [mpf(1)] * k
    for j in range(k - 1):
        b = exp(theta[j])
        x[j + 1] = x[j] * (exp(Sa[j]) + b) / (exp(Sa[j + 1]) + b)
    c = [x[j] * exp(Sa[j]) for j in range(k)]
    return c, x


def phi_val(c, x, beta):
    return max(c[j] + x[j] * beta for j in range(len(c)))


def dmat(theta, S):
    data = [lines(theta, S[a]) for a in range(4)]
    pts = [exp(t) for t in theta]
    D = [[mpf(0)] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(i + 1, 4):
            ci, xi = data[i]
            cj, xj = data[j]
            vals = [log(phi_val(cj, xj, b)) - log(phi_val(ci, xi, b))
                    for b in pts]
            vals.append(log(max(cj)) - log(max(ci)))          # beta = 0
            vals.append(log(max(xj)) - log(max(xi)))          # beta = infinity
            D[i][j] = D[j][i] = max(vals) - min(vals)
    return D, data


def main():
    theta, S, s, margins = build(SEED["z"], SEED["perm"])
    mp.dps = 40
    D, data = dmat(theta, S)
    err = max(abs(D[i][j] - s * C4[i][j]) for i in range(4) for j in range(4))
    print("=== C_4 witness in the cone, certified at 40 digits ===")
    print(f"  scale s = {mp.nstr(s, 30)}")
    print("  breakpoints theta = " + ", ".join(mp.nstr(t, 20) for t in theta))
    for a in range(4):
        print(f"   S_{a} = " + ", ".join(mp.nstr(v, 20) for v in S[a]))
    print("\n  d matrix (40 digits):")
    for i in range(4):
        print("   " + "  ".join(mp.nstr(D[i][j], 22) for j in range(4)))
    print(f"\n  max_ij | d_ij - s * (C_4)_ij |  =  {mp.nstr(err, 6)}")

    print("\n  strict-inequality margins that the construction uses:")
    mono = []
    for a in range(4):
        for k in range(4):
            mono.append(float(S[a][k] - S[a][k + 1]))
    print(f"   min over a,k of  S_(a,k) - S_(a,k+1)   = {min(mono):.6e}"
          "   (>= 0 required; 0 = two lines merge, still legal)")
    print(f"   min over a,k of the nondegenerate ones = "
          f"{min(v for v in mono if v > 1e-12):.6e}")
    print(f"   min over a,k of  min(u, L-u)           = "
          f"{min(m[1] for m in margins):.6e}")
    print(f"   min node gap                           = "
          f"{float(min(theta[i+1]-theta[i] for i in range(3))):.6e}")
    print(f"   min pairwise distance d_ij             = "
          f"{float(min(D[i][j] for i in range(4) for j in range(4) if i != j)):.6e}")

    out = {
        "scale": mp.nstr(s, 34),
        "theta": [mp.nstr(t, 34) for t in theta],
        "S": [[mp.nstr(v, 34) for v in row] for row in S],
        "c": [[mp.nstr(v, 34) for v in data[a][0]] for a in range(4)],
        "x": [[mp.nstr(v, 34) for v in data[a][1]] for a in range(4)],
        "d": [[mp.nstr(D[i][j], 34) for j in range(4)] for i in range(4)],
        "max_abs_error_vs_sC4": mp.nstr(err, 6),
    }
    (HERE / "i_certify.json").write_text(json.dumps(out, indent=1))
    print("\n  written: i_certify.json")

    # independent cross-check: a dense grid in log beta, double precision,
    # sharing no code with the exact evaluation above
    import numpy as np
    cc = np.array([[float(v) for v in data[a][0]] for a in range(4)])
    xx = np.array([[float(v) for v in data[a][1]] for a in range(4)])
    th = np.linspace(-40.0, 40.0, 4000001)
    bb = np.exp(th)
    Y = np.array([np.log(np.max(cc[a][None, :] + np.outer(bb, xx[a]), axis=1))
                  for a in range(4)])
    Dg = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            v = Y[j] - Y[i]
            e0 = math.log(cc[j].max()) - math.log(cc[i].max())
            e1 = math.log(xx[j].max()) - math.log(xx[i].max())
            Dg[i, j] = Dg[j, i] = (max(v.max(), e0, e1)
                                   - min(v.min(), e0, e1))
    C4a = np.array(C4, dtype=float)
    iu = np.triu_indices(4, 1)
    rr = Dg[iu] / C4a[iu]
    print("\n  independent double-precision cross-check on a 4e6-point grid:")
    print(f"   max |d - s C_4| = {np.abs(Dg - float(s)*C4a).max():.3e}"
          "   (grid resolution, not an error in the witness)")
    print(f"   distortion      = {rr.max()/rr.min():.10f}")


if __name__ == "__main__":
    main()
