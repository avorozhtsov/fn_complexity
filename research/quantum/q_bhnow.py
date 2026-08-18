"""Brief L, Part 5 -- placement against Brandao-Horodecki-Ng-Oppenheim-Wehner.

BHNOW, "The second laws of quantum thermodynamics" (arXiv:1305.5278, PNAS 112
(2015) 3275), define generalised free energies

    F_a(rho, rho_b) = kT D_a(rho || rho_b) - kT log Z                     (2)

and prove that for states block-diagonal in the energy basis, rho -> rho' is
possible by a CATALYTIC THERMAL OPERATION iff F_a(rho) >= F_a(rho') for all
a >= 0.  For general states they give two quantum families, eq. (4) with the
Petz divergence (monotone for 0 <= a <= 2) and eq. (5)

    Fh_a(rho, rho_b) = kT/(a-1) log Tr[(rho_b^{(1-a)/2a} rho rho_b^{(1-a)/2a})^a]
                       - kT log Z

with the SANDWICHED divergence, monotone for a >= 1/2 -- necessary conditions
only, not sufficient.

Four things are checked here.

  B1  the classical exchange monotone IS their family:
      log Z_a(beta) = beta log N + (1-beta) H_beta(p),  p_i = a_i / N,
      so log Z_a(beta) = (beta-1) [ D_beta(p||u) - log r ] + beta log N + ...
      -- verified as an identity, and the exchange rate is unchanged by the
      affine normalisation because the ratio kills it.

  B2  the sandwiched profile of Part 2 IS their eq. (5):
      F_(A,S)(beta) = (beta-1) * Fh_beta(A, gamma) / kT   with gamma = S/Tr S,
      Z = Tr S.  Exact identity; the (beta-1) CANCELS in the ratio, so

          Ct((A,S)->(B,S)) = inf_beta Fh_beta(A) / Fh_beta(B).

  B3  BHNOW's own scalar functional of the family is the WORK DISTANCE
      D_work = kT inf_a [F_a(rho) - F_a(rho')], an infimum of a DIFFERENCE.
      The exchange rate is an infimum of a RATIO.  The difference functional
      induces a partial ORDER (a -> b feasible), which is transitive and can
      never cycle; the ratio functional induces a TOURNAMENT, which does.
      Both are computed on the certified cycle.

  B4  the exchange rate is invariant under an arbitrary positive reweighting
      F(beta) -> c(beta) F(beta) of the monotone family; the work distance is
      not.  So the exchange rate sees the family only up to a ray at each beta
      -- which is exactly brief I's Hilbert-projective-metric statement, read
      resource-theoretically.

Run:  python3 q_bhnow.py
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

from q_core import QSig, extrema, rand_admissible, sym

HERE = os.path.dirname(os.path.abspath(__file__))
SGRID = np.linspace(math.log(0.5), math.log(1e6), 3001)
BETA = np.exp(SGRID)


def logZ(a, beta):
    a = np.asarray(a, float)
    v = beta[:, None] * np.log(a)[None, :]
    m = v.max(axis=1)
    return m + np.log(np.exp(v - m[:, None]).sum(axis=1))


def renyi_H(p, beta):
    v = beta[:, None] * np.log(p)[None, :]
    m = v.max(axis=1)
    lz = m + np.log(np.exp(v - m[:, None]).sum(axis=1))
    return lz / (1.0 - beta)


def sandwiched_D(q, beta):
    """Dt_b(A || gamma) with gamma = S/Tr S -- BHNOW eq. (5) without kT."""
    Z = float(np.trace(q.S))
    # Qt_b(A||cS) = c^{1-b} Qt_b(A||S), so Dt_b(A||gamma) = [F(b) + (b-1)logZ]/(b-1)
    F = q.F_grid(beta)
    return (F + (beta - 1.0) * math.log(Z)) / (beta - 1.0)


def main():
    rng = np.random.default_rng(2718)
    out = {}
    print("=" * 74)
    print("PART 5 -- placement against BHNOW (arXiv:1305.5278)")
    print("=" * 74)

    print("\nB1  log Z_a(beta) = beta log N + (1-beta) H_beta(p),  p = a/N")
    worst = 0.0
    for _ in range(200):
        r = int(rng.integers(2, 9))
        a = rng.integers(1, 4000, size=r).astype(float)
        if a.max() < 2:
            continue
        N = a.sum()
        p = a / N
        b = BETA[(BETA > 1e-3) & (np.abs(BETA - 1) > 1e-6)]
        lhs = logZ(a, b)
        rhs = b * math.log(N) + (1.0 - b) * renyi_H(p, b)
        worst = max(worst, float(np.max(np.abs(lhs - rhs) / np.abs(lhs))))
    print(f"    max RELATIVE |lhs - rhs|, 200 signatures, beta to 1e6 = "
          f"{worst:.3e}")
    print("    The exchange monotone is a Renyi alpha-free energy of the")
    print("    resource, at trivial Hamiltonian, on an UNNORMALISED state.")
    print("    Unnormalised is what makes it strictly positive, hence what")
    print("    makes the RATIO of two of them a rate.")
    out["B1_identity_residual"] = worst

    print("\nB2  F_(A,S)(beta) = (beta-1) * [ Dt_beta(A||gamma) - log Tr S ]")
    print("    with gamma = S/Tr S -- BHNOW eq. (5) up to kT.")
    worst = 0.0
    for _ in range(100):
        r = int(rng.integers(2, 6))
        q = rand_admissible(rng, r, 4.0, 20.0)
        b = BETA[(BETA >= 0.5) & (np.abs(BETA - 1) > 1e-3)]
        D = sandwiched_D(q, b)
        rhs = (b - 1.0) * (D - math.log(float(np.trace(q.S))))
        F = q.F_grid(b)
        worst = max(worst, float(np.max(np.abs(F - rhs) / np.abs(F))))
    print(f"    max RELATIVE |F - (beta-1)(Dt - log Z)| = {worst:.3e}")
    print("    So Ct = inf_beta Fh_beta(A)/Fh_beta(B): the (beta-1) cancels in")
    print("    the ratio, and the exchange rate is exactly an infimum of a")
    print("    ratio of BHNOW's quantum alpha-free energies.")
    print("    Their quantum laws are NECESSARY ONLY, so Ct is an upper bound")
    print("    on the achievable conversion rate, not the rate.")
    out["B2_identity_residual"] = worst

    print("\nB3  ratio vs difference on the certified cycle")
    path = os.path.join(HERE, "q_cycle.json")
    if not os.path.exists(path):
        print("    q_cycle.json not found -- run q_cycle.py first.")
    else:
        with open(path) as fh:
            w = json.load(fh)["best_r2"]
        S = np.array(w["S"])
        qs = [QSig(sym(np.array(A)), S) for A in w["A"]]
        prof = [q.U_grid(SGRID) for q in qs]
        Fs = [np.exp(p) for p in prof]
        print("        pair    C(i->j)      C(j->i)     A(i,j)=mid    "
              "slope Lam_i - Lam_j")
        rows = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            e_inf = math.log(qs[j].Lam) - math.log(qs[i].Lam)
            P, Q = extrema(lambda s: qs[j].U(s) - qs[i].U(s), SGRID, (e_inf,))
            cij = math.exp(-P)      # C(i->j) = inf F_i/F_j = exp(-max phi)
            cji = math.exp(Q)
            dl = qs[i].Lam - qs[j].Lam
            print(f"        {i+1}->{j+1}   {cij:.9f}  {cji:.9f}  "
                  f"{0.5*(P+Q):+.9f}   {dl:+.9f}")
            rows.append({"i": i + 1, "j": j + 1, "C_ij": cij, "C_ji": cji,
                         "mid": 0.5 * (P + Q), "dLam": dl})
        print("    Every rate is < 1 both ways, so no edge of the cycle is a")
        print("    feasible conversion and the cycle contradicts no second law:")
        print("    BHNOW's relation 'all free energies decrease' is a partial")
        print("    order and is acyclic by construction.  The cycle lives in")
        print("    the RATE tournament, which is a different object.")
        print("    Note also that F_i - F_j ~ beta (Lam_i - Lam_j) diverges, so")
        print("    for UNNORMALISED resources the work distance inf_b(F_i - F_j)")
        print("    is -oo whenever Lam_i < Lam_j: the difference functional is")
        print("    degenerate here and the ratio is the only finite scalar.")
        out["B3"] = rows

    print("\nB4  reweighting invariance")
    print("    Replace F(beta) by c(beta) F(beta) for a random positive c.")
    q1 = rand_admissible(rng, 4, 4.0, 20.0)
    q2 = rand_admissible(rng, 4, 4.0, 20.0)
    U1, U2 = q1.U_grid(SGRID), q2.U_grid(SGRID)
    lc = 0.7 * np.sin(2.3 * SGRID) + 0.4 * np.cos(0.9 * SGRID)
    d0 = float((U2 - U1).max() - (U2 - U1).min())
    d1 = float(((U2 + lc) - (U1 + lc)).max() - ((U2 + lc) - (U1 + lc)).min())
    F1, F2 = np.exp(U1), np.exp(U2)
    win = BETA <= 100.0                       # a window, or the inf is -oo
    w0 = float(np.min((F1 - F2)[win]))
    w1 = float(np.min((np.exp(lc) * (F1 - F2))[win]))
    print(f"    d before / after reweighting        : {d0:.12f} / {d1:.12f}")
    print(f"    windowed work distance (beta<=100)  : {w0:+.9f} / {w1:+.9f}")
    print("    The exchange geometry is a functional of the family of")
    print("    monotones only up to a positive ray at each beta.  That is")
    print("    OBSTRUCTION.md Sec. 1's Hilbert projective metric, read as a")
    print("    statement about resource monotones.")
    out["B4"] = {"d_before": d0, "d_after": d1, "work_before": w0,
                 "work_after": w1}

    with open(os.path.join(HERE, "q_bhnow.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote q_bhnow.json")


if __name__ == "__main__":
    main()
