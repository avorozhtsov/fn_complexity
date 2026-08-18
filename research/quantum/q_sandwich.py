"""Brief L, Part 2 -- the sandwiched exchange rate and its structure theorem.

A quantum signature is a pair (A, S) of positive definite operators with

        A >= S >= I                                     (admissibility)

-- the operator form of FINDINGS Sec. 1.1's 'multiplicities >= 1, fiber sizes
>= 1'.  Its profile is

        F(b) = log Tr[ (S^{(1-b)/2b} A S^{(1-b)/2b})^b ],   R = F(0) = log Tr S,
        Lam  = log lambda_max(S^{-1/2} A S^{-1/2}),         sigma = log(R/Lam).

Verified here on random admissible ensembles:

  Q0  F(0+) = log Tr S  and  F(b)/b -> Lam                       (endpoints)
  Q1  max(R, b Lam) <= F(b) <= R + b Lam                 (PROVED in FINDINGS.md)
  Q2  F increasing, F convex, F - b*Lam nonincreasing (F' <= Lam),
      U' in (0,1), U convex, w unimodal with peak at sigma, height <= log 2
  Q3  pinching: F_(A,S) - F_(P_S A, S) has the sign of (b-1) on [1/2, oo),
      and vanishes at b = 1 and b = 0                            (PROVED; check)
  Q4  non-spectrality: isospectral (A,S) pairs with different profiles and
      different exchange rates -- so Ct does NOT factor through the spectra.
  Q5  the pinched (commuting) profile is exactly a classical cone point.

Run:  python3 q_sandwich.py
"""
from __future__ import annotations

import math

import numpy as np

from q_core import QSig, osc_mid_fast, rand_admissible, sym

LABEL = {
    "low_R": "F >= R",
    "low_bLam": "F >= beta*Lam",
    "up": "F <= R + beta*Lam",
    "F_incr": "F nondecreasing",
    "G_decr": "F - beta*Lam nonincreasing  (F' <= Lam)",
    "F_conv": "F convex in beta",
    "U_slope_lo": "U' >= 0",
    "U_slope_hi": "U' <= 1",
    "U_conv": "U convex in s",
    "w_lo": "w >= 0",
    "w_env": "w <= log(1+e^{-|s-sigma|})",
}
ORDER = list(LABEL)


def widest_grid(q, n=3001, bmax=3000.0):
    """Log-beta grid down to the smallest beta this (A,S) can do in double."""
    bmin = q.beta_safe * 1.02
    return np.linspace(math.log(bmin), math.log(bmax), n)


def profile_checks(q, sgrid):
    """Worst violation of each structural claim (negative = claim holds)."""
    beta = np.exp(sgrid)
    F = q.F_grid(beta)
    R, Lam = q.R, q.Lam
    out = {}
    out["low_R"] = float(np.max(R - F))
    out["low_bLam"] = float(np.max(beta * Lam - F))
    out["up"] = float(np.max(F - (R + beta * Lam)))
    out["F_incr"] = float(-np.min(np.diff(F)))
    out["G_decr"] = float(np.max(np.diff(F - beta * Lam)))
    sl = np.diff(F) / np.diff(beta)
    out["F_conv"] = float(-np.min(np.diff(sl)))
    U = np.log(F)
    slu = np.diff(U) / np.diff(sgrid)
    out["U_slope_lo"] = float(-np.min(slu))
    out["U_slope_hi"] = float(np.max(slu) - 1.0)
    out["U_conv"] = float(-np.min(np.diff(slu)))
    sigma = math.log(R / Lam)
    w = U - math.log(Lam) - np.maximum(sigma, sgrid)
    env = np.log1p(np.exp(-np.abs(sgrid - sigma)))
    out["w_lo"] = float(-np.min(w))
    out["w_env"] = float(np.max(w - env))
    out["_peak"] = float(abs(sgrid[int(np.argmax(w))] - sigma))
    out["_height"] = float(np.max(w))
    out["_sigma_in"] = float(sgrid[0] <= sigma <= sgrid[-1])
    return out


def ensemble(rng, n, pinch=False):
    for _ in range(n):
        r = int(rng.integers(2, 7))
        q = rand_admissible(rng, r, rng.uniform(1.2, 8.0), rng.uniform(1.5, 60.0))
        yield q.pinched() if pinch else q


def main():
    rng = np.random.default_rng(20260818)

    print("=" * 74)
    print("PART 2 -- the sandwiched (background-relative) exchange rate")
    print("=" * 74)
    print("Admissible ensemble  A >= S >= I,  r = 2..6,  600 draws.")
    print("Each draw is checked on its own widest double-precision beta grid,")
    print("[beta_safe, 3000], 3001 log-spaced points.\n")

    worst = {k: -math.inf for k in ORDER + ["_peak", "_height"]}
    coh, sig_in = [], 0
    for q in ensemble(rng, 600):
        g = widest_grid(q)
        res = profile_checks(q, g)
        sig_in += int(res["_sigma_in"])
        coh.append(q.coherence())
        for k in worst:
            worst[k] = max(worst[k], res[k])

    print("Q1-Q2  worst violation over the ensemble (negative = claim holds)")
    for k in ORDER:
        print(f"    {LABEL[k]:<42s} {worst[k]:+.3e}")
    print(f"    max bump height w                          "
          f"{worst['_height']:.6f}   (log 2 = {math.log(2):.6f})")
    print(f"    max |argmax w - sigma|                     {worst['_peak']:.3e}")
    print(f"    draws whose sigma lay inside the grid      {sig_in}/600")
    print(f"    coherence sampled: median {np.median(coh):.4f}, max {max(coh):.4f}")

    print("\n    Noise floor for the two second-difference statistics: the same")
    print("    numbers on the PINCHED profiles, where convexity is a theorem.")
    rng2 = np.random.default_rng(20260818)
    wc = wu = -math.inf
    for q in ensemble(rng2, 200, pinch=True):
        res = profile_checks(q, widest_grid(q))
        wc = max(wc, res["F_conv"])
        wu = max(wu, res["U_conv"])
    print(f"    F convex (pinched)  {wc:+.3e}      U convex (pinched)  {wu:+.3e}")

    print("\nQ0  the two endpoints, on 5 fresh draws")
    print("      log Tr S      F(beta_safe)     F(2000)/2000        Lam")
    for _ in range(5):
        q = rand_admissible(rng, 4, 5.0, 25.0)
        bs = q.beta_safe * 1.02
        print(f"    {q.R:12.8f}  {q.F(bs):12.8f}  "
              f"{q.F(2000.0)/2000.0:14.10f}  {q.Lam:12.10f}   (beta_safe "
              f"{bs:.4f})")

    print("\nQ3  pinching: sign(F_(A,S) - F_(P_S A,S)) = sign(beta - 1)?")
    bad_lo = bad_hi = at_one = 0.0
    ncheck = 0
    for q in ensemble(rng, 300):
        p = q.pinched()
        b = np.exp(widest_grid(q))
        b = b[b >= 0.5]
        D = q.F_grid(b) - p.F_grid(b)
        lo, hi = b < 1.0, b > 1.0
        if lo.any():
            bad_lo = max(bad_lo, float(np.max(D[lo])))
        bad_hi = max(bad_hi, float(-np.min(D[hi])))
        at_one = max(at_one, abs(q.F(1.0) - p.F(1.0)))
        ncheck += 1
    print(f"    worst positive value of F-F_pinched on 1/2<=beta<1 : {bad_lo:+.3e}")
    print(f"    worst negative value of F-F_pinched on beta>1      : {-bad_hi:+.3e}")
    print(f"    max |F(1) - F_pinched(1)|                          : {at_one:.3e}")
    print("    (both are log Tr A exactly; the shear is pinned at beta = 1)")

    print("\nQ4  non-spectrality.  Fix S and rotate A by a random O(r):")
    print("    spec(A) and spec(S) unchanged; the profile is not.")
    q0 = rand_admissible(rng, 4, 4.0, 20.0)
    S, A = q0.S, q0.A
    base = QSig(A, S)
    g = widest_grid(base)
    Ub = base.U_grid(g)
    print("      max|dU|       d(base,rot)   Ct(base->rot)  Ct(rot->base)")
    for _ in range(6):
        O, _ = np.linalg.qr(rng.normal(size=(4, 4)))
        qr_ = QSig(sym(O @ A @ O.T), S)
        Ur = qr_.U_grid(g)
        dd, _ = osc_mid_fast(base, qr_, g, Ub, Ur, include_zero=False)
        c_ab = math.exp(min(float((Ub - Ur).min()),
                            math.log(base.Lam) - math.log(qr_.Lam)))
        c_ba = math.exp(min(float((Ur - Ub).min()),
                            math.log(qr_.Lam) - math.log(base.Lam)))
        print(f"      {np.abs(Ub-Ur).max():.6e}  {dd:.6e}  "
              f"{c_ab:.10f}  {c_ba:.10f}")
    print("    Every row has the same two spectra.  A nonzero d is a quantity")
    print("    the spectral rate of Part 1 cannot see.")

    print("\nQ5  the commuting case really is classical.")
    print("    For [A,S] = 0, F(b) = log sum_i m_i x_i^b with m_i = s_i,")
    print("    x_i = A_ii/s_i -- a cone point of OBSTRUCTION.md Theorem 1.")
    worst5 = 0.0
    for q in ensemble(rng, 200):
        p = q.pinched()
        m, x = p.classical_atoms()
        b = np.exp(widest_grid(p))
        lm, lx = np.log(m), np.log(x)
        E = lm[None, :] + b[:, None] * lx[None, :]
        mx = E.max(axis=1)
        cls = mx + np.log(np.exp(E - mx[:, None]).sum(axis=1))
        worst5 = max(worst5, float(np.max(np.abs(cls - p.F_grid(b)))))
    print(f"    max |F_pinched - classical log-sum-exp| = {worst5:.3e}")


if __name__ == "__main__":
    main()
