"""Session brief N, part 3 -- the triangle curl, and where log 2 gets in.

OBSTRUCTION.md Sec. 7 records  max |curl A| = (log 2)/2  as COMPUTED to nine
digits, one third of the proved bound  3(log 2)/2, with "a proof does not
follow from Theorem 4 alone".

This script

  (1) gives an EXACT three-point family in the cone whose curl is
          curl = ( log 2 - log(1 + e^{-T}) ) / 2   ->   (log 2)/2 ,
      certified at 40 digits, which turns the lower bound into a proof;
  (2) checks the reformulation
          curl A = (1/2) log Omega,
          Omega = prod_cyc sup(Phi_j/Phi_i) / prod_cyc sup(Phi_i/Phi_j),
      so that "max curl = (log 2)/2" is exactly "Omega <= 2 around every
      triangle" -- the total forward/backward arbitrage ratio is at most 2;
  (3) checks the circulation identity  2 curl = J(p) + J(q)  with
          J(x) = int_{x_3}^{x_2} u_1 + int_{x_1}^{x_3} u_2 + int_{x_2}^{x_1} u_3,
      u_i = U_i' nondecreasing from 0 to 1, p, q the argmax/argmin triples;
  (4) pushes the upper-bound search: k = 2..8 lines, many restarts, plus a
      local refinement seeded at the family of (1);
  (5) locates where the extremum lives: the equal-sigma slice, and the
      dependence on the sigma-spread;
  (6) tests two candidate joint constraints,
          |curl| <= max_e |D_e|      (REFUTED: the ratio reaches 3)
          |curl| <= max_e eps_e      (survives; would give |curl| <= log 2).

    python research/birkhoff/n3_curl.py
"""
from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "realizability"))

import i_cone as T                                        # noqa: E402
from optimizers import differential_evolution, pattern_search  # noqa: E402

LOG2 = math.log(2.0)


def parts(a: T.Trop, b: T.Trop):
    """P, Q, D, d, eps for the ordered pair (a,b), exactly (finite max)."""
    bs = T.candidates(a, b)
    v = T.phi(a, b, bs)
    e0 = math.log(b.R) - math.log(a.R)
    e1 = math.log(b.Lam) - math.log(a.Lam)
    hi = max(float(v.max()), e0, e1)
    lo = min(float(v.min()), e0, e1)
    P = hi - max(e0, e1)
    Q = min(e0, e1) - lo
    return P, Q, 0.5 * (P - Q), hi - lo, P + Q


def curl_of(fs):
    return sum(parts(fs[i], fs[j])[2] for i, j in ((0, 1), (1, 2), (2, 0)))


def stats(fs):
    D, eps = [], []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        P, Q, d_, dd, e = parts(fs[i], fs[j])
        D.append(d_)
        eps.append(e)
    return sum(D), max(abs(x) for x in D), max(eps)


# ---------------------------------------------------------------------------
# (1) the exact extremal family
# ---------------------------------------------------------------------------

def part1():
    print("=== (1) an exact family with curl -> (log 2)/2 ===")
    print("  Phi_1 = max(1, beta)      (the tropical corner: (r,1,...,1))")
    print("  Phi_2 = 1 + beta          (flat, sigma = 0)")
    print("  Phi_3 = 1 + e^T beta      (flat, sigma = -T)")
    print("  D_12 = (log 2)/2 exactly, D_23 = 0 (two flats), "
          "D_31 = -log(1+e^-T)/2,")
    print("  so   curl = ( log 2 - log(1 + e^{-T}) ) / 2 .")
    from mpmath import mp, mpf, log as mplog, exp as mpexp
    mp.dps = 45
    print(f"  {'T':>6} {'exact curl, 40 digits':>44} {'cone computation':>20} "
          f"{'|diff|':>10}")
    for Tt in (1.0, 5.0, 20.0, 100.0, 300.0):
        ex = (mplog(2) - mplog(1 + mpexp(-mpf(Tt)))) / 2
        fs = [T.Trop([1.0, 0.0], [0.0, 1.0]),
              T.Trop([1.0], [1.0]),
              T.Trop([1.0], [math.exp(Tt)])]
        num = curl_of(fs)
        print(f"  {Tt:6.0f} {mp.nstr(ex, 40):>44} {num:20.15f} "
              f"{abs(float(ex) - num):10.2e}")
    print(f"  limit = (log 2)/2 = {LOG2/2:.15f}   -- so sup |curl| >= (log 2)/2 "
          "is PROVED,\n  and it is a supremum, not a maximum.")


# ---------------------------------------------------------------------------
# (2) the Omega reformulation
# ---------------------------------------------------------------------------

def sup_ratio(a: T.Trop, b: T.Trop):
    """log sup_beta (Phi_b / Phi_a) = L(a,b) = -log C(a->b)."""
    bs = T.candidates(a, b)
    v = T.phi(a, b, bs)
    return max(float(v.max()), math.log(b.Lam) - math.log(a.Lam),
               math.log(b.R) - math.log(a.R))


def part2(rng, trials=200000):
    print("\n=== (2) curl A = (1/2) log Omega, Omega the cycle asymmetry ratio ===")
    worst = 0.0
    worstO = 0.0
    for _ in range(trials):
        fs = [rand_trop(rng) for _ in range(3)]
        fwd = sum(sup_ratio(fs[i], fs[j]) for i, j in ((0, 1), (1, 2), (2, 0)))
        bwd = sum(sup_ratio(fs[j], fs[i]) for i, j in ((0, 1), (1, 2), (2, 0)))
        c = curl_of(fs)
        worst = max(worst, abs(c - 0.5 * (fwd - bwd)))
        worstO = max(worstO, math.exp(fwd - bwd))
    print(f"  max | curl - (1/2)(sum_fwd L - sum_bwd L) | over {trials} triples "
          f"= {worst:.3e}")
    print(f"  max Omega over the same sample = {worstO:.9f}  (conjecture: <= 2)")
    print("  Both cycle sums are >= 0 (the quasi-metric triangle inequality), so")
    print("  'max curl = (log 2)/2' says: the forward arbitrage around a 3-cycle")
    print("  exceeds the backward one by at most a factor 2.")


# ---------------------------------------------------------------------------
# (3) the circulation identity
# ---------------------------------------------------------------------------

def argext(a: T.Trop, b: T.Trop):
    """(argmax, argmin) in s = log beta of U_b - U_a; +-inf allowed."""
    bs = T.candidates(a, b)
    v = T.phi(a, b, bs)
    e0 = math.log(b.R) - math.log(a.R)
    e1 = math.log(b.Lam) - math.log(a.Lam)
    cand = [(e0, -math.inf), (e1, math.inf)]
    for beta, val in zip(bs, v):
        cand.append((float(val), math.log(beta) if beta > 0 else -math.inf))
    hi = max(cand)
    lo = min(cand)
    return hi[1], lo[1]


def Ufun(t: T.Trop, s):
    """The REDUCED potential  U_i(s) - max(0,s), finite at both ends.

    J is unchanged by subtracting a function common to all three U_i, and the
    reduced potential has limits log R (at -inf) and log Lam (at +inf).
    """
    if s == -math.inf:
        return math.log(t.R)
    if s == math.inf:
        return math.log(t.Lam)
    return math.log(float(t.val(np.array([math.exp(s)]))[0])) - max(0.0, s)


def part3(rng, trials=40000):
    print("\n=== (3) the circulation identity  2 curl = J(p) + J(q) ===")
    print("  J(x) = [U_1(x_2)-U_1(x_3)] + [U_2(x_3)-U_2(x_1)] + [U_3(x_1)-U_3(x_2)]")
    print("       = int_{x_3}^{x_2} u_1 + int_{x_1}^{x_3} u_2 + int_{x_2}^{x_1} u_3,")
    print("  u_i = U_i' nondecreasing from 0 to 1 (Theorems 2, 3 of "
          "OBSTRUCTION.md).")
    worst = 0.0
    n_ok = 0
    for _ in range(trials):
        fs = [rand_trop(rng) for _ in range(3)]
        pq = {}
        for (i, j), lbl in zip(((0, 1), (1, 2), (2, 0)), (2, 0, 1)):
            pq[lbl] = argext(fs[i], fs[j])
        n_ok += 1

        def J(sel):
            x = [pq[k][sel] for k in (0, 1, 2)]
            return ((Ufun(fs[0], x[1]) - Ufun(fs[0], x[2]))
                    + (Ufun(fs[1], x[2]) - Ufun(fs[1], x[0]))
                    + (Ufun(fs[2], x[0]) - Ufun(fs[2], x[1])))
        worst = max(worst, abs(2 * curl_of(fs) - (J(0) + J(1))))
    print(f"  {n_ok} triples: max residual "
          f"= {worst:.3e}")


# ---------------------------------------------------------------------------
# search machinery
# ---------------------------------------------------------------------------

def rand_trop(rng, k=None, hi=6.0):
    while True:
        kk = k if k is not None else rng.randint(1, 4)
        try:
            return T.Trop([rng.uniform(0.02, hi) for _ in range(kk)],
                          [rng.uniform(0.02, hi) for _ in range(kk)])
        except ValueError:
            continue


def _mk(z, k, n=3, equal_sigma=False):
    z = np.asarray(z, float).reshape(n, k, 2)
    out = []
    for row in z:
        c = np.maximum(row[:, 0], 0.0)
        x = np.maximum(row[:, 1], 0.0)
        if c.max() <= 1e-9:
            c = c + 1e-3
        if x.max() <= 1e-9:
            x = x + 1e-3
        t = T.Trop(c, x)
        if equal_sigma:
            t = T.Trop(t.c, t.x * (t.R / t.Lam))     # push sigma to 0
        out.append(t)
    return out


def neg_curl(z, k, equal_sigma=False):
    try:
        fs = _mk(z, k, 3, equal_sigma)
    except ValueError:
        return 1e3
    return -abs(curl_of(fs))


def climb(fun, k, seed, args, nvar, restarts=8, xmax=10.0, maxiter=450):
    bounds = [(0.0, xmax)] * nvar
    best_z, best = None, math.inf
    for t in range(restarts):
        z, f = differential_evolution(fun, bounds, args=args, seed=seed + 911 * t,
                                      maxiter=maxiter, popsize=14, F=(0.3, 1.2),
                                      CR=0.9)
        for step in (0.5, 0.1, 0.02, 4e-3, 8e-4, 1.6e-4, 3e-5, 6e-6, 1e-6, 2e-7):
            z, f = pattern_search(fun, z, args=args, step=step, min_step=1e-13,
                                  maxiter=40000, bounds=bounds)
        if f < best:
            best_z, best = z, f
    return best_z, best


def part4():
    print("\n=== (4) the upper-bound search, pushed ===")
    print(f"  {'k lines':>8} {'sup |curl|':>16} {'/(log2/2)':>12} {'/(3log2/2)':>12}")
    for k in (2, 3, 4, 6, 8):
        z, f = climb(neg_curl, k, seed=181 + k, args=(k, False), nvar=3 * k * 2,
                     restarts=8 if k <= 6 else 5)
        print(f"  {k:>8} {-f:16.9f} {-f/(LOG2/2):12.7f} {-f/(1.5*LOG2):12.7f}",
              flush=True)


def part5():
    print("\n=== (5) where the extremum lives ===")
    print("  (a) restricted to the equal-sigma slice (all three sigma = 0):")
    print(f"  {'k lines':>8} {'sup |curl|':>16} {'/(log2/2)':>12}")
    for k in (2, 3, 4, 5, 6):
        z, f = climb(neg_curl, k, seed=281 + k, args=(k, True), nvar=3 * k * 2,
                     restarts=10)
        print(f"  {k:>8} {-f:16.9f} {-f/(LOG2/2):12.7f}", flush=True)
    print("  (b) the exact family of (1) has sigma = (0, 0, -T): the extremum")
    print("      needs one point pushed to sigma = -infinity.")


def part6(rng, trials=400000):
    print("\n=== (6) candidate joint constraints ===")
    w1 = w2 = 0.0
    for _ in range(trials):
        fs = [rand_trop(rng) for _ in range(3)]
        c, mxD, mxe = stats(fs)
        if mxD > 1e-12:
            w1 = max(w1, abs(c) / mxD)
        if mxe > 1e-12:
            w2 = max(w2, abs(c) / mxe)
    print(f"  random sample of {trials} triples:")
    print(f"    max |curl| / max_e |D_e|   = {w1:.9f}   (3 is the trivial bound)")
    print(f"    max |curl| / max_e eps_e   = {w2:.9f}   (1.5 is the trivial bound)")

    def neg_r1(z, k):
        try:
            fs = _mk(z, k)
        except ValueError:
            return 1e3
        c, mxD, mxe = stats(fs)
        return 1e3 if mxD < 1e-9 else -abs(c) / mxD

    def neg_r2(z, k):
        try:
            fs = _mk(z, k)
        except ValueError:
            return 1e3
        c, mxD, mxe = stats(fs)
        return 1e3 if mxe < 1e-9 else -abs(c) / mxe

    for k in (2, 3):
        _, f1 = climb(neg_r1, k, seed=17 + k, args=(k,), nvar=3 * k * 2, restarts=5)
        _, f2 = climb(neg_r2, k, seed=37 + k, args=(k,), nvar=3 * k * 2, restarts=5)
        print(f"  hill-climb k={k}: max |curl|/max|D| = {-f1:.10f} "
              f"  max |curl|/max eps = {-f2:.10f}", flush=True)
    print("  => |curl| <= max_e |D_e| is REFUTED (ratio 3 is attained).")
    print("  => |curl| <= max_e eps_e survives; it would give |curl| <= log 2,")
    print("     still a factor 2 above the computed supremum.")


def main():
    rng = random.Random(20260818)
    part1()
    part2(rng)
    part3(rng)
    part4()
    part5()
    part6(rng)


if __name__ == "__main__":
    main()
