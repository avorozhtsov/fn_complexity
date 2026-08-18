"""Session brief N, part 2 -- log 2 IS a projective diameter.

Claim (FINDINGS Sec. 4):  the sigma-fibre

    C_0 = { Phi in C : Phi(0) = Lam_Phi }        (equivalently sigma_Phi = 0)

is a subcone of C, and its diameter in the Hilbert (= exchange) metric is
EXACTLY log 2, attained between the two extreme rays

    Phi_min(beta) = max(1, beta)      (the tropical corner, "(r,1,...,1)")
    Phi_max(beta) = 1 + beta          (the flat resource)

so brief I's sharp defect  d <= |Dsigma| + log(1+e^{-|Dsigma|})  is the
statement that C fibres over the sigma-line with fibres of diameter log 2, and
the Birkhoff ratio attached to that diameter is  tanh((log 2)/4).

This script checks, all with the exact tropical machinery of i_cone:

  (1) C_0 is a subcone; the diameter is log 2 and nothing larger;
  (2) the naive order-interval bound is 2 log 2 -- the sharpening is real;
  (3) the two-fibre law  sup{ d : Dsigma = D } = D + log(1+e^{-D}), at 40 digits;
  (4) the Birkhoff constants at 40 digits;
  (5) a COUNTEREXAMPLE to the domain-restricted form of Birkhoff's theorem
      (the step it is easiest to wave at);
  (6) an explicit positive linear map with image in C_0, and the verification
      that it does contract by tanh((log 2)/4);
  (7) the dilation action: sigma is equivariant, orbits are geodesic lines.

    python research/birkhoff/n2_diameter.py
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
PHI_MIN = T.Trop([1.0, 0.0], [0.0, 1.0])      # max(1, beta)
PHI_MAX = T.Trop([1.0], [1.0])                # 1 + beta


def to_fibre(t: T.Trop) -> T.Trop:
    """Dilate beta so that R = Lam, i.e. sigma = 0.  Dilation is an isometry."""
    return T.Trop(t.c, t.x * (t.R / t.Lam))


def rand_trop(rng, k=None, hi=5.0):
    while True:
        kk = k if k is not None else rng.randint(1, 5)
        try:
            return T.Trop([rng.uniform(0.02, hi) for _ in range(kk)],
                          [rng.uniform(0.02, hi) for _ in range(kk)])
        except ValueError:
            continue


def tensor(a, b):
    return T.Trop((a.c[:, None] + b.c[None, :]).ravel(),
                  (a.x[:, None] + b.x[None, :]).ravel())


# ---------------------------------------------------------------------------
def part1(rng, trials=200000):
    print("=== (1) C_0 is a subcone; its Hilbert diameter is log 2 ===")
    bad = 0.0
    for _ in range(2000):
        p, q = to_fibre(rand_trop(rng)), to_fibre(rand_trop(rng))
        s = tensor(p, q)
        bad = max(bad, abs(s.sigma))
    print(f"  Phi, Psi in C_0  =>  Phi + Psi in C_0:  max |sigma| = {bad:.3e}")
    print(f"  extreme pair  d(max(1,beta), 1+beta) = "
          f"{T.hilbert(PHI_MIN, PHI_MAX):.15f}   log 2 = {LOG2:.15f}")
    worst = 0.0
    arg = None
    for _ in range(trials):
        p, q = to_fibre(rand_trop(rng)), to_fibre(rand_trop(rng))
        d = T.hilbert(p, q)
        if d > worst:
            worst, arg = d, (p, q)
    print(f"  {trials} random pairs in C_0: max d = {worst:.12f} "
          f"(= {worst/LOG2:.9f} log 2)")
    # hill-climb, to be sure nothing beats log 2
    def neg_d(z, k):
        z = np.asarray(z, float).reshape(2, k, 2)
        try:
            ts = [to_fibre(T.Trop(np.maximum(r[:, 0], 1e-9),
                                  np.maximum(r[:, 1], 1e-9))) for r in z]
        except ValueError:
            return 1e3
        return -T.hilbert(*ts)
    for k in (2, 3, 5):
        bounds = [(0.0, 9.0)] * (2 * k * 2)
        best = math.inf
        for t in range(4):
            z, f = differential_evolution(neg_d, bounds, args=(k,), seed=11 + 77 * t + k,
                                          maxiter=350, popsize=14, F=(0.3, 1.2), CR=0.9)
            for step in (0.5, 0.1, 0.02, 4e-3, 8e-4, 1.6e-4, 3e-5, 6e-6, 1e-6):
                z, f = pattern_search(neg_d, z, args=(k,), step=step, min_step=1e-13,
                                      maxiter=30000, bounds=bounds)
            best = min(best, f)
        print(f"  hill-climb, k = {k} lines: sup d over C_0 = {-best:.12f}  "
              f"({-best/LOG2:.10f} log 2)")


def part2(rng, trials=50000):
    print("\n=== (2) the naive order-interval bound is 2 log 2 ===")
    print("  Phi in C_0  =>  max(1,beta) <= Phi <= 1+beta pointwise, and the "
          "generic\n  order-interval bound gives diam <= 2 log(sup ratio) = "
          f"{2*LOG2:.9f}.")
    worst = 0.0
    for _ in range(trials):
        p = to_fibre(rand_trop(rng))
        worst = max(worst, T.hilbert(PHI_MIN, p) + T.hilbert(p, PHI_MAX))
    print(f"  max over samples of d(Phi_min,Phi) + d(Phi,Phi_max) = {worst:.9f}"
          f"  (a triangle-inequality witness that the fibre is 'thin')")


def part3():
    print("\n=== (3) the two-fibre law  sup{d : Dsigma = D} = D + log(1+e^-D) ===")
    from mpmath import mp, mpf, log as mplog, exp as mpexp, tanh as mptanh
    mp.dps = 45
    print(f"  {'D':>6} {'exact D + log(1+e^-D)':>44} {'hill-climb in C':>18}")

    def neg_d_at(z, k, D):
        z = np.asarray(z, float).reshape(2, k, 2)
        try:
            a = to_fibre(T.Trop(np.maximum(z[0][:, 0], 1e-9),
                                np.maximum(z[0][:, 1], 1e-9)))
            b = to_fibre(T.Trop(np.maximum(z[1][:, 0], 1e-9),
                                np.maximum(z[1][:, 1], 1e-9)))
        except ValueError:
            return 1e3
        b = T.Trop(b.c, b.x * math.exp(-D))     # push b to sigma = D
        return -T.hilbert(a, b)

    for D in (0.0, 0.25, 0.5, 1.0, 2.0):
        ex = mpf(D) + mplog(1 + mpexp(-mpf(D)))
        best = math.inf
        for k in (2, 3):
            bounds = [(0.0, 9.0)] * (2 * k * 2)
            for t in range(3):
                z, f = differential_evolution(neg_d_at, bounds, args=(k, D),
                                              seed=101 + 313 * t + k,
                                              maxiter=300, popsize=14,
                                              F=(0.3, 1.2), CR=0.9)
                for step in (0.5, 0.1, 0.02, 4e-3, 8e-4, 1.6e-4, 3e-5, 6e-6, 1e-6):
                    z, f = pattern_search(neg_d_at, z, args=(k, D), step=step,
                                          min_step=1e-13, maxiter=30000,
                                          bounds=bounds)
                best = min(best, f)
        print(f"  {D:6.2f} {mp.nstr(ex, 40):>44} {-best:18.12f}")


def part4():
    print("\n=== (4) the Birkhoff constants, 40 digits ===")
    from mpmath import mp, mpf, log as mplog, exp as mpexp, tanh as mptanh
    mp.dps = 45
    l2 = mplog(2)
    print(f"  log 2                     = {mp.nstr(l2, 40)}")
    print(f"  (log 2)/2                 = {mp.nstr(l2/2, 40)}")
    print(f"  tanh((log 2)/4)           = {mp.nstr(mptanh(l2/4), 40)}")
    print(f"  closed form 3 - 2 sqrt 2  = {mp.nstr(3 - 2*mp.sqrt(2), 40)}")
    print(f"  difference                = "
          f"{mp.nstr(mptanh(l2/4) - (3 - 2*mp.sqrt(2)), 10)}")
    print(f"  {'D':>6} {'diam of the pair of fibres':>44} {'tanh(./4)':>44}")
    for D in (0.0, 0.5, 1.0, 2.0, 5.0):
        v = mpf(D) + mplog(1 + mpexp(-mpf(D)))
        print(f"  {D:6.2f} {mp.nstr(v, 40):>44} {mp.nstr(mptanh(v/4), 40):>44}")


def part5(rng):
    print("\n=== (5) why Delta must be the diameter of T(K), not of T(domain) ===")
    print("  If the domain-restricted form were valid one could take K' = C_W =")
    print("  {Phi in C : |sigma| <= W}, a subcone of finite d-diameter")
    print("  <= 2W + log 2, and T = identity, whose Lipschitz constant is 1 but")
    print("  whose 'tanh(diam/4)' would be < 1.  Explicitly, with W = 0:")
    d = T.hilbert(PHI_MIN, PHI_MAX)
    print(f"    diam C_0 = {d:.12f},  tanh(diam/4) = {math.tanh(d/4):.12f},")
    print(f"    but d(Id Phi_min, Id Phi_max)/d(Phi_min,Phi_max) = 1.")
    print("  So Birkhoff's Delta is the diameter of the image of the WHOLE cone")
    print("  whose order defines the metric.  (This is the step brief N warned")
    print("  was easy to wave at.)")


def part6(rng, trials=100000):
    print("\n=== (6) a positive linear map with image in C_0 really does "
          "contract by tanh((log 2)/4) = 3 - 2 sqrt 2 = 0.17157... ===")
    print("  T Phi = Phi(b1) * max(1,beta) + Phi(b2) * (1+beta):  R(T Phi) = "
          "Lam(T Phi)\n  for every Phi, so T(K) is contained in C_0 and "
          "Delta(T) <= log 2.")
    b1, b2 = 0.7, 4.0

    def Tmap(p: T.Trop) -> T.Trop:
        al = float(p.val(np.array([b1]))[0])
        ga = float(p.val(np.array([b2]))[0])
        return T.Trop([al + ga, ga], [ga, al + ga])      # al*max(1,b) + ga*(1+b)

    # Delta(T): diameter of the image = sup over (al,ga), (al',ga') of d
    def neg_img(z):
        al, ga, al2, ga2 = [max(v, 1e-12) for v in z]
        u = T.Trop([al + ga, ga], [ga, al + ga])
        v = T.Trop([al2 + ga2, ga2], [ga2, al2 + ga2])
        return -T.hilbert(u, v)
    bounds = [(0.0, 40.0)] * 4
    best = math.inf
    for t in range(6):
        z, f = differential_evolution(neg_img, bounds, seed=7 + 91 * t,
                                      maxiter=300, popsize=14, F=(0.3, 1.2), CR=0.9)
        for step in (1.0, 0.2, 0.04, 8e-3, 1.6e-3, 3e-4, 6e-5, 1e-5, 2e-6):
            z, f = pattern_search(neg_img, z, step=step, min_step=1e-13,
                                  maxiter=30000, bounds=bounds)
        best = min(best, f)
    Delta = -best
    print(f"  Delta(T) computed = {Delta:.12f}   (log 2 = {LOG2:.12f})")
    print(f"  Birkhoff ratio tanh(Delta/4) = {math.tanh(Delta/4):.12f}")
    worst = 0.0
    for _ in range(trials):
        p, q = rand_trop(rng), rand_trop(rng)
        d0 = T.hilbert(p, q)
        if d0 < 1e-9:
            continue
        worst = max(worst, T.hilbert(Tmap(p), Tmap(q)) / d0)
    print(f"  empirical Lipschitz ratio over {trials} random pairs of C: "
          f"{worst:.12f}   (<= tanh(log2/4) = 0.171573 required)")


def part7(rng, trials=20000):
    print("\n=== (7) the beta-dilation action: sigma is equivariant, orbits are "
          "geodesic lines ===")
    worst_eq = 0.0
    worst_orb = 0.0
    for _ in range(trials):
        p = rand_trop(rng)
        c = math.exp(rng.uniform(-5, 5))
        q = T.Trop(p.c, p.x * c)
        worst_eq = max(worst_eq, abs((q.sigma - p.sigma) + math.log(c)))
        worst_orb = max(worst_orb, abs(T.hilbert(p, q) - abs(math.log(c))))
    print(f"  max |sigma(Phi(c.)) - sigma(Phi) + log c| = {worst_eq:.3e}")
    print(f"  max |d(Phi, Phi(c.)) - |log c||           = {worst_orb:.3e}")
    print("  (so every dilation orbit is an isometric copy of the line, and")
    print("   sigma : (C,d) -> (R,|.|) is 1-Lipschitz with fibres of diameter log 2)")
    worst_qi = (0.0, 0.0)
    for _ in range(trials):
        p, q = rand_trop(rng), rand_trop(rng)
        d = T.hilbert(p, q)
        ds = abs(p.sigma - q.sigma)
        worst_qi = (max(worst_qi[0], ds - d), max(worst_qi[1], d - ds - LOG2))
    print(f"  quasi-isometry check: max(|Dsigma| - d) = {worst_qi[0]:.3e}, "
          f"max(d - |Dsigma| - log 2) = {worst_qi[1]:.3e}")


def main():
    rng = random.Random(20260818)
    part1(rng)
    part2(rng)
    part3()
    part4()
    part5(rng)
    part6(rng)
    part7(rng)


if __name__ == "__main__":
    main()
