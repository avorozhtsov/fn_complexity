"""G1 -- every tournament, by separating the potential from the defect.

The Cartesian power a -> a^{(x)k} multiplies both r and max by exponent k, so

    sigma unchanged,  w unchanged,  log Lam -> log Lam + log k,
    psi -> psi + log k,   d unchanged,   D unchanged,

hence   A(a^{(x)j}, b^{(x)k}) = A(a,b) + log(k/j)   exactly.  (This is FINDINGS
T1.5 Theorem B refined: the power acts on the *potential only*.)

So the realisable tournaments are exactly

    { sign(c_b - c_a + D_ab) : D a realisable defect matrix, c in R^n },

because log(k/j) over positive integers is dense in R.  Realising a tournament
therefore splits into two independent steps:

  1. find a subset of the pool whose defect matrix D admits a potential c with
     sign(c_b - c_a + D_ab) = T_ab   -- a linear feasibility problem in c;
  2. realise c by Cartesian powers k_a ~ e^{c_a} N.

Two-fiber signatures are used as the base, because (u,v)^{(x)k} has only k+1
distinct values (binomial multiplicities), so the witnesses stay small enough
to certify exactly at 40 digits.

    python research/realizability/g1_potential.py [n]
"""
from __future__ import annotations

import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from g1_atoms import A_matrix, two_fiber_grid  # noqa: E402
from g1_mine import all_classes, perm_index, canon_batch  # noqa: E402

HERE = Path(__file__).resolve().parent


def min_mean_cycle(W):
    """Minimum mean weight over directed cycles of the digraph with weights W
    (W[i,j] finite when the edge i->j is present).  Karp's algorithm.

    The realisability criterion: the constraints  c_j - c_i >= m - D_ij  along
    the edges of the tournament T are feasible iff every directed cycle of T
    has  sum D_e >= m * (cycle length).  So the largest achievable margin is
    exactly the minimum mean cycle weight of D over T's directed cycles, and
    T is realisable from D iff that minimum is > 0.
    """
    n = W.shape[0]
    INF = np.inf
    F = np.full((n + 1, n), INF)
    F[0, :] = 0.0                      # allow starting anywhere
    for k in range(1, n + 1):
        prev = F[k - 1]
        M = prev[:, None] + W
        F[k] = M.min(axis=0)
    best = INF
    for v in range(n):
        if not np.isfinite(F[n, v]):
            continue
        worst = -INF
        for k in range(n):
            if np.isfinite(F[k, v]):
                worst = max(worst, (F[n, v] - F[k, v]) / (n - k))
        best = min(best, worst)
    return best


def lp_margin(Dsub, S, iters=0):
    """Largest achievable margin, exactly, plus a potential attaining it.

    T is read off S:  edge i->j present iff S[i,j] > 0.  Weight of i->j is
    D_ij (== -D_ji).  See min_mean_cycle.
    """
    n = Dsub.shape[0]
    W = np.where(S > 0, Dsub, np.inf)
    np.fill_diagonal(W, np.inf)
    m = min_mean_cycle(W)
    if not np.isfinite(m):
        m = float(np.min(np.abs(Dsub[np.triu_indices(n, 1)])))
    if m <= 0:
        return float(m), np.zeros(n)
    # Bellman-Ford potential for the feasible system c_j - c_i >= 0.5*m - D_ij
    mm = 0.5 * m
    c = np.zeros(n)
    for _ in range(4 * n + 40):
        upd = False
        for i in range(n):
            for j in range(n):
                if i != j and S[i, j] > 0 and c[j] < c[i] + mm - Dsub[i, j]:
                    c[j] = c[i] + mm - Dsub[i, j]
                    upd = True
        if not upd:
            break
    c -= c.mean()
    M = S * (c[None, :] - c[:, None] + Dsub)
    return float(M[np.triu_indices(n, 1)].min()), c


def anneal(D, T, n, rng, rounds=25, steps=1500):
    N = D.shape[0]
    S = np.where(T, 1.0, -1.0)
    S = np.triu(S, 1) - np.triu(S, 1).T
    best = (-math.inf, None, None)
    for _ in range(rounds):
        idx = list(rng.choice(N, n, replace=False))
        m, c = lp_margin(D[np.ix_(idx, idx)], S)
        for t in range(steps):
            k = int(rng.integers(0, n))
            new = int(rng.integers(0, N))
            if new in idx:
                continue
            cand = list(idx)
            cand[k] = new
            m2, c2 = lp_margin(D[np.ix_(cand, cand)], S)
            temp = 0.02 * (1 - t / steps) + 1e-4
            if m2 > m or rng.random() < math.exp((m2 - m) / temp):
                idx, m, c = cand, m2, c2
            if m > best[0]:
                best = (m, list(idx), c.copy())
        if best[0] > 5e-3:
            break
    return best


def rationalise(c, kmax=250):
    """Integers k_a with log(k_a/k_b) ~ c_a - c_b."""
    best = None
    for N in range(1, 1200):
        ks = [max(1, int(round(N * math.exp(t - max(c))))) for t in c]
        if max(ks) > kmax or min(ks) < 1:
            continue
        lg = np.array([math.log(k) for k in ks])
        err = np.max(np.abs((lg - lg.mean()) - (c - c.mean())))
        if best is None or err < best[0]:
            best = (err, ks)
    return best


def power(sig, k):
    """Cartesian power of a two-value signature: k+1 values, binomial mults."""
    if k == 1:
        return sig
    xs, ms = list(sig.xs), list(sig.mults)
    vals, mults = {0.0: 1.0}, None
    cur = {0.0: 1.0}
    for _ in range(k):
        nxt = {}
        for x0, m0 in cur.items():
            for x1, m1 in zip(xs, ms):
                key = round(x0 + x1, 12)
                nxt[key] = nxt.get(key, 0.0) + m0 * m1
        cur = nxt
    ks = sorted(cur, reverse=True)
    return C.Sig.from_logs(ks, [cur[t] for t in ks])


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    sigs, params = two_fiber_grid(55, 55, 0.25, 4.0)
    t0 = time.time()
    A = A_matrix(sigs, C.make_grid(-11.0, 11.0, 0.01))
    ps = np.array([s.psi for s in sigs])
    D = A - (ps[None, :] - ps[:, None])
    print(f"two-fiber pool: {len(sigs)} signatures, D computed in "
          f"{time.time()-t0:.0f}s;  max |D| = {np.abs(D).max():.6f}")

    classes = all_classes(n)
    perms = perm_index(n)
    iu = np.triu_indices(n, 1)
    rng = np.random.default_rng(4242)
    rows, found = [], 0
    print(f"\n=== n = {n}: {len(classes)} isomorphism classes, two-fiber bases "
          f"+ Cartesian powers ===")
    for ci, (code, T) in enumerate(sorted(classes.items())):
        t1 = time.time()
        m, idx, c = anneal(D, T, n, rng)
        scores = "-".join(str(int(T[i].sum())) for i in range(n))
        if m <= 0 or idx is None:
            print(f"  class {ci:2d} code {code:<6} scores {scores}  "
                  f"NOT FOUND (best LP margin {m:.2e})  {time.time()-t1:.0f}s")
            rows.append([code, scores, "no", "", ""])
            continue
        err, ks = rationalise(c)
        fam = [power(sigs[i], k) for i, k in zip(idx, ks)]
        Ac = np.zeros((n, n))
        e = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                _, _, A_mp, _, ee = C.certified_A_d(fam[i], fam[j])
                Ac[i, j], Ac[j, i] = A_mp, -A_mp
                e = max(e, ee)
        code2 = int(canon_batch((Ac > 0)[None, :, :], perms)[0])
        mm = float(np.min(np.abs(Ac[iu])))
        ok = code2 == code
        found += ok
        base = [tuple(round(math.exp(x)) for x in sigs[i].xs) for i in idx]
        print(f"  class {ci:2d} code {code:<6} scores {scores}  "
              f"{'OK ' if ok else 'MISS'}  LP margin={m:.2e}  "
              f"certified min|A|={mm:.3e}  rational err={err:.1e}  "
              f"powers={ks}  {time.time()-t1:.0f}s")
        rows.append([code, scores, "yes" if ok else "no", f"{mm:.3e}",
                     json.dumps({"bases": [list(b) for b in base],
                                 "powers": list(map(int, ks))})])
    print(f"\n{found}/{len(classes)} isomorphism classes realised and certified")

    with (HERE / f"g1_potential_n{n}.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["iso_code", "score_sequence", "realised",
                     "certified_min_absA", "witness"])
        for row in rows:
            wr.writerow(row)
    print(f"wrote {HERE/f'g1_potential_n{n}.csv'}")


if __name__ == "__main__":
    main()
