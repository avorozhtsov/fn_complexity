"""Brief L, Part 3 -- a 3-cycle created by coherence alone.

Three admissible quantum signatures (A_1,S), (A_2,S), (A_3,S) that share one
background S.  Comparison is brief D's  a < b  <=>  A(a,b) > 0,  where A is the
MIDRANGE of  phi = U_b - U_a  and  U = log F.  The triple is a strict 3-cycle
when A(1,2), A(2,3), A(3,1) all have the same sign.

We look for a triple which

    * is a strict 3-cycle, and
    * whose PINCHED shadow (A_i -> sum_k P_k A_i P_k, the decoherence of A_i in
      the eigenbasis of S) is a strict *transitive* tournament.

By Theorem Q3 the pinching moves each profile by an amount pinned to zero at
beta = 1 and signed by (beta - 1), so the coherence acts on the triple as three
independent shears about beta = 1.  The question is whether that shear can flip
one edge of the tournament.

The infimum defining the rate is taken over beta in [1/2, oo): that is the
range on which the sandwiched Renyi divergence is a monotone (BHNOW's quantum
second laws hold for alpha >= 1/2), and it is also the range on which the
double-precision pencil (A, S^{(b-1)/b}) is perfectly conditioned.

Run:  python3 q_cycle.py            (default r = 2, then r = 3)
      python3 q_cycle.py --quick    (fewer restarts)
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "realizability"))

from optimizers import differential_evolution, pattern_search  # noqa: E402

from q_core import QSig, osc_mid_fast, sym  # noqa: E402

EDGES = ((0, 1), (1, 2), (2, 0))
SEARCH_S = np.linspace(math.log(0.5), math.log(1000.0), 321)
FINE_S = np.linspace(math.log(0.5), math.log(3000.0), 4001)


def cyc_margin(v):
    """max(min v, min -v):  > 0 iff the three midranges have a common sign."""
    return max(min(v), min(-x for x in v))


def unpack(z, r):
    """z -> (S, [A_1, A_2, A_3]);  S = diag(1 + e^y) >= I,  A_i = S^.5 K_i S^.5."""
    y = z[:r]
    s = 1.0 + np.exp(y)
    S = np.diag(s)
    Sh = np.diag(np.sqrt(s))
    ops = []
    off = r
    for _ in range(3):
        G = z[off:off + r * r].reshape(r, r)
        off += r * r
        K = np.eye(r) + G @ G.T
        ops.append(sym(Sh @ K @ Sh))
    return S, ops


def midranges(qs, sgrid, refine=False):
    prof = [q.U_grid(sgrid) for q in qs]
    return [osc_mid_fast(qs[i], qs[j], sgrid, prof[i], prof[j],
                         include_zero=False, refine=refine)[1]
            for i, j in EDGES]


def evaluate(z, r, sgrid, refine=False):
    """(quantum cycle margin, shadow acyclicity margin)."""
    S, ops = unpack(z, r)
    try:
        qs = [QSig(A, S) for A in ops]
    except (ValueError, np.linalg.LinAlgError):
        return None
    ps = [q.pinched() for q in qs]
    try:
        Aq = midranges(qs, sgrid, refine)
        As = midranges(ps, sgrid, refine)
    except np.linalg.LinAlgError:
        return None
    return cyc_margin(Aq), -cyc_margin(As), Aq, As


def objective(z, r, sgrid):
    got = evaluate(z, r, sgrid)
    if got is None:
        return 10.0
    return -min(got[0], got[1])


def search(r, seed, maxiter, popsize):
    dim = r + 3 * r * r
    bounds = [(-3.0, 2.0)] * r + [(-2.5, 2.5)] * (3 * r * r)
    z, val = differential_evolution(objective, bounds, args=(r, SEARCH_S),
                                    seed=seed, maxiter=maxiter, popsize=popsize)
    z, val = pattern_search(objective, z, args=(r, SEARCH_S), step=0.2,
                            min_step=1e-9, bounds=bounds)
    return z, -val


def report(z, r, tag, out):
    got = evaluate(z, r, FINE_S, refine=True)
    S, ops = unpack(z, r)
    qs = [QSig(A, S) for A in ops]
    mq, ms, Aq, As = got
    print(f"\n--- {tag}:  r = {r} ---")
    print(f"    S = diag({np.diag(S)})")
    print(f"    coherence of A_i in S's basis: "
          f"{[round(q.coherence(), 5) for q in qs]}")
    print(f"    quantum  midranges A(1,2), A(2,3), A(3,1): "
          f"{[f'{x:+.9f}' for x in Aq]}")
    print(f"    shadow   midranges A(1,2), A(2,3), A(3,1): "
          f"{[f'{x:+.9f}' for x in As]}")
    print(f"    quantum cycle margin      = {mq:+.10f}")
    print(f"    shadow acyclicity margin  = {ms:+.10f}")
    verdict = ("QUANTUM 3-CYCLE WITH TRANSITIVE SHADOW"
               if mq > 1e-10 and ms > 1e-10 else "no separation")
    print(f"    verdict: {verdict}")
    out[tag] = {
        "r": r,
        "z": list(map(float, z)),
        "S": [list(map(float, row)) for row in S],
        "A": [[list(map(float, row)) for row in M] for M in ops],
        "midranges_quantum": list(map(float, Aq)),
        "midranges_shadow": list(map(float, As)),
        "quantum_cycle_margin": float(mq),
        "shadow_acyclicity_margin": float(ms),
        "verdict": verdict,
    }
    return mq, ms


def main():
    quick = "--quick" in sys.argv
    maxiter = 120 if quick else 400
    popsize = 12 if quick else 25
    restarts = 3 if quick else 10
    out = {}
    for r in (2, 3):
        best, bestval = None, -math.inf
        for k in range(restarts):
            z, val = search(r, 1000 * r + k, maxiter, popsize)
            print(f"  r={r} restart {k}: margin {val:+.8f}")
            sys.stdout.flush()
            if val > bestval:
                best, bestval = z, val
        report(best, r, f"best_r{r}", out)
    here = os.path.dirname(os.path.abspath(__file__))
    name = "q_cycle.json"
    if "--tag" in sys.argv:
        name = f"q_cycle_{sys.argv[sys.argv.index('--tag') + 1]}.json"
    with open(os.path.join(here, name), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {name}")


if __name__ == "__main__":
    main()
