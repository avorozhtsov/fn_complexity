"""G1 -- which tournaments does the exchange comparison realise?

For every isomorphism class of tournament on n = 3, 4, 5 vertices (2, 4 and 12
classes) and for a sample on n = 6, 7, search for signatures realising it, at
each atom count r.  Records the smallest r that works, the margin, and an
integer witness certified against 40-digit mpmath.

    python research/realizability/g1_tournaments.py [nmax] [rmax]
"""
from __future__ import annotations

import csv
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
import realize as R  # noqa: E402

HERE = Path(__file__).resolve().parent


# --------------------------------------------------------------------------
# tournament isomorphism classes
# --------------------------------------------------------------------------

def canon(M):
    """Lexicographically least adjacency bitstring over vertex permutations."""
    n = M.shape[0]
    best = None
    for p in itertools.permutations(range(n)):
        bits = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    bits = (bits << 1) | int(M[p[i], p[j]])
        if best is None or bits < best:
            best = bits
    return best


def iso_classes(n):
    """One representative per isomorphism class of tournament on n vertices."""
    pairs = list(itertools.combinations(range(n), 2))
    reps = {}
    for bits in range(1 << len(pairs)):
        M = np.zeros((n, n), dtype=bool)
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                M[i, j] = True
            else:
                M[j, i] = True
        c = canon(M)
        if c not in reps:
            reps[c] = M.copy()
    return list(reps.values())


def score_seq(M):
    return tuple(sorted(int(M[i].sum()) for i in range(M.shape[0])))


def quadratic_residue_tournament(q):
    """The doubly regular tournament on F_q, q = 3 mod 4: i -> j iff j-i is a QR."""
    qr = {(k * k) % q for k in range(1, q)}
    M = np.zeros((q, q), dtype=bool)
    for i in range(q):
        for j in range(q):
            if i != j and (j - i) % q in qr:
                M[i, j] = True
    return M


# --------------------------------------------------------------------------

def attempt(T, r, seed, warm=None, maxiter=180, restarts=2):
    n = T.shape[0]
    x, m = R.realise(T, r=r, seed=seed, maxiter=maxiter, restarts=restarts)
    if warm is not None and warm.size == n * (r - 1):
        # re-run once warm-started from the (r-1)-atom witness, padded
        pad = np.concatenate([warm.reshape(n, r - 1),
                              warm.reshape(n, r - 1)[:, -1:]], axis=1).ravel()
        S = np.where(T, 1.0, -1.0)
        S = np.triu(S, 1) - np.triu(S, 1).T
        bnds = [(0.0, 5.0)] * (n * r)
        from optimizers import pattern_search
        y, _ = pattern_search(R.soft, pad, args=(n, r, S), step=0.2,
                              min_step=1e-6, maxiter=6000, bounds=bnds)
        for step in (0.1, 0.01, 1e-3):
            y, _ = pattern_search(R.margin, y, args=(n, r, S), step=step,
                                  min_step=1e-8, maxiter=6000, bounds=bnds)
        m2 = -R.margin(y, n, r, S)
        if m2 > m:
            x, m = y, m2
    return x, m


def certify(T, x, n, r):
    """Integerise and re-verify the witness at 40 digits."""
    sigs = R._sigs_from(x, n, r)
    ints, p = R.to_integers(sigs)
    A, err = R.certified_matrix(ints)
    ok = all((A[i, j] > 0) == bool(T[i, j])
             for i in range(n) for j in range(n) if i != j)
    iu = np.triu_indices(n, 1)
    cert_margin = float(np.min(np.abs(A[iu])))
    return ints, A, ok, cert_margin, err


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    rmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    rows = []
    witnesses = {}
    for n in range(3, nmax + 1):
        classes = iso_classes(n)
        print(f"\n=== n = {n}: {len(classes)} isomorphism classes, r = 2..{rmax} ===")
        for ci, T in enumerate(classes):
            best = None
            warm = None
            for r in range(2, rmax + 1):
                t0 = time.time()
                x, m = attempt(T, r, seed=1000 * n + ci, warm=warm)
                if m > 0:
                    warm = np.asarray(x)
                got = ""
                if m > 1e-6:
                    ints, A, ok, cm, err = certify(T, x, n, r)
                    got = (f" -> integer witness ok={ok} margin={cm:.3e} "
                           f"mp_err={err:.1e}")
                    if ok and best is None:
                        best = (r, m, cm, ints)
                        witnesses[f"n{n}_c{ci}"] = {
                            "score_seq": list(score_seq(T)),
                            "r": r, "real_margin": m, "integer_margin": cm,
                            "signatures": [list(map(str, t)) for t in ints],
                            "T": T.astype(int).tolist(),
                        }
                print(f"  class {ci:2d} scores {score_seq(T)}  r={r}  "
                      f"margin={m:+.3e}  {time.time()-t0:.0f}s{got}")
                if best is not None:
                    break
            rows.append({
                "n": n, "class": ci, "scores": score_seq(T),
                "min_r": best[0] if best else None,
                "margin": best[1] if best else 0.0,
                "certified_margin": best[2] if best else 0.0,
            })
        done = sum(1 for x in rows if x["n"] == n and x["min_r"])
        print(f"  --- n={n}: {done}/{len(classes)} classes realised")

    with (HERE / "g1_tournaments.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["n", "class", "score_sequence", "min_atoms_r",
                     "search_margin", "certified_margin"])
        for row in rows:
            wr.writerow([row["n"], row["class"], "-".join(map(str, row["scores"])),
                         row["min_r"] or "", f"{row['margin']:.6e}",
                         f"{row['certified_margin']:.6e}"])
    (HERE / "g1_witnesses.json").write_text(json.dumps(witnesses, indent=1))
    print(f"\nwrote {HERE/'g1_tournaments.csv'} and {HERE/'g1_witnesses.json'}")


if __name__ == "__main__":
    main()
