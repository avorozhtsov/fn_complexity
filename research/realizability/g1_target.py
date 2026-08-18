"""G1 -- targeted search for every isomorphism class of tournament.

Random sampling of subsets finds the common classes; the rare ones need a
directed search.  Here a merged pool of signatures is built, its A matrix is
computed once, and for each target class an annealed subset search minimises
the number of mismatched edges (minimised over all relabelings).

    python research/realizability/g1_target.py [n] [pool_per_bucket]
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
from g1_mine import all_classes, perm_index, canon_batch, bucket  # noqa: E402

HERE = Path(__file__).resolve().parent
GRID = C.make_grid(-13.0, 13.0, 0.005)


def merged_pool(per=450):
    specs = [(6, 12), (6, 13), (5, 12), (7, 12), (6, 20), (4, 12), (8, 12),
             (6, 11), (5, 20), (7, 20)]
    sigs = []
    for (r, M) in specs:
        sigs.extend(bucket(r, M, limit=per, seed=r * 991 + M))
    # deduplicate
    seen, out = set(), []
    for s in sigs:
        key = tuple(s.atoms)
        if key not in seen:
            seen.add(key)
            out.append(s)
    return out


def A_matrix(sigs, grid=GRID):
    tab = np.vstack([s.U(grid) for s in sigs])
    eR = np.array([math.log(s.R) for s in sigs])
    eL = np.array([math.log(s.Lam) for s in sigs])
    n = len(sigs)
    A = np.zeros((n, n))
    for i in range(n):
        diff = tab - tab[i]
        hi = np.maximum(diff.max(axis=1), np.maximum(eR - eR[i], eL - eL[i]))
        lo = np.minimum(diff.min(axis=1), np.minimum(eR - eR[i], eL - eL[i]))
        A[i] = 0.5 * (hi + lo)
    np.fill_diagonal(A, 0.0)
    return A


def mismatch(S, idx, T, perms, iu):
    """min over relabelings of the number of edges where S[idx] differs from T."""
    sub = S[np.ix_(idx, idx)]
    best = len(iu[0]) + 1
    for p in perms:
        sp = sub[p][:, p]
        m = int((sp[iu] != T[iu]).sum())
        if m < best:
            best = m
            if m == 0:
                return 0
    return best


def anneal(S, amp, T, n, perms, iu, rng, rounds=60, steps=900, tol=1e-6):
    N = S.shape[0]
    best_idx, best_m = None, len(iu[0]) + 1
    for _ in range(rounds):
        idx = list(rng.choice(N, n, replace=False))
        m = mismatch(S, idx, T, perms, iu)
        temp0 = 2.0
        for t in range(steps):
            temp = temp0 * (1.0 - t / steps) + 1e-3
            k = int(rng.integers(0, n))
            new = int(rng.integers(0, N))
            if new in idx:
                continue
            cand = list(idx)
            cand[k] = new
            sub = amp[np.ix_(cand, cand)]
            if sub[iu].min() <= tol:
                continue
            m2 = mismatch(S, cand, T, perms, iu)
            if m2 <= m or rng.random() < math.exp(-(m2 - m) / temp):
                idx, m = cand, m2
            if m == 0:
                break
        if m < best_m:
            best_idx, best_m = list(idx), m
        if best_m == 0:
            break
    return best_idx, best_m


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 450
    t0 = time.time()
    sigs = merged_pool(per)
    print(f"merged pool: {len(sigs)} signatures")
    A = A_matrix(sigs)
    print(f"A matrix computed in {time.time()-t0:.0f}s")

    S = A > 0
    amp = np.abs(A)
    classes = all_classes(n)
    perms = perm_index(n)
    iu = np.triu_indices(n, 1)
    rng = np.random.default_rng(2024)

    print(f"\n=== n = {n}: {len(classes)} isomorphism classes ===")
    found, failed = {}, []
    for ci, (code, T) in enumerate(sorted(classes.items())):
        t1 = time.time()
        idx, m = anneal(S, amp, T, n, perms, iu, rng)
        scores = tuple(sorted(int(T[i].sum()) for i in range(n)))
        if m == 0:
            fam = [sigs[i] for i in idx]
            Ac = np.zeros((n, n))
            err = 0.0
            for i in range(n):
                for j in range(i + 1, n):
                    _, _, A_mp, _, e = C.certified_A_d(fam[i], fam[j])
                    Ac[i, j], Ac[j, i] = A_mp, -A_mp
                    err = max(err, e)
            code2 = int(canon_batch((Ac > 0)[None, :, :], perms)[0])
            mm = float(np.min(np.abs(Ac[iu])))
            ok = code2 == code
            found[code] = {"scores": list(scores),
                           "signatures": [list(map(int, s.atoms)) for s in fam],
                           "certified_min_absA": mm, "confirmed": ok,
                           "mp_err": err}
            print(f"  class {ci:2d} code {code:<6} scores {scores}  FOUND  "
                  f"min|A|={mm:.3e}  confirmed={ok}  {time.time()-t1:.0f}s")
        else:
            failed.append((code, scores, m))
            print(f"  class {ci:2d} code {code:<6} scores {scores}  "
                  f"not found (best mismatch {m}/{len(iu[0])})  {time.time()-t1:.0f}s")
    print(f"\n{len(found)}/{len(classes)} classes realised and certified")
    if failed:
        print("  not realised in this pool: "
              + ", ".join(f"{c}({'-'.join(map(str,s))}, miss {m})" for c, s, m in failed))

    out = HERE / f"g1_target_n{n}.json"
    out.write_text(json.dumps(
        {"pool_size": len(sigs), "n": n, "classes": len(classes),
         "found": {str(k): v for k, v in found.items()},
         "failed": [{"code": c, "scores": list(s), "best_mismatch": m}
                    for c, s, m in failed]}, indent=1))
    with (HERE / f"g1_target_n{n}.csv").open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["iso_code", "score_sequence", "realised",
                     "certified_min_absA", "signatures"])
        for code, T in sorted(classes.items()):
            sc = "-".join(str(int(T[i].sum())) for i in range(n))
            if code in found:
                f = found[code]
                wr.writerow([code, sc, "yes", f"{f['certified_min_absA']:.3e}",
                             " | ".join(str(tuple(t)) for t in f["signatures"])])
            else:
                wr.writerow([code, sc, "no", "", ""])
    print(f"wrote {out} and {HERE/f'g1_target_n{n}.csv'}")


if __name__ == "__main__":
    main()
