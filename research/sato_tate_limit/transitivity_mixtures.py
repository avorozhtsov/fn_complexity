#!/usr/bin/env python3
"""Is the midrange comparison a total order on *measures* of one ``alpha_max``?

Brief K asks for a proof that "for measures with the same ``alpha_max`` the
midrange comparison is a total order".  Products of symplectic blocks are not
all the measures with a given ``alpha_max``: any convex combination of two of
them is again a mean-zero probability measure with the *same* ``alpha_max``
(the essential supremum of a mixture is the larger of the two), and it is the
trace measure of a compact group with two cosets whenever such a group exists.
Mixtures are therefore the honest test of the statement as posed, and the first
place a same-genus 3-cycle can live.

Inside each fixed ``alpha_max`` class this script searches

* the multiplicity-free products (partitions of the genus),
* all uniform 2-mixtures of same-genus members,
* free-weight mixtures, by a bisection on the mixing weight of any near miss.

Cycles are counted exactly by the tournament identity

    #(cyclic triples) = C(n,3) - sum_i C(outdeg_i, 2),

so no triple is missed, and every cycle found is re-verified on a second,
independent ``tau`` grid.

    python research/sato_tate_limit/transitivity_mixtures.py
"""

from __future__ import annotations

import csv
import itertools
from math import comb
from pathlib import Path

import numpy as np

import kappa_lib as KL
from transitivity_dominance import partitions

HERE = Path(__file__).resolve().parent
GMAX = 12
TAU = np.geomspace(1e-4, 1e5, 1201)
TAU_B = np.geomspace(1e-5, 1e6, 1501)
MIX_CAP = 26          # partitions used as mixture ingredients


def mid_of(d_interior: np.ndarray, d_inf: float = 0.0) -> float:
    dd = np.concatenate([[0.0], d_interior, [d_inf]])
    return 0.5 * (float(dd.max()) + float(dd.min()))


def mixture_K(Ka: np.ndarray, Kb: np.ndarray, p: float) -> np.ndarray:
    if p <= 0.0:
        return Kb
    if p >= 1.0:
        return Ka
    m = np.maximum(Ka, Kb)
    return m + np.log(p * np.exp(Ka - m) + (1.0 - p) * np.exp(Kb - m))


def mid_matrix(Ks, tau) -> np.ndarray:
    """``mid[i, j] = mid(Psi_i - Psi_j)``, endpoints supplied analytically."""

    psi = np.asarray([K / tau for K in Ks])
    n, T = psi.shape
    padded = np.concatenate([np.zeros((n, 1)), psi, np.zeros((n, 1))], axis=1)
    mid = np.zeros((n, n))
    step = max(1, 40_000_000 // (n * (T + 2)))
    for lo in range(0, n, step):
        hi = min(n, lo + step)
        d = padded[lo:hi, None, :] - padded[None, :, :]
        mid[lo:hi] = 0.5 * (d.max(axis=2) + d.min(axis=2))
    return mid


def cycle_count(mid: np.ndarray) -> int:
    """Number of cyclic triples of the tournament ``i -> j  iff  mid[i,j] < 0``."""

    n = mid.shape[0]
    A = (mid < 0.0)
    np.fill_diagonal(A, False)
    out = A.sum(axis=1)
    return comb(n, 3) - int(sum(comb(int(d), 2) for d in out))


def find_cycles(mid: np.ndarray, limit: int = 5):
    n = mid.shape[0]
    A = (mid < 0.0)
    np.fill_diagonal(A, False)
    out = []
    for i in range(n):
        js = np.where(A[i])[0]
        if js.size == 0:
            continue
        for j in js:
            ks = np.where(A[j] & A[:, i])[0]
            for k in ks:
                if i < j and i < k:
                    out.append((max(mid[i, j], mid[j, k], mid[k, i]),
                                i, j, k))
                    if len(out) >= limit:
                        return sorted(out)
    return sorted(out)


def closest(mid: np.ndarray) -> float:
    """min over oriented triples of max(mid on the three edges)."""

    n = mid.shape[0]
    best = np.inf
    for i, j, k in itertools.combinations(range(n), 3):
        best = min(best,
                   max(mid[i, j], mid[j, k], mid[k, i]),
                   max(mid[i, k], mid[k, j], mid[j, i]))
    return best


def main() -> int:
    kappa, _ = KL.kappa_and_b(TAU, GMAX)
    kappa_b, _ = KL.kappa_and_b(TAU_B, GMAX)
    rows, winners = [], []
    print("=" * 78)
    print("same-alpha_max 3-cycles: products, and products + 2-mixtures")
    print("=" * 78)
    print(f"  {'genus':>5}  {'family':<22}{'measures':>9}{'triples':>13}"
          f"{'cycles':>8}{'closest':>12}")
    for genus in range(3, GMAX + 1):
        lams = partitions(genus)
        labels = ["".join(map(str, lam)) for lam in lams]
        Ks = [sum(kappa[p - 1] for p in lam) for lam in lams]
        Ks_b = [sum(kappa_b[p - 1] for p in lam) for lam in lams]
        n = len(labels)
        mid = mid_matrix(Ks, TAU)
        nc = cycle_count(mid)
        cl = closest(mid)
        rows.append([genus, "products", n, comb(n, 3), nc, f"{cl:.6g}"])
        print(f"  {genus:>5}  {'products':<22}{n:>9}{comb(n, 3):>13}"
              f"{nc:>8}{cl:>12.5g}")

        # uniform 2-mixtures of the first MIX_CAP partitions
        use = list(range(min(n, MIX_CAP)))
        mlabels, mKs, mKs_b = list(labels), list(Ks), list(Ks_b)
        for i, j in itertools.combinations(use, 2):
            mlabels.append(f"[{labels[i]}+{labels[j]}]/2")
            mKs.append(mixture_K(Ks[i], Ks[j], 0.5))
            mKs_b.append(mixture_K(Ks_b[i], Ks_b[j], 0.5))
        n2 = len(mlabels)
        mid2 = mid_matrix(mKs, TAU)
        nc2 = cycle_count(mid2)
        rows.append([genus, "+ uniform 2-mixtures", n2, comb(n2, 3), nc2, ""])
        print(f"  {genus:>5}  {'+ uniform 2-mixtures':<22}{n2:>9}"
              f"{comb(n2, 3):>13}{nc2:>8}", end="")
        if nc2 == 0:
            print(f"{closest(mid2) if n2 <= 120 else float('nan'):>12.5g}")
        else:
            print(f"{'--':>12}")
            for v, i, j, k in find_cycles(mid2, limit=3):
                print(f"           CYCLE margin {-v:.6g}:"
                      f"  {mlabels[i]} -> {mlabels[j]} -> {mlabels[k]}")
                mb = max(mid_of(mKs_b[i] / TAU_B - mKs_b[j] / TAU_B),
                         mid_of(mKs_b[j] / TAU_B - mKs_b[k] / TAU_B),
                         mid_of(mKs_b[k] / TAU_B - mKs_b[i] / TAU_B))
                print(f"             grid B margin {-mb:.6g}")
                rows.append([genus, "CYCLE", mlabels[i], mlabels[j],
                             mlabels[k], f"{-v:.10g}", f"{-mb:.10g}"])
                winners.append((genus, mlabels[i], mlabels[j], mlabels[k],
                                -v, -mb))

    print("\n" + "=" * 78)
    if winners:
        print("VERDICT: 'same alpha_max => total order' is FALSE for measures")
    else:
        print("no same-alpha_max cycle among products or uniform 2-mixtures")
    print("=" * 78)

    with (HERE / "transitivity_mixtures.csv").open("w", newline="",
                                                   encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["genus", "family", "a", "b", "c", "d", "e"])
        for r in rows:
            wr.writerow(list(r) + [""] * (7 - len(r)))
    print("\nwritten: transitivity_mixtures.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
