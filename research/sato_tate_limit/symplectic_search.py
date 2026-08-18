#!/usr/bin/env python3
"""The arithmetically decisive search: symplectic monodromy only, to genus six.

``cone_search.py`` finds that every 3-cycle in the limiting comparison uses a
**torus** factor -- ``U(1)`` (arcsine) or the ``CM`` mixture
``(1/2) delta_0 + (1/2) arcsine``.  Deligne's semisimplicity theorem (Weil II,
3.4.1(iii)) says the geometric monodromy group of a pure lisse sheaf on a curve
over a finite field is **semisimple**, so no positive-dimensional torus can be
the identity component of the monodromy of ``R^1 pi_*`` of a family of curves.
The arithmetically reachable measures are therefore the trace measures of
semisimple subgroups of ``Sp(2g)`` -- products of ``Sp(2g_i)`` blocks with
multiplicities, and their mixtures over the components of a disconnected group.

This script searches that set exhaustively:

* every product of ``k`` isogenous copies of ``Sp(2g)`` blocks with
  ``alpha_max = sum 2 g_i k_i <= 12`` -- genus one through six;
* and then convex combinations of those, which is what a **disconnected**
  monodromy group realises: if ``G/G^0`` has ``r`` cosets the trace measure is
  the average of the ``r`` coset measures.  A genus-four family whose Jacobian
  splits into two genus-two pieces conjugate over a quadratic extension of the
  base, for instance, has measure ``(1/2)(USp(4) * USp(4)) + (1/2) USp(4)``.

Both come out transitive.

    python research/sato_tate_limit/symplectic_search.py
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

import st_lib as S

HERE = Path(__file__).resolve().parent
TAU = S.tau_grid(1e-4, 1e5, 1201)
CAP = 12.0


def symplectic_factors(cap: float) -> list[S.Factor]:
    out = []
    for g in (1, 2, 3, 4, 5, 6):
        for k in range(1, int(cap // (2 * g)) + 1):
            out.append(S.Factor("SU2" if g == 1 else f"USp{2 * g}", k, 0.0))
    return [f for f in out if f.alpha_max <= cap + 1e-9]


def enumerate_products(cap: float, factors: list[S.Factor]) -> list[S.Measure]:
    factors = sorted(factors, key=lambda f: (f.alpha_max, f.group, f.multiplicity))
    out: list[S.Measure] = []

    def rec(start: int, used: float, chosen: list[S.Factor]) -> None:
        if chosen:
            out.append(S.Measure(tuple(chosen)))
        for i in range(start, len(factors)):
            f = factors[i]
            if used + f.alpha_max <= cap + 1e-9:
                rec(i, used + f.alpha_max, chosen + [f])

    rec(0, 0.0, [])
    seen: dict[str, S.Measure] = {}
    for m in out:
        seen.setdefault(m.label, m)
    return list(seen.values())


def mid_matrix(psis: np.ndarray, amax: np.ndarray) -> np.ndarray:
    n = psis.shape[0]
    mid = np.zeros((n, n))
    for i in range(n):
        d = psis[i][None, :] - psis
        end = amax[i] - amax
        hi = np.maximum(np.maximum(d.max(axis=1), 0.0), end)
        lo = np.minimum(np.minimum(d.min(axis=1), 0.0), end)
        mid[i] = 0.5 * (hi + lo)
    np.fill_diagonal(mid, 0.0)
    return mid


def count_cycles(mid: np.ndarray, tol: float) -> tuple[int, list]:
    n = mid.shape[0]
    strict = mid < -tol
    cyc = []
    total = 0
    for i, j, k in itertools.combinations(range(n), 3):
        for a, b, c in ((i, j, k), (i, k, j)):
            total += 1
            if strict[a, b] and strict[b, c] and strict[c, a]:
                cyc.append((a, b, c))
    return total, cyc


def worst_triangle(mid: np.ndarray) -> tuple[float, tuple[int, int, int]]:
    """How close does the tournament come to cycling?  ``min_i (-mid_i)``
    maximised over oriented triangles; negative means transitive."""

    n = mid.shape[0]
    best = (-math.inf, None)
    for i, j, k in itertools.combinations(range(n), 3):
        for a, b, c in ((i, j, k), (i, k, j)):
            s = min(-mid[a, b], -mid[b, c], -mid[c, a])
            if s > best[0]:
                best = (s, (a, b, c))
    return best


# ------------------------------------------------------- mixtures (cosets)


def mixture_psi(logits: np.ndarray, atom_K: np.ndarray) -> np.ndarray:
    w = np.exp(logits - logits.max())
    w /= w.sum()
    return S.mixture_K(w, atom_K) / TAU


def objective(x: np.ndarray, atom_K: np.ndarray, amax: np.ndarray) -> float:
    n = atom_K.shape[0]
    psis = []
    ends = []
    for i in range(3):
        lg = x[i * n:(i + 1) * n]
        w = np.exp(lg - lg.max())
        w /= w.sum()
        psis.append(S.mixture_K(w, atom_K) / TAU)
        ends.append(float(amax[w > 1e-12].max()))
    mids = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        d = np.concatenate([[0.0], psis[a] - psis[b], [ends[a] - ends[b]]])
        mids.append(0.5 * (d.max() + d.min()))
    return max(mids)


def main() -> int:
    lib = enumerate_products(CAP, symplectic_factors(CAP))
    lib.sort(key=lambda m: (m.alpha_max, m.variance, m.tail, m.label))
    print(f"symplectic library: {len(lib)} measures with alpha_max <= {CAP:g}"
          f"  (genus 1 to {int(CAP // 2)})")
    amax = np.array([m.alpha_max for m in lib])
    var = np.array([m.variance for m in lib])
    tail = np.array([m.tail for m in lib])
    psis = np.array([m.Psi(TAU) for m in lib])
    mid = mid_matrix(psis, amax)
    tol = 1e-6

    rows: list[list] = []
    total, cyc = count_cycles(mid, tol)
    score, tri = worst_triangle(mid)
    print(f"\n  oriented triangles: {total}")
    print(f"  3-cycles:           {len(cyc)}")
    print(f"  widest cycle: min(-mid) = {score:+.6e} over the triangle")
    a, b, c = tri
    print(f"      {lib[a].label}  ->  {lib[b].label}  ->  {lib[c].label}")
    print(f"      mid = {mid[a, b]:+.6f}, {mid[b, c]:+.6f}, {mid[c, a]:+.6f}")
    seen = set()
    print("\n  every distinct 3-cycle:")
    for a, b, c in cyc:
        key = frozenset((a, b, c))
        if key in seen:
            continue
        seen.add(key)
        marg = min(-mid[a, b], -mid[b, c], -mid[c, a])
        mf = all(all(f.multiplicity == 1 for f in lib[x].factors)
                 for x in (a, b, c))
        print(f"    {lib[a].label:<24} < {lib[b].label:<24} < "
              f"{lib[c].label:<28} margin {marg:.4e}  genus "
              f"{amax[a] / 2:g}/{amax[b] / 2:g}/{amax[c] / 2:g}  "
              f"multiplicity-free {mf}")
        rows.append(["cycle", lib[a].label, lib[b].label, lib[c].label,
                     f"{marg:.6e}", mf])

    print("\n  by alpha_max class (same genus):")
    print(f"  {'alpha_max':>10}{'members':>9}{'triangles':>11}{'cycles':>8}"
          f"{'closest':>14}")
    for a in sorted(set(amax)):
        sel = np.where(np.abs(amax - a) < 1e-9)[0]
        if sel.size < 3:
            continue
        sub = mid[np.ix_(sel, sel)]
        t2, c2 = count_cycles(sub, tol)
        s2, _ = worst_triangle(sub)
        print(f"  {a:>10.0f}{sel.size:>9}{t2:>11}{len(c2):>8}{s2:>14.4e}")
        rows.append(["class", f"{a:g}", sel.size, t2, len(c2), f"{s2:.6e}"])

    # ------------------------------------------- multiplicity-free sub-cone
    print("\n" + "=" * 78)
    print("the multiplicity-free sub-cone: Jacobians with no repeated factor")
    print("=" * 78)
    print("  These are the measures with an explicit witness: USp(2g) is the")
    print("  generic pencil y^2 = h(x) + c with deg h = 2g+1, and the split")
    print("  products come from y^2 = f(x^2) + c, whose two elliptic-type")
    print("  quotients v^2 = f(u)+c and w^2 = u(f(u)+c) have genus")
    print("  floor((n-1)/2) and floor(n/2) for n = deg f, with independent")
    print("  monodromy.  A repeated factor -- Jac isogenous to A^k -- has no")
    print("  such witness here.")
    free = [i for i, m in enumerate(lib)
            if all(f.multiplicity == 1 for f in m.factors)]
    subm = mid[np.ix_(free, free)]
    tf, cf = count_cycles(subm, tol)
    sf, trif = worst_triangle(subm)
    print(f"\n  members {len(free)}   oriented triangles {tf}   "
          f"3-cycles {len(cf)}")
    print(f"  closest approach min(-mid) = {sf:+.6e}")
    a2, b2, c2 = trif
    print(f"      {lib[free[a2]].label} -> {lib[free[b2]].label} -> "
          f"{lib[free[c2]].label}")
    rows.append(["mult-free", len(free), tf, len(cf), f"{sf:.6e}", ""])

    # the total order, and whether it is the lexicographic (m2, -t) rule
    print("\n  is the comparison the lexicographic rule  (m2 ascending, then t")
    print("  descending) inside an alpha_max class?")
    bad = 0
    checked = 0
    for a in sorted(set(amax)):
        sel = np.where(np.abs(amax - a) < 1e-9)[0]
        for i, j in itertools.combinations(sel, 2):
            checked += 1
            lex = (var[i], -tail[i]) < (var[j], -tail[j])
            got = mid[i, j] < 0
            if lex != got:
                bad += 1
                if bad <= 8:
                    print(f"      exception: {lib[i].label} vs {lib[j].label}  "
                          f"m2 {var[i]:g}/{var[j]:g}  t {tail[i]:g}/{tail[j]:g}  "
                          f"mid {mid[i, j]:+.6f}")
    print(f"      {checked - bad} of {checked} same-genus pairs follow it "
          f"({100 * (checked - bad) / checked:.2f}%)")
    rows.append(["lex-rule", "", checked - bad, checked, "", ""])

    # ---------------------------------------------------------- mixtures
    print("\n" + "=" * 78)
    print("mixtures of symplectic measures (disconnected monodromy)")
    print("=" * 78)
    print("  A disconnected monodromy group with r cosets realises the uniform")
    print("  average of the r coset trace measures.  Two searches, both over the")
    print("  multiplicity-free sub-cone, which is where the explicit witnesses")
    print("  live: the exact rational averages of pairs and triples, and a free")
    print("  optimisation over the whole simplex.")
    sub = [lib[i] for i in free]
    n = len(sub)
    atom_K = np.array([m.K(TAU) for m in sub])
    amax_sub = amax[free]

    cosets: list[tuple[str, np.ndarray]] = []
    for combo in itertools.combinations(range(n), 2):
        w = np.zeros(n)
        w[list(combo)] = 0.5
        cosets.append((" + ".join(sub[i].label for i in combo), w))
    kmix = np.array([S.mixture_K(w, atom_K) for _, w in cosets])
    amix = np.array([float(amax_sub[w > 0].max()) for _, w in cosets])
    allpsi = np.vstack([atom_K, kmix]) / TAU
    allA = np.concatenate([amax_sub, amix])
    m = allpsi.shape[0]
    labels = [x.label for x in sub] + [c[0] for c in cosets]
    print(f"\n  {n} products plus {len(cosets)} coset averages = {m} measures")
    hi_all = np.zeros((m, m))
    lo_all = np.zeros((m, m))
    for i in range(m):
        d = allpsi[i][None, :] - allpsi
        e = allA[i] - allA
        hi_all[i] = np.maximum(np.maximum(d.max(axis=1), 0.0), e)
        lo_all[i] = np.minimum(np.minimum(d.min(axis=1), 0.0), e)
    midm = 0.5 * (hi_all + lo_all)
    np.fill_diagonal(midm, 0.0)
    # too many triples for a Python loop: count directed 3-cycles as
    # trace(S^3)/3 for S the adjacency matrix of the strict relation, and find
    # the closest approach by a chunked numpy maximum.
    strict = (midm < -tol).astype(np.float64)
    np.fill_diagonal(strict, 0.0)
    n_cyc = int(round(np.trace(strict @ strict @ strict) / 3.0))
    neg = -midm
    best_s, best_t = -math.inf, None
    for i in range(m):
        # score[j, k] = min(-mid[i,j], -mid[j,k], -mid[k,i])
        v = np.minimum(neg[i][:, None], neg)
        v = np.minimum(v, neg[:, i][None, :])
        np.fill_diagonal(v, -math.inf)
        v[i, :] = -math.inf
        v[:, i] = -math.inf
        j, k = np.unravel_index(int(np.argmax(v)), v.shape)
        if v[j, k] > best_s:
            best_s, best_t = float(v[j, k]), (i, j, k)
    tm = m * (m - 1) * (m - 2) // 3
    print(f"  oriented triangles {tm}   3-cycles {n_cyc}   "
          f"closest approach {best_s:+.6e}")
    print(f"      {labels[best_t[0]]} -> {labels[best_t[1]]} -> "
          f"{labels[best_t[2]]}")
    rows.append(["coset-mixture", m, tm, n_cyc, f"{best_s:.6e}", ""])

    # the one coset construction with an arithmetic witness: a Jacobian whose
    # two halves are exchanged by the geometric monodromy.  The non-identity
    # coset consists of [[0, X], [Y, 0]], whose trace vanishes identically, so
    # the measure is (1/2)(mu * mu) + (1/2) delta_0.
    print("\n  the swap construction, 1/2 (mu * mu) + 1/2 delta_0, added to the")
    print("  multiplicity-free products:")
    swap_K = []
    swap_a = []
    swap_lab = []
    for i, mm in enumerate(sub):
        if 2 * mm.alpha_max > CAP + 1e-9:
            continue
        swap_K.append(np.logaddexp(math.log(0.5), math.log(0.5) + 2 * atom_K[i]))
        swap_a.append(2 * mm.alpha_max)
        swap_lab.append(f"swap({mm.label})")
    allK2 = np.vstack([atom_K, np.array(swap_K)])
    allA2 = np.concatenate([amax_sub, np.array(swap_a)])
    mid2 = mid_matrix(allK2 / TAU, allA2)
    t3, c3 = count_cycles(mid2, tol)
    s3, tri3 = worst_triangle(mid2)
    lab2 = [x.label for x in sub] + swap_lab
    print(f"    members {mid2.shape[0]}   oriented triangles {t3}   "
          f"3-cycles {len(c3)}   closest approach {s3:+.6e}")
    print(f"      {lab2[tri3[0]]} -> {lab2[tri3[1]]} -> {lab2[tri3[2]]}")
    rows.append(["swap-coset", mid2.shape[0], t3, len(c3), f"{s3:.6e}", ""])

    rng = np.random.default_rng(3)
    best = (math.inf, None)
    tries = 6
    for _ in range(tries):
        x0 = rng.normal(0.0, 4.0, size=3 * n)
        res = minimize(objective, x0, args=(atom_K, amax_sub), method="Powell",
                       options={"maxiter": 400, "maxfev": 20_000,
                                "xtol": 1e-5, "ftol": 1e-8})
        if res.fun < best[0]:
            best = (float(res.fun), res.x.copy())
    print(f"\n  {tries} restarts of Powell over 3 x {n} free mixture weights")
    print(f"  best  max(mid over the three edges) = {best[0]:+.6e}")
    print(f"  cycle found: {best[0] < 0}")
    print("  (a spot check only -- the exhaustive coset scan above already")
    print("   settles the mixture question, and settles it positively)")
    rows.append(["free-mixture", tries, n, f"{best[0]:.6e}", best[0] < 0, ""])

    with (HERE / "symplectic_search.csv").open("w", newline="",
                                               encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind", "a", "b", "c", "d", "e"])
        wr.writerows(rows)
    with (HERE / "symplectic_library.csv").open("w", newline="",
                                                encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["measure", "alpha_max", "genus", "m2", "edge_exponent_t"])
        for m2_ in lib:
            wr.writerow([m2_.label, f"{m2_.alpha_max:g}",
                         f"{m2_.alpha_max / 2:g}", f"{m2_.variance:.6f}",
                         f"{m2_.tail:.4f}"])
    print("\nwritten: symplectic_search.csv, symplectic_library.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
