#!/usr/bin/env python3
"""The multiplicity-free symplectic cone: the order is the dominance order.

A multiplicity-free symplectic measure is a product ``prod_i USp(2 lambda_i)``
indexed by a partition ``lambda`` of the genus ``G``; ``alpha_max = 2G``,
``m2 = length(lambda)``, ``t = sum lambda_i^2 + G/2``.  Because ``Psi`` is
additive over independent factors,

    Psi_lambda(tau) = ( sum_i kappa_{lambda_i}(tau) ) / tau,
    kappa_g := K_{USp(2g)},   kappa_0 = 0.

**Reduction (proved in TRANSITIVITY.md).**  If the sequence ``g -> kappa_g(tau)``
is concave for every ``tau > 0``, then ``lambda`` dominating ``mu`` implies
``Psi_lambda <= Psi_mu`` pointwise, hence ``lambda`` precedes ``mu`` with no
crossing.  Concavity of ``kappa`` is equivalent to

    M_{g+1} M_{g-1} / M_g^2 = 4 b_g(tau)^2 <= 1,

``b_g`` the Jacobi off-diagonal coefficient of the tilted Chebyshev weight
``e^{2 tau x} sqrt(1-x^2)`` (free value ``1/2``).  This script checks that
inequality over eleven decades of ``tau`` and every rank in reach, then walks
the dominance order genus by genus, counts crossings and looks for 3-cycles.

    python research/sato_tate_limit/transitivity_dominance.py
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

import numpy as np
from mpmath import mp, mpf

import kappa_lib as KL

HERE = Path(__file__).resolve().parent
GMAX = int(__import__("os").environ.get("GMAX", "15"))
TAU = np.geomspace(1e-4, 1e5, 1201)


# ------------------------------------------------------------- partitions


def partitions(n: int, cap: int | None = None) -> list[tuple[int, ...]]:
    if cap is None:
        cap = n
    if n == 0:
        return [()]
    out = []
    for first in range(min(n, cap), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return out


def dominates(a: tuple[int, ...], b: tuple[int, ...]) -> bool:
    """``a >= b`` in the dominance order (same total)."""

    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0
        sb += b[i] if i < len(b) else 0
        if sa < sb:
            return False
    return True


def stats(lam: tuple[int, ...]) -> tuple[float, float]:
    """``(m2, t)`` of ``prod USp(2 lambda_i)``."""

    return float(len(lam)), float(sum(l * l for l in lam) + sum(lam) / 2)


# ---------------------------------------------------------------- the run


def main() -> int:
    print("=" * 78)
    print("1.  kappa_g is concave in g:  M_{g+1} M_{g-1} <= M_g^2")
    print("=" * 78)
    kappa, bsq = KL.kappa_and_b(TAU, GMAX)
    worst = []
    for g in range(1, GMAX):
        v = 4.0 * bsq[g - 1]
        i = int(v.argmax())
        worst.append((g, float(v[i]), float(TAU[i])))
    print(f"  {'g':>3}{'max_tau 4 b_g^2':>20}{'at tau':>14}"
          f"{'   (must be <= 1)':<20}")
    for g, v, t in worst:
        print(f"  {g:>3}{v:>20.15f}{t:>14.5g}")
    print(f"\n  overall max = {max(v for _, v, _ in worst):.15f}")
    print("  (the value 1 is the tau -> 0 limit, where b_g = 1/2 exactly)")

    # 40-digit spot checks of the two ends
    mp.dps = 60
    print("\n  40-digit endpoint checks of 4 b_g(tau)^2 - 1:")
    for g in (1, 3, 6, 9):
        for t in (mpf("1e-3"), mpf(1), mpf(10), mpf("1e4")):
            h = KL.hankels(t, g + 1, KL.working_dps(g + 1, float(t)) + 60)
            val = 4 * h[g + 1] * h[g - 1] / h[g] ** 2 - 1
            print(f"      g = {g:>2}  tau = {mp.nstr(t, 5):>8}   "
                  f"{mp.nstr(val, 40)}")

    # higher-order structure of the increment sequence
    print("\n  finite differences of Delta_g = kappa_g - kappa_{g-1} "
          "(complete monotonicity):")
    kap = np.vstack([np.zeros_like(TAU), kappa])          # kappa_0 = 0
    delta = np.diff(kap, axis=0)                          # Delta_1..Delta_GMAX
    row = delta.copy()
    print(f"      {'order j':>8}{'max sign-violation of (-1)^j D^j Delta':>44}")
    for j in range(0, GMAX):
        sgn = (-1) ** j
        viol = float((-sgn * row).max()) if row.size else 0.0
        print(f"      {j:>8}{viol:>44.3e}")
        if row.shape[0] <= 1:
            break
        row = np.diff(row, axis=0)

    # ---------------------------------------------------------------- order
    print("\n" + "=" * 78)
    print("2.  the multiplicity-free order, genus by genus")
    print("=" * 78)
    rows, incomp_rows = [], []
    print(f"  {'G':>3}{'parts':>7}{'pairs':>7}{'dom-comp':>10}{'incomp':>8}"
          f"{'crossing':>10}{'cycles':>8}{'closest':>12}")
    for G in range(1, GMAX + 1):
        ps = partitions(G)
        n = len(ps)
        psi = np.array([sum(kappa[p - 1] for p in lam) / TAU for lam in ps])
        mid = np.zeros((n, n))
        cross = np.zeros((n, n), dtype=int)
        ncomp = ninc = ncross = 0
        for i, j in itertools.combinations(range(n), 2):
            d = psi[i] - psi[j]
            dd = np.concatenate([[0.0], d, [0.0]])
            m = 0.5 * (float(dd.max()) + float(dd.min()))
            mid[i, j], mid[j, i] = m, -m
            sc = _sign_changes(d)
            cross[i, j] = cross[j, i] = sc
            ncross += sc > 0
            comp = dominates(ps[i], ps[j]) or dominates(ps[j], ps[i])
            ncomp += comp
            ninc += not comp
            if not comp:
                m2i, ti = stats(ps[i])
                m2j, tj = stats(ps[j])
                incomp_rows.append([G, "".join(map(str, ps[i])),
                                    "".join(map(str, ps[j])),
                                    m2i - m2j, ti - tj, f"{m:.10f}", sc])
            if comp:
                hi = ps[i] if dominates(ps[i], ps[j]) else ps[j]
                # dominating partition must precede: mid(dominating - other) < 0
                sgn = m if hi == ps[i] else -m
                if sgn > 0 or sc > 0:
                    print(f"      !! dominance violated: {ps[i]} vs {ps[j]}"
                          f"  mid={m:+.6f} crossings={sc}")
        # 3-cycles
        cycles, closest = 0, -np.inf
        for i, j, k in itertools.combinations(range(n), 3):
            for a, b, c in ((i, j, k), (i, k, j)):
                trio = (mid[a, b], mid[b, c], mid[c, a])
                if max(trio) < 0:
                    cycles += 1
                    closest = max(closest, max(trio))
                else:
                    closest = max(closest, -max(trio)) if False else closest
        # closest approach: the best (least positive) max over oriented triples
        best = np.inf
        for i, j, k in itertools.combinations(range(n), 3):
            for a, b, c in ((i, j, k), (i, k, j)):
                best = min(best, max(mid[a, b], mid[b, c], mid[c, a]))
        rows.append([G, n, n * (n - 1) // 2, ncomp, ninc, ncross, cycles,
                     f"{best:.6g}"])
        print(f"  {G:>3}{n:>7}{n * (n - 1) // 2:>7}{ncomp:>10}{ninc:>8}"
              f"{ncross:>10}{cycles:>8}{best:>12.4g}")

    print("\n  incomparable pairs (dominance) and what decides them:")
    print(f"  {'G':>3}  {'lambda':<14}{'mu':<14}{'dm2':>6}{'dt':>8}"
          f"{'mid':>14}{'crossings':>11}")
    for r in incomp_rows:
        if r[0] <= 8:
            print(f"  {r[0]:>3}  {r[1]:<14}{r[2]:<14}{r[3]:>6g}{r[4]:>8g}"
                  f"{float(r[5]):>14.6f}{r[6]:>11}")
    print(f"  ... {len(incomp_rows)} incomparable pairs in all to genus "
          f"{GMAX}")
    ncomono = sum(1 for r in incomp_rows if r[3] * r[4] > 0)
    print(f"  of these, {ncomono} are (m2, t)-comonotone, hence must cross;"
          f"  {sum(1 for r in incomp_rows if r[6] > 0)} do cross")

    with (HERE / "transitivity_dominance.csv").open("w", newline="",
                                                    encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["genus", "partitions", "pairs", "dominance_comparable",
                     "incomparable", "crossing_pairs", "cycles",
                     "closest_approach"])
        wr.writerows(rows)
        wr.writerow([])
        wr.writerow(["genus", "lambda", "mu", "dm2", "dt", "mid",
                     "sign_changes"])
        wr.writerows(incomp_rows)
    print("\nwritten: transitivity_dominance.csv")
    return 0


def _sign_changes(d: np.ndarray) -> int:
    tol = 1e-9 * max(float(np.abs(d).max()), 1e-300)
    s = np.sign(np.where(np.abs(d) < tol, 0.0, d))
    s = s[s != 0]
    if s.size == 0:
        return 0
    return int(np.count_nonzero(np.diff(s) != 0))


if __name__ == "__main__":
    raise SystemExit(main())
