#!/usr/bin/env python3
"""Which multiplicity-free symplectic measures actually have a curve family?

``TRANSITIVITY.md`` finds the first same-genus 3-cycle at genus 13, among
multiplicity-free products of ``USp(2g)`` blocks, and observes that no vertex
needs a repeated isogeny factor -- so the arithmetic input it needs is a
Jacobian splitting into four or five **distinct** independent factors of
prescribed dimensions.  This script asks which such splittings the one
construction that reliably produces them -- an elementary abelian 2-cover of
``P^1`` -- can supply, and whether a same-genus 3-cycle survives inside that
realisable set.

**The construction.**  Let ``G = (Z/2)^r`` act on a curve ``C`` over a field of
characteristic ``!= 2`` with ``C/G = P^1``.  Inertia is cyclic, hence of order
two, so Riemann-Hurwitz gives

    2 g(C) - 2 = 2^r (-2) + B 2^{r-1},        B = (2 g + 2^{r+1} - 2) / 2^{r-1}

with ``B`` the number of branch points, and ``H^1(C) = sum_{chi != 1} H^1(C)_chi``
with ``H^1(C)_chi = H^1(C / ker chi)``, the double cover branched at
``B_chi = {P : chi(v_P) = -1}``, of genus ``|B_chi|/2 - 1``.  The data is the
multiset of inertia vectors ``{v_P} in (F_2^r \\ 0)^B`` with ``sum v_P = 0``
(the cover exists) and ``<v_P> = G`` (it is connected).  Enumerating those
multisets enumerates the achievable partitions.

**The consequence, immediate from ``B``.**  A ``(Z/2)^r``-cover has at most
``2^r - 1`` blocks and every block has genus at most ``B/2 - 1``, and ``B``
*shrinks* as ``r`` grows.  Many blocks and one large block are incompatible.
At genus 13:

    r = 2   B = 16   at most  3 blocks, each of genus <= 7
    r = 3   B = 10   at most  7 blocks, each of genus <= 4
    r = 4   B =  7   at most 15 blocks, each of genus <= 3

so ``(7,2,2,1,1)``, ``(5,4,2,2)`` and ``(6,4,1,1,1)`` -- the three vertices of
``TRANSITIVITY.md``'s genus-13 cycle -- are all **out of reach**.

    python research/sato_tate_limit/realisable_partitions.py
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

import numpy as np

import kappa_lib as KL

HERE = Path(__file__).resolve().parent
GMAX_GENUS = 11
TAU = np.geomspace(1e-4, 1e5, 151)


# ------------------------------------------------------------- realisability


def branch_count(g: int, r: int):
    """``B`` for a ``(Z/2)^r``-cover of ``P^1`` of genus ``g``, or None."""
    num = 2 * g + 2 ** (r + 1) - 2
    den = 2 ** (r - 1)
    if num % den:
        return None
    B = num // den
    return B if B > 0 else None


def par(u: int, v: int) -> int:
    return bin(u & v).count("1") % 2


def profiles_r(g: int, r: int) -> set[tuple[int, ...]]:
    """the partitions of ``g`` realised by ``(Z/2)^r``-covers of ``P^1``."""
    B = branch_count(g, r)
    out: set[tuple[int, ...]] = set()
    if B is None:
        return out
    nz = list(range(1, 2 ** r))
    if B > 24:                      # enumeration would blow up; r = 2 is closed
        return out
    for v in itertools.combinations_with_replacement(nz, B):
        s = 0
        for z in v:
            s ^= z
        if s:
            continue
        span = {0}
        for z in v:
            span |= {y ^ z for y in span}
        if len(span) != 2 ** r:
            continue
        gs = []
        ok = True
        for chi in nz:
            b = sum(1 for z in v if par(chi, z))
            if b == 0:
                ok = False
                break
            gs.append(b // 2 - 1)
        if not ok:
            continue
        lam = tuple(sorted((x for x in gs if x > 0), reverse=True))
        if sum(lam) == g:
            out.add(lam)
    return out


def profiles_r2(g: int) -> set[tuple[int, ...]]:
    """``r = 2`` in closed form: the data is ``(d1,d2,d3)`` with
    ``sum d_i = B = g+3`` and ``|B_i| = d_j + d_k``."""
    B = g + 3
    out = set()
    for d1 in range(B + 1):
        for d2 in range(B - d1 + 1):
            d3 = B - d1 - d2
            if sum(1 for d in (d1, d2, d3) if d > 0) < 2:
                continue
            sizes = (d2 + d3, d3 + d1, d1 + d2)
            if any(s < 2 for s in sizes):
                continue
            lam = tuple(sorted((s // 2 - 1 for s in sizes if s // 2 - 1 > 0),
                               reverse=True))
            if sum(lam) == g:
                out.add(lam)
    return out


_ACH: dict[int, dict[tuple[int, ...], str]] = {}


def realisable(g: int) -> dict[tuple[int, ...], str]:
    """Partition -> the construction that realises it.

    Three rules, applied to closure:

    * ``r=1``  the generic hyperelliptic pencil ``y^2 = h(x) + c``, one block;
    * ``r=2,3,4``  an elementary abelian 2-cover of ``P^1`` (above);
    * ``prym``  a double cover ``C -> X`` with ``X`` of genus ``h`` already
      realised: ``Jac(C) ~ Jac(X) x Prym``, ``dim Prym = g - h``, and
      Riemann-Hurwitz forces ``2g - 2 = 2(2h-2) + (branch points) >= 4h - 4``,
      i.e. ``g >= 2h - 1``: **the new block must be at least ``h - 1``, the
      genus of the base minus one**.  So a partition is reachable this way only
      if its largest part is at least the sum of the others minus one, and that
      condition then recurses.
    """
    if g in _ACH:
        return _ACH[g]
    out: dict[tuple[int, ...], str] = {(g,): "r=1"}
    for lam in profiles_r2(g):
        out.setdefault(lam, "r=2")
    for r in (3, 4):
        for lam in profiles_r(g, r):
            out.setdefault(lam, f"r={r}")
    for h in range(1, g):                      # Prym of a double cover
        p = g - h
        if p < h - 1:
            continue                           # Riemann-Hurwitz
        for mu in realisable(h):
            lam = tuple(sorted(mu + (p,), reverse=True))
            out.setdefault(lam, "prym")
    _ACH[g] = out
    return out


# --------------------------------------------------------------- the search


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


def cycles_among(lams, kappa):
    psis = np.array([KL.psi_from_kappa(kappa, TAU, l) for l in lams])
    n = len(lams)
    mid = np.zeros((n, n))
    for i in range(n):
        d = psis[i][None, :] - psis
        hi = np.maximum(d.max(axis=1), 0.0)
        lo = np.minimum(d.min(axis=1), 0.0)
        mid[i] = 0.5 * (hi + lo)          # same alpha_max: no endpoint term
    np.fill_diagonal(mid, 0.0)
    strict = mid < -1e-9
    out = []
    for i, j, k in itertools.combinations(range(n), 3):
        for a, b, c in ((i, j, k), (i, k, j)):
            if strict[a, b] and strict[b, c] and strict[c, a]:
                out.append((a, b, c,
                            min(-mid[a, b], -mid[b, c], -mid[c, a])))
    seen, uniq = set(), []
    for a, b, c, m in out:
        key = frozenset((a, b, c))
        if key not in seen:
            seen.add(key)
            uniq.append((a, b, c, m))
    return uniq


def label(lam):
    return " x ".join(f"USp{2 * p}" if p > 1 else "SU2" for p in lam)


def main() -> int:
    rows: list[list] = []

    print("=" * 96)
    print("1.  what an elementary abelian 2-cover of P^1 can split a Jacobian"
          " into")
    print("=" * 96)
    print(f"  {'genus':>6}{'r':>4}{'B':>5}{'max blocks':>12}"
          f"{'max block genus':>17}")
    for g in (13, 14, 15):
        for r in (2, 3, 4, 5):
            B = branch_count(g, r)
            if B is None:
                print(f"  {g:>6}{r:>4}{'-':>5}{'-':>12}{'(no such cover)':>17}")
                continue
            print(f"  {g:>6}{r:>4}{B:>5}{2 ** r - 1:>12}{B // 2 - 1:>17}")
            rows.append(["shape", g, r, B, 2 ** r - 1, B // 2 - 1])
        print()

    R = realisable(13)
    print("  the vertices of TRANSITIVITY.md's two genus-13 cycles:")
    for tag, tri in (("cycle A, margin 8.0e-3",
                      [(7, 2, 2, 1, 1), (5, 4, 2, 2), (6, 4, 1, 1, 1)]),
                     ("cycle B, margin 4.4e-4",
                      [(8, 1, 1, 1, 1, 1), (4, 4, 4, 1), (6, 3, 2, 2)])):
        print(f"    {tag}")
        for lam in tri:
            print(f"      {str(lam):<20} blocks {len(lam)}  max {max(lam)}   "
                  f"realisable: {R.get(lam, 'NO')}")
            rows.append(["cycle13-vertex", tag, str(lam), len(lam), max(lam),
                         R.get(lam, "NO")])
    print("    at most ONE character can have |B_chi| = B: if two did, their")
    print("    product would be trivial on every v_P, hence trivial, since the")
    print("    v_P span.  That alone kills (4,4,4,1), which would need three.")

    print()
    print("=" * 96)
    print("2.  same-genus 3-cycles restricted to the realisable partitions")
    print("    (genus 12 and 13 are settled by the Corollary instead:"
          " TRANSITIVITY.md")
    print("     is exhaustive over all partitions and finds 0 cycles at genus"
          " <= 12 and")
    print("     exactly 2 at genus 13, and part 1 above shows each of those 2"
          " has an")
    print("     unrealisable vertex -- so 0 realisable same-genus cycles at"
          " genus <= 13)")
    print("=" * 96)
    kappa, _ = KL.kappa_and_b(TAU, GMAX_GENUS)
    print(f"  {'genus':>6}{'partitions':>12}{'realisable':>12}"
          f"{'cycles (all)':>14}{'cycles (realisable)':>21}")
    for g in range(2, GMAX_GENUS + 1):
        allp = partitions(g)
        R = realisable(g)
        rp = [l for l in allp if l in R]
        ca = cycles_among(allp, kappa) if len(allp) <= 200 else None
        cr = cycles_among(rp, kappa)
        print(f"  {g:>6}{len(allp):>12}{len(rp):>12}"
              f"{(len(ca) if ca is not None else -1):>14}{len(cr):>21}")
        rows.append(["search", g, len(allp), len(rp),
                     len(ca) if ca is not None else "", len(cr)])
        for a, b, c, m in sorted(cr, key=lambda t: -t[3])[:5]:
            print(f"        {label(rp[a])} < {label(rp[b])} < {label(rp[c])}"
                  f"   margin {m:.4e}")
            rows.append(["realisable-cycle", g, str(rp[a]), str(rp[b]),
                         str(rp[c]), f"{m:.6e}"])
        if ca:
            for a, b, c, m in sorted(ca, key=lambda t: -t[3])[:5]:
                tag = " ".join(R.get(l, "NO") for l in (allp[a], allp[b],
                                                        allp[c]))
                print(f"     [all] {label(allp[a])} < {label(allp[b])} < "
                      f"{label(allp[c])}   margin {m:.4e}   realisable: {tag}")
                rows.append(["all-cycle", g, str(allp[a]), str(allp[b]),
                             str(allp[c]), f"{m:.6e}", tag])

    with (HERE / "realisable_partitions.csv").open("w", newline="",
                                                   encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind"] + [f"c{i}" for i in range(7)])
        for r in rows:
            wr.writerow(list(r) + [""] * (8 - len(r)))
    print("\nwritten: realisable_partitions.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
