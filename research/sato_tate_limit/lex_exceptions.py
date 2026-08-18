#!/usr/bin/env python3
"""The 20 exceptions to the same-genus lexicographic rule of ``FINDINGS.md``.

``FINDINGS.md`` reports that inside a fixed ``alpha_max`` class the limiting
comparison agrees with "``m2`` ascending, ties broken by ``t`` descending" on
745 of 765 pairs of the symplectic library, and asks (Open, item 4) whether
same-genus transitivity is a theorem.  This script isolates the 20 exceptions
and asks what separates them from the 745.

Result (computed): every one of the 20 exceptions has ``dm2 < 0`` **and**
``dt < 0`` -- they are exactly the comonotone (hence crossing) pairs on which
the edge term outweighs a small ``m2`` gap.  Replacing the lexicographic rule by
the one-parameter trade-off

    sign mid(Psi_mu - Psi_nu)  =  sign( dm2 - kappa dt ),     kappa ~ 0.668

raises the agreement from 715 to 731 of the 735 pairs with ``dm2 != 0``, and all
30 pairs tied in ``m2`` are decided by ``t`` (larger ``t`` precedes), 30 of 30.
Four pairs still escape, so this is a sharper description and not a proof.

    python research/sato_tate_limit/lex_exceptions.py
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path

import numpy as np

import st_lib as S

HERE = Path(__file__).resolve().parent
TAU = S.tau_grid(1e-4, 1e5, 1201)
CAP = 12.0


def factors(cap: float) -> list[S.Factor]:
    out = []
    for g in (1, 2, 3, 4, 5, 6):
        for k in range(1, int(cap // (2 * g)) + 1):
            out.append(S.Factor("SU2" if g == 1 else f"USp{2 * g}", k, 0.0))
    return [f for f in out if f.alpha_max <= cap + 1e-9]


def products(cap: float) -> list[S.Measure]:
    facs = sorted(factors(cap),
                  key=lambda f: (f.alpha_max, f.group, f.multiplicity))
    out: list[S.Measure] = []

    def rec(start: int, used: float, chosen: list[S.Factor]) -> None:
        if chosen:
            out.append(S.Measure(tuple(chosen)))
        for i in range(start, len(facs)):
            f = facs[i]
            if used + f.alpha_max <= cap + 1e-9:
                rec(i, used + f.alpha_max, chosen + [f])

    rec(0, 0.0, [])
    seen: dict[str, S.Measure] = {}
    for m in out:
        seen.setdefault(m.label, m)
    return list(seen.values())


def main() -> int:
    lib = products(CAP)
    lib.sort(key=lambda m: (m.alpha_max, m.variance, m.tail, m.label))
    amax = np.array([m.alpha_max for m in lib])
    var = np.array([m.variance for m in lib])
    tail = np.array([m.tail for m in lib])
    psis = np.array([m.Psi(TAU) for m in lib])

    rows = []
    exc, agree = [], []
    for a in sorted(set(amax)):
        sel = np.where(np.abs(amax - a) < 1e-9)[0]
        for i, j in itertools.combinations(sel, 2):
            d = psis[i] - psis[j]
            dd = np.concatenate([[0.0], d, [0.0]])
            mid = 0.5 * (float(dd.max()) + float(dd.min()))
            lex = (var[i], -tail[i]) < (var[j], -tail[j])
            got = mid < 0
            rec = (lib[i].label, lib[j].label, a / 2, var[i] - var[j],
                   tail[i] - tail[j], mid)
            (agree if lex == got else exc).append(rec)

    print(f"same-genus pairs: {len(agree) + len(exc)};  "
          f"lexicographic rule holds on {len(agree)}, fails on {len(exc)}")
    print("\nthe exceptions:")
    print(f"  {'mu':<28}{'nu':<28}{'g':>3}{'dm2':>6}{'dt':>8}{'mid':>12}")
    for la, lb, g, dm, dt, mid in exc:
        print(f"  {la:<28}{lb:<28}{g:>3g}{dm:>6g}{dt:>8.1f}{mid:>12.5f}")
        rows.append(["exception", la, lb, g, dm, dt, f"{mid:.8f}"])

    # What separates them?  Scan the single threshold kappa in
    #     sign mid  =  sign( dm2 - kappa dt )
    # over all pairs with dm2 != 0; kappa = 0 is the plain m2 rule.
    ks = np.concatenate([np.linspace(0, 0.2, 4001), np.linspace(0.2, 5, 2001)])
    best = (-1, None)
    pairs = [(dm, dt, mid) for _, _, _, dm, dt, mid in agree + exc if dm != 0]
    for k in ks:
        ok = sum(1 for dm, dt, mid in pairs
                 if (dm - k * dt > 0) == (mid > 0))
        if ok > best[0]:
            best = (ok, k)
    print(f"\n  pairs with dm2 != 0: {len(pairs)}")
    print(f"  best single-threshold rule  sign(mid) = sign(dm2 - kappa dt):")
    print(f"      kappa = {best[1]:.5f}   correct on {best[0]} of {len(pairs)}"
          f"  ({100 * best[0] / len(pairs):.2f}%)")
    print(f"  (kappa = 0 is the plain m2 rule: "
          f"{sum(1 for dm, dt, mid in pairs if (dm > 0) == (mid > 0))} of "
          f"{len(pairs)})")
    rows.append(["threshold", f"{best[1]:.6f}", best[0], len(pairs), "", "", ""])

    ties = [(la, lb, g, dm, dt, mid)
            for la, lb, g, dm, dt, mid in agree + exc if dm == 0]
    okt = sum(1 for _, _, _, _, dt, mid in ties if (dt > 0) == (mid < 0))
    print(f"\n  pairs tied in m2: {len(ties)};  larger t precedes on {okt} "
          f"({100 * okt / max(len(ties), 1):.2f}%)")
    rows.append(["m2-ties", len(ties), okt, "", "", "", ""])

    with (HERE / "lex_exceptions.csv").open("w", newline="",
                                            encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind", "a", "b", "c", "d", "e", "f"])
        for r in rows:
            wr.writerow(list(r) + [""] * (7 - len(r)))
    print("\nwritten: lex_exceptions.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
