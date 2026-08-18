#!/usr/bin/env python3
"""Route 1 of session brief K: is the same-genus comparison a trade-off in
``(dm2, dt)``?

``lex_exceptions.py`` fitted ``sign mid = sign(dm2 - kappa dt)`` with
``kappa ~ 0.668`` and brief K asks whether the true constant is ``2/3``.  This
script settles the question by tabulating **every** same-genus pair of the
symplectic library together with its ``(dm2, dt)`` and asking whether
``sign(mid)`` is a function of that pair at all.

It also records, for each pair, where the sup and the inf of
``D = Psi_mu - Psi_nu`` sit, which is the input to the level lemma.

    python research/sato_tate_limit/transitivity_pairs.py
"""

from __future__ import annotations

import csv
import itertools
from collections import defaultdict
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
    buckets: dict[tuple[float, float], list] = defaultdict(list)
    for a in sorted(set(amax)):
        sel = np.where(np.abs(amax - a) < 1e-9)[0]
        for i, j in itertools.combinations(sel, 2):
            d = psis[i] - psis[j]
            dd = np.concatenate([[0.0], d, [0.0]])
            hi, lo = float(dd.max()), float(dd.min())
            mid = 0.5 * (hi + lo)
            # where the extrema sit: index 0 is tau=0, index -1 is tau=inf
            ihi, ilo = int(dd.argmax()), int(dd.argmin())
            thi = 0.0 if ihi == 0 else (np.inf if ihi == len(dd) - 1
                                        else float(TAU[ihi - 1]))
            tlo = 0.0 if ilo == 0 else (np.inf if ilo == len(dd) - 1
                                        else float(TAU[ilo - 1]))
            nsc = S.sign_changes(d)
            dm2, dt = var[i] - var[j], tail[i] - tail[j]
            rows.append([lib[i].label, lib[j].label, a / 2, dm2, dt,
                         f"{mid:.10f}", f"{hi:.10f}", f"{lo:.10f}",
                         f"{thi:g}", f"{tlo:g}", nsc])
            buckets[(dm2, dt)].append((mid, lib[i].label, lib[j].label))

    print("=" * 78)
    print("same-genus pairs of the symplectic library (alpha_max <= 12)")
    print("=" * 78)
    print(f"  measures {len(lib)},  same-genus pairs {len(rows)}")

    # --- Is sign(mid) a function of (dm2, dt) at all? --------------------
    bad = []
    for key, entries in sorted(buckets.items()):
        signs = {np.sign(m) for m, _, _ in entries}
        if len(signs) > 1:
            bad.append((key, entries))
    print(f"\n  distinct (dm2, dt) values: {len(buckets)}")
    print(f"  values carrying BOTH signs of mid: {len(bad)}")
    for key, entries in bad:
        pos = [e for e in entries if e[0] > 0]
        neg = [e for e in entries if e[0] < 0]
        print(f"\n    (dm2, dt) = ({key[0]:g}, {key[1]:g}):"
              f"  {len(pos)} with mid > 0, {len(neg)} with mid < 0")
        for m, a, b in sorted(entries)[:1] + sorted(entries)[-1:]:
            print(f"        mid = {m:+.8f}   {a}   vs   {b}")

    # --- the exact feasible interval for kappa ---------------------------
    pairs = [(r[3], r[4], float(r[5])) for r in rows if r[3] != 0]
    lows, highs = [], []
    for dm, dt, mid in pairs:
        # want sign(dm - k dt) == sign(mid)
        if dt == 0:
            continue
        thr = dm / dt                       # dm - k dt = 0 at k = dm/dt
        if dt > 0:                          # dm - k dt > 0  <=>  k < thr
            (highs if mid > 0 else lows).append(thr)
        else:                               # dm - k dt > 0  <=>  k > thr
            (lows if mid > 0 else highs).append(thr)
    lo_k = max(lows) if lows else -np.inf
    hi_k = min(highs) if highs else np.inf
    print(f"\n  linear rule sign(mid) = sign(dm2 - kappa dt):")
    print(f"      constraints force  kappa > {lo_k:.10f}  and  "
          f"kappa < {hi_k:.10f}")
    print(f"      feasible: {'YES' if lo_k < hi_k else 'NO'}")
    best = (-1, None)
    ks = np.unique(np.concatenate([np.array([dm / dt for dm, dt, _ in pairs
                                             if dt != 0]),
                                   np.linspace(0, 5, 20001)]))
    ks = np.concatenate([ks, ks + 1e-9])
    for k in ks:
        ok = sum(1 for dm, dt, mid in pairs if (dm - k * dt > 0) == (mid > 0))
        if ok > best[0]:
            best = (ok, float(k))
    print(f"      best kappa = {best[1]:.10f} correct on {best[0]} of "
          f"{len(pairs)}")
    # the binding pair
    binding = [(dm, dt, mid) for dm, dt, mid in pairs
               if dt != 0 and abs(dm / dt - lo_k) < 1e-12]
    from fractions import Fraction
    frac = Fraction(lo_k).limit_denominator(1000)
    print(f"      kappa > {lo_k:.10f} = {frac} comes from {len(binding)} "
          f"pair(s) with (dm2, dt) = "
          + (str(binding[0][:2]) if binding else "None"))

    with (HERE / "transitivity_pairs.csv").open("w", newline="",
                                                encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["mu", "nu", "genus", "dm2", "dt", "mid", "sup", "inf",
                     "tau_sup", "tau_inf", "sign_changes"])
        wr.writerows(rows)
    print("\nwritten: transitivity_pairs.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
