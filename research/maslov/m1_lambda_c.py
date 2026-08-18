#!/usr/bin/env python3
"""Brief M item 1 -- is the critical lambda prior-free?

Four experiments.

**1a Extremum census.**  For a pair ``(a,b)`` each extremum of
``f = u_a - u_b`` is either *interior* (nondegenerate, at a finite ``s``) or a
*plateau* (the extremum is the limit at ``s = +-inf``, approached
exponentially).  On the known 3-cycle two of the three edges carry a plateau
extremum; inside a ``phi``-class (equal ``(r, M)``) ``f(+-inf) = 0`` and both
extrema are always interior.

**1b Laplace expansion.**  With two interior nondegenerate extrema,

    A_lambda = A_inf + (1/(2 lambda)) log[ (rho(s_+)/rho(s_-))
                                            sqrt(|f''(s_-)|/|f''(s_+)|) ] + O(1/lambda^2)

-- the ``sqrt(2 pi / lambda)`` factors cancel between softmax and softmin, and
so does any *constant* value of ``rho``.  So for a uniform prior the correction
is prior-free.  A plateau extremum at the right instead gives
``-(1/lambda) log(S_+ - log lambda - log c - gamma)`` where ``S_+`` is the right
end of the support: the *width* still cancels (rho is constant) but the
*distance to the endpoint* does not.

**1c The prediction, tested.**  lambda_c of the known 3-cycle must then depend
on the uniform support only through its right endpoint ``S_+``.  It does.

**1d The general statement.**  lambda_c measured for every 3-cycle of the
``F_11`` pool under five uniform priors of very different support.

Writes ``m1_lambda_c.csv``.
"""

from __future__ import annotations

import csv
import math
import sys

import numpy as np

import common as C

DS = 1.0e-3  # common grid spacing for every prior, so grids are comparable
CENSUS_GRID = np.linspace(-40.0, 40.0, 800001)


def matched_uniform(lo, hi, ds=DS):
    n = int(round((hi - lo) / ds)) + 1
    return C.uniform_prior(lo, hi, n)


# ------------------------------------------------------------------ extrema


def extremum_census(a, b, s=CENSUS_GRID):
    f = C.u_of(a, s) - C.u_of(b, s)
    fm = math.log(math.log(len(a))) - math.log(math.log(len(b)))
    fp = math.log(math.log(max(a))) - math.log(math.log(max(b)))
    out = {"f_minus_inf": fm, "f_plus_inf": fp}
    for tag, k in (("max", int(f.argmax())), ("min", int(f.argmin()))):
        val = float(f[k])
        interior = (val > max(fm, fp) + 1e-13) if tag == "max" else (val < min(fm, fp) - 1e-13)
        loc, curv = float(s[k]), float("nan")
        if interior and 0 < k < s.size - 1:
            h = float(s[1] - s[0])
            d2 = f[k + 1] - 2 * f[k] + f[k - 1]
            curv = float(d2 / h**2)
            if d2 != 0:
                loc = float(s[k] - 0.5 * h * (f[k + 1] - f[k - 1]) / d2)
        end = None
        if not interior:
            end = "+inf" if abs(val - fp) <= abs(val - fm) else "-inf"
        out[tag] = {"value": val, "s": loc, "interior": interior,
                    "plateau_end": end, "curvature": curv}
    out["mid"] = 0.5 * (out["max"]["value"] + out["min"]["value"])
    out["osc"] = out["max"]["value"] - out["min"]["value"]
    return out


def chi_prior_free(cen):
    if not (cen["max"]["interior"] and cen["min"]["interior"]):
        return float("nan")
    return 0.25 * math.log(abs(cen["min"]["curvature"]) / abs(cen["max"]["curvature"]))


# ------------------------------------------------------------- critical lambda


def critical_lambda(Us, w, lo=1e-2, hi=1e9, iters=80):
    """Smallest lambda above which the three ``u``-rows form a directed cycle."""

    fs = [Us[i] - Us[(i + 1) % 3] for i in range(3)]

    def is_cycle(L):
        return all(C.soft_mid(f, w, L) > C.TIE for f in fs)

    if not is_cycle(hi):
        return None
    if is_cycle(lo):
        return lo
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        if is_cycle(mid):
            hi = mid
        else:
            lo = mid
    return hi


def main():
    rows, out = [], []

    def say(t=""):
        print(t)
        out.append(t)

    say("=" * 78)
    say("1a.  Extremum census, known 3-cycle {(6,3,3),(7,2,1),(6,5,1)}")
    say("=" * 78)
    for i in range(3):
        a, b = C.CYCLE[i], C.CYCLE[(i + 1) % 3]
        cen = extremum_census(a, b)
        say(f"  {a} -> {b}   A_inf = {cen['mid']:+.9f}   osc = {cen['osc']:.9f}")
        for tag in ("max", "min"):
            e = cen[tag]
            kind = "interior" if e["interior"] else f"PLATEAU at s={e['plateau_end']}"
            extra = f"  f''={e['curvature']:+.6e}" if e["interior"] else ""
            say(f"     {tag}: f={e['value']:+.9f}  s={e['s']:+.5f}  {kind}{extra}")
        say(f"     f(-inf)={cen['f_minus_inf']:+.9f}   f(+inf)={cen['f_plus_inf']:+.9f}"
            f"   chi={chi_prior_free(cen)}")
        rows.append(dict(block="census", a=str(a), b=str(b), A_inf=cen["mid"], osc=cen["osc"],
                         max_interior=cen["max"]["interior"], min_interior=cen["min"]["interior"],
                         s_max=cen["max"]["s"], s_min=cen["min"]["s"],
                         chi=chi_prior_free(cen)))

    say()
    say("=" * 78)
    say("1b.  lambda_c of the known cycle vs the uniform support (matched ds=1e-3)")
    say("     PREDICTION: depends on the right endpoint S_+ only.")
    say("=" * 78)
    say(f"  {'support':>18s} {'width':>7s} {'S_+':>6s} {'lambda_c':>12s}")
    supports = [(-8, 8), (-12, 12), (-16, 16), (-24, 24), (-32, 32),
                (-6, 20), (-20, 6), (-6, 6), (-4, 4),
                (-40, 6), (-8, 6), (-40, 12), (-8, 12), (-40, 20), (-8, 20)]
    for lo, hi in supports:
        s, w = matched_uniform(lo, hi)
        Us = [C.u_of(a, s) for a in C.CYCLE]
        lc = critical_lambda(Us, w)
        say(f"  {f'[{lo},{hi}]':>18s} {hi-lo:7d} {hi:6d} {lc:12.2f}")
        rows.append(dict(block="known_support", lo=lo, hi=hi, width=hi - lo, lambda_c=lc))

    say()
    say("     non-flat priors (rho(s_+) != rho(s_-) -> genuine prior dependence):")
    for name, s, w in [("gaussian sd=2", *C.gaussian_prior(0.0, 2.0, 40001)),
                       ("gaussian sd=4", *C.gaussian_prior(0.0, 4.0, 40001)),
                       ("gaussian sd=8", *C.gaussian_prior(0.0, 8.0, 40001)),
                       ("logistic sc=2", *C.logistic_prior(0.0, 2.0, 40001)),
                       ("logistic sc=4", *C.logistic_prior(0.0, 4.0, 40001))]:
        Us = [C.u_of(a, s) for a in C.CYCLE]
        lc = critical_lambda(Us, w)
        say(f"  {name:>18s} {'-':>7s} {'-':>6s} {lc:12.2f}")
        rows.append(dict(block="known_support", prior=name, lambda_c=lc))

    say()
    say("=" * 78)
    say("1c.  The F_11 pool: cycle census and where the cycles sit")
    say("=" * 78)
    pool = C.f11_pool()
    sT, wT = matched_uniform(-14, 14, 2e-3)
    U, um, up = C.u_matrix(pool, sT)
    A = C.tropical_matrix(U, um, up)
    T = A > C.TIE
    n = len(pool)
    tri = []
    for i in range(n):
        for j in np.nonzero(T[i])[0]:
            if j < i:
                continue
            for k in np.nonzero(T[j] & T[:, i])[0]:
                if k > i:
                    tri.append((i, int(j), int(k)))
    tags = [(len(a), max(a)) for a in pool]
    same = [t for t in tri if tags[t[0]] == tags[t[1]] == tags[t[2]]]
    allint = 0
    for t in tri:
        ok = True
        for x in range(3):
            cen = extremum_census(pool[t[x]], pool[t[(x + 1) % 3]])
            if not (cen["max"]["interior"] and cen["min"]["interior"]):
                ok = False
                break
        allint += ok
    say(f"  {n} signatures, {len(tri)} directed 3-cycles")
    say(f"  inside one (r,M) bucket : {len(same)}")
    say(f"  all six extrema interior: {allint}")
    rows.append(dict(block="f11_census", n_sig=n, n_cycles=len(tri),
                     n_same_bucket=len(same), n_all_interior=allint))

    say()
    say("=" * 78)
    say("1d.  lambda_c of every F_11 3-cycle under five uniform priors")
    say("=" * 78)
    test = [("[-8,8]", *matched_uniform(-8, 8)),
            ("[-14,14]", *matched_uniform(-14, 14)),
            ("[-24,24]", *matched_uniform(-24, 24)),
            ("[-6,20]", *matched_uniform(-6, 20)),
            ("[-20,6]", *matched_uniform(-20, 6))]
    Ucache = [{} for _ in test]
    spreads, lcs_all = [], []
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else len(tri)
    for t in tri[:limit]:
        lcs = []
        for m, (_, s, w) in enumerate(test):
            for idx in t:
                if idx not in Ucache[m]:
                    Ucache[m][idx] = C.u_of(pool[idx], s)
            lcs.append(critical_lambda([Ucache[m][i] for i in t], w))
        good = [x for x in lcs if x]
        sp = max(good) / min(good) if len(good) == len(lcs) else float("nan")
        spreads.append(sp)
        lcs_all.append(lcs)
        rows.append(dict(block="f11_lambda_c", triple=str(t),
                         same_bucket=tags[t[0]] == tags[t[1]] == tags[t[2]],
                         **{f"lc_{nm}": x for (nm, _, _), x in zip(test, lcs)}, spread=sp))
    sp = np.array([x for x in spreads if x == x])
    say(f"  {'triple':>16s} " + " ".join(f"{nm:>12s}" for nm, _, _ in test) + f" {'spread':>9s}")
    for t, lcs, s_ in list(zip(tri, lcs_all, spreads))[:12]:
        say(f"  {str(t):>16s} " + " ".join(f"{(x if x else float('nan')):12.3f}" for x in lcs)
            + f" {s_:9.6f}")
    say(f"  ... {len(spreads)} cycles in all")
    say(f"  max_lambda_c / min_lambda_c over the five priors:  "
        f"median {np.median(sp):.6f}   max {sp.max():.6f}")
    say(f"  known 3-cycle, same five priors, for contrast:")
    kn = []
    for nm, s, w in test:
        Us = [C.u_of(a, s) for a in C.CYCLE]
        kn.append(critical_lambda(Us, w))
    say("      " + "  ".join(f"{nm}={x:.2f}" for (nm, _, _), x in zip(test, kn))
        + f"   spread {max(kn)/min(kn):.4f}")
    rows.append(dict(block="summary", median_spread=float(np.median(sp)),
                     max_spread=float(sp.max()), known_spread=max(kn) / min(kn)))

    with open("m1_lambda_c.csv", "w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    say("\nwrote m1_lambda_c.csv")


if __name__ == "__main__":
    main()
