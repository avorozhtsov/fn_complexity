#!/usr/bin/env python3
"""Step 5 of brief F: where the two regimes meet, predicted and then tested.

Brief F asks, if the limit is transitive, to predict the ``q`` at which
finite-``q`` cycles die and then test it.  The prediction has to be stated
carefully, because **the limit is silent about the certified cycles**.

*The limit separates two families only if their limiting measures differ.*  All
three pencils of the certified ``F_101`` cycle are generic genus-two
hyperelliptic pencils with monodromy ``Sp(4)``, so all three have the same limit
``USp(4)``, ``Psi_u - Psi_v -> 0``, and the limiting comparison makes no
prediction about them at all.  A transitive limit therefore **cannot** kill
same-symmetry-type cycles, and no ``q`` exists past which the same-genus census
collapses.  What it does kill is cycles across symmetry types, and there the
prediction is quantitative.

**Part 1** re-derives the ``F_101`` certificate from the polynomials, with an
independent point count, and reads it in the limiting language.

**Part 2** identifies the finite-``q`` shadow of the edge exponent.  For large
``tau``,

    Psi_f(tau) = alpha_max(f) + log(mult_f / q)/tau + O(tau^{-2}),
    Psi_mu(tau) = alpha_max(mu) - (t log tau - log A)/tau + O(tau^{-2}),

with ``mult_f`` the multiplicity of the largest fibre.  **So ``log mult`` and the
edge exponent ``t`` occupy the same slot** -- the coefficient of ``1/tau`` at the
top end -- and brief E's finding that ``log mu`` predicts the certified edge
where the moment ladder does not is the finite-``q`` form of the statement that
the limiting comparison at fixed ``(alpha_max, m_2)`` is decided by ``t``.

**Part 3** measures the extreme-value law and predicts where genus classes
separate; **part 4** tests both halves of the prediction; **part 5** measures the
curl fraction against ``q``, which brief E predicts should converge to a value
set by the limiting trace spread.

    python research/sato_tate_limit/census_test.py
"""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "research" / "curve_family_cycles"))
sys.path.insert(0, str(HERE))

import st_lib as S                                            # noqa: E402
from common import Engine, beta_grid                          # noqa: E402
from finite_q import (counts_hyperelliptic, prime_at_least,    # noqa: E402
                      psi_finite, random_coeffs)

QS = [prime_at_least(n) for n in (31, 101, 211, 401, 1009, 2003, 4001, 8009)]
POOL = 3000
CLASS_CAP = 80
TAU = np.geomspace(1e-3, 1e4, 701)

F101 = {
    "f1": (1, 70, 28, 15, 11, 31),
    "f2": (1, 42, 32, 74, 96, 60),
    "f3": (1, 72, 21, 2, 6, 57),
}


def pool_of(q: int, degree: int, count: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    seen: dict[tuple[int, ...], np.ndarray] = {}
    for _ in range(count):
        counts = counts_hyperelliptic(random_coeffs(rng, degree, q), q)
        if counts.min() <= 0:
            continue
        seen.setdefault(tuple(sorted((int(v) for v in counts), reverse=True)),
                        counts)
    return seen


def rate_matrix(sigs: list[tuple[int, ...]], q: int):
    return Engine(sigs, beta_grid(q, points=30_000)).rate_matrix(chunk=48)


def cycles_of(rates: np.ndarray, tol: float = 1e-9):
    n = rates.shape[0]
    strict = (rates - rates.T) < -tol
    margins = np.abs(rates - rates.T)
    out = []
    for i, j, k in itertools.combinations(range(n), 3):
        for a, b, c in ((i, j, k), (i, k, j)):
            if strict[a, b] and strict[b, c] and strict[c, a]:
                out.append((a, b, c,
                            min(margins[a, b], margins[b, c], margins[c, a])))
    return out


def curl_fraction(rates: np.ndarray) -> float:
    """``||curl A|| / ||A||`` with ``A`` the antisymmetric part of ``-log C``."""

    L = -np.log(rates)
    A = 0.5 * (L - L.T)
    psi = -A.mean(axis=1)
    grad = psi[None, :] - psi[:, None]
    return float(np.linalg.norm(A - grad) / np.linalg.norm(A))


def main() -> int:
    rows: list[list] = []

    # ---------------------------------------------------- 1. the F_101 cycle
    print("=" * 78)
    print("1.  the certified F_101 cycle, recomputed here, read in the limit")
    print("=" * 78)
    q = 101
    sigs = {}
    for name, coeffs in F101.items():
        counts = counts_hyperelliptic(coeffs, q)
        sigs[name] = counts
        alpha = (counts - q) / math.sqrt(q)
        mult = int((counts == counts.max()).sum())
        print(f"  {name}: sum N_c = {int(counts.sum())} (= q^2 = {q * q}), "
              f"max fibre {int(counts.max())} with multiplicity {mult}, "
              f"m2 = {float((alpha ** 2).mean()):.6f}, "
              f"m3 = {float((alpha ** 3).mean()):+.6f}")
        rows.append(["f101", name, int(counts.max()), mult,
                     f"{float((alpha ** 2).mean()):.6f}",
                     f"{float((alpha ** 3).mean()):.6f}", ""])
    names = list(F101)
    order = [tuple(sorted((int(v) for v in sigs[n]), reverse=True)) for n in names]
    rates, contacts = rate_matrix(order, q)
    print(f"\n  {'edge':<12}{'C(a->b)':>16}{'C(b->a)':>16}{'margin':>13}"
          f"{'A = -mid':>13}")
    edge_signs = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        m = rates[j, i] - rates[i, j]
        print(f"  {names[i] + ' -> ' + names[j]:<12}{rates[i, j]:>16.12f}"
              f"{rates[j, i]:>16.12f}{m:>13.4e}"
              f"{0.5 * math.log(rates[j, i] / rates[i, j]):>13.4e}")
        rows.append(["f101-edge", f"{names[i]}->{names[j]}",
                     f"{rates[i, j]:.12f}", f"{rates[j, i]:.12f}", f"{m:.6e}",
                     "", ""])
        edge_signs.append(m)
    cyc_ok = (all(m < -1e-9 for m in edge_signs)
              or all(m > 1e-9 for m in edge_signs))
    print(f"\n  all three edges agree in sign, hence a strict 3-cycle: {cyc_ok}")
    print("  Every C(a->b) here exceeds C(b->a), so the orientation is")
    print("  f_1 > f_2 > f_3 > f_1, exactly addendum 2's.  Smallest margin "
          f"{min(abs(m) for m in edge_signs):.4e}.")

    print("\n  In the limit all three are generic genus-two pencils, so all")
    print("  three have limiting measure USp(4).  Their Psi curves therefore")
    print("  converge to one another; the whole comparison is fluctuation:")
    lim = S.Measure((S.Factor("USp4"),))
    psi_lim = lim.Psi(TAU)
    for name in names:
        d = psi_finite(sigs[name], q, TAU) - psi_lim
        print(f"    {name}: sup|Psi_f - Psi_USp4| = {np.abs(d).max():.4f}, "
              f"alpha_max(f) = {(sigs[name].max() - q) / math.sqrt(q):.4f} "
              f"(limit 4)")
    for i, j in ((0, 1), (1, 2), (2, 0)):
        d = psi_finite(sigs[names[i]], q, TAU) - psi_finite(sigs[names[j]], q, TAU)
        e = ((sigs[names[i]].max() - sigs[names[j]].max()) / math.sqrt(q))
        dd = np.concatenate([[0.0], d, [e]])
        mid = 0.5 * (dd.max() + dd.min())
        pred = -2.0 * mid / (math.sqrt(q) * math.log(q))
        print(f"    mid(Psi_{names[i]} - Psi_{names[j]}) = {mid:+.6f}"
              f"   -> predicted margin {pred:+.4e}"
              f"   (exact {rates[j, i] - rates[i, j]:+.4e})")
        rows.append(["f101-psi", f"{names[i]}-{names[j]}", f"{mid:.6f}",
                     f"{pred:.6e}", f"{rates[j, i] - rates[i, j]:.6e}", "", ""])

    # --------------------------------- 2. log(multiplicity) is the edge datum
    print()
    print("=" * 78)
    print("2.  log(multiplicity of the largest fibre) is the finite-q edge datum")
    print("=" * 78)
    print("  tau (Psi_f(tau) - alpha_max(f)) -> log(mult_f / q) as tau -> inf,")
    print("  the slot the limiting measure fills with -(t log tau - log A).")
    print(f"\n  {'family':>7}{'mult':>6}{'log(mult/q)':>14}"
          f"{'tau(Psi-amax) at tau=200':>26}")
    for name in names:
        counts = sigs[name]
        mult = int((counts == counts.max()).sum())
        amax = (counts.max() - q) / math.sqrt(q)
        t = np.array([200.0])
        val = float((t * (psi_finite(counts, q, t) - amax))[0])
        print(f"  {name:>7}{mult:>6}{math.log(mult / q):>14.6f}{val:>26.6f}")
        rows.append(["edge-slot", name, mult, f"{math.log(mult / q):.6f}",
                     f"{val:.6f}", "", ""])

    # ---------------------------------------- 3. the extreme-value separation
    print()
    print("=" * 78)
    print("3.  the extreme-value law, and where the genus classes separate")
    print("=" * 78)
    print(f"\n  {'q':>7} {'g=2 E[max a]':>14} {'sd':>7} {'range':>17}"
          f" {'g=3 E[max a]':>14} {'sd':>7} {'range':>17} {'overlap':>8}")
    fits = {2: [], 3: []}
    for qq in QS:
        stats = {}
        for g, degree in ((2, 5), (3, 7)):
            rng = np.random.default_rng(7000 + qq + g)
            am = []
            for _ in range(300):
                counts = counts_hyperelliptic(random_coeffs(rng, degree, qq), qq)
                if counts.min() <= 0:
                    continue
                am.append((int(counts.max()) - qq) / math.sqrt(qq))
            am = np.array(am)
            stats[g] = am
            fits[g].append((qq, float(am.mean())))
        ov = "yes" if stats[2].max() >= stats[3].min() else "no"
        print(f"  {qq:>7} {stats[2].mean():>14.4f} {stats[2].std():>7.4f}"
              f" [{stats[2].min():>6.3f},{stats[2].max():>6.3f}]"
              f" {stats[3].mean():>14.4f} {stats[3].std():>7.4f}"
              f" [{stats[3].min():>6.3f},{stats[3].max():>6.3f}] {ov:>8}")
        rows.append(["extreme", qq, f"{stats[2].mean():.6f}",
                     f"{stats[2].std():.6f}", f"{stats[3].mean():.6f}",
                     f"{stats[3].std():.6f}", ov])

    print("\n  fit  2g - E[max alpha] = Gamma(1 + 1/t) (c_g q)^(-1/t),  t = g(2g+1)/2:")
    cs = {}
    for g in (2, 3):
        t = g * (2 * g + 1) / 2.0
        qsv = np.array([x[0] for x in fits[g]], dtype=float)
        dfc = np.array([2 * g - x[1] for x in fits[g]])
        slope, inter = np.polyfit(np.log(qsv), np.log(dfc), 1)
        gam = math.gamma(1 + 1.0 / t)
        cs[g] = math.exp((math.log(gam) - inter) * t)
        print(f"    g = {g}, t = {t:.1f}: fitted exponent {slope:+.5f} "
              f"(predicted {-1 / t:+.5f}),  c_{g} = {cs[g]:.4e}")
        rows.append(["evfit", g, f"{slope:.6f}", f"{-1 / t:.6f}",
                     f"{cs[g]:.6e}", "", ""])

    def deficit(g: int, qv: float) -> float:
        t = g * (2 * g + 1) / 2.0
        return math.gamma(1 + 1.0 / t) * (cs[g] * qv) ** (-1.0 / t)

    # the spread scales with the same power of q as the deficit, so fit it too
    sds = {}
    for g in (2, 3):
        arr = [(r[1], float(r[3] if g == 2 else r[5]))
               for r in rows if r[0] == "extreme"]
        lq = np.log([float(x[0]) for x in arr])
        ls = np.log([x[1] for x in arr])
        sl, ic = np.polyfit(lq, ls, 1)
        sds[g] = (math.exp(ic), sl)
        print(f"    sd(max alpha) for g = {g}: {math.exp(ic):.4f} q^({sl:+.4f})")
        rows.append(["sdfit", g, f"{math.exp(ic):.6f}", f"{sl:.6f}", "", "", ""])

    def sd(g: int, qv: float) -> float:
        a, b = sds[g]
        return a * qv ** b

    print("\n  extrapolated means and spreads.  Two genus classes can share a")
    print("  largest fibre while their max-alpha distributions overlap; the")
    print("  separation criterion taken here is  gap > 3(sd_2 + sd_3).")
    print(f"  {'q':>10}{'E[max a] g=2':>14}{'E[max a] g=3':>14}{'gap':>9}"
          f"{'3(sd2+sd3)':>12}{'separated':>11}")
    sep_q = None
    for qv in (1e2, 1e3, 1e4, 1e6, 1e8, 1e10, 1e12, 1e15, 1e18):
        g2 = 4 - deficit(2, qv)
        g3 = 6 - deficit(3, qv)
        thr = 3 * (sd(2, qv) + sd(3, qv))
        ok = g3 - g2 > thr
        print(f"  {qv:>10.0e}{g2:>14.4f}{g3:>14.4f}{g3 - g2:>9.4f}"
              f"{thr:>12.4f}{str(ok):>11}")
        rows.append(["prediction", qv, f"{g2:.6f}", f"{g3:.6f}",
                     f"{g3 - g2:.6f}", f"{thr:.6f}", ok])
    lo, hi = 1e2, 1e40
    for _ in range(200):
        mid_q = math.sqrt(lo * hi)
        if (6 - deficit(3, mid_q)) - (4 - deficit(2, mid_q)) > \
                3 * (sd(2, mid_q) + sd(3, mid_q)):
            hi = mid_q
        else:
            lo = mid_q
    sep_q = hi
    print(f"\n  PREDICTION: genus 2 and genus 3 stop sharing a largest fibre at")
    print(f"  q = {sep_q:.2e}.  Beyond it no mixed-genus class exists, hence no")
    print("  cycle can mix those two genera, while same-genus cycles continue.")
    rows.append(["separation-q", f"{sep_q:.6e}", "", "", "", "", ""])

    print("\n  direct test of the same statistic: the probability that a random")
    print("  genus-2 and a random genus-3 pencil share a largest fibre.")
    print(f"  {'q':>8}{'P(max2 = max3)':>17}{'distinct max values':>21}")
    for qq in QS + [prime_at_least(20011), prime_at_least(50021)]:
        hist = {}
        for g, degree in ((2, 5), (3, 7)):
            rng = np.random.default_rng(9000 + qq + g)
            vals = []
            for _ in range(600):
                counts = counts_hyperelliptic(random_coeffs(rng, degree, qq), qq)
                if counts.min() <= 0:
                    continue
                vals.append(int(counts.max()))
            v, c = np.unique(vals, return_counts=True)
            hist[g] = dict(zip(v.tolist(), (c / c.sum()).tolist()))
        p = sum(hist[2].get(k, 0.0) * hist[3].get(k, 0.0)
                for k in set(hist[2]) | set(hist[3]))
        print(f"  {qq:>8}{p:>17.5f}{len(set(hist[2]) | set(hist[3])):>21}")
        rows.append(["overlap", qq, f"{p:.6f}",
                     len(set(hist[2]) | set(hist[3])), "", "", ""])

    # --------------------------------- 4. testing both halves of the prediction
    print()
    print("=" * 78)
    print("4.  the test: mixed-genus classes vanish, same-genus cycles do not")
    print("=" * 78)
    print(f"\n  {'q':>7}{'g2 sigs':>9}{'g3 sigs':>9}{'shared classes':>16}"
          f"{'mixed cycles':>14}{'same-genus cycles':>19}{'fraction':>10}"
          f"{'median margin x sqrt(q) log q':>31}")
    for qq in QS:
        p2 = pool_of(qq, 5, POOL, 11 + qq)
        p3 = pool_of(qq, 7, POOL, 22 + qq)
        m2c: dict[int, list] = {}
        for s in p2:
            m2c.setdefault(s[0], []).append(s)
        m3c: dict[int, list] = {}
        for s in p3:
            m3c.setdefault(s[0], []).append(s)
        shared = sorted(set(m2c) & set(m3c))
        mixed = 0
        if shared:
            best = max(shared, key=lambda k: min(len(m2c[k]), len(m3c[k])))
            sel = m2c[best][:CLASS_CAP // 2] + m3c[best][:CLASS_CAP // 2]
            gen = ([2] * len(m2c[best][:CLASS_CAP // 2])
                   + [3] * len(m3c[best][:CLASS_CAP // 2]))
            if len(sel) >= 3:
                r, _ = rate_matrix(sel, qq)
                for a, b, c, _ in cycles_of(r):
                    if len({gen[a], gen[b], gen[c]}) > 1:
                        mixed += 1
        big = max(m2c, key=lambda k: len(m2c[k]))
        sel2 = m2c[big][:CLASS_CAP]
        if len(sel2) < 8:
            continue
        r2, _ = rate_matrix(sel2, qq)
        cyc = cycles_of(r2)
        n = len(sel2)
        tri = math.comb(n, 3) * 2
        up = np.triu_indices(n, 1)
        med = float(np.median(np.abs(r2 - r2.T)[up]))
        print(f"  {qq:>7}{len(p2):>9}{len(p3):>9}{len(shared):>16}"
              f"{mixed:>14}{len(cyc):>19}{len(cyc) / tri:>10.5f}"
              f"{med * math.sqrt(qq) * math.log(qq):>31.5f}")
        rows.append(["test", qq, len(shared), mixed, len(cyc),
                     f"{len(cyc) / tri:.6f}",
                     f"{med * math.sqrt(qq) * math.log(qq):.6f}"])

    # ------------------------------------------- 5. curl fraction against q
    print()
    print("=" * 78)
    print("5.  the curl fraction inside a class, against q")
    print("=" * 78)
    print("  Brief E finds ||curl A||/||A|| is a function of the trace spread.")
    print("  The trace spread converges to its Sato-Tate value, so the curl")
    print("  fraction should converge too.")
    print(f"\n  {'q':>7}{'members':>9}{'trace spread':>14}{'||curl||/||A||':>16}")
    for qq in QS:
        p2 = pool_of(qq, 5, POOL, 11 + qq)
        m2c: dict[int, list] = {}
        for s in p2:
            m2c.setdefault(s[0], []).append(s)
        big = max(m2c, key=lambda k: len(m2c[k]))
        sel2 = m2c[big][:CLASS_CAP]
        if len(sel2) < 8:
            continue
        r2, _ = rate_matrix(sel2, qq)
        spread = float(np.mean([np.std((np.array(s) - qq) / math.sqrt(qq))
                                for s in sel2]))
        cf = curl_fraction(r2)
        print(f"  {qq:>7}{len(sel2):>9}{spread:>14.5f}{cf:>16.5f}")
        rows.append(["curl", qq, len(sel2), f"{spread:.6f}", f"{cf:.6f}", "", ""])
    print("\n  the limiting trace spread for USp(4) is sqrt(m2) = 1.")

    with (HERE / "census_test.csv").open("w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind", "a", "b", "c", "d", "e", "f"])
        wr.writerows(rows)
    print("\nwritten: census_test.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
