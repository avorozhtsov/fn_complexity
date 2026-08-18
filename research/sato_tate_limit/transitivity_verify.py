#!/usr/bin/env python3
"""The headline numbers of ``TRANSITIVITY.md``, recomputed at 40 digits.

Nothing here uses a grid for its final answer: every extremum is polished by
golden section on ``log tau`` and every value is an ``mpmath`` evaluation of the
Weyl/Andreief determinant at working precision.  Both endpoints are supplied
analytically -- ``Psi(0) = 0`` because every limit measure has mean zero, and
``Psi(inf) = alpha_max`` -- so the two contacts a grid can never see are exact.

Verified here:

* the same-``alpha_max`` 3-cycle among mixtures at genus 8 (the counterexample
  to "same alpha_max => total order"),
* the genus-11 cycle with a pure-product vertex,
* the tightest same-genus triple of the full symplectic library, whose
  certificate is the number ``0.0805...``,
* the concavity margin ``1 - 4 b_g(tau)^2`` at the ranks and arguments that
  carry the dominance theorem.

    python research/sato_tate_limit/transitivity_verify.py
"""

from __future__ import annotations

import csv
from pathlib import Path

from mpmath import mp, mpf

import kappa_lib as KL

HERE = Path(__file__).resolve().parent
mp.dps = 60
EXTRA = 30


class Prod:
    """``prod_i USp(2 lambda_i)`` as an exact ``K``."""

    def __init__(self, lam):
        self.lam = tuple(lam)
        self.alpha_max = mpf(2 * sum(lam))
        self._z = {}
        for g in set(lam):
            dps = KL.working_dps(g + 1, 1.0) + EXTRA
            h = KL.hankels(0.0, g, dps)
            with mp.workdps(dps):
                self._z[g] = mp.log(h[g])

    def K(self, tau):
        t = mpf(tau)
        out = mpf(0)
        for g in self.lam:
            dps = KL.working_dps(g + 1, float(t)) + EXTRA
            h = KL.hankels(t, g, dps)
            with mp.workdps(dps):
                out += mp.log(h[g]) - self._z[g]
        return out

    @property
    def label(self):
        return "".join(str(p) for p in self.lam)


class Mix:
    """``p mu + (1-p) nu`` for two products."""

    def __init__(self, a: Prod, b: Prod, p="0.5"):
        self.a, self.b, self.p = a, b, mpf(p)
        self.alpha_max = max(a.alpha_max, b.alpha_max)

    def K(self, tau):
        ka, kb = self.a.K(tau), self.b.K(tau)
        m = max(ka, kb)
        return m + mp.log(self.p * mp.e ** (ka - m)
                          + (1 - self.p) * mp.e ** (kb - m))

    @property
    def label(self):
        p = mp.nstr(self.p, 6)
        return f"[{p}*{self.a.label} + {mp.nstr(1 - self.p, 6)}*{self.b.label}]"


def D(u, v, tau):
    t = mpf(tau)
    return (u.K(t) - v.K(t)) / t


def golden(f, lo, hi, iters=90, maximise=True):
    """golden-section search on ``log tau``; returns ``(tau, value)``."""

    phi = (mp.sqrt(5) - 1) / 2
    a, b = mp.log(mpf(lo)), mp.log(mpf(hi))
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = f(mp.e ** c), f(mp.e ** d)
    for _ in range(iters):
        better = fc > fd if maximise else fc < fd
        if better:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = f(mp.e ** c)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = f(mp.e ** d)
        if b - a < mpf(10) ** (-30):
            break
    x = mp.e ** ((a + b) / 2)
    return x, f(x)


def midrange(u, v, brackets_hi, brackets_lo):
    """``mid``, with both endpoints analytic and the interior polished."""

    ends = [mpf(0), u.alpha_max - v.alpha_max]
    hi = max(ends)
    lo = min(ends)
    at_hi, at_lo = "endpoint", "endpoint"
    for br in brackets_hi:
        t, val = golden(lambda s: D(u, v, s), br[0], br[1], maximise=True)
        if val > hi:
            hi, at_hi = val, t
    for br in brackets_lo:
        t, val = golden(lambda s: D(u, v, s), br[0], br[1], maximise=False)
        if val < lo:
            lo, at_lo = val, t
    return (hi + lo) / 2, hi, lo, at_hi, at_lo


BRACKETS_HI = [(mpf("0.02"), mpf(2)), (mpf(2), mpf(40)),
               (mpf(40), mpf("2e3")), (mpf("2e3"), mpf("1e6"))]
BRACKETS_LO = BRACKETS_HI


def report(name, cyc, rows):
    print(f"\n  {name}")
    worst = None
    for u, v in zip(cyc, cyc[1:] + cyc[:1]):
        m, hi, lo, thi, tlo = midrange(u, v, BRACKETS_HI, BRACKETS_LO)
        print(f"    mid({u.label} - {v.label})")
        print(f"        = {mp.nstr(m, 40)}")
        print(f"        sup {mp.nstr(hi, 20)} at tau = "
              f"{thi if isinstance(thi, str) else mp.nstr(thi, 10)}")
        print(f"        inf {mp.nstr(lo, 20)} at tau = "
              f"{tlo if isinstance(tlo, str) else mp.nstr(tlo, 10)}")
        rows.append([name, u.label, v.label, mp.nstr(m, 40),
                     mp.nstr(hi, 20), mp.nstr(lo, 20),
                     thi if isinstance(thi, str) else mp.nstr(thi, 12),
                     tlo if isinstance(tlo, str) else mp.nstr(tlo, 12)])
        worst = m if worst is None else max(worst, m)
    print(f"    margin = {mp.nstr(-worst, 25)}   "
          f"({'3-CYCLE' if worst < 0 else 'transitive'})")
    return worst


def main() -> int:
    rows = []
    print("=" * 78)
    print("headline numbers at 40 digits")
    print("=" * 78)

    # --- genus 8: the counterexample --------------------------------------
    p8 = Prod((8,))
    a = Mix(p8, Prod((4, 1, 1, 1, 1)))
    b = Mix(p8, Prod((2, 2, 2, 2)))
    c = Mix(Prod((5, 1, 1, 1)), Prod((3, 2, 2, 1)))
    report("genus 8 mixture 3-cycle (alpha_max = 16)", [a, b, c], rows)

    # --- genus 11: a cycle with a pure-product vertex ----------------------
    a = Prod((8, 1, 1, 1))
    b = Mix(Prod((9, 1, 1)), Prod((6, 4, 1)))
    c = Mix(Prod((8, 3)), Prod((7, 2, 2)))
    report("genus 11 3-cycle with a pure vertex (alpha_max = 22)",
           [a, b, c], rows)

    # --- genus 11 near-miss among pure products ---------------------------
    report("genus 11 tightest pure-product triple (transitive)",
           [Prod((6, 2, 1, 1, 1)), Prod((4, 3, 3, 1)), Prod((5, 2, 2, 2))],
           rows)

    # --- the first crossing pair of the multiplicity-free cone ------------
    print("\n  genus 7: the first multiplicity-free crossing pair")
    m, hi, lo, thi, tlo = midrange(Prod((5, 1, 1)), Prod((4, 3)),
                                   BRACKETS_HI, BRACKETS_LO)
    print(f"    mid(511 - 43) = {mp.nstr(m, 40)}")
    print(f"        sup {mp.nstr(hi, 20)} at tau = "
          f"{thi if isinstance(thi, str) else mp.nstr(thi, 10)}")
    print(f"        inf {mp.nstr(lo, 20)} at tau = "
          f"{tlo if isinstance(tlo, str) else mp.nstr(tlo, 10)}")
    rows.append(["genus 7 first crossing", "511", "43", mp.nstr(m, 40),
                 mp.nstr(hi, 20), mp.nstr(lo, 20),
                 thi if isinstance(thi, str) else mp.nstr(thi, 12),
                 tlo if isinstance(tlo, str) else mp.nstr(tlo, 12)])

    # --- the concavity margin ---------------------------------------------
    print("\n  concavity margin  1 - 4 b_g(tau)^2  (must be > 0 for tau > 0)")
    for g in (1, 2, 3, 4, 5, 6, 8, 11):
        line = []
        for t in ("0.1", "1", "5", "30"):
            tt = mpf(t)
            dps = KL.working_dps(g + 2, float(tt)) + 60
            h = KL.hankels(tt, g + 1, dps)
            with mp.workdps(dps):
                val = 1 - 4 * h[g + 1] * h[g - 1] / h[g] ** 2
            line.append(mp.nstr(val, 12))
            rows.append(["concavity", f"g={g}", f"tau={t}",
                         mp.nstr(val, 40), "", "", "", ""])
        print(f"    g = {g:>2}:  " + "   ".join(f"{s:>16}" for s in line))

    with (HERE / "transitivity_verify.csv").open("w", newline="",
                                                 encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["item", "mu", "nu", "mid", "sup", "inf", "tau_sup",
                     "tau_inf"])
        wr.writerows(rows)
    print("\nwritten: transitivity_verify.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
