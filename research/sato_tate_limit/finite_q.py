#!/usr/bin/env python3
"""Does the finite-``q`` comparison converge to the limiting one, and how fast?

Brief F insists that ``Psi_f -> Psi_mu`` be *demonstrated* before anything is
concluded from ``Psi_mu``, because the leading-order midrange law gets the sign
of a certified ``F_11`` edge wrong.  This script does that, and then measures the
two deviation exponents that decide when the limit takes over.

Exactly, for a map ``f : A^2 -> A^1`` over ``F_q`` with ``alpha_c = (N_c-q)/sqrt q``,

    Psi_f(tau) = log( (1/q) sum_c e^{tau alpha_c} ) / tau,
    Psi_f(0) = 0,   Psi_f(inf) = alpha_max(f) = (max_c N_c - q)/sqrt q,

and by Deligne equidistribution ``Psi_f -> Psi_mu`` pointwise, ``mu`` the trace
measure of the geometric monodromy group.

**Two deviation scales, with different exponents.**

*Bulk.*  ``K_f(tau) - K_mu(tau)`` is the error of a ``q``-sample empirical mean,
so it is ``O(q^{-1/2})`` at fixed ``tau``.

*Edge.*  ``Psi_f(inf) - Psi_mu(inf) = max_c alpha_c - 2g`` is an extreme-value
deficit.  The ``USp(2g)`` lower tail ``P(2g - tr < eps) ~ K_g eps^{d/2}`` with
``d = dim USp(2g) = g(2g+1)`` gives, over ``q`` samples,

    2g - max_c alpha_c  ~  Gamma(1 + 2/d) (K_g q)^{-2/d},

i.e. ``q^{-1/5}`` at genus 2 and ``q^{-2/21}`` at genus 3 -- recorded already as
T2.1's fitted exponents ``-0.1969`` and ``-0.0931``.  **The edge deviation
dominates the bulk one by ``q^{3/10}`` at genus two**, so the finite-``q``
``Psi_f`` differs from its limit mainly near ``tau = inf``, which is exactly
where ``phi`` reads.

Traces are computed exactly by FFT from
``N_c = sum_x (1 + chi(h(x)+c)) = q + (m_h * chi)[c]`` (the method of
``research/m_and_e_and_a_c/t2_1_genus_scaling.py``), and cross-checked against a
brute-force count at small ``q``.

    python research/sato_tate_limit/finite_q.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

import st_lib as S

HERE = Path(__file__).resolve().parent

TAU = np.geomspace(1e-3, 1e4, 701)


# ------------------------------------------------------------- finite fields


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_at_least(n: int) -> int:
    while not is_prime(n):
        n += 1
    return n


def legendre_table(q: int) -> np.ndarray:
    chi = -np.ones(q, dtype=np.int64)
    chi[(np.arange(1, q, dtype=np.int64) ** 2) % q] = 1
    chi[0] = 0
    return chi


def poly_values(coeffs: tuple[int, ...], q: int) -> np.ndarray:
    """``h(x) mod q`` for all ``x``; ``coeffs`` from the leading term down."""

    x = np.arange(q, dtype=np.int64)
    v = np.zeros(q, dtype=np.int64)
    for c in coeffs:
        v = (v * x + int(c)) % q
    return v


def counts_hyperelliptic(coeffs: tuple[int, ...], q: int) -> np.ndarray:
    """``N_c = #{(x,y) : y^2 = h(x) + c}`` for all ``c``, by FFT."""

    mult = np.bincount(poly_values(coeffs, q), minlength=q).astype(float)
    chi = legendre_table(q).astype(float)
    corr = np.fft.irfft(np.conj(np.fft.rfft(mult)) * np.fft.rfft(chi), n=q)
    return q + np.rint(corr).astype(np.int64)


def counts_hyperelliptic_brute(coeffs: tuple[int, ...], q: int) -> np.ndarray:
    h = poly_values(coeffs, q)
    y = np.arange(q, dtype=np.int64)
    vals = (y[:, None] ** 2 - h[None, :]) % q
    return np.bincount(vals.ravel(), minlength=q)


def counts_split(coeffs: tuple[int, ...], q: int) -> np.ndarray:
    """``y^2 = h(x^2) + c``: genus two with a split Jacobian ``E_1 x E_2``.

    The involution ``x -> -x`` gives the two elliptic quotients
    ``v^2 = h(u)+c`` and ``w^2 = u(h(u)+c)`` with ``u = x^2``, so the monodromy
    is ``SL_2 x SL_2`` and the limit measure is ``SU(2) x SU(2)``.
    """

    x = np.arange(q, dtype=np.int64)
    xs = (x * x) % q
    v = np.zeros(q, dtype=np.int64)
    for c in coeffs:
        v = (v * xs + int(c)) % q
    mult = np.bincount(v, minlength=q).astype(float)
    chi = legendre_table(q).astype(float)
    corr = np.fft.irfft(np.conj(np.fft.rfft(mult)) * np.fft.rfft(chi), n=q)
    return q + np.rint(corr).astype(np.int64)


def counts_isotrivial(q: int) -> np.ndarray:
    """``y^2 = x^3 + c``: the CM pencil, finite monodromy."""

    return counts_hyperelliptic((1, 0, 0, 0), q)


# ----------------------------------------------------------------- Psi and mid


def psi_finite(counts: np.ndarray, q: int, tau: np.ndarray) -> np.ndarray:
    """``Psi_f`` from a fibre signature, using the ``O(g sqrt q)`` distinct values."""

    alpha = (counts.astype(np.float64) - q) / math.sqrt(q)
    vals, mult = np.unique(alpha, return_counts=True)
    w = mult.astype(np.float64) / counts.size
    z = np.outer(tau, vals)
    m = z.max(axis=1, keepdims=True)
    return (np.log((w[None, :] * np.exp(z - m)).sum(axis=1)) + m[:, 0]) / tau


def mid_of(diff: np.ndarray, end: float) -> float:
    d = np.concatenate([[0.0], diff, [end]])
    return 0.5 * (float(d.max()) + float(d.min()))


def exact_mid_delta(sig_a: np.ndarray, sig_b: np.ndarray, q: int) -> float:
    """``mid_beta (log log Z_a - log log Z_b)`` -- the exact comparison.

    ``a < b`` iff this is negative (brief B).  Computed on a beta grid running
    to ``360 q``, the horizon that isolates the largest fibre.
    """

    betas = np.geomspace(1e-3, 360.0 * q, 40_000)

    def loglogz(sig: np.ndarray) -> np.ndarray:
        vals, mult = np.unique(sig, return_counts=True)
        ell = np.log(vals.astype(float))
        top = ell.max()
        lz = betas * top + np.log(mult.astype(float) @ np.exp(np.outer(ell - top,
                                                                      betas)))
        return np.log(lz)

    d = loglogz(sig_a) - loglogz(sig_b)
    end0 = math.log(math.log(len(sig_a))) - math.log(math.log(len(sig_b)))
    endinf = math.log(math.log(sig_a.max())) - math.log(math.log(sig_b.max()))
    dd = np.concatenate([[end0], d, [endinf]])
    return 0.5 * (float(dd.max()) + float(dd.min()))


# ------------------------------------------------------------------- families


def random_coeffs(rng, degree: int, q: int) -> tuple[int, ...]:
    return tuple([1] + [int(v) for v in rng.integers(0, q, size=degree)])


def main() -> int:
    rows: list[list] = []

    # ------------------------------------------------------------------ check
    print("=" * 78)
    print("0.  FFT fibre counts against brute force")
    print("=" * 78)
    ok = True
    for q in (101, 211):
        for coeffs in ((1, 0, 0, 1, 0, 0), (1, 2, 3, 4, 5, 6, 7, 8)):
            a = counts_hyperelliptic(coeffs, q)
            b = counts_hyperelliptic_brute(coeffs, q)
            same = np.array_equal(a, b)
            ok &= same and int(a.sum()) == q * q
            print(f"  q={q:>4} deg={len(coeffs)-1}: identical={same}  "
                  f"sum N_c = {int(a.sum())} (= q^2 = {q*q})")
    print(f"  all checks passed: {ok}\n")

    # -------------------------------------------------- 1. Psi_f -> Psi_mu
    print("=" * 78)
    print("1.  Psi_f converges to Psi_mu, and at what rate")
    print("=" * 78)
    limits = {
        "g1 generic": ("SU2", S.Measure((S.Factor("SU2"),))),
        "g2 generic": ("USp4", S.Measure((S.Factor("USp4"),))),
        "g3 generic": ("USp6", S.Measure((S.Factor("USp6"),))),
        "g2 split": ("SU2xSU2", S.Measure((S.Factor("SU2"), S.Factor("SU2")))),
    }
    builders = {
        "g1 generic": (lambda rng, q: counts_hyperelliptic(random_coeffs(rng, 3, q), q)),
        "g2 generic": (lambda rng, q: counts_hyperelliptic(random_coeffs(rng, 5, q), q)),
        "g3 generic": (lambda rng, q: counts_hyperelliptic(random_coeffs(rng, 7, q), q)),
        "g2 split": (lambda rng, q: counts_split(random_coeffs(rng, 3, q), q)),
    }
    qs = [prime_at_least(n) for n in (101, 401, 1601, 6421, 25601, 102407, 409601)]
    nfam = 12
    print(f"\n  {'family':<12}{'q':>8}{'m2':>9}{'limit m2':>10}"
          f"{'sup|dPsi|':>11}{'|dPsi| at tau=2':>16}{'edge deficit':>14}")
    for name, (lab, meas) in limits.items():
        psi_lim = meas.Psi(TAU)
        for q in qs:
            if q > 30000 and name == "g3 generic":
                pass
            rng = np.random.default_rng(1000 + q)
            sup = []
            bulk = []
            edge = []
            m2s = []
            for _ in range(nfam):
                counts = builders[name](rng, q)
                if counts.min() <= 0:
                    continue
                psi = psi_finite(counts, q, TAU)
                alpha = (counts - q) / math.sqrt(q)
                d = psi - psi_lim
                amax_end = float(alpha.max()) - meas.alpha_max
                sup.append(max(float(np.abs(d).max()), abs(amax_end)))
                bulk.append(abs(float(np.interp(2.0, TAU, d))))
                edge.append(amax_end)
                m2s.append(float((alpha ** 2).mean()))
            print(f"  {name:<12}{q:>8}{np.mean(m2s):>9.4f}"
                  f"{meas.variance:>10.4f}{np.mean(sup):>11.5f}"
                  f"{np.mean(bulk):>16.6f}{np.mean(edge):>14.5f}")
            rows.append(["convergence", name, lab, q, f"{np.mean(m2s):.6f}",
                         f"{meas.variance:.6f}", f"{np.mean(sup):.6f}",
                         f"{np.mean(bulk):.6f}", f"{np.mean(edge):.6f}"])
        print()

    print("  fitted exponents  d ~ q^(-e), against the general edge law")
    print("  E[alpha_max(mu) - max_c alpha_c] ~ Gamma(1 + 1/t) (c q)^(-1/t):")
    print(f"  {'family':<12}{'bulk e':>10}{'edge e':>10}{'t':>7}{'1/t':>10}")
    for name, (lab, meas) in limits.items():
        sel = [r for r in rows if r[0] == "convergence" and r[1] == name]
        lq = np.log([float(r[3]) for r in sel])
        eb = -np.polyfit(lq, np.log([float(r[7]) for r in sel]), 1)[0]
        ee = -np.polyfit(lq, np.log([abs(float(r[8])) for r in sel]), 1)[0]
        pred = 1.0 / meas.tail
        print(f"  {name:<12}{eb:>10.4f}{ee:>10.4f}{meas.tail:>7.1f}{pred:>10.4f}")
        rows.append(["exponent", name, lab, "", f"{eb:.6f}", f"{ee:.6f}",
                     f"{pred:.6f}", f"{meas.tail:.4f}", ""])
    print("\n  T2.1 quotes 2/dim USp(2g) = 2/(g(2g+1)); that is the special case")
    print("  t = dim/2 of the general law, which also covers the split family.")

    with (HERE / "convergence.csv").open("w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["kind", "family", "limit", "q", "m2", "limit_m2",
                     "sup_dPsi", "bulk_dPsi_tau2", "edge_deficit"])
        wr.writerows(rows)
    print("\nwritten: convergence.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
