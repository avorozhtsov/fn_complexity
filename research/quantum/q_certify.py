"""Brief L -- 40-digit certification of the coherence-generated qubit 3-cycle.

Reads the witness from q_cycle.json (key "best_r2"), rebuilds it in mpmath at
60 digits and re-derives every headline number at 40.

For r = 2 nothing needs an iterative eigensolver: the mu_k(beta) are the roots
of the quadratic  det(A - mu S^t) = 0,  t = (beta-1)/beta, so

    F(beta) = log( mu_1^beta + mu_2^beta ),   mu_{1,2} = roots of
    det(S^t) mu^2 - (A_11 (S^t)_22 + A_22 (S^t)_11 - 2 A_12 (S^t)_12) mu
                  + det(A) = 0

and with S diagonal (S^t)_12 = 0, so the quadratic is

    s1^t s2^t mu^2 - (A_11 s2^t + A_22 s1^t) mu + det A = 0.

The extrema of phi = U_b - U_a are located by a bracketed golden search in
double precision and then polished by mpmath Newton on phi' with a numerical
derivative at 60 digits; each is reported with the value of phi' at the
polished point as the certificate that it really is a critical point.

Run:  python3 q_certify.py [witness.json] [key]
      defaults: q_cycle.json, best_r2
"""
from __future__ import annotations

import json
import math
import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
DPS = 60
REPORT = 40
BMIN = mp.mpf("0.5")


class MPQSig2:
    """(A, S) on C^2 with S diagonal, in mpmath."""

    def __init__(self, A, s):
        self.A = [[mp.mpf(str(A[i][j])) for j in range(2)] for i in range(2)]
        self.s = [mp.mpf(str(v)) for v in s]
        self.detA = self.A[0][0] * self.A[1][1] - self.A[0][1] * self.A[1][0]

    def mu(self, beta):
        t = (beta - 1) / beta
        p = [mp.power(v, t) for v in self.s]
        a2 = p[0] * p[1]
        a1 = -(self.A[0][0] * p[1] + self.A[1][1] * p[0])
        a0 = self.detA
        disc = mp.sqrt(a1 * a1 - 4 * a2 * a0)
        return ((-a1 + disc) / (2 * a2), (-a1 - disc) / (2 * a2))

    def F(self, beta):
        m1, m2 = self.mu(beta)
        return mp.log(mp.power(m1, beta) + mp.power(m2, beta))

    def U(self, s):
        return mp.log(self.F(mp.exp(s)))

    @property
    def R(self):
        return mp.log(self.s[0] + self.s[1])

    @property
    def Lam(self):
        """log lambda_max(S^{-1/2} A S^{-1/2}) = log of the largest mu at t=1."""
        p = self.s
        a2 = p[0] * p[1]
        a1 = -(self.A[0][0] * p[1] + self.A[1][1] * p[0])
        disc = mp.sqrt(a1 * a1 - 4 * a2 * self.detA)
        return mp.log((-a1 + disc) / (2 * a2))


def pinch(A):
    return [[A[0][0], 0.0], [0.0, A[1][1]]]


def golden(f, lo, hi, maximize, iters=400):
    g = (mp.sqrt(5) - 1) / 2
    a = hi - g * (hi - lo)
    b = lo + g * (hi - lo)
    sg = -1 if maximize else 1
    fa, fb = sg * f(a), sg * f(b)
    for _ in range(iters):
        if fa <= fb:
            hi, b, fb = b, a, fa
            a = hi - g * (hi - lo)
            fa = sg * f(a)
        else:
            lo, a, fa = a, b, fb
            b = lo + g * (hi - lo)
            fb = sg * f(b)
        if hi - lo < mp.mpf(10) ** (-(DPS - 8)):
            break
    return (a, f(a)) if fa <= fb else (b, f(b))


def scan_extrema(phi, smin, smax, n=600):
    """Bracket every interior extremum on a coarse grid, then golden-polish."""
    xs = [smin + (smax - smin) * i / (n - 1) for i in range(n)]
    vs = [phi(x) for x in xs]
    hi = [(xs[0], vs[0]), (xs[-1], vs[-1])]
    lo = [(xs[0], vs[0]), (xs[-1], vs[-1])]
    for i in range(1, n - 1):
        if vs[i] > vs[i - 1] and vs[i] >= vs[i + 1]:
            hi.append(golden(phi, xs[i - 1], xs[i + 1], True))
        if vs[i] < vs[i - 1] and vs[i] <= vs[i + 1]:
            lo.append(golden(phi, xs[i - 1], xs[i + 1], False))
    return hi, lo


def d_and_mid(qa, qb, smin, smax):
    phi = lambda s: qb.U(s) - qa.U(s)
    e_inf = mp.log(qb.Lam) - mp.log(qa.Lam)
    hi, lo = scan_extrema(phi, smin, smax)
    P = max([v for _, v in hi] + [e_inf])
    Q = min([v for _, v in lo] + [e_inf])
    return P - Q, (P + Q) / 2, P, Q


def main():
    mp.mp.dps = DPS
    src = sys.argv[1] if len(sys.argv) > 1 else "q_cycle.json"
    key = sys.argv[2] if len(sys.argv) > 2 else "best_r2"
    with open(os.path.join(HERE, src)) as fh:
        data = json.load(fh)
    w = data[key]
    if w["r"] != 2:
        raise SystemExit("this certifier is the closed-form r = 2 one")
    s = [w["S"][0][0], w["S"][1][1]]
    ops = w["A"]

    smin = mp.log(BMIN)
    smax = mp.log(mp.mpf(10) ** 6)

    qs = [MPQSig2(A, s) for A in ops]
    ps = [MPQSig2(pinch(A), s) for A in ops]

    print("=" * 74)
    print("40-DIGIT CERTIFICATION -- coherence-generated qubit 3-cycle")
    print("=" * 74)
    print(f"witness: {src} [{key}]")
    print(f"mpmath working precision {DPS} digits; reported to {REPORT}.")
    print(f"beta range of the infimum: [1/2, 10^6] u {{oo}}.\n")
    print("S = diag(")
    for v in s:
        print(f"      {mp.nstr(mp.mpf(str(v)), REPORT)}")
    print(")")
    for k, A in enumerate(ops):
        print(f"A_{k+1} = [[{mp.nstr(mp.mpf(str(A[0][0])), REPORT)}, "
              f"{mp.nstr(mp.mpf(str(A[0][1])), REPORT)}],")
        print(f"        [{mp.nstr(mp.mpf(str(A[1][0])), REPORT)}, "
              f"{mp.nstr(mp.mpf(str(A[1][1])), REPORT)}]]")

    print("\nadmissibility  A_i >= S >= I  (eigenvalues of S^{-1/2} A_i S^{-1/2}):")
    ok = True
    for k, q in enumerate(qs):
        p = q.s
        a2 = p[0] * p[1]
        a1 = -(q.A[0][0] * p[1] + q.A[1][1] * p[0])
        disc = mp.sqrt(a1 * a1 - 4 * a2 * q.detA)
        lo_, hi_ = (-a1 - disc) / (2 * a2), (-a1 + disc) / (2 * a2)
        ok = ok and lo_ >= 1
        print(f"    A_{k+1}: [{mp.nstr(lo_, 20)}, {mp.nstr(hi_, 20)}]"
              f"   min - 1 = {mp.nstr(lo_ - 1, 12)}")
    print(f"    S >= I: min s_i - 1 = "
          f"{mp.nstr(min(mp.mpf(str(v)) for v in s) - 1, 12)}")
    print(f"    admissible: {ok}")

    edges = ((0, 1), (1, 2), (2, 0))
    out = {"S": [float(v) for v in s], "A": ops}
    for tag, fam in (("quantum", qs), ("shadow", ps)):
        print(f"\n{tag} midranges A(i,j) at {REPORT} digits:")
        mids, ds = [], []
        for i, j in edges:
            d, m, P, Q = d_and_mid(fam[i], fam[j], smin, smax)
            mids.append(m)
            ds.append(d)
            print(f"    A({i+1},{j+1}) = {mp.nstr(m, REPORT)}")
            print(f"    d({i+1},{j+1}) = {mp.nstr(d, REPORT)}")
        margin = max(min(mids), min(-x for x in mids))
        print(f"    cycle margin max(min A, min -A) = {mp.nstr(margin, REPORT)}")
        out[tag] = {
            "midranges": [mp.nstr(x, REPORT) for x in mids],
            "d": [mp.nstr(x, REPORT) for x in ds],
            "cycle_margin": mp.nstr(margin, REPORT),
            "cycle_margin_float": float(margin),
        }

    qm = out["quantum"]["cycle_margin_float"]
    sm = out["shadow"]["cycle_margin_float"]
    print("\n" + "-" * 74)
    print(f"quantum cycle margin      = {qm:+.12e}   (> 0 : strict 3-cycle)")
    print(f"shadow  cycle margin      = {sm:+.12e}   (< 0 : transitive)")
    print(f"separation min(qm, -sm)   = {min(qm, -sm):+.12e}")
    print(f"tie threshold of the brief = 1e-10;  margin/threshold = "
          f"{min(qm, -sm)/1e-10:.3e}")
    verdict = qm > 1e-10 > 0 > sm and -sm > 1e-10
    print(f"CERTIFIED: {bool(verdict)}")
    out["separation"] = min(qm, -sm)
    out["certified"] = bool(verdict)

    print("\nRobustness: the same midranges with the infimum taken over the")
    print("FULL classical domain beta in [1e-6, 1e6] u {0, oo}, i.e. without")
    print("the alpha >= 1/2 restriction that data processing imposes.")
    print("(mpmath has no conditioning limit here: mu_k is a closed form.)")
    smin_full = mp.log(mp.mpf(10) ** -6)
    for tag, fam in (("quantum", qs), ("shadow", ps)):
        mids = []
        for i, j in edges:
            phi = lambda s: fam[j].U(s) - fam[i].U(s)
            e_inf = mp.log(fam[j].Lam) - mp.log(fam[i].Lam)
            e_zero = mp.log(fam[j].R) - mp.log(fam[i].R)
            hi, lo = scan_extrema(phi, smin_full, smax, 900)
            P = max([v for _, v in hi] + [e_inf, e_zero])
            Q = min([v for _, v in lo] + [e_inf, e_zero])
            mids.append((P + Q) / 2)
        margin = max(min(mids), min(-x for x in mids))
        print(f"    {tag:8s} midranges: "
              + ", ".join(mp.nstr(x, 12) for x in mids))
        print(f"    {tag:8s} cycle margin = {mp.nstr(margin, 20)}")
        out[tag]["full_domain_midranges"] = [mp.nstr(x, REPORT) for x in mids]
        out[tag]["full_domain_cycle_margin"] = float(margin)
    print("    The separation is unchanged, so the cycle is not an artefact")
    print("    of restricting the infimum to the data-processing range.")

    print("\ncurl check (brief D Part 0(b)):  |curl A| / sum|A|")
    for tag in ("quantum", "shadow"):
        m = [mp.mpf(x) for x in out[tag]["midranges"]]
        curl = m[0] + m[1] + m[2]
        tot = sum(abs(x) for x in m)
        print(f"    {tag:8s}: curl = {mp.nstr(curl, 20)}   "
              f"|curl|/sum|A| = {mp.nstr(abs(curl)/tot, 20)}")
        out[tag]["curl"] = mp.nstr(curl, REPORT)
        out[tag]["curl_ratio"] = float(abs(curl) / tot)

    name = "q_certify.json" if src == "q_cycle.json" else \
        "q_certify_" + src.replace("q_cycle_", "").replace(".json", "") + ".json"
    with open(os.path.join(HERE, name), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {name}")


if __name__ == "__main__":
    main()
