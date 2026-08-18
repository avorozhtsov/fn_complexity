"""Session brief I -- is distortion exactly 1 attained by a signature family?

The cone witness (i_certify.py) realises d = s.C_4 exactly, and the signature
families it induces have distortion 1 + O(1/log r) (i_witness.py).  Whether an
actual signature family attains 1 is a different question: the cone point is a
tropical (piecewise-linear) boundary point of the achievable set.

This script tests it directly.  Fix the multiplicities at the rung-B values and
solve the five scale-free equations

    d_ij / (C_4)_ij  =  const                      (5 equations)

by damped Gauss-Newton, first over the 20 real log-atoms and then over all 40
parameters (log-atoms and log-multiplicities together).  d is evaluated by a
grid plus a golden-section refinement inside the bracketing cell, which is
rigorous because |(U_b - U_a)'| <= 1.

Result (recorded in i_solve_sig_output.txt): the residual does NOT go to zero.
It converges to a stationary point at ~1.8e-3 (B = 100) and ~6.2e-4 (B = 300),
improving the distortion from 1.0547 to 1.01586 and from 1.0148 to 1.00594.  So
the five-line ansatz is too rigid; more distinct atoms are the obvious next
move.

    python research/realizability/i_solve_sig.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import i_logsig as L  # noqa: E402
import i_witness as W  # noqa: E402

C4 = W.C4
IU = W.IU
GR = (math.sqrt(5) - 1) / 2


def sharp(a, b, lo=-25.0, hi=25.0, step=2e-3):
    g = np.arange(lo, hi + 0.5 * step, step)
    phi = b.U(g) - a.U(g)
    e0 = math.log(b.R) - math.log(a.R)
    e1 = math.log(b.Lam) - math.log(a.Lam)

    def f(t):
        return float(b.U(np.array([t]))[0] - a.U(np.array([t]))[0])

    def golden(l, r, sign):
        x1, x2 = r - GR * (r - l), l + GR * (r - l)
        y1, y2 = sign * f(x1), sign * f(x2)
        for _ in range(120):
            if r - l < 1e-14:
                break
            if y1 <= y2:
                r, x2, y2 = x2, x1, y1
                x1 = r - GR * (r - l)
                y1 = sign * f(x1)
            else:
                l, x1, y1 = x1, x2, y2
                x2 = l + GR * (r - l)
                y2 = sign * f(x2)
        return sign * min(y1, y2)

    i = int(np.argmax(phi))
    j = int(np.argmin(phi))
    hv = max(e0, e1, golden(g[max(i - 1, 0)], g[min(i + 1, len(g) - 1)], -1.0))
    lv = min(e0, e1, golden(g[max(j - 1, 0)], g[min(j + 1, len(g) - 1)], +1.0))
    return hv - lv


def dm(sigs):
    D = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            D[i, j] = D[j, i] = sharp(sigs[i], sigs[j])
    return D


def make_all(p):
    p = p.reshape(4, 10)
    return [L.LogSig(np.abs(p[a, :5]), np.abs(p[a, 5:])) for a in range(4)]


def make_x(p, cfix):
    xs = p.reshape(4, 5)
    return [L.LogSig(cfix[a], np.abs(xs[a])) for a in range(4)]


def gn(p, mk, iters=60, label=""):
    def res(q):
        D = dm(mk(q))
        r = D[IU] / C4[IU]
        return r[1:] - r[0]

    def dist(q):
        D = dm(mk(q))
        r = D[IU] / C4[IU]
        return float(r.max() / r.min())

    f = res(p)
    print(f"   {label} start  max|res| = {np.max(np.abs(f)):.3e}  "
          f"dist = {dist(p):.12f}", flush=True)
    for it in range(iters):
        J = np.empty((5, len(p)))
        for i in range(len(p)):
            h = 1e-7 * max(1.0, abs(p[i]))
            q = p.copy()
            q[i] += h
            J[:, i] = (res(q) - f) / h
        dx, *_ = np.linalg.lstsq(J, -f, rcond=1e-10)
        ls, best = 1.0, None
        for _ in range(45):
            q = p + ls * dx
            fq = res(q)
            if np.max(np.abs(fq)) < np.max(np.abs(f)):
                best = (q, fq)
                break
            ls *= 0.5
        if best is None:
            break
        p, f = best
        if np.max(np.abs(f)) < 1e-13:
            break
    print(f"   {label} final  max|res| = {np.max(np.abs(f)):.3e}  "
          f"dist = {dist(p):.15f}", flush=True)
    return p, f


def main():
    cs, xs, s = W.cone_lines()
    print("=== Gauss-Newton on the five scale-free equations ===")
    for B in (100.0, 300.0):
        lam = B / cs.max()
        print(f"  --- rung B = {B:.0f} ---")
        cfix = lam * cs
        p = (lam * xs).reshape(-1).copy()
        gn(p, lambda q: make_x(q, cfix), label="atoms only ")
        p2 = np.hstack([np.hstack([lam * cs[a], lam * xs[a]])
                        for a in range(4)])
        gn(p2, make_all, label="all 40 par ")


if __name__ == "__main__":
    main()
