#!/usr/bin/env python3
"""Cycles are a zero-temperature phenomenon in the Maslov parameter.

A(a,b) = mid_s(u_a - u_b) = (max + min)/2 is the TROPICAL (lambda -> infinity)
member of the family

    softmax_L(f) = (1/L) log INT e^{L f} rho ,   softmin_L(f) = -(1/L) log INT e^{-L f} rho
    A_L(a,b)     = (softmax_L + softmin_L)/2  of  f = u_a - u_b

As L -> infinity, A_L -> mid.  As L -> 0, BOTH soft-extrema tend to INT f rho, so

    A_0(a,b) = INT (u_a - u_b) rho = Psi(a) - Psi(b),   Psi(a) = INT u_a rho,

an exact potential difference -- hence a total order with no cycles, for any rho.
So every cycle in this framework lives above a critical L, and averaging over
temperature restores a scalar complexity.
"""
import numpy as np

CYCLE = [(6, 3, 3), (7, 2, 1), (6, 5, 1)]


def u_of(sig, s):
    """u_a(s) = log log Z_a(e^s), computed stably."""
    beta = np.exp(s)
    logs = np.log(np.array(sig, float))
    z = np.outer(beta, logs)
    m = z.max(axis=1)
    return np.log(m + np.log(np.exp(z - m[:, None]).sum(axis=1)))


def soft_mid(f, rho, lam):
    """(softmax + softmin)/2; lam=None gives the tropical (max+min)/2."""
    if lam is None:
        return 0.5 * (f.max() + f.min())
    mx, mn = f.max(), f.min()
    smax = (np.log((rho * np.exp(lam * (f - mx))).sum()) + lam * mx) / lam
    smin = -(np.log((rho * np.exp(-lam * (f - mn))).sum()) - lam * mn) / lam
    return 0.5 * (smax + smin)


def critical_lambda(smin=-12.0, smax=12.0, n=24001, cycle=CYCLE):
    s = np.linspace(smin, smax, n)
    rho = np.ones_like(s) / n
    U = {a: u_of(a, s) for a in cycle}
    A = lambda a, b, L: soft_mid(U[a] - U[b], rho, L)
    is_cycle = lambda L: all(A(cycle[i], cycle[(i + 1) % 3], L) > 0 for i in range(3))
    if not is_cycle(1e6):
        return None
    lo, hi = 1.0, 1e6
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        if is_cycle(mid):
            hi = mid
        else:
            lo = mid
    return hi


if __name__ == "__main__":
    s = np.linspace(-12, 12, 24001)
    rho = np.ones_like(s) / len(s)
    U = {a: u_of(a, s) for a in CYCLE}
    A = lambda a, b, L: soft_mid(U[a] - U[b], rho, L)

    print("tropical (lambda = infinity):")
    v = [A(CYCLE[i], CYCLE[(i + 1) % 3], None) for i in range(3)]
    print("   A =", [f"{x:+.6f}" for x in v], " cycle:", all(x > 0 for x in v))

    print("\nlambda -> 0: A_0 is a potential difference, Psi(a) = INT u_a rho")
    print("  ", {a: round(float((rho * U[a]).sum()), 6) for a in CYCLE})

    print("\n  lambda    cycle?   A(1,2)      A(2,3)      A(3,1)")
    for L in (1, 10, 100, 300, 1000, 10000):
        v = [A(CYCLE[i], CYCLE[(i + 1) % 3], L) for i in range(3)]
        print(f"  {L:<9} {str(all(x>0 for x in v)):6}  " + "  ".join(f"{x:+.6f}" for x in v))

    print("\ncritical lambda vs the prior's support (uniform in s = log beta):")
    for w in [(-8, 8), (-12, 12), (-16, 16), (-6, 20), (-20, 6)]:
        print(f"   s in {str(w):10s} -> lambda_c = {critical_lambda(*w):,.1f}")
