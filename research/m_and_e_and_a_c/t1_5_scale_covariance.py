#!/usr/bin/env python3
"""Is dilation covariance of the Weil pairing equivalent to RH?  No -- it is the
functional equation, and holds unconditionally.

FINDINGS T1.5 and section 3b(ii) of `exchange_positivity_and_weil.md` state
Theorem A as `E_{lambda a, lambda b} = lambda^{2 Re rho} E_ab`, flagged "this
step uses RH", with a counterfactual that moves 60 of 1200 zeros to Re = 0.7 and
observes the invariance break from 6e-14 to 1.3e-1.  That experiment measures the
wrong matrix.  There are two matrices here and they differ exactly when RH fails:

    W_ab = sum_rho Z_a(rho) Z_b(1 - rho)          the Weil pairing
    G_ab = sum_rho Z_a(rho) conj(Z_b(rho))        the Gram surrogate

`riemann_hypothesis_exchange_matrices.md` defines the pairing as W; G is what W
BECOMES under RH, since rho on the critical line gives 1 - rho = conj(rho).  The
research thread computed G throughout, which is harmless while the zeros used are
genuinely on the line, and wrong the moment a zero is moved off it.

The scaling of each is one line, since Z_{lambda a}(s) = lambda^s Z_a(s):

    W_{lambda a, lambda b} = sum_rho lambda^{rho} lambda^{1 - rho} Z_a Z_b
                           = lambda * W_ab                     ALWAYS
    G_{lambda a, lambda b} = sum_rho lambda^{2 Re rho} Z_a conj(Z_b)

W is homogeneous of degree exactly 1 for ANY multiset of zeros, because the
pairing puts rho against 1 - rho and rho + (1 - rho) = 1.  That is the functional
equation, not the Riemann hypothesis.  G is homogeneous only when all the Re rho
agree.  So:

    scale covariance  <=>  the functional equation      (unconditional)
    positivity        <=>  the Riemann hypothesis       (Weil's criterion)

and the two must not be conflated.  Theorem A stands as stated in
`T1_5_multiplicative_design.md`, where it is derived from the Landau form and is
correctly unconditional; what must go is the "this step uses RH" callout and the
counterfactual attached to it.

This script verifies all of it.  Zero multisets are built closed under both
involutions rho -> 1 - rho and rho -> conj(rho), as Weil's zero set is; the
earlier counterfactual moved zeros singly, which breaks that closure and is the
second reason it cannot be read as a test of anything.

    python research/m_and_e_and_a_c/t1_5_scale_covariance.py
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

FAMILY = [
    (2, 2), (3, 1), (4, 2), (5, 3), (6, 1), (3, 1, 1), (8, 4), (9, 3),
    (4, 4, 1), (7, 5), (5, 5), (7, 1), (9, 1), (6, 3), (10, 5), (11, 7),
]


def z(signature, s: np.ndarray) -> np.ndarray:
    """Z_a(s) = sum_i a_i^s, vectorised over s."""
    return np.sum(np.array(signature, dtype=float)[:, None] ** s[None, :], axis=0)


def pairings(family, rhos: np.ndarray, scale: float = 1.0):
    """Return (W, G) for the family scaled by `scale`."""
    scaled = [tuple(scale * x for x in a) for a in family]
    za = np.array([z(a, rhos) for a in scaled])
    zb = np.array([z(a, 1.0 - rhos) for a in scaled])
    w = (za @ zb.T)
    g = (za @ za.conj().T)
    return w, g


def closed_zero_set(gammas: np.ndarray, delta: float) -> np.ndarray:
    """Zeros at Re = 1/2 +- delta, closed under rho -> 1-rho and rho -> conj(rho)."""
    if delta == 0.0:
        return np.concatenate([0.5 + 1j * gammas, 0.5 - 1j * gammas])
    offsets = np.array([delta, -delta])
    return np.concatenate([
        (0.5 + o) + 1j * s * gammas for o in offsets for s in (1.0, -1.0)
    ])


def report(name: str, rhos: np.ndarray, lambdas=(2.0, 7.0, 210.0)) -> None:
    w1, g1 = pairings(FAMILY, rhos)
    print(f"\n{name}   ({len(rhos)} zeros)")
    print(f"  ||W - G||_max                     = {np.abs(w1 - g1).max():.3e}")
    print(f"  max |Im W|                        = {np.abs(w1.imag).max():.3e}"
          "   (W is real symmetric)")
    print(f"  ||W - W^T||_max                   = {np.abs(w1 - w1.T).max():.3e}")
    for lam in lambdas:
        wl, gl = pairings(FAMILY, rhos, scale=lam)
        rel_w = np.abs(wl - lam * w1).max() / np.abs(w1).max()
        rel_g = np.abs(gl - lam * g1).max() / np.abs(g1).max()
        print(f"  lambda = {lam:6.1f}   W: rel dev = {rel_w:.3e}"
              f"      G: rel dev = {rel_g:.3e}")
    print(f"  lambda_min(W) = {np.linalg.eigvalsh(w1.real):.6e}"
          if w1.shape[0] == 1 else
          f"  lambda_min(W) = {np.linalg.eigvalsh(w1.real)[0]: .6e}"
          f"    lambda_min(G) = {np.linalg.eigvalsh(g1)[0]: .6e}")


def positivity_threshold(family, gammas: np.ndarray, deltas) -> None:
    """Does the truncated Weil form detect off-line zeros at all?"""
    print(f"\n  family of {len(family)} signatures")
    for delta in deltas:
        rhos = closed_zero_set(gammas, delta)
        za = np.array([z(a, rhos) for a in family])
        zb = np.array([z(a, 1.0 - rhos) for a in family])
        w = (za @ zb.T).real
        eigenvalues = np.linalg.eigvalsh(0.5 * (w + w.T))
        scale = np.abs(w).max()
        print(f"    Re = 1/2 +- {delta:<5} lambda_min/scale = {eigenvalues[0] / scale: .3e}"
              f"   {'NEGATIVE' if eigenvalues[0] < -1e-9 * scale else 'positive semidefinite'}")


def random_family(count: int = 40, seed: int = 3):
    rng = random.Random(seed)
    return sorted({
        tuple(sorted((rng.randint(1, 40) for _ in range(rng.randint(2, 6))), reverse=True))
        for _ in range(count)
    })


def main() -> None:
    gammas = np.load(HERE / "zeta_zeros_1200.npy")[:300]
    print("W_ab = sum_rho Z_a(rho) Z_b(1-rho)   vs   G_ab = sum_rho Z_a(rho) conj(Z_b(rho))")
    print("scaling claim:  W is homogeneous of degree 1 for ANY zero multiset;")
    print("                G is homogeneous only if every Re rho agrees.")

    report("all zeros on the critical line (RH holds)", closed_zero_set(gammas, 0.0))
    for delta in (0.05, 0.2):
        report(f"quadruples off the line, Re = 1/2 +- {delta}",
               closed_zero_set(gammas, delta))

    print("\n\nWhat IS sensitive to the critical line is positivity of W, not its")
    print("covariance -- that is Weil's criterion.  How sensitive, on atomic measures:")
    deltas = (0.0, 0.1, 0.2, 0.3, 0.45, 0.49)
    positivity_threshold(FAMILY, gammas, deltas)
    positivity_threshold(random_family(), gammas, deltas)
    print("\n  So the truncated form has a DETECTION THRESHOLD.  The 16-signature family")
    print("  of T1.4 stays positive semidefinite even at Re = 0.99; the 40-signature one")
    print("  first goes negative between Re = 0.7 and Re = 0.8, i.e. the Re = 0.7 of the")
    print("  FINDINGS T1.5 counterfactual sits just below detection, and the quantity that")
    print("  broke there (covariance of G) was not the one that carries the information.")
    print("  Atomic measures are not admissible Weil test functions, so this is a")
    print("  property of the truncation and not evidence about zeta either way.")

    print("\nverdict: W is covariant to machine precision in every case, including the")
    print("off-line ones; G is not.  Dilation covariance is the functional equation.")
    print("It is not equivalent to RH, so it cannot serve as a criterion for it.")


if __name__ == "__main__":
    main()
