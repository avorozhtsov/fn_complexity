#!/usr/bin/env python3
"""Szegedy quantisation of the exchange market, and where it fails.

Szegedy's construction turns any stochastic matrix ``P`` on ``n`` states into a
unitary ``W(P)`` on ``C^n (x) C^n``.  Its spectrum is controlled by the
DISCRIMINANT

    D(P)_{ab} = sqrt(P_ab * P_ba),

the entrywise geometric mean of the forward and backward transition rates; the
eigenvalues of ``W`` are ``exp(+- i arccos(sigma))`` for ``sigma`` a singular
value of ``D``, plus a trivial part.  Reversibility is NOT required, which is why
this is the right quantisation for a market that does not obey detailed balance.

TWO CHOICES OF P, AND WHY.

(1) Direct.  Set ``P_ab prop. C(a->b)^theta`` off the diagonal.  This is
    self-contained but arbitrary: the exponent ``theta`` and the diagonal mass
    are free, and the discriminant depends on both.

(2) Metropolis on the exchange metric.  Propose uniformly and accept with
    ``min(1, exp(-t (d(a,b) - 0)))``.  This is the standard reversible chain, so
    ``D(P) = P`` symmetrised and the walk carries no information about the
    ANTISYMMETRIC part A -- which is precisely the part that carries the
    comparison and the cycles.  It quantises the wrong half.

The point of the file is that (2) is what the textbook route (Rokhsar-Kivelson,
stochastic matrix form) forces on you, and it throws away A.  (1) keeps A, and
what it exposes is that the discriminant is

    D(P)_{ab} prop. sqrt(C(a->b) C(b->a))^theta = exp(-theta d(a,b)/2),

i.e. the ``t = theta/2`` member of the Schoenberg family -- INDEPENDENT of the
antisymmetric part after all, because the geometric mean of forward and backward
rates is exactly what A cancels out of.  So:

    * the discriminant of the exchange market is the Gibbs kernel at t = theta/2;
    * at the natural normalisation theta = 1 that is t = 1/2, and whether the
      discriminant is a Gram matrix is decided by t* vs 1/2 (T1.2);
    * A survives only in the walk's non-Hermitian part, never in its spectrum.

The last line is the result.  **Szegedy quantisation of the exchange matrix is
blind to the comparison**, for the same structural reason Rokhsar-Kivelson is:
both are built from sqrt(P_ab P_ba), and A is the odd part of log P.  There is no
choice of theta or of diagonal that repairs this.

    python research/m_and_e_and_a_c/szegedy_walk.py
"""
from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))

import t1_2_common as T  # noqa: E402
from fn_complexity import exchange_rate  # noqa: E402
from gauge_decomposition import CERT5, CYCLE  # noqa: E402


def rate_matrix(family) -> np.ndarray:
    n = len(family)
    m = np.eye(n)
    for i, j in itertools.permutations(range(n), 2):
        m[i, j] = exchange_rate(family[i], family[j])
    return m


def stochastic(m: np.ndarray, theta: float) -> np.ndarray:
    """The LAZY chain ``P_ab prop. C(a->b)^theta``, diagonal included.

    Keeping ``a = b`` matters: ``C(a->a) = 1``, so the diagonal of the
    discriminant comes out as ``1/R_a``, which is exactly what makes the identity
    of `main` hold on the diagonal as well as off it.
    """
    w = m ** theta
    return w / w.sum(axis=1, keepdims=True)


def stationary(p: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eig(p.T)
    pi = np.real(vectors[:, int(np.argmin(np.abs(values - 1.0)))])
    return pi / pi.sum()


def time_reversal(p: np.ndarray) -> np.ndarray:
    """``P*_ba = pi_a P_ab / pi_b``; equals ``P`` iff the chain is reversible."""
    pi = stationary(p)
    return (p * pi[:, None]).T / pi[:, None]


def discriminant(p: np.ndarray) -> np.ndarray:
    return np.sqrt(p * p.T)


def szegedy_unitary(p: np.ndarray) -> np.ndarray:
    """W = R2 R1 with R_k = 2 Pi_k - I, the standard two-reflection walk."""
    n = p.shape[0]
    psi = np.zeros((n, n * n))
    for a in range(n):
        psi[a] = np.kron(np.eye(n)[a], np.sqrt(p[a]))
    proj1 = psi.T @ psi
    star = time_reversal(p)
    phi = np.zeros((n, n * n))
    for b in range(n):
        phi[b] = np.kron(np.sqrt(star[b]), np.eye(n)[b])
    proj2 = phi.T @ phi
    eye = np.eye(n * n)
    return (2 * proj2 - eye) @ (2 * proj1 - eye)


def main() -> None:
    print("Szegedy discriminant of the exchange market")
    print("claim:  D(P) = Delta * exp(-theta d / 2) * Delta   exactly, "
          "Delta = diag(R_a^-1/2)\n")
    for name, family in (("known 3-cycle", CYCLE), ("minimal certificate", CERT5)):
        print(f"=== {name} ({len(family)} signatures) ===")
        m = rate_matrix(family)
        d = T.distance_matrix(family)
        for theta in (0.5, 1.0, 4.0):
            p = stochastic(m, theta)
            disc = discriminant(p)
            delta = 1.0 / np.sqrt((m ** theta).sum(axis=1))
            model = np.outer(delta, delta) * np.exp(-0.5 * theta * d)
            rel = np.abs(disc - model).max() / np.abs(disc).max()
            ev_disc = np.linalg.eigvalsh(disc)
            ev_gibbs = np.linalg.eigvalsh(np.exp(-0.5 * theta * d))
            same = np.array_equal(np.sign(np.where(np.abs(ev_disc) < 1e-12, 0, ev_disc)),
                                  np.sign(np.where(np.abs(ev_gibbs) < 1e-12, 0, ev_gibbs)))
            print(f"  theta = {theta:>4} (t = {theta / 2})   "
                  f"max rel dev = {rel:.3e}   "
                  f"lambda_min(D) = {ev_disc[0]: .4e}   "
                  f"inertia matches e^(-td): {same}")

        star = T.psd_threshold(d)
        print(f"  t* = {'-' if star is None else f'{star:.6f}'}"
              f"   -> at theta = 1 the discriminant is "
              f"{'NOT ' if star is not None and star > 0.5 else ''}"
              f"positive semidefinite")

        p = stochastic(m, 1.0)
        w = szegedy_unitary(p)
        print(f"  W unitary: ||W W^* - I|| = "
              f"{np.abs(w @ w.conj().T - np.eye(w.shape[0])).max():.2e}"
              f"    chain reversible: "
              f"{np.allclose(time_reversal(p), p, atol=1e-12)}")

        # A enters only through the row sums, i.e. only as one number per node.
        p_rev = stochastic(m.T, 1.0)
        ev_a = np.linalg.eigvalsh(discriminant(p))
        ev_b = np.linalg.eigvalsh(discriminant(p_rev))
        print(f"  flipping A -> -A moves the discriminant spectrum by "
              f"{np.abs(ev_a - ev_b).max():.3e}, and its inertia by "
              f"{int(np.sum(np.sign(ev_a) != np.sign(ev_b)))} signs\n")

    print("verdict.  Off the diagonal sqrt(P_ab P_ba) = Delta_a Delta_b e^{-theta d/2},")
    print("because the geometric mean of the forward and backward rates is exactly what")
    print("A cancels out of; the lazy diagonal makes it hold everywhere.  So D(P) is a")
    print("POSITIVE DIAGONAL CONGRUENCE of the Gibbs kernel, and by Sylvester's law of")
    print("inertia its signature -- in particular whether it is positive semidefinite --")
    print("is the Schoenberg question at t = theta/2 and does not involve A at all.")
    print("A reaches the walk only through the row sums R_a: one number per node, i.e.")
    print("exactly the kind of data a potential is, and exactly the kind that cannot")
    print("produce a cycle.  Szegedy and Rokhsar-Kivelson quantise the same half of the")
    print("exchange matrix, and it is the half that carries no arithmetic.")


if __name__ == "__main__":
    main()
