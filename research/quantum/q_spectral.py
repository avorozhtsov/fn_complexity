"""Brief L, Part 1 -- the naive quantisation is empty.

    Z_A(beta) = Tr A^beta      C(A -> B) = inf_beta log Tr A^b / log Tr B^b

Tr A^b = sum_k lambda_k(A)^b depends on A only through its spectrum, so

    C(A -> B) = C(spec A -> spec B)

with the right-hand side the *classical* exchange rate of the two positive
r-tuples.  This script is a smoke test of that one-line proof and of its two
corollaries:

  (a) conjugation invariance -- C(UAU*, B) = C(A, B) for every unitary U,
      so no arrangement of non-commuting operators can be told from a
      commuting one;
  (b) no quantum 3-cycle without a classical shadow -- the tournament, the
      metric d, the midrange A and every curl of a family of positive
      operators equal those of the family of their spectra, exactly.

Run:  python3 q_spectral.py
"""
from __future__ import annotations

import math

import numpy as np

from q_core import osc_mid_fast, rand_spd, spectral_QSig

GRID = np.linspace(math.log(1e-3), math.log(1e6), 3001)


def classical_from_spectrum(A):
    """The commuting operator with the same spectrum: diag(eigenvalues)."""
    return np.diag(np.linalg.eigvalsh(A))


def main():
    rng = np.random.default_rng(20260818)

    print("=" * 74)
    print("PART 1  --  Tr A^beta sees only the spectrum")
    print("=" * 74)
    print(f"beta grid: 1e-3 .. 1e6, {len(GRID)} log-spaced points, "
          "plus the analytic endpoints beta = 0, oo")

    print("\n(a) C(A->B) vs C(diag spec A -> diag spec B), 300 random pairs")
    worst_d = worst_mid = worst_comm = 0.0
    for _ in range(300):
        r = int(rng.integers(2, 8))
        A = rand_spd(rng, r, 1.05, 40.0)
        B = rand_spd(rng, r, 1.05, 40.0)
        qa, qb = spectral_QSig(A), spectral_QSig(B)
        ca = spectral_QSig(classical_from_spectrum(A))
        cb = spectral_QSig(classical_from_spectrum(B))
        dq, mq = osc_mid_fast(qa, qb, GRID, qa.U_grid(GRID), qb.U_grid(GRID))
        dc, mc = osc_mid_fast(ca, cb, GRID, ca.U_grid(GRID), cb.U_grid(GRID))
        worst_d = max(worst_d, abs(dq - dc))
        worst_mid = max(worst_mid, abs(mq - mc))
        worst_comm = max(worst_comm, float(np.linalg.norm(A @ B - B @ A)))
    print(f"    max |d_quantum - d_classical|   = {worst_d:.3e}")
    print(f"    max |A_quantum - A_classical|   = {worst_mid:.3e}")
    print(f"    largest commutator norm sampled = {worst_comm:.4f}")
    print("    (zero up to LAPACK eigenvalue error -- the two computations")
    print("     are the same sum written twice)")

    print("\n(b) conjugation invariance: d, A for U A U^T against A")
    A = rand_spd(rng, 5, 1.05, 40.0)
    B = rand_spd(rng, 5, 1.05, 40.0)
    qa, qb = spectral_QSig(A), spectral_QSig(B)
    base = osc_mid_fast(qa, qb, GRID, qa.U_grid(GRID), qb.U_grid(GRID))
    worst = 0.0
    for _ in range(50):
        U, _ = np.linalg.qr(rng.normal(size=(5, 5)))
        qu = spectral_QSig(U @ A @ U.T)
        got = osc_mid_fast(qu, qb, GRID, qu.U_grid(GRID), qb.U_grid(GRID))
        worst = max(worst, abs(got[0] - base[0]), abs(got[1] - base[1]))
    print(f"    max deviation over 50 random unitaries = {worst:.3e}")

    print("\n(c) no quantum 3-cycle without a classical shadow")
    print("    2000 random non-commuting triples; compare the tournament of")
    print("    the operators with the tournament of their spectra.")
    diffs = cyc_q = cyc_c = 0
    for _ in range(2000):
        r = int(rng.integers(2, 6))
        ops = [rand_spd(rng, r, 1.05, 30.0) for _ in range(3)]
        qs = [spectral_QSig(M) for M in ops]
        cs = [spectral_QSig(classical_from_spectrum(M)) for M in ops]
        pq = [q.U_grid(GRID) for q in qs]
        pc = [q.U_grid(GRID) for q in cs]
        sq = [osc_mid_fast(qs[i], qs[j], GRID, pq[i], pq[j], refine=False)[1]
              for i, j in ((0, 1), (1, 2), (2, 0))]
        sc = [osc_mid_fast(cs[i], cs[j], GRID, pc[i], pc[j], refine=False)[1]
              for i, j in ((0, 1), (1, 2), (2, 0))]
        if any(np.sign(x) != np.sign(y) for x, y in zip(sq, sc)):
            diffs += 1
        cyc_q += int(all(x > 1e-10 for x in sq) or all(x < -1e-10 for x in sq))
        cyc_c += int(all(x > 1e-10 for x in sc) or all(x < -1e-10 for x in sc))
    print("    triples scanned                       : 2000")
    print(f"    tournaments differing from the shadow : {diffs}")
    print(f"    3-cycles quantum / classical          : {cyc_q} / {cyc_c}")

    print("\nVERDICT: the extension of C to positive operators factors through")
    print("the pair of spectra.  It is the classical exchange rate composed")
    print("with A |-> spec A, which is onto the positive r-tuples, so it is")
    print("neither a restriction nor an enlargement.  There is no quantum")
    print("cycle without a classical shadow, and there cannot be one.")


if __name__ == "__main__":
    main()
