#!/usr/bin/env python3
"""Does the within-class potential of E2.c predict the ``f_1`` vs ``f_2`` edge?

The ``F_101`` witness of brief B's addendum 2 has two pencils with the same
largest fiber, the same ``m_2`` and the same ``nu(P)`` -- so ``phi`` and
``phi~`` are exactly tied and the comparison is decided by the interior.  E2.c
says the within-class potential is the moment ladder ``m2, m3, m4, log mu``.
This script fits that potential on a ``q = 101`` pool, restricted to each class
of equal largest fiber, and then asks whether it gets the ``f_1``/``f_2`` edge
right -- a genuine out-of-sample prediction, since the fit never sees that pair.

Writes ``e2_f101_class.csv``.
"""

from __future__ import annotations

import csv

import numpy as np

import curves
import flux
import pools
from e2_potential import COARSE, FINE, r2, statistics_of

Q = 101
N = 700
PENCILS = {
    "f_1": {5: 1, 4: 70, 3: 28, 2: 15, 1: 11},
    "f_2": {5: 1, 4: 42, 3: 32, 2: 74, 1: 96},
}
MODELS = [
    ("m2", ["m2"]),
    ("log mult", ["log_mult"]),
    ("m3", ["m3"]),
    ("m2, m3", ["m2", "m3"]),
    ("m2, m3, m4", ["m2", "m3", "m4"]),
    ("m2, m3, m4, log mult", ["m2", "m3", "m4", "log_mult"]),
]


def main() -> None:
    full = pools.sampled_pool(Q, draws=4000)
    rng = np.random.default_rng(101)
    S = full[np.sort(rng.choice(len(full), min(N, len(full)), replace=False))]
    witness = {k: np.array(sorted((f.n_affine for f in curves.pencil_fibers(v, Q, 5)),
                                  reverse=True))
               for k, v in PENCILS.items()}
    # make sure the two witnesses are in the pool exactly once each
    keep = [row for row in S if not any((row == w).all() for w in witness.values())]
    S = np.array(keep + [witness["f_1"], witness["f_2"]], dtype=np.int64)
    i1, i2 = len(S) - 2, len(S) - 1

    A = flux.flux_matrix(S, flux.beta_grid(Q, COARSE, FINE))
    st = statistics_of(S, Q)
    maxN = st["maxN"]
    print(f"pool n = {len(S)}, {len(np.unique(maxN))} classes of equal largest fiber")
    print(f"whole-pool curl/|A| = {flux.hodge(A)['curl_frac']:.6f}")
    print(f"the witness edge: A(f_1, f_2) = {A[i1, i2]:+.6e}  "
          f"(negative: f_2 precedes f_1)")

    rows = []
    fits = {}
    num = den = 0.0
    curl_num = curl_den = 0.0
    used = 0
    agg = {name: [0.0, 0.0] for name, _ in MODELS}
    for value in np.unique(maxN):
        idx = np.flatnonzero(maxN == value)
        if len(idx) < 3:
            continue
        used += 1
        B = A[np.ix_(idx, idx)]
        h = flux.hodge(B)
        psi = h["psi"]
        curl_num += np.linalg.norm(B - (psi[None, :] - psi[:, None])) ** 2
        curl_den += np.linalg.norm(B) ** 2
        sub = {k: v[idx] for k, v in st.items()}
        for name, cols in MODELS:
            X = np.column_stack([np.ones(len(idx))] + [sub[c] for c in cols])
            beta, *_ = np.linalg.lstsq(X, psi, rcond=None)
            res = psi - X @ beta
            agg[name][0] += float((res**2).sum())
            agg[name][1] += float(((psi - psi.mean()) ** 2).sum())
            if value == maxN[i1]:
                fits[name] = (beta, cols, idx)
        x = sub["m2"] - sub["m2"].mean()
        num += float(x @ (psi - psi.mean()))
        den += float(x @ x)

    print(f"{used} classes with 3 or more members; "
          f"within-class curl/|A| = {np.sqrt(curl_num/curl_den):.6f}")
    print(f"pooled within-class slope of psi on m2: {num/den:+.8f}")
    print("\nwithin-class regressions")
    for name, _ in MODELS:
        value = 1 - agg[name][0] / agg[name][1]
        rows.append({"q": Q, "model": name, "R2": value})
        print(f"  R2 = {value:.6f}   {name}")

    print(f"\nout-of-sample check on the witness pair "
          f"(class max fiber = {int(maxN[i1])}, "
          f"{int((maxN == maxN[i1]).sum())} members)")
    truth = np.sign(A[i1, i2])
    for name, (beta, cols, idx) in fits.items():
        X = np.column_stack([np.ones(2)] + [np.array([st[c][i1], st[c][i2]]) for c in cols])
        pred = X @ beta
        got = np.sign(pred[1] - pred[0])
        rows.append({"q": Q, "model": f"[witness edge] {name}",
                     "R2": float(got == truth)})
        print(f"  {name:24s} predicts psi(f_2) - psi(f_1) = {pred[1]-pred[0]:+.3e}  "
              f"-> {'RIGHT' if got == truth else 'wrong'}")

    with open("e2_f101_class.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["q", "model", "R2"])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
