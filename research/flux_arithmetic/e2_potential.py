#!/usr/bin/env python3
"""E2 -- what is the least-squares potential ``psi_opt`` of the flux?

Three things, and the third is the one that matters.

1. Nested least-squares regressions of ``psi_opt`` on the statistics the rate is
   known to see.  **Read the table with the caveat printed alongside it:** at a
   fixed ``q`` every pool member has exactly ``q`` fibers, so ``M = max_c N_c -
   q``, ``log max_c N_c`` and ``(1/2) log phi = (1/2) log(log q * log max_c
   N_c)`` are three *link functions of the same integer* ``max_c N_c``.  The
   contest between them is a contest of shapes, not of statistics, and the
   honest ceiling for the family is the categorical fit -- the best function of
   ``max_c N_c`` whatever it is.

2. The unfitted test.  ``psi_end(a) = (1/2) log(log q * log max_a)`` is not a
   regressor here, it is the exact potential the flux has wherever both infima
   sit at an endpoint (brief D Part 0(c)).  Compare ``||A - grad psi_end||``
   with ``||A - grad psi_opt||`` directly, with no fitting at all, and read off
   the slope of ``psi_opt`` against ``psi_end``.

3. Brief B's addendum scalar ``phi~ = M - ((3-2 sqrt 2)/2) m2``, which the
   addendum derived in the ``beta = O(1)`` regime and which brief B's own
   measurement places ``sqrt q`` away from the operative scale.

Writes ``e2_regressions.csv`` and ``e2_potential.csv``.
"""

from __future__ import annotations

import csv
import math

import numpy as np

import flux
import pools

COARSE, FINE = 2400, 13600


def statistics_of(S: np.ndarray, q: int) -> dict:
    N = S.astype(float)
    maxN = N.max(axis=1)
    minN = N.min(axis=1)
    a = q - N
    return {
        "maxN": maxN,
        "M": maxN - q,
        "log_maxN": np.log(maxN),
        "half_log_phi": 0.5 * np.log(math.log(q) * np.log(maxN)),
        "m2": (a**2).sum(axis=1) / q**2,
        "m3": (a**3).sum(axis=1) / q**3,
        "m4": (a**4).sum(axis=1) / q**4,
        "log_mult": np.log((N == maxN[:, None]).sum(axis=1).astype(float)),
        "minN": minN,
        "log_minN": np.log(minN),
    }


def r2(y: np.ndarray, cols: list[np.ndarray]) -> float:
    X = np.column_stack([np.ones(len(y))] + cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return 1.0 - resid.var() / y.var()


def r2_categorical(y: np.ndarray, key: np.ndarray) -> float:
    """The best possible fit by *any* function of ``key`` -- the group means."""

    total = 0.0
    for value in np.unique(key):
        mask = key == value
        total += ((y[mask] - y[mask].mean()) ** 2).sum()
    return 1.0 - total / len(y) / y.var()


MODELS = [
    ("M", ["M"]),
    ("log maxN", ["log_maxN"]),
    ("(1/2) log phi", ["half_log_phi"]),
    ("m2", ["m2"]),
    ("M, m2  (addendum phi~)", ["M", "m2"]),
    ("M, m2, log mult", ["M", "m2", "log_mult"]),
    ("M, m2, m3, m4, log mult", ["M", "m2", "m3", "m4", "log_mult"]),
    ("(1/2) log phi, m2", ["half_log_phi", "m2"]),
    ("(1/2) log phi, log mult", ["half_log_phi", "log_mult"]),
    ("(1/2) log phi, m2, m3, m4, log mult",
     ["half_log_phi", "m2", "m3", "m4", "log_mult"]),
    ("(1/2) log phi, log minN", ["half_log_phi", "log_minN"]),
]


def order_agreement(A: np.ndarray, scalar: np.ndarray) -> float:
    """Fraction of ordered pairs whose direction a scalar gets right.

    A tie in the scalar is a failure, not a half-success: the comparison ``A``
    is strict on every pair of these pools.  This is the metric that matters --
    two scalars that are monotone functions of each other score identically,
    however different their ``R^2`` against ``psi_opt``.
    """

    G = scalar[None, :] - scalar[:, None]
    mask = np.abs(A) > flux.TIE
    return float((np.sign(A[mask]) == np.sign(G[mask])).mean())


def fit(y: np.ndarray, cols: list[np.ndarray]) -> np.ndarray:
    X = np.column_stack([np.ones(len(y))] + cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


def run(q: int, S: np.ndarray, tag: str, rows: list[dict], pot_rows: list[dict]) -> None:
    betas = flux.beta_grid(q, COARSE, FINE)
    A = flux.flux_matrix(S, betas)
    h = flux.hodge(A)
    psi = h["psi"]
    st = statistics_of(S, q)
    nA = np.linalg.norm(A)

    print(f"\n=== q={q} ({tag}), n={len(S)} ===")
    print(f"  |grad|/|A| = {h['grad_frac']:.6f}   |curl|/|A| = {h['curl_frac']:.6f}")
    print("  psi_opt regressions")
    for name, cols in MODELS:
        value = r2(psi, [st[c] for c in cols])
        rows.append({"q": q, "pool": tag, "n": len(S), "model": name, "R2": value})
        print(f"    R2 = {value:.6f}   {name}")
    cat_max = r2_categorical(psi, st["maxN"])
    cat_max_mult = r2_categorical(
        psi, st["maxN"] * 1000 + np.exp(st["log_mult"]).round()
    )
    for name, value in (
        ("[ceiling] any function of maxN", cat_max),
        ("[ceiling] any function of (maxN, mult)", cat_max_mult),
    ):
        rows.append({"q": q, "pool": tag, "n": len(S), "model": name, "R2": value})
        print(f"    R2 = {value:.6f}   {name}")

    # --- the unfitted comparison
    psi_end = st["half_log_phi"]
    psi_end = psi_end - psi_end.mean()
    G_opt = psi[None, :] - psi[:, None]
    G_end = psi_end[None, :] - psi_end[:, None]
    slope = float(np.polyfit(psi_end, psi, 1)[0])
    scaled = psi_end * float((psi @ psi_end) / (psi_end @ psi_end))
    G_scaled = scaled[None, :] - scaled[:, None]
    rec = {
        "q": q,
        "pool": tag,
        "n": len(S),
        "grad_frac": h["grad_frac"],
        "curl_frac": h["curl_frac"],
        "order_agreement": h["order_agreement"],
        "resid_opt": float(np.linalg.norm(A - G_opt) / nA),
        "resid_end_unfitted": float(np.linalg.norm(A - G_end) / nA),
        "resid_end_scaled": float(np.linalg.norm(A - G_scaled) / nA),
        "slope_psiopt_vs_psiend": slope,
        "R2_half_log_phi": r2(psi, [psi_end]),
        "ceiling_maxN": cat_max,
        "order_agreement_end": float(
            (np.sign(A[np.abs(A) > flux.TIE]) == np.sign(G_end[np.abs(A) > flux.TIE])).mean()
        ),
    }
    print(f"  unfitted endpoint potential: ||A - grad psi_end||/||A|| = "
          f"{rec['resid_end_unfitted']:.6f}  (best scalar multiple "
          f"{rec['resid_end_scaled']:.6f}; least squares {rec['resid_opt']:.6f})")
    print(f"  slope of psi_opt on psi_end = {slope:.6f}")

    # --- the addendum's fixed coefficient against the fitted one
    kappa = (3 - 2 * math.sqrt(2)) / 2
    phi_tilde = st["M"] - kappa * st["m2"]
    beta_free = fit(psi, [st["M"], st["m2"]])
    ratio = beta_free[2] / beta_free[1] if beta_free[1] else float("nan")
    rec["R2_phi_tilde_fixed"] = r2(psi, [phi_tilde])
    rec["fitted_m2_over_M"] = float(ratio)
    rec["addendum_m2_over_M"] = -kappa
    print(f"  addendum: phi~ with its own coefficient -{kappa:.8f} scores "
          f"R2 = {rec['R2_phi_tilde_fixed']:.6f}; the free fit of (M, m2) puts "
          f"the ratio at {ratio:+.6f}")

    # --- what actually matters: the induced order
    best_two = fit(psi, [psi_end, st["m2"]])
    scalars = {
        "M  (= any function of max fiber)": st["M"],
        "log maxN": st["log_maxN"],
        "(1/2) log phi": psi_end,
        "phi~ = M - 0.0858 m2": phi_tilde,
        "(1/2) log phi + fitted m2": best_two[1] * psi_end + best_two[2] * st["m2"],
        "psi_opt (least squares)": psi,
    }
    print("  order agreement (fraction of ordered pairs given the right sign)")
    for name, vector in scalars.items():
        value = order_agreement(A, vector)
        rec[f"order::{name}"] = value
        print(f"    {value:.6f}   {name}")
    # how the misses split: ties of phi against genuine violations
    Gend = psi_end[None, :] - psi_end[:, None]
    mask = np.abs(A) > flux.TIE
    tied = int((mask & (np.abs(Gend) < 1e-14)).sum())
    violating = int((mask & (np.abs(Gend) >= 1e-14) & (np.sign(A) != np.sign(Gend))).sum())
    rec["phi_blind_ordered_pairs"] = tied
    rec["phi_violating_ordered_pairs"] = violating
    print(f"    of the {tied + violating} ordered pairs phi gets wrong, "
          f"{tied} are exact phi-ties and {violating} are genuine violations")
    pot_rows.append(rec)


def main() -> None:
    rows: list[dict] = []
    pot_rows: list[dict] = []
    for q in (11, 13):
        S, _ = pools.arithmetic_pool(q)
        run(q, S, "exhaustive genus 2", rows, pot_rows)
    for q in (17, 19, 23):
        S = pools.sampled_pool(q, draws=12000)
        if len(S) > 600:
            rng = np.random.default_rng(q)
            S = S[np.sort(rng.choice(len(S), 600, replace=False))]
        run(q, S, "sampled genus 2", rows, pot_rows)
    # a mixed-genus pool, to see whether genus enters at all
    for q in (13,):
        S = pools.sampled_pool(q, degrees=(5, 6, 7, 8, 9, 10), draws=8000)
        if len(S) > 600:
            rng = np.random.default_rng(q + 1)
            S = S[np.sort(rng.choice(len(S), 600, replace=False))]
        run(q, S, "sampled genus 2-4", rows, pot_rows)

    with open("e2_regressions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open("e2_potential.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pot_rows[0]))
        w.writeheader()
        w.writerows(pot_rows)


if __name__ == "__main__":
    main()
