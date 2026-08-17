#!/usr/bin/env python3
"""T2.2 survey: rate vectors versus trace moments over a large pool of maps.

Builds a pool of maps ``f : A^2 -> A^1`` over ``F_q`` (structured families plus
random polynomials), records the signature, the trace moments and the vector of
exchange rates against a small fixed reference family, then asks how well the
rate vector predicts ``m2``, ``m3`` and the extreme traces.

Two regimes are reported separately:

*   **Weil regime** -- every fiber non-empty and ``max_c |a_c| <= 6 sqrt(q)``,
    i.e. the fibers really are curves of small genus.  Here the ``1/q``
    expansion of ``log Z_f`` is valid and the questions have sharp answers.
*   **All maps** -- includes pushforwards like ``f = x^d`` whose fibers are
    empty or of size ``dq``; the expansion is meaningless there and the rates
    are dominated by the leading-order mismatch.

Run:  python research/m_and_e_and_a_c/t2_2_survey.py [q ...]
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from t2_2_common import (  # noqa: E402
    BETA_STAR, KAPPA, Sig, classify, moments, rate, references, sig_from_counts,
)
import ffmaps as F  # noqa: E402

PROBES = [
    ("C(L->f)", "L", "impl"),
    ("C(f->L)", "L", "impd"),
    ("C(Xsplit->f)", "Xsplit", "impl"),
    ("C(f->Xsplit)", "Xsplit", "impd"),
    ("C(Xaniso->f)", "Xaniso", "impl"),
    ("C(f->Xaniso)", "Xaniso", "impd"),
    ("C(Sq->f)", "Sq", "impl"),
    ("C(f->Sq)", "Sq", "impd"),
]
PROBE_NAMES = [p[0] for p in PROBES]


def rate_vector(f: Sig, refs: dict[str, Sig]) -> tuple[list[float], list[float]]:
    vals, betas = [], []
    for _, ref, role in PROBES:
        g = refs[ref]
        v, b = rate(g, f) if role == "impl" else rate(f, g)
        vals.append(v)
        betas.append(b)
    return vals, betas


def build_pool(q: int, rng: np.random.Generator, n_random: int = 250) -> list[dict]:
    pool: list[dict] = []

    def add(tag: str, detail: str, counts: np.ndarray) -> None:
        assert counts.sum() == q * q
        pool.append({"tag": tag, "detail": detail, "counts": counts})

    for deg, label in ((3, "ell"), (5, "g2"), (7, "g3"), (9, "g4"), (11, "g5")):
        for _ in range(60):
            c = [int(v) for v in rng.integers(0, q, size=deg)] + [1]
            add(label, f"y^2=P{deg}+c {c}", F.hyperelliptic(q, c))
    for d in (3, 4, 5, 6, 7, 8, 9):
        add("monomial", f"y^2=x^{d}+c", F.hyperelliptic(q, [0] * d + [1]))
    for r in (2, 3, 4, 5):
        for d in (2, 3, 4, 5, 6, 7):
            add("superell", f"y^{r}=x^{d}+c", F.superelliptic(q, r, [0] * d + [1]))
    for _ in range(80):
        dp, dq = int(rng.integers(2, 8)), int(rng.integers(2, 8))
        pc = [int(v) for v in rng.integers(0, q, size=dp)] + [1]
        qc = [int(v) for v in rng.integers(0, q, size=dq)] + [1]
        add("additive", f"P{dp}(x)+Q{dq}(y)", F.additive(q, pc, qc))
    for d in (2, 3, 4, 5, 6):
        A = np.zeros((d + 1, 1), dtype=np.int64)
        A[d, 0] = 1
        add("pushfwd", f"f=x^{d}", F.bilinear_family(q, A))
    for _ in range(n_random):
        dx, dy = int(rng.integers(2, 6)), int(rng.integers(2, 6))
        add("random", f"bidegree<({dx},{dy})", F.bilinear_family(q, rng.integers(0, q, size=(dx, dy))))
    return pool


def analyse(q: int, rng: np.random.Generator) -> dict:
    refs = references(q)
    rows = []
    seen: dict[tuple, list[str]] = {}
    for item in build_pool(q, rng):
        counts = item["counts"]
        sig = sig_from_counts(counts)
        key = (tuple(int(v) for v in sig.values), tuple(int(m) for m in sig.mults))
        seen.setdefault(key, []).append(item["tag"])
        mom = moments(counts, q)
        vals, betas = rate_vector(sig, refs)
        rows.append({"tag": item["tag"], "detail": item["detail"], "sig": sig,
                     "key": key, "mom": mom, "rates": vals, "betas": betas,
                     "logz": None})
    return {"q": q, "refs": refs, "rows": rows, "seen": seen}


# ------------------------------------------------------------------ analysis


def regress(X: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    A = np.hstack([np.ones((X.shape[0], 1)), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    return r2, float(np.abs(resid).max())


def nearest_neighbour(X: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Leave-one-out 1-NN prediction of y from the (standardised) rate vector."""
    Z = (X - X.mean(0)) / np.maximum(X.std(0), 1e-15)
    d = ((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d, np.inf)
    j = d.argmin(1)
    err = np.abs(y - y[j])
    return float(np.median(err)), float(err.max())


def two_probe_moments(row: dict, q: int, refs: dict[str, Sig]) -> tuple[float, float] | None:
    """Recover (m2, m3) from two interior tangencies via the 1/q expansion.

    Each interior probe gives ``log Z_f(beta0)`` exactly; subtracting the flat
    part leaves ``O(beta) = q(log Z_f(beta) - (beta+1) log q)
                          = C(beta,2) m2 - C(beta,3) m3 / sqrt(q) + O(1/q)``.
    """
    picks = []
    for name, ref, role in PROBES:
        i = PROBE_NAMES.index(name)
        b = row["betas"][i]
        if math.isinf(b) or b <= 1e-6 or b > 8.0:
            continue
        if classify(refs[ref], row["sig"]) != "interior" and role == "impl":
            continue
        if classify(row["sig"], refs[ref]) != "interior" and role == "impd":
            continue
        picks.append(b)
    picks = sorted(set(round(b, 9) for b in picks))
    if len(picks) < 2:
        return None
    b1, b2 = picks[0], picks[-1]
    if abs(b1 - b2) < 0.2:
        return None
    obs = []
    A = []
    for b in (b1, b2):
        lz = float(row["sig"].log_z(np.array([b]))[0])
        obs.append(q * (lz - (b + 1) * math.log(q)))
        A.append([b * (b - 1) / 2, -b * (b - 1) * (b - 2) / 6 / math.sqrt(q)])
    sol = np.linalg.solve(np.array(A), np.array(obs))
    return float(sol[0]), float(sol[1])


def report(q: int, res: dict) -> dict:
    rows, refs, seen = res["rows"], res["refs"], res["seen"]
    print(f"\n{'=' * 80}\nq = {q}   pool = {len(rows)} maps   "
          f"distinct signatures = {len(seen)}")
    by_tag: dict[str, list[dict]] = {}
    for r in rows:
        by_tag.setdefault(r["tag"], []).append(r)
    print("  signature collisions by family:")
    for tag, group in sorted(by_tag.items()):
        keys = {r["key"] for r in group}
        print(f"    {tag:10s} {len(group):4d} maps -> {len(keys):4d} distinct signatures")

    weil = [r for r in rows
            if r["mom"]["n_image"] == q and r["mom"]["a_absmax"] <= 6 * math.sqrt(q)]
    print(f"\n  Weil regime (all fibers non-empty, |a_c| <= 6 sqrt q): {len(weil)} maps")

    # which probes are interior at all
    print("  probe tangency type on the Weil subpool (fraction interior):")
    for name, ref, role in PROBES:
        kinds = [classify(refs[ref], r["sig"]) if role == "impl"
                 else classify(r["sig"], refs[ref]) for r in weil]
        frac = sum(k == "interior" for k in kinds) / max(len(kinds), 1)
        i = PROBE_NAMES.index(name)
        bs = [r["betas"][i] for r, k in zip(weil, kinds) if k == "interior"]
        rng_s = f"beta in [{min(bs):.4f}, {max(bs):.4f}]" if bs else "-"
        print(f"    {name:14s} interior {frac:6.1%}   {rng_s}")

    m2 = np.array([r["mom"]["m2"] for r in weil])
    m3 = np.array([r["mom"]["m3"] for r in weil])
    c_fL = np.array([r["rates"][PROBE_NAMES.index("C(f->L)")] for r in weil])
    b_fL = np.array([r["betas"][PROBE_NAMES.index("C(f->L)")] for r in weil])
    ok = np.isfinite(b_fL) & (b_fL > 0)
    print(f"\n  argmin beta of C(f->L): [{b_fL[ok].min():.6f}, {b_fL[ok].max():.6f}]  "
          f"(sqrt2-1 = {BETA_STAR:.6f}); {int((~ok).sum())} degenerate (m2 = 0)")

    m2_hat = 2 * q * math.log(q) * (1.0 - c_fL) / KAPPA
    rel = np.abs(m2_hat - m2) / np.maximum(np.abs(m2), 1e-12)
    good = m2 > 1e-9
    print(f"  m2 from C(f->L) alone: median rel err {np.median(rel[good]):.3e}, "
          f"max {rel[good].max():.3e}   (q^-1/2 = {q ** -0.5:.3e})")

    # two-probe recovery of (m2, m3)
    e2, e3 = [], []
    for r, mm2, mm3 in zip(weil, m2, m3):
        got = two_probe_moments(r, q, refs)
        if got is None:
            continue
        e2.append(abs(got[0] - mm2))
        e3.append(abs(got[1] - mm3))
    if e2:
        print(f"  two-probe (m2, m3) recovery on {len(e2)} maps: "
              f"median |dm2| = {np.median(e2):.3e} (max {max(e2):.3e}), "
              f"median |dm3| = {np.median(e3):.3e} (max {max(e3):.3e})")
        print(f"    spread of true m3 over the subpool: "
              f"{m3.min():.4f} .. {m3.max():.4f}")

    X = np.array([r["rates"] for r in weil])
    print("\n  8-rate vector -> moment, on the Weil subpool:")
    for target in ("m2", "m3", "m4", "a_absmax", "amin", "amax"):
        y = np.array([r["mom"][target] for r in weil])
        r2, mx = regress(X, y)
        nn_med, nn_max = nearest_neighbour(X, y)
        print(f"    {target:9s} spread {y.max() - y.min():10.4f}  "
              f"lin R^2 = {r2:10.7f}  max|resid| = {mx:9.3e}   "
              f"1-NN median err = {nn_med:9.3e}")

    Xall = np.array([r["rates"] for r in rows])
    print("\n  8-rate vector -> moment, on the FULL pool (expansion invalid):")
    for target in ("m2", "m3", "a_absmax"):
        y = np.array([r["mom"][target] for r in rows])
        r2, mx = regress(Xall, y)
        nn_med, nn_max = nearest_neighbour(Xall, y)
        print(f"    {target:9s} spread {y.max() - y.min():10.4f}  "
              f"lin R^2 = {r2:10.7f}  max|resid| = {mx:9.3e}   "
              f"1-NN median err = {nn_med:9.3e}")

    # exact integer readouts from the endpoints
    nmax_from_rate = [round(math.exp(math.log(q) / r["rates"][PROBE_NAMES.index("C(L->f)")]))
                      for r in weil]
    exact = sum(a == r["mom"]["max_fiber"] for a, r in zip(nmax_from_rate, weil))
    print(f"\n  largest fiber recovered from C(L->f) = log q / log N_max: "
          f"{exact}/{len(weil)} exact")
    nfib = [round(math.exp(math.log(q) / r["rates"][PROBE_NAMES.index("C(Sq->f)")]
                           * math.log((q + 1) // 2) / math.log(q)))
            for r in rows]
    return {"q": q, "pool": len(rows), "distinct": len(seen), "weil": len(weil)}


def main(qs: list[int]) -> None:
    rng = np.random.default_rng(20260817)
    out = {}
    for q in qs:
        out[q] = report(q, analyse(q, rng))
    (HERE / "t2_2_survey_summary.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [101, 211, 503])
