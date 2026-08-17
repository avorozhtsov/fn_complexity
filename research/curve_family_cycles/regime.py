#!/usr/bin/env python3
"""Why the flat regime permits cycles, quantitatively.

The brief for this search expected the opposite result: fiber signatures of
curve families are nearly flat, ``C(L->f)`` sits at the ``beta = infinity``
endpoint, and the endpoint-regime theorem forbids cycles there.  What that
argument misses is that flatness makes the endpoint data *degenerate* rather
than decisive.  This script measures the three quantities that decide it.

Write ``a_c = q - N_c``, ``alpha_c = -a_c/sqrt(q)``, and put

    Lambda_f(beta) = log( (1/q) sum_c (N_c/q)^beta ),
    G_uv(beta)     = (Lambda_u(beta) - Lambda_v(beta)) / (1 + beta),

so that, exactly,

    log Z_f(beta) = (1 + beta) log q + Lambda_f(beta),
    C(u->v)       = 1 + inf_beta G_uv(beta) / log q  + O(1/(q log q)).

``G_uv(0) = 0`` always, and ``G_uv(infinity) = log(max_u / max_v)`` is the only
thing the index ``phi`` sees.  Since ``C(v->u) = 1 + inf(-G_uv)/log q``,

    u < v   iff   min_beta G_uv + max_beta G_uv < 0,

a *midrange* comparison, and the margin is that midrange over ``log q``.

Substituting ``beta = tau sqrt(q)`` turns ``Lambda`` into the cumulant
generating function of the normalised traces and gives

    G_uv(beta) = (Psi_u(tau) - Psi_v(tau)) / sqrt(q) + O(1/q),
    Psi_f(tau) = log( (1/q) sum_c exp(tau alpha_c) ) / tau,

an increasing function from ``Psi(0) = 0`` --- the first moment vanishes
identically --- to ``Psi(infinity) = alpha_max = (max_c N_c - q)/sqrt(q)``.

Three consequences are measured below.

1.  Interior structure lives at ``beta ~ sqrt(q)`` and has size ``1/sqrt(q)``,
    not the ``1/q`` of the ``beta = O(1)`` bottleneck of T2.2.
2.  The endpoint gap that ``phi`` reads between neighbouring classes is
    ``log((N+1)/N) ~ 1/q``, a factor ``sqrt(q)`` smaller.  So an interior
    tangency can overturn a ``phi``-gap of up to ``~ sqrt(q)`` units of largest
    fiber, and ``phi``-violations are generic rather than marginal.
3.  The midrange law predicts the observed margins.

    python research/curve_family_cycles/regime.py
"""

from __future__ import annotations

import collections
import csv
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from common import Engine, beta_grid, build_pool  # noqa: E402

QS = (31, 61, 101, 211, 401)
CLASS_CAP = 120


def reduced_curves(signatures: list[tuple[int, ...]], q: int, betas: np.ndarray) -> np.ndarray:
    """``Lambda_f(beta)`` for every signature, on a shared grid."""

    out = np.empty((len(signatures), betas.size))
    for index, signature in enumerate(signatures):
        ell = np.log(np.asarray(signature, dtype=float) / q)
        top = ell.max()
        out[index] = betas * top + np.log(np.exp(np.outer(ell - top, betas)).sum(axis=0) / len(ell))
    return out


def main() -> int:
    rows = []
    print("Interior structure of the comparison on nearly flat signatures\n")
    header = (
        f"{'q':>5} {'pairs':>8} {'median beta*':>13} {'beta*/sqrt q':>13} "
        f"{'median |1-C|':>13} {'x sqrt(q) log q':>16} {'midrange law':>14}"
    )
    print(header)
    print("-" * len(header))
    for q in QS:
        pool = build_pool(
            q,
            budget=dict(hyperelliptic=1500, superelliptic=600, twist=600, additive=600, dense=600),
        )
        classes = collections.defaultdict(list)
        for entry in pool:
            classes[entry.max_fiber].append(entry)
        members = max(classes.values(), key=len)[:CLASS_CAP]
        signatures = [entry.signature for entry in members]

        betas = beta_grid(q, points=30_000)
        engine = Engine(signatures, betas)
        rates, contacts = engine.rate_matrix(chunk=64)
        margins = rates - rates.T
        upper = np.triu_indices(len(signatures), 1)

        # the midrange law, evaluated on the same grid
        reduced = reduced_curves(signatures, q, betas)
        weight = 1.0 / (1.0 + betas)
        predicted = np.empty_like(margins)
        for i in range(len(signatures)):
            g = (reduced[i][None, :] - reduced) * weight
            predicted[i] = (g.min(axis=1) + g.max(axis=1)) / math.log(q)
        error = np.abs(predicted[upper] - margins[upper]) / np.abs(margins[upper])

        interior = contacts[np.isfinite(contacts) & (contacts > 0)]
        deviation = np.abs(1.0 - rates[upper])
        median_beta = float(np.median(interior))
        median_dev = float(np.median(deviation))
        print(
            f"{q:>5} {len(upper[0]):>8} {median_beta:>13.3f} "
            f"{median_beta / math.sqrt(q):>13.3f} {median_dev:>13.3e} "
            f"{median_dev * math.sqrt(q) * math.log(q):>16.4f} "
            f"{float(np.median(error)) * 100:>13.2f}%"
        )
        rows.append(
            [
                q,
                len(signatures),
                f"{median_beta:.6f}",
                f"{median_beta / math.sqrt(q):.6f}",
                f"{median_dev:.9e}",
                f"{median_dev * math.sqrt(q) * math.log(q):.6f}",
                f"{float(np.median(error)):.6f}",
            ]
        )

    print(
        "\n  beta*/sqrt(q) is flat in q and |1-C| sqrt(q) log q is flat in q: the\n"
        "  contact sits at tau = beta/sqrt(q) of order one, and the interior\n"
        "  correction is Theta(1/(sqrt(q) log q)).  The last column is the median\n"
        "  relative error of the midrange law against the computed margins."
    )

    print("\n\nHow large a phi-gap an interior tangency overturns\n")
    header = (
        f"{'q':>5} {'signatures':>11} {'phi-violating':>14} "
        f"{'max overturned':>15} {'sqrt(q)':>9}"
    )
    print(header)
    print("-" * len(header))
    for q in QS:
        pool = build_pool(
            q,
            budget=dict(hyperelliptic=1500, superelliptic=600, twist=600, additive=600, dense=600),
        )
        # one representative per largest-fiber class keeps the matrix small while
        # covering the whole phi range, which is what a violation search needs.
        by_max: dict[int, list] = collections.defaultdict(list)
        for entry in pool:
            by_max[entry.max_fiber].append(entry)
        chosen = [entry for members in by_max.values() for entry in members[:12]]
        signatures = [entry.signature for entry in chosen]
        rates, _ = comparison_matrix(signatures, q)
        difference = rates - rates.T
        strict = difference < -1e-9
        tops = np.array([max(s) for s in signatures])
        gap = tops[:, None] - tops[None, :]
        violating = strict & (gap > 0)  # u wins although u has the larger max fiber
        overturned = int(gap[violating].max()) if violating.any() else 0
        print(
            f"{q:>5} {len(signatures):>11} {int(violating.sum()):>14} "
            f"{overturned:>15} {math.sqrt(q):>9.2f}"
        )
        rows.append([q, len(signatures), "", "", "", "", f"violations={int(violating.sum())}"])

    print(
        "\n  The largest overturned gap tracks sqrt(q), as the two scales predict:\n"
        "  the endpoint gap per unit of largest fiber is log((N+1)/N) ~ 1/q while\n"
        "  the interior dip is ~ 1/sqrt(q)."
    )

    with (HERE / "regime.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "q",
                "signatures",
                "median_contact_beta",
                "median_contact_beta_over_sqrt_q",
                "median_abs_1_minus_C",
                "rescaled",
                "midrange_law_relative_error",
            ]
        )
        writer.writerows(rows)
    print("\nwritten: regime.csv")
    return 0


def comparison_matrix(signatures: list[tuple[int, ...]], q: int):
    return Engine(signatures, beta_grid(q, points=20_000)).rate_matrix(chunk=64)


if __name__ == "__main__":
    raise SystemExit(main())
