#!/usr/bin/env python3
"""``C_sig(Q -> x)`` at 40 digits for the two binary quadratic classes over F_q.

Fibre signatures of a map ``F_q^2 -> F_q``:

* ``x``                 -- ``q`` fibres of size ``q``;
* ``Q`` anisotropic     -- the norm form of ``F_{q^2}``: one fibre of size ``1``
                           over ``0`` and ``q-1`` fibres of size ``q+1``;
* ``Q`` isotropic       -- one fibre of size ``2q-1`` over ``0`` and ``q-1``
                           fibres of size ``q-1``.

``x^2 + y^2`` is anisotropic exactly when ``-1`` is a non-residue, i.e.
``q = 3 mod 4``; ``x^2 - n y^2`` with ``n`` a non-residue is anisotropic for
every odd ``q``.

``C_sig(g -> f) = inf_{beta in [0, inf]} log Z_g(beta) / log Z_f(beta)`` with
``Z(beta) = sum_i n_i^beta``.  Every signature here sums to ``q^2``, so the
ratio is exactly ``1`` at ``beta = 1``; the endpoints ``beta = 0`` and
``beta = inf`` have the closed forms ``log(#fibres_g)/log(#fibres_f)`` and
``log(max_g)/log(max_f)``, which are printed for the check.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from mpmath import mp, mpf, log, exp, nstr

mp.dps = 80
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "isotropy_csig.csv"
GAP_OUTPUT = HERE / "isotropy_csig_gap.csv"

ODD_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)


def anisotropic_signature(q: int):
    return tuple([q + 1] * (q - 1) + [1])


def isotropic_signature(q: int):
    return tuple([2 * q - 1] + [q - 1] * (q - 1))


def linear_signature(q: int):
    return tuple([q] * q)


def log_partition(signature, beta):
    return log(sum(exp(beta * log(mpf(n))) for n in signature))


def ratio(source, target, beta):
    return log_partition(source, beta) / log_partition(target, beta)


def csig(source, target):
    """Infimum over ``beta in [0, inf]``, with the location of the contact."""

    zero = log(mpf(len(source))) / log(mpf(len(target)))
    infinity = log(mpf(max(source))) / log(mpf(max(target)))
    grid = [mpf(10) ** (mpf(e) / 400) for e in range(-1600, 1601)]
    values = sorted(((ratio(source, target, b), b) for b in grid),
                    key=lambda pair: pair[0])
    best, beta = values[0]
    ordered = sorted(grid)
    position = ordered.index(beta)
    low = ordered[max(position - 1, 0)]
    high = ordered[min(position + 1, len(ordered) - 1)]
    for _ in range(400):
        c = low + (high - low) / 3
        d = high - (high - low) / 3
        if ratio(source, target, c) <= ratio(source, target, d):
            high = d
        else:
            low = c
    interior = ratio(source, target, (low + high) / 2)
    candidates = [(zero, "beta = 0"), (infinity, "beta = inf"),
                  (interior, f"beta = {nstr((low + high) / 2, 12)}")]
    candidates.sort(key=lambda pair: pair[0])
    return candidates[0][0], candidates[0][1], zero, infinity


def main() -> None:
    rows = []
    print(f"{'q':>4} {'family':<14} {'aniso?':>7}  C_sig(Q -> x) at 40 digits")
    for q in ODD_PRIMES:
        target = linear_signature(q)
        for name, anisotropic in (("x^2 - n y^2", True),
                                  ("x y", False),
                                  ("x^2 + y^2", q % 4 == 3)):
            source = (anisotropic_signature(q) if anisotropic
                      else isotropic_signature(q))
            value, where, zero, infinity = csig(source, target)
            print(f"{q:>4} {name:<14} {str(anisotropic):>7}  {nstr(value, 40)}"
                  f"   ({where})")
            rows.append({
                "q": q, "family": name, "anisotropic": int(anisotropic),
                "signature": " ".join(str(n) for n in sorted(source, reverse=True)),
                "c_sig_40": nstr(value, 40),
                "contact": where,
                "endpoint_beta_0": nstr(zero, 40),
                "endpoint_beta_inf": nstr(infinity, 40),
                "c_aff": "2/3" if anisotropic else "1",
            })
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # how much more sharply C_aff separates isotropic from anisotropic than C_sig
    lookup = {(row["q"], row["family"]): mpf(row["c_sig_40"]) for row in rows}
    gaps = []
    print()
    print(f"{'q':>4} {'C_sig(xy->x)':>14} {'C_sig(aniso->x)':>17}"
          f" {'C_sig gap':>12} {'(1/3) / gap':>13}")
    for q in ODD_PRIMES:
        isotropic = lookup[(q, "x y")]
        anisotropic = lookup[(q, "x^2 - n y^2")]
        gap = isotropic - anisotropic
        amplification = mpf(1) / 3 / gap
        gaps.append({"q": q,
                     "c_sig_isotropic_40": nstr(isotropic, 40),
                     "c_sig_anisotropic_40": nstr(anisotropic, 40),
                     "c_sig_gap_40": nstr(gap, 40),
                     "c_aff_gap": "1/3",
                     "amplification": nstr(amplification, 20)})
        print(f"{q:>4} {float(isotropic):>14.10f} {float(anisotropic):>17.10f}"
              f" {float(gap):>12.8f} {float(amplification):>13.2f}")
    with GAP_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gaps[0]))
        writer.writeheader()
        writer.writerows(gaps)
    print()
    print(f"wrote {OUTPUT} and {GAP_OUTPUT}")


if __name__ == "__main__":
    main()
