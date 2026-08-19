#!/usr/bin/env python3
"""``N_1`` and ``N_2`` for ``Q -> x`` over ``F_q``: witnesses plus exhaustive checks.

Two families of binary quadratic forms:

* ``x^2 - n y^2`` with ``n`` the least non-residue -- anisotropic for every odd
  ``q`` (the norm form of ``F_{q^2}``);
* ``x^2 + y^2`` -- anisotropic exactly when ``q = 3 mod 4``.

For the anisotropic case the script builds the witnesses of Proposition A2

    N_1 <= 2 :  alpha_1(x,y) = (x, y),
                alpha_2(x,y) = (x + s, y + t)      with  M (s,t) = (1,0),
    N_2 <= 3 :  alpha_1 = (x_1, x_2),
                alpha_2 = (x_1 + s_2, x_2 + t_2)   with  M (s_2,t_2) = (1,0),
                alpha_3 = (x_1 + s_3, x_2 + t_3)   with  M (s_3,t_3) = (0,1),

where ``M = [[2a, b], [b, 2c]]`` is the polarisation matrix of
``Q = a u^2 + b u v + c v^2``, and re-verifies each on *all* ``q^2`` resp.
``q^4`` points by evaluating the actual functions.  For the isotropic case it
builds the one-step witness ``N_1 = 1``.

The exhaustive columns run the jet search of ``isotropy_atoms.py``, which is
feasible for ``k = 1`` up to ``q = 13`` and for ``k = 2`` up to ``q = 5``.
"""

from __future__ import annotations

import csv
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from isotropy_atoms import (  # noqa: E402
    JetSpace, atom_codes, is_anisotropic, least_non_residue,
    n_k_le_hull_plus_one, n_k_small, target_jets,
)

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "isotropy_family.csv"

PRIMES = (3, 5, 7, 11, 13)
K1_EXHAUSTIVE = (3, 5, 7, 11, 13)
K2_EXHAUSTIVE = (3, 5)


def polarisation_inverse(a, b, c, q):
    """Solve ``M (s,t) = target`` for ``M = [[2a,b],[b,2c]]``."""

    det = (4 * a * c - b * b) % q
    assert det % q != 0, "degenerate form"
    inverse_det = pow(det, q - 2, q)

    def solve(target):
        u, v = target
        s = (2 * c * u - b * v) * inverse_det % q
        t = (-b * u + 2 * a * v) * inverse_det % q
        return s, t

    return solve


def evaluate(coefficients, u, v, q):
    a, b, c, d, e = coefficients
    return (a * u * u + b * u * v + c * v * v + d * u + e * v) % q


def verify_k1(coefficients, q):
    """Check ``x = h_2 - h_1 - Q(s,t)`` on all ``q^2`` points of ``F_q^2``."""

    a, b, c = coefficients[:3]
    solve = polarisation_inverse(a, b, c, q)
    s, t = solve((1, 0))
    constant = evaluate(coefficients, s, t, q)
    for x, y in product(range(q), repeat=2):
        h1 = evaluate(coefficients, x, y, q)
        h2 = evaluate(coefficients, (x + s) % q, (y + t) % q, q)
        if (h2 - h1 - constant) % q != x % q:
            return False, (s, t)
    return True, (s, t)


def verify_k2(coefficients, q):
    """Check ``x_1`` and ``x_2`` from three atoms on all ``q^4`` points."""

    a, b, c = coefficients[:3]
    solve = polarisation_inverse(a, b, c, q)
    s2, t2 = solve((1, 0))
    s3, t3 = solve((0, 1))
    constant2 = evaluate(coefficients, s2, t2, q)
    constant3 = evaluate(coefficients, s3, t3, q)
    for x1, y1, x2, y2 in product(range(q), repeat=4):
        h1 = evaluate(coefficients, x1, x2, q)
        h2 = evaluate(coefficients, (x1 + s2) % q, (x2 + t2) % q, q)
        h3 = evaluate(coefficients, (x1 + s3) % q, (x2 + t3) % q, q)
        if (h2 - h1 - constant2) % q != x1 % q:
            return False, None
        if (h3 - h1 - constant3) % q != x2 % q:
            return False, None
    return True, ((s2, t2), (s3, t3))


def verify_isotropic_one_step(coefficients, q):
    """For isotropic ``Q`` find an affine ``alpha`` with ``Q o alpha`` linear in x."""

    a, b, c = coefficients[:3]
    zeros = [(u, v) for u in range(q) for v in range(q)
             if (u, v) != (0, 0) and (a * u * u + b * u * v + c * v * v) % q == 0]
    for (u0, v0) in zeros:
        for u1, v1 in product(range(q), repeat=2):
            values = []
            for x in range(q):
                u = (x * u0 + u1) % q
                v = (x * v0 + v1) % q
                values.append((a * u * u + b * u * v + c * v * v) % q)
            differences = {(values[(i + 1) % q] - values[i]) % q for i in range(q)}
            if len(set(values)) == q and len(differences) == 1:
                # Q o alpha is x |-> lambda x + mu with lambda invertible
                return True, ((u0, v0), (u1, v1))
    return False, None


def main() -> None:
    rows = []
    for q in PRIMES:
        n = least_non_residue(q)
        families = {
            f"x^2 - {n} y^2": (1, 0, (-n) % q, 0, 0),
            "x^2 + y^2": (1, 0, 1, 0, 0),
        }
        for name, coefficients in families.items():
            anisotropic = is_anisotropic(coefficients, q)
            record = {"q": q, "form": name, "anisotropic": int(anisotropic),
                      "q_mod_4": q % 4}
            if anisotropic:
                ok1, shift1 = verify_k1(coefficients, q)
                ok2, shifts2 = verify_k2(coefficients, q)
                record["witness_N1_le_2"] = int(ok1)
                record["witness_N2_le_3"] = int(ok2)
                record["shift_k1"] = str(shift1)
                record["shifts_k2"] = str(shifts2)
                record["witness_N1_eq_1"] = ""
            else:
                ok, data = verify_isotropic_one_step(coefficients, q)
                record["witness_N1_le_2"] = ""
                record["witness_N2_le_3"] = ""
                record["shift_k1"] = ""
                record["shifts_k2"] = ""
                record["witness_N1_eq_1"] = f"{int(ok)} {data}"

            if q in K1_EXHAUSTIVE:
                start = time.time()
                space = JetSpace.make(q, 2)
                atoms = atom_codes(space, coefficients)
                got = n_k_small(space, atoms, target_jets(space, (0, 0, 0, 1, 0), 1), 4)
                record["N1_exhaustive"] = str(got)
                record["N1_seconds"] = f"{time.time()-start:.1f}"
            else:
                record["N1_exhaustive"] = ""
                record["N1_seconds"] = ""

            if q in K2_EXHAUSTIVE:
                start = time.time()
                space = JetSpace.make(q, 4)
                atoms = atom_codes(space, coefficients)
                got = n_k_le_hull_plus_one(space, atoms,
                                           target_jets(space, (0, 0, 0, 1, 0), 2))
                record["N2_le3_exhaustive"] = str(got) if got is not None else ">3"
                record["N2_seconds"] = f"{time.time()-start:.1f}"
            else:
                record["N2_le3_exhaustive"] = ""
                record["N2_seconds"] = ""

            record["C_aff_proved"] = "2/3" if anisotropic else "1"
            rows.append(record)
            print(f"q={q:>3} {name:<12} anisotropic={anisotropic!s:<5} "
                  f"N1={record['N1_exhaustive']:<5} "
                  f"N2<=3? {record['N2_le3_exhaustive']:<3} "
                  f"witnesses {record['witness_N1_le_2']}{record['witness_N2_le_3']}"
                  f"{record['witness_N1_eq_1']}", flush=True)

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print()
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
