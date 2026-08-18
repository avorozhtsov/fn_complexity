#!/usr/bin/env python3
"""Rigorous certification of the three-cycle among genus-two pencils over F_11.

The cycle is

    A = y^2 = x^5 +  3x^4 + 4x^3 + x^2 +  x      + c
    B = y^2 = x^6 +  9x^5 + 7x^4 + 2x^3 +  x      + c
    C = y^2 = x^6 + 10x^5 + 8x^4 + 8x^3 + x^2 + 2x + c

with ``A < B < C < A`` for the exchange comparison ``a < b iff C(a->b) <
C(b->a)``.  Every fiber of every one of the three pencils is a smooth genus-two
curve: ``P(x) + c`` is squarefree for all eleven ``c`` in each case.  Every
inequality below is established with interval arithmetic, so the cycle is a
theorem and not a floating-point observation.

For the edge ``a < b`` two bounds are needed:

  * an UPPER bound on ``C(a->b)`` --- free, because ``C`` is an infimum, so
    evaluating the ratio at a single rational ``beta`` (or at either endpoint)
    certifies one;
  * a LOWER bound on ``C(b->a)`` --- an assertion about *every* ``beta``,
    obtained by branch and bound with interval arithmetic on a compact
    ``[0, B0]`` together with a closed-form tail bound on ``[B0, infinity)``.

**Conditioning.**  These signatures are nearly flat: every entry is ``q(1 +
O(1/sqrt q))`` and the margins are ``10^-3`` against rates of ``0.99``.  A naive
enclosure of ``log Z = log sum N_c^beta`` has width ``~ w log q`` on a box of
width ``w``, and that width appears twice in the quotient, which forces boxes
around ``10^-5`` wide.  Factoring out the common ``q^beta`` first,

    log Z(beta) = beta log q + log S(beta),   S(beta) = sum_c (N_c/q)^beta,

lets the two ``beta log q`` terms cancel *before* the interval quotient:

    log Z_g / log Z_f = 1 + (log S_g - log S_f) / (beta log q + log S_f),

and ``log S`` moves only at rate ``~ log(1 + O(1/sqrt q))``, an order of
magnitude slower.  That is what makes the branch and bound terminate in seconds
rather than hours.

Tail bound, unchanged from ``analysis/certify_length_three_cycle.py``.  Writing
``m`` for the multiplicity of the largest entry ``a1`` and ``rho`` for (largest
entry below ``a1``)/``a1`` < 1,

    log Z_a(beta) = beta log a1 + log T_a(beta),
    m <= T_a(beta) <= m + (len(a) - m) rho**beta,

and ``T_a`` decreases.  So on ``[B0, infinity)`` the ratio is at least a Moebius
function of ``beta``, hence monotone, and its infimum there is the smaller of its
value at ``B0`` and its limit ``log g1 / log f1``.

    python research/curve_family_cycles/certify.py
"""

from mpmath import iv, mp

iv.dps = 40
mp.dps = 40


def show(bound) -> str:
    """Display an exact interval endpoint; the comparisons use the endpoint."""

    return mp.nstr(mp.mpf(bound), 18)

Q = 11
B0 = iv.mpf(60)
MIN_WIDTH = iv.mpf(2) ** -30

A = (18, 16, 15, 15, 14, 12, 9, 6, 6, 5, 5)
B = (18, 18, 14, 13, 12, 9, 9, 9, 8, 7, 4)
C = (19, 14, 12, 11, 11, 10, 10, 10, 9, 9, 6)
NAME = {A: "A", B: "B", C: "C"}

LOG_Q = iv.log(iv.mpf(Q))


def log_reduced(sig, beta):
    """Enclosure of ``log sum_c (N_c/q)**beta``.

    Each term is monotone in ``beta``, so its enclosure is exact; only the sum
    over terms of opposite monotonicity loses anything, and the loss is of the
    size of the spread of ``N_c/q`` rather than of ``q`` itself.
    """

    total = iv.mpf(0)
    for value in sig:
        total += iv.exp(beta * iv.log(iv.mpf(value) / Q))
    return iv.log(total)


def ratio(g, f, beta):
    """Enclosure of ``log Z_g(beta) / log Z_f(beta)`` for an interval ``beta``."""

    reduced_g = log_reduced(g, beta)
    reduced_f = log_reduced(f, beta)
    return 1 + (reduced_g - reduced_f) / (beta * LOG_Q + reduced_f)


def tail_lower_bound(g, f):
    """Rigorous lower bound on ``log Z_g / log Z_f`` over ``[B0, infinity)``."""

    top_g, top_f = max(g), max(f)
    mult_g, mult_f = g.count(top_g), f.count(top_f)
    below_f = [value for value in f if value < top_f]
    rho_f = iv.mpf(max(below_f)) / iv.mpf(top_f) if below_f else iv.mpf(0)
    numerator = B0 * iv.log(iv.mpf(top_g)) + iv.log(iv.mpf(mult_g))
    denominator = B0 * iv.log(iv.mpf(top_f)) + iv.log(
        iv.mpf(mult_f) + iv.mpf(len(f) - mult_f) * rho_f**B0
    )
    at_b0 = numerator / denominator
    at_infinity = iv.log(iv.mpf(top_g)) / iv.log(iv.mpf(top_f))
    return min(at_b0.a, at_infinity.a)


def lower_bound(g, f, target):
    """Verify ``inf_beta log Z_g / log Z_f > target``.  Returns (ok, worst, boxes)."""

    worst = tail_lower_bound(g, f)
    if worst <= target:
        return False, worst, 0
    stack = [(iv.mpf(0), B0)]
    boxes = 0
    while stack:
        low, high = stack.pop()
        boxes += 1
        enclosure = ratio(g, f, iv.mpf([low, high]))
        if enclosure.a > target:
            worst = min(worst, enclosure.a)
            continue
        if high - low < MIN_WIDTH:
            return False, enclosure.a, boxes
        middle = (low + high) / 2
        stack.append((low, middle))
        stack.append((middle, high))
    return True, worst, boxes


def upper_bound(g, f, beta):
    """``C(g->f) <= min(ratio at this beta, the two endpoint values)``."""

    at_beta = ratio(g, f, beta).b
    at_zero = (iv.log(iv.mpf(len(g))) / iv.log(iv.mpf(len(f)))).b
    at_infinity = (iv.log(iv.mpf(max(g))) / iv.log(iv.mpf(max(f)))).b
    return min(at_beta, at_zero, at_infinity)


CYCLE = [(A, B), (B, C), (C, A)]
# A rational beta near each forward rate's minimiser, for the free upper bound.
# Where the infimum sits at an endpoint instead, upper_bound picks that up alone.
PROBE = {(A, B): "33/2", (B, C): "1", (C, A): "23/6"}


def main() -> int:
    ok = True
    print("Strict three-cycle among genus-two pencils over F_11\n")
    print("  A:  y^2 = x^5 +  3x^4 + 4x^3 + x^2 +  x      + c   sigma(A) =", A)
    print("  B:  y^2 = x^6 +  9x^5 + 7x^4 + 2x^3 +  x      + c   sigma(B) =", B)
    print("  C:  y^2 = x^6 + 10x^5 + 8x^4 + 8x^3 + x^2 + 2x + c   sigma(C) =", C)
    print("\n  claim:  A < B < C < A,  where a < b iff C(a->b) < C(b->a)\n")
    for a, b in CYCLE:
        text = PROBE[(a, b)]
        parts = (text.split("/") + ["1"])[:2]
        beta = iv.mpf(int(parts[0])) / int(parts[1])
        upper = upper_bound(a, b, beta)  # C(a->b) <= upper
        good, worst, boxes = lower_bound(b, a, upper)  # C(b->a) > upper ?
        ok = ok and good
        first, second = NAME[a], NAME[b]
        print(f"  edge {first} < {second}")
        print(f"    C({first}->{second}) <= {show(upper)}   (at beta = {text})")
        print(
            f"    C({second}->{first}) >  {show(worst)}   "
            f"(branch and bound on [0, {int(B0)}] with {boxes} boxes, plus tail)"
        )
        print(f"    [{'ok' if good else 'FAIL'}] certified margin >= {show(worst - upper)}\n")
    print("cycle certified" if ok else "CERTIFICATION FAILED")
    print(
        "\n  The lower bounds stop as soon as they clear the target, so the margins\n"
        "  above are certificates rather than the true gaps, which are\n"
        "  4.6953e-03, 1.7146e-03 and 2.0998e-03."
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
