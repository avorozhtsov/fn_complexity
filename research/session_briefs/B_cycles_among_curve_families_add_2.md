# Addendum 2 to brief B — the cycle exists

**Status: the question of brief B is answered, affirmatively and with a
certificate.** Arithmetic families of curves do run in circles.

Read with `B_cycles_among_curve_families.md` and its addendum 1. Reproduce with
`analysis/certify_curve_family_cycle.py`.

---

## The certificate

Three hyperelliptic pencils `f_i : A² → A¹`, `f_i(x,y) = y² − P_i(x)`, over
`F_101`, whose fibers over `c ∈ F_101` are the genus-2 curves `y² = P_i(x) + c`:

```
P_1 = x⁵ + 70x⁴ + 28x³ + 15x² + 11x + 31
P_2 = x⁵ + 42x⁴ + 32x³ + 74x² + 96x + 60
P_3 = x⁵ + 72x⁴ + 21x³ +  2x² +  6x + 57
```

The comparison `a ≺ b ⟺ C(a→b) < C(b→a)` runs

```
f_1 ≻ f_2 ≻ f_3 ≻ f_1
```

| edge | `A = ½ log(C(b→a)/C(a→b))` | argmin `β`, both directions |
|---|---:|---|
| `f_1 → f_2` | `−1.619719·10⁻⁴` | 20.99 / 93.17 |
| `f_2 → f_3` | `−1.301503·10⁻⁴` | 38.10 / ∞ |
| `f_3 → f_1` | `−1.164720·10⁻⁴` | ∞ / 24.78 |

`|curl A| / Σ|A| = 1.000000000000`, which for an antisymmetric edge function is
equivalent to all three edges agreeing in sign, i.e. to a strict cycle.

**Verification, three independent ways.** The package solver, the grid solver of
`research/m_and_e_and_a_c/t2_2_common.py`, and a 40-digit `mpmath` scan on
`[0, 2000]` agree in sign on all three edges, with spreads `1.6·10⁻⁹`,
`3.6·10⁻⁹`, `8.0·10⁻⁸`. **The smallest margin is `1.165·10⁻⁴`, six orders above
the `10⁻¹⁰` tie floor.** Signatures are recomputed by point count inside the
script, so nothing depends on the search RNG; each totals `q² = 10201`.

No arbitrage holds, as it must: the cycle products are `0.996192` forward and
`0.995378` reverse, both `< 1`. **The comparison cycles; the market does not.**

## Why it was allowed to exist

The endpoint theorem forbids cycles when both rates of a pair are attained at an
endpoint. Here **four of the six rates are attained at interior `β`** — 20.99,
93.17, 38.10, 24.78 — and two at `β = ∞`. The hypothesis fails on exactly the
edges it has to.

Addendum 1's refinement `φ̃ = M − ((3−2√2)/2)·m₂` is also violated, and in the
precise way it predicted. Addendum 1 argued a cycle would need exact `(M, m₂)`
degeneracy, since `m₂ = ν(P)/q − 1` makes equality an integer condition on
`ν(P)`. That is what happened, and it was not designed for:

| | max fiber | `m₂` | `ν(P)` | `φ` | `φ̃` |
|---|---:|---:|---:|---:|---:|
| `f_1` | 123 | 0.851485 | 187 | 22.208810749 | 21.926954 |
| `f_2` | 123 | 0.851485 | 187 | 22.208810749 | 21.926954 |
| `f_3` | 122 | 0.990099 | 201 | 22.171136087 | 20.915063 |

So `f_1` and `f_2` are **exactly tied** by `φ`, by `m₂` and hence by `φ̃`, and
their comparison is decided by the interior alone. Of the three edges: one
(`f_2 ≻ f_3`) agrees with `φ̃`, one (`f_3 ≻ f_1`) contradicts it, one is a `φ̃`
tie. That is the "one or two `φ`-violating edges — never zero, never three"
pattern of the parent brief, now with the tie case exhibited.

## How it was found: `curl A` instead of `φ`-disagreements

The parent brief says to search for pairs where the comparison disagrees with a
scalar; addendum 1 shows the scalar has to be refined and the refinement is
*still* a total order. That is a ladder with no top. The gauge decomposition of
`research/m_and_e_and_a_c/gauge_decomposition.py` ends it. Write

```
L(a,b) = −log C(a→b) = S(a,b) + A(a,b),    S = (L+Lᵀ)/2 = d/2,   A = (L−Lᵀ)/2
```

so `a ≺ b ⟺ A(a,b) > 0`. For an antisymmetric function on a complete graph the
triangle sums `curl A = A(a,b) + A(b,c) + A(c,a)` all vanish iff `A = dψ` for a
potential `ψ`. Hence:

* **every scalar invariant, at every order of any expansion, contributes exactly
  zero to `curl A`** — `φ`, `φ̃`, and whatever comes next are annihilated without
  being computed;
* a strict cycle forces `|curl A| = Σ|A|`, so `r = |curl A|/Σ|A| ∈ [0,1]` equals
  `1` exactly on cycles: a smooth objective where sign-hunting gives none.

Scanning `r` over pools of 90 hyperelliptic signatures found the cycle
immediately. `research/m_and_e_and_a_c/curl_on_curve_families.py`:

| `q` | genus | `n` | triangles | max `r` | max `|curl A|` | median `|A|` | margin |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 101 | 2 | 90 | 117480 | **1.000000** | 3.051·10⁻³ | 2.835·10⁻³ | 1.165·10⁻⁴ |
| 211 | 2 | 90 | 117480 | **1.000000** | 1.811·10⁻³ | 1.527·10⁻³ | 7.499·10⁻⁵ |

At `q = 101` there are **9 strict 3-cycles among the 117480 triangles**, with
margins from `1.165·10⁻⁴` down. Cycles are not rare and not marginal.

## What this does and does not establish

**Does.** The success criterion of the parent brief is met: three explicit
families over an explicit `F_q`, with rates, contact temperatures, margins and
independent verification. The endpoint regime does not swallow curve families,
and the `√q` dominance argument of addendum 1 is defeated exactly where that
addendum said it could be — at `(M, m₂)` degeneracies, at small `q`.

**Does not.** This is a cycle among **signatures**. The signature merges families
the geometry separates, so it is not yet a cycle among the underlying pencils in
any stronger sense. That distinction is the crux of brief C and must be settled
before any claim of the form "no scalar invariant of arithmetic families exists".
Specifically, open:

1. Are `f_1, f_2, f_3` inequivalent as fibrations, or do two of them differ by a
   coordinate change that the signature cannot see? `f_1` and `f_2` share `M`,
   `m₂` and `ν(P)`; check whether they share more.
2. Does the cycle survive base change — is there a cycle at `q = 211, 503, 1009`
   among *the same three pencils reduced mod different primes*, or only among
   families chosen afresh at each `q`? The scan answers the second question
   affirmatively; the first is untouched and is the arithmetically meaningful one.
3. The `q = 211` cycle should be certified to the same standard as the `q = 101`
   one before both are quoted.

## Next

* Certify the `q = 211` cycle (same script, new constants).
* Run the scan at `q = 503, 1009` — the margin should fall like `1/(q log q)`,
  and the point at which it crosses `10⁻¹⁰` is where the phenomenon becomes
  unobservable rather than absent. That crossover is worth reporting.
* Hand items 1–3 above to brief C, which now has content.
