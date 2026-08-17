# Session brief B — do arithmetic families of curves ever cycle?

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (branch; commit or
stash first).

**Read first:** `research/m_and_e_and_a_c/FINDINGS.md` (Notation, T2.1, T2.2) and
`paper_finite_fields_maps/docs/finite_field_exchange_matrix.md` (the
endpoint-regime theorem and the cycle sections).

## The question

The exchange comparison is `a ≺ b ⟺ C(a→b) < C(b→a)`. It is known to run in
circles for integer signatures and for quadratic and cubic homogeneous maps
`F_3³ → F_3³` (7 and 586 distinct 3-cycles, margins ~1e−2). It does **not**
cycle among the quadratic or cubic map classes over `F_q`.

**Does it cycle among families of curves?** Take many maps `f : A² → A¹` over a
fixed `F_q` — elliptic and hyperelliptic pencils `y² = P(x) + c`, quadratic and
quartic twist families, superelliptic `y^r = P(x) + c`, families of different
genus and monodromy — compute the exchange matrix on their fiber signatures, and
search for strict 3-cycles.

## Why it is the decisive experiment

There is a proved theorem: if both rates of a pair are attained at an endpoint,
then `a ≺ b ⟺ φ(a) < φ(b)` with `φ = log(#fibers)·log(max fiber)`, so the
comparison is a total preorder and **no cycle is possible**. A cycle therefore
requires a pair on which an interior tangency overturns `φ`. Every cycle found so
far has exactly one or two `φ`-violating edges — never zero, never three.

So the search is directed: compute `φ` for every family, find pairs where the
computed comparison disagrees with `φ`, and look for cycles through those pairs.

## The obstacle to expect

For curve families `C(L→f)` is always the `β = ∞` endpoint and the signatures are
close to flat, which is exactly the regime the theorem forbids cycles in. If no
cycle is found, that is itself a result and should be explained quantitatively:
how large would an interior correction have to be to overturn a typical
`φ`-gap, and how far is the actual correction (`O(1/(q log q))`) from that?

## Traps, learned the hard way

* `−½JDJ` always has the constant vector in its kernel, so its smallest
  eigenvalue is capped at 0 and gives no search gradient. Work in an orthonormal
  basis of `{Σx = 0}`.
* Grids truncated below `β ≈ 500` can hide the phenomena entirely; verify any
  headline number on a dense grid out to `β ~ 10³`, and against 30+ digit mpmath.
* Signature collisions are structural: 400 random elliptic fibrations give only
  `gcd(4,q−1)+1` distinct signatures. Use genus ≥ 2 for variety.
* `exchange_rate` is accurate to ~1e−13; treat differences below 1e−10 as ties.

## Success criterion

Three explicit families over an explicit `F_q` with `C` values, contact
temperatures, margins, and independent verification. Or a quantitative
explanation of why the regime forbids it.
