# Session brief G — what can an exchange geometry be?

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (branch; commit or
stash first).

**Read first:** `research/m_and_e_and_a_c/FINDINGS.md` (T1.1, T1.2, T1.3 — all
about the symmetric part), `research/m_and_e_and_a_c/PLAN.md` (Conventions),
`research/session_briefs/D_quantum_mechanics_of_the_exchange_matrix.md` Part 0.

## The question nobody has asked

Every result so far is of the form "the exchange geometry has property P" —
it is a metric, it is not of negative type, its minimal violation has five
points, its `l2`-distortion is `1.3375`, its comparison has cycles. **None of
them says which geometries are possible.** That is the structural question, and
it decides how much of the programme is about arithmetic and how much is about
the framework.

Set up. Every signature gives a function on the `β`-line,

```
Z_a(β) = Σ_i a_i^β,     u_a(β) = log log Z_a(β),     β ∈ [0,∞]
```

and both halves of the exchange matrix are functionals of differences:

```
d(a,b) = osc_β(u_a − u_b)          the metric        (symmetric part, ×2)
A(a,b) = mid_β(u_a − u_b)          the comparison    (antisymmetric part)
```

So the exchange space is a subset of `C[0,∞]/constants` with the oscillation
seminorm — a **sup-norm geometry**. Every finite metric space embeds
isometrically in `ℓ∞` (Fréchet), so if `u` were arbitrary, *every* finite metric
would be realisable and the geometry would carry no information at all.

**But `u` is not arbitrary.** Writing `ν_a = Σ_i δ_{log a_i}` for the measure on
log-fiber-sizes,

```
log Z_a(β) = K_{ν_a}(β) + log r,      r = #atoms
```

is exactly a cumulant generating function of a positive measure, shifted. So
`u_a = log(K_ν(β) + log r)`: log of a *convex, increasing, CGF-shaped* function.
That is a strong constraint, and it is presumably why the observed
`l2`-distortion is `≈ 1.1` instead of the `O(log n)` a general metric may need.

## The three questions

### G1 — which tournaments are realisable?

`a ≺ b ⟺ A(a,b) > 0`. The comparison on `n` signatures is a tournament.

**Seed measurement, already done, reproduce it first.** 1500 random 4-subsets of
a pool of 298 random integer signatures (2–7 entries, values 1–40) give the score
sequences

```
(0,1,2,3)  transitive   1499
(0,2,2,2)  one 3-cycle     1
(1,1,1,3)                  0
(1,1,2,2)                  0
```

(`research/realizability/tournament_seed.py`, `random.seed(11)`.)

So two of the four isomorphism classes of 4-tournament were not seen at all in a
random pool, and cycles occur at rate `≈ 7·10⁻⁴`. Yet brief B found **132
three-cycles in an exhaustive 296-signature arithmetic pool at `F_11`** and 1475
at `F_13`, and the integer-signature literature in this repo has 7 and 586
cycles among `F_3` tensor classes. **The rate is not the point; realisability
is.** Determine:

* Is every tournament on `n` vertices realisable by `n` signatures, for every
  `n`? (Try `n = 4` exhaustively by targeted search, then `n = 5`.) A single
  non-realisable tournament is a theorem and a genuine obstruction.
* If everything is realisable, then **the comparison structure carries no
  information and all the arithmetic lives in the map `family ↦ signature`.**
  That is a decisive, publishable answer to brief A's one-way-flow objection —
  in the negative, and it should be written that way.
* The natural intermediate: which tournaments are realisable with *bounded*
  signature size `r`, or by signatures with a fixed number of atoms? A dimension
  count is available — `r` atoms is `r` real parameters, `n` signatures give
  `nr` parameters against `n(n−1)/2` sign conditions — so expect a threshold
  `r ≳ (n−1)/2`. **Test it.**

### G2 — which metrics are realisable?

Same question for `d`. Known: `d` is a metric, is not always of negative type
(the five-point certificate), fails the pentagonal inequality, so exchange
metrics leave `HYP` at the first opportunity. Unknown: whether they fill `MET`.

* Is every finite metric realisable up to a scale factor? Given the `ℓ∞`
  picture the answer is plausibly yes with enough atoms — **but the CGF
  constraint might block it, and that is exactly what the observed
  `c2 ≈ 1.1 ≪ O(log n)` hints at.** Settle which.
* Concrete test: the metrics with maximal `l2`-distortion are expander-like.
  Try to realise `K_{m,n}` (whose `t*` is known exactly, `½ log((m−1)(n−1))`,
  T1.2), the 5-cycle, the Petersen graph metric, and a small expander. Report
  the best achievable distortion to each.
* If some metric is *not* realisable, characterise the obstruction. The first
  place to look: `u_a − u_b` is a difference of logs of CGFs, so it has bounded
  variation and at most a controlled number of sign changes in `β` — which caps
  how many "independent directions" the sup-norm can use. Count the sign changes
  empirically first (**this is cheap and might immediately give the theorem**).

### G3 — how much curl is possible?

Brief D showed a strict 3-cycle is exactly `|curl A| = Σ|A|` on the triangle.
Random integer pools give a global Hodge split of `‖grad‖/‖A‖ ≈ 0.9987–0.9992`,
residual `0.041–0.051` at `n = 8, 16, 24` with **zero** 3-cycles in 2640
triangles (`research/realizability/tournament_seed.py`; reproduce first — the
residual is mildly pool-dependent, an earlier ad-hoc pool gave `0.065–0.088`, so
quote the script's output and not this sentence). Questions:

* Is `‖curl‖/‖A‖` bounded away from 1 for exchange geometries, and if so by
  what? A bound would be a rigidity theorem: **the exchange comparison is
  necessarily almost-scalar**, which is the strongest possible form of brief A's
  objection and worth knowing.
* Does the bound depend on `n`, on the number of atoms, or on the *spread* of
  the signatures? Hill-climb `‖curl‖/‖A‖` and report the maximum found together
  with the signatures attaining it.
* Cross-check against the `no-arbitrage` inequality of brief D(d): `Σ S ≥ |Σ A|`
  around every loop. That already bounds the flux by the metric. Is it tight?
  On the known 3-cycle it is `0.201398 ≥ 0.015814` — a factor 13 of slack, so
  probably not tight, and the true bound is the thing to find.

## Why this matters more than it looks

The programme has two possible self-assessments and this brief chooses between
them.

* **If exchange geometries are universal** (every tournament, every metric),
  then the framework is a rich but structureless container; every theorem proved
  about `d` or `A` is a theorem about the *particular* signatures arithmetic
  supplies, and the honest paper is about the map from arithmetic to signatures.
* **If they are rigid** (some tournament or metric is unreachable, or `‖curl‖`
  is bounded), then the framework itself has content, the observed
  `c2 ≈ 1.1`, the pentagonal failure and the `≈ 8%` curl are all instances of
  one rigidity theorem, and *that* theorem is the paper.

Either answer is publishable and the current papers cannot be finished without
knowing which holds.

## Traps

* `−½JDJ` has the constant vector in its kernel — no search gradient. Work in an
  orthonormal basis of `{Σx = 0}` (T1.1). An earlier session reported "no 5-point
  violation" purely because of this.
* Any windowed `β` computation needs `β ~ 10³`; `[0.001, 200]` destroyed the
  13-signature violation with truncation error `3.5e−3` against a violation of
  `5.45e−4` (T1.2).
* `exchange_rate` is good to `~1e−13`; treat differences below `1e−10` as ties.
  A realisability search that produces margins at `1e−12` has produced nothing —
  demand `1e−6` or better and verify with 40-digit mpmath, as T1.1 did.
* Signatures are non-empty collections of positive integers, and `(1,)` is the
  zero-complexity identity with infinite rate from every target — exclude it from
  every pool or it will dominate every search.
* Realisability searches over *real-valued* atoms are legitimate for the
  structural question and much easier than integer ones; if you use them, say so,
  and check whether the witness can be rounded to integers.

## Success criterion

For each of G1, G2, G3: either a construction showing realisability with an
explicit recipe, or an obstruction with a proof. "Searched and did not find" is
acceptable only with the search space stated and a dimension count explaining
why the search was adequate.

## Reproduce / build on

`research/realizability/tournament_seed.py` (the seed numbers above),
`analysis/exchange_positivity.py`, `analysis/negative_type_certificate.py`,
`analysis/l2_distortion.py`, `research/m_and_e_and_a_c/gauge_decomposition.py`,
`research/m_and_e_and_a_c/t1_2_part1_psd_in_t.py`.
