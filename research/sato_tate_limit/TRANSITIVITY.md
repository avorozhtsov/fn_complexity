# Findings — is same-genus transitivity a theorem?

Answer to session brief K. **No. It is a theorem up to genus 6 and false from
genus 13 — and the counterexample needs no repeated isogeny factor, which
removes the arithmetic bottleneck `FINDINGS.md` reduced everything to.**

## Summary

* **Route 1 is dead.** `sign mid(Ψ_μ − Ψ_ν)` is **not a function of
  `(Δm₂, Δt)` at all**, so no trade-off `Δm₂ − c·Δt` — linear or otherwise —
  exists. Two `(Δm₂, Δt)` values of the symplectic library carry both signs of
  `mid`, and the linear rule is infeasible by an exact certificate: the library
  forces `κ > 3/2` and `κ < 4/9` at once. The suspicious `0.668` is the grid
  point just above `2/3`, the decision boundary of one family of pairs with
  `(Δm₂, Δt) = (−1, −3/2)`; `2/3 = 1/(3/2)` is that family's `Δm₂/Δt` and
  nothing more.
* **Route 2 is dead too, and now measurably so.** The level lemma is proved
  outright below (Theorems 4 and 5: `≤ 2` levels give `0`, `≥ 3` give exactly
  `1/4`), but its hypothesis is badly false on the library: on all 119 triples
  that survive the crossing reduction the three interior minima sit at **three
  different `τ`**, spread by a median factor of `57`, and the three maxima by a
  median factor of `4.5·10⁴`. The library uses six levels, not two. Whatever
  makes genus `≤ 12` transitive, it is not a level count.
* **Route 3 succeeded — twice.**
  * **Genus 13, multiplicity-free products.** The first same-genus 3-cycle among
    honest products of `USp(2g)` blocks:

    ```
    USp14×USp4²×SU2²   ≺   USp10×USp8×USp4²   ≺   USp12×USp8×SU2³   ≺   USp14×USp4²×SU2²
       (7,2,2,1,1)              (5,4,2,2)              (6,4,1,1,1)
    ```

    all of genus 13, `α_max = 26`, midranges (polished, independent grid)
    `−0.0133900010, −0.0079979243, −0.0344360664`, margin `7.998·10⁻³`, all six
    contacts interior at `τ ∈ [1.6, 21]`. **No vertex has a repeated isogeny
    factor.** Genus 13 is the first genus at which a cycle exists; genus 14 has
    3 and genus 15 has 6, and the search over partitions is exhaustive at every
    genus to 15.
  * **Genus 4, mixtures.** Allowing the convex combinations that a
    two-component monodromy group would give, a cycle appears already at
    genus 4 (`α_max = 8`), margin `8.674·10⁻³`, one of whose vertices is the
    pure product `USp4 × SU2 × SU2`. So the brief's statement — *for measures
    with the same `α_max` the midrange comparison is a total order* — is false,
    and false at a genus the existing libraries already reach.
* **What is a theorem.** Because `Ψ` is additive over independent factors, a
  multiplicity-free symplectic measure is a partition `λ` of the genus and
  `Ψ_λ = (Σ_i κ_{λ_i})/τ`, `κ_g = K_{USp(2g)}`. If `g ↦ κ_g(τ)` is concave then
  `λ ⊵ μ` in the **dominance order** forces `Ψ_λ ≤ Ψ_μ` *pointwise*. Dominance
  is a total order on partitions of `G` exactly when `G ≤ 5`, and the two
  incomparable pairs at `G = 6` are settled by convexity of the increments.

  > **Theorem.** For genus `G ≤ 6` the multiplicity-free symplectic measures of
  > genus `G` are totally ordered by pointwise domination of `Ψ`; not one pair
  > of them crosses.

  Concavity of `κ` in `g` is exactly `4 b_g(τ)² ≤ 1` for the Jacobi coefficients
  of the tilted Chebyshev weight `e^{2τx}√(1−x²)` (free value `b_g(0) = 1/2`).
  Verified over eleven decades of `τ` for all `g ≤ 14`, and at 40 digits; both
  ends proved. That one inequality is the session's only analytic gap.
* **Where the proof breaks is exactly where the cycles start to become
  possible.** `(4,3)` and `(5,1,1)` at genus 7 are the first
  dominance-incomparable pair with `(m₂, t)` comonotone, so they must cross —
  and the multiplicity-free closest approach then falls
  `0.170 → 0.0790 → 0.0505 → 0.0022 → −0.0080` from genus 7 to genus 13.
* **The structure that does most of the work.** Inside a genus, **718 of the 765
  pairs of the symplectic library never cross**, and no pair crosses more than
  once. Pointwise domination is transitive, so at most one edge of a 3-cycle can
  be a non-crossing edge; the 6840 same-genus triples of that library collapse
  to **119**, all transitive, tightest certificate `+0.08050704`.
* **Two theorems of independent interest.** (i) The comparison is
  **Legendre-invariant**: `sup`, `inf`, `mid` and `osc` of `Ψ_μ − Ψ_ν` equal
  those of `J_μ − J_ν`, `J = I^{-1}` the inverse large-deviation rate function —
  a second, independent way to compute everything, in which the `Ψ(∞) = α_max`
  trap cannot occur. (ii) The level lemma, with its exact constants.

Four earlier claims are corrected in [Corrections](#corrections).

---

## Notation

Conventions of `FINDINGS.md`. For a limit measure `μ`,

```
K_μ(τ) = log E_μ[e^{τα}],   Ψ_μ(τ) = K_μ(τ)/τ,   Ψ_μ(0) = 0,  Ψ_μ(∞) = α_max
μ ≺ ν  ⟺  mid_τ(Ψ_μ − Ψ_ν) < 0,   mid = ½(sup + inf) over [0, ∞]
```

`m₂ = E[α²]`, `t` the edge exponent, `κ_g := K_{USp(2g)}` with `κ_0 = 0`.
A multiplicity-free symplectic measure is a partition `λ` of the genus `G`:
`α_max = 2G`, `m₂ = ℓ(λ)`, `t = Σλ_i² + G/2`, and

```
Ψ_λ(τ) = ( Σ_i κ_{λ_i}(τ) ) / τ .
```

Two consequences used throughout. **`Ψ_μ − Ψ_ν` is unchanged when both `μ` and
`ν` are multiplied by a common independent factor**, so the comparison depends
only on the reduced difference — which is why the same `mid` values (`0.04026`,
`0.23322`, …) recur verbatim across genera in `FINDINGS.md`'s exception table,
and why cycles propagate upward in the genus by appending a common block. And
the differences form a lattice: `Ψ_μ − Ψ_ν = Σ_B c_B Ψ_B` with `c` integral and
`Σ_B c_B α_max(B) = 0` inside a genus.

All `τ` grids below run to at least `10⁵`, all endpoints are analytic
(`Ψ(0) = 0`, `Ψ(∞) = α_max`), and every measure used is validated against its
cumulants (`κ₂ = 1` for all `USp(2g)`; `κ₄ = −1` at `g = 1` and `0` for
`g ≥ 2`) before it enters a comparison.

---

## 1. Route 1: there is no `(Δm₂, Δt)` rule *(computed, exact certificate)*

`transitivity_pairs.py` tabulates all **765** same-genus pairs of the
71-measure symplectic library (`α_max ≤ 12`) with their `(Δm₂, Δt)`.

* The 765 pairs realise **387 distinct `(Δm₂, Δt)` values**.
* **Two of those values carry both signs of `mid`:**

| `(Δm₂, Δt)` | pair | `mid` |
|---|---|---:|
| `(−2, −9/2)` | `USp4 × USp4 × USp4` vs `2·SU2 × USp8` | `−0.06388822` |
| | `SU2⁴ × 2·SU2` vs `3·SU2 × USp6` | `+0.29139540` |
| `(−2, −3)` | `2·USp6` vs `SU2 × 2·SU2 × USp6` | `−0.11002471` |
| | `3·USp4` vs `SU2 × USp4 × 3·SU2` | `+0.26598888` |

  So `sign mid` is not a function of `(Δm₂, Δt)`, let alone a monotone or linear
  one. **Route 1 as posed has no answer.**
* The linear rule `sign mid = sign(Δm₂ − κ·Δt)` is infeasible with an exact
  certificate: over the 735 pairs with `Δm₂ ≠ 0` the constraints force
  `κ > 3/2` (from `(Δm₂, Δt) = (−3, −2)`: `SU2 × 2·SU2 × 3·SU2` vs
  `USp4 × 4·SU2`) **and** `κ < 4/9` (from the `(−2, −9/2)` pair above).
* **Where `0.668` comes from.** The best single threshold is `κ = 2/3`, correct
  on 731 of 735 — and `2/3` is `Δm₂/Δt` for the four pairs with
  `(Δm₂, Δt) = (−1, −3/2)` (`2·USp4` vs `2·SU2 × USp4` and its common-factor
  translates), the exceptions the plain `m₂` rule gets wrong. `lex_exceptions.py`
  scans `κ` on a grid of spacing `0.0024` above `0.2` and returns the grid point
  `0.668`. **The constant is a fitted decision boundary set by one family of
  pairs; `2/3` has no derivation behind it and the full system has no solution.**

---

## 2. The comparison is Legendre-invariant *(proved; verified)*

Let `I_μ = K_μ*` be the large-deviation rate function and
`J_μ = (I_μ|_{[0,α_max)})^{-1} : [0,∞) → [0, α_max)`, concave and increasing,
`J_μ(0) = 0`, `J_μ(∞) = α_max`.

**Lemma.** With `u = 1/τ`,  `Ψ_μ(1/u) = sup_{v ≥ 0} ( J_μ(v) − u v )`.

*Proof.* `K_μ(τ) = sup_x(τx − I_μ(x))`, so `Ψ_μ = sup_x(x − u I_μ(x))`; the
supremum is attained at some `x ≥ 0` since the value at `x = 0` is `0` and the
expression is negative for `x < 0`. Substitute `v = I_μ(x)`, `x = J_μ(v)`. ∎

**Theorem 1.** For any two compactly supported mean-zero `μ, ν`,

```
sup_τ (Ψ_μ − Ψ_ν) = sup_v (J_μ − J_ν),      inf_τ (Ψ_μ − Ψ_ν) = inf_v (J_μ − J_ν),
```

so `mid` and `osc` may be computed on either side.

*Proof.* Fix `u > 0` and let `v*` attain `Ψ_μ(1/u)`; then
`Ψ_ν(1/u) ≥ J_ν(v*) − u v*`, so `Ψ_μ − Ψ_ν ≤ J_μ(v*) − J_ν(v*) ≤ sup_v(J_μ−J_ν)`.
Conversely fix `v₀ > 0` and take `u` in the superdifferential of `J_ν` at `v₀`;
then `v₀` maximises `J_ν(v) − uv`, so `Ψ_ν(1/u) = J_ν(v₀) − u v₀` while
`Ψ_μ(1/u) ≥ J_μ(v₀) − u v₀`, giving `(Ψ_μ − Ψ_ν)(1/u) ≥ J_μ(v₀) − J_ν(v₀)`.
Apply the result to `(ν, μ)` for the `inf`. ∎

Both endpoints are automatic on the `J` side: `δ(0) = 0`,
`δ(∞) = α_max(μ) − α_max(ν)`. **A grid in `v` cannot underestimate `Ψ(∞)`,
because `α_max` is where `J` ends** — the trap the brief warns about does not
exist in this coordinate. The critical points translate too: `δ = J_μ − J_ν` has
`δ'(v) = 0` exactly where `J_μ'(v) = J_ν'(v)`, so the interior extrema of
`Ψ_μ − Ψ_ν` are the **parallel tangencies of two concave curves**. The two
asymptotics dualise to

```
J(v) ≈ √(2 m₂ v)   (v → 0),        α_max − J(v) ≈ C e^{−v/t}   (v → ∞).
```

*Verified* (`transitivity_duality.py`) on six pairs including all three edges of
the genus-11 near miss, with `J` built parametrically from `x = K'(τ)`,
`v = τK'(τ) − K(τ)` at 40 digits. The residual between the two sides falls
`4.687·10⁻⁴ → 5.499·10⁻⁵` when the parametric sample is refined `220 → 880`
points: it is interpolation error, and the identity holds.

---

## 3. Inside a genus, the comparison is mostly pointwise domination

### 3.1 At most one sign change *(computed; mechanism proved)*

Sign changes of `Ψ_μ − Ψ_ν` are sign changes of
`M_μ − M_ν = ∫ e^{τα} d(μ−ν)`. The kernel `e^{τα}` is totally positive, so by
Schoenberg's variation-diminishing property the transform has at most as many
sign changes as `μ − ν` does in `α`; both measures being symmetric the transform
is even, so **the count on `(0, ∞)` is at most half** the count on `ℝ`. If the
densities cross once on `(0, α_max)`, `Ψ_μ − Ψ_ν` crosses at most once.

*Computed:* over the **765** same-genus pairs of the symplectic library the
number of sign changes on a 1201-point grid spanning eleven decades is **0 on
718 pairs and 1 on 47 — never 2 or more**; and over the **35683** same-genus
pairs of the multiplicity-free cone at every genus from 2 to 15, **0 on 32535
and 1 on 3148 — again never 2 or more**. That is 36448 pairs and not one
double crossing.

### 3.2 The crossing criterion, made exact *(proved, given 3.1)*

`sign(Ψ_μ − Ψ_ν) = sign(Δm₂)` near `τ = 0` and `= −sign(Δt)` near `τ = ∞`, so
with `Δm₂` and `Δt` both non-zero the number of sign changes is odd iff
`Δm₂·Δt > 0`. With at most one sign change available:

> **`Δm₂·Δt > 0` ⟹ exactly one crossing. `Δm₂·Δt < 0` ⟹ no crossing, and then
> `Ψ_μ − Ψ_ν` has constant sign: the pair is ordered by pointwise domination,
> and the `mid` comparison is decided without any further computation.**

*Computed:* the criterion has **zero genuine failures** — on the symplectic
library, 765 of 765 with no vacuous case at all; on the multiplicity-free cone
to genus 15, 35496 of 35683, the 187 remaining being exactly the cases where the
criterion says nothing (185 with `Δm₂ = 0`, 2 with `Δt = 0`). Note that an
`m₂`-tie *can* cross — 185 of them do — so the tie cases genuinely need the
next cumulant, and `FINDINGS.md`'s "`t` is the decider on a matched pair but not
an infallible one" is the right description of them.

### 3.3 The reduction *(proved + computed)*

**Proposition (P1).** Pointwise domination of `Ψ` is a partial order and fixes
the sign of `mid`. Hence **at most one edge of a 3-cycle can be a non-crossing
edge**: two non-crossing edges of a triangle are consecutive, and chaining
`Ψ_1 ≤ Ψ_2 ≤ Ψ_3` gives `Ψ_1 ≤ Ψ_3`, contradicting the third edge.

*Computed* (`transitivity_certificate.py`):

| genus | measures | same-genus triples | with `≥ 2` crossing edges |
|---:|---:|---:|---:|
| 2 | 3 | 1 | 0 |
| 3 | 5 | 10 | 0 |
| 4 | 11 | 165 | 1 |
| 5 | 17 | 680 | 3 |
| 6 | 34 | 5984 | 115 |
| **total** | **71** | **6840** | **119** |

106 of the 119 have exactly two crossing edges, 13 have three. **All 119 are
transitive**, tightest certificate

```
SU2² × 2·USp4   /   SU2² × 2·SU2 × USp4   /   3·SU2 × USp6        genus 6
max_i mid_i = +0.08050704
```

which recovers `FINDINGS.md`'s closest same-genus approach `−8.05·10⁻²`. The
exhaustive search over 6840 triples is thereby replaced by a proved structural
step plus a check on 119.

---

## 4. The multiplicity-free cone is the dominance order *(proved, modulo one inequality)*

### 4.1 The reduction to a one-variable inequality

**Theorem 2.** Suppose `g ↦ κ_g(τ)` is concave for every `τ > 0`. If `λ ⊵ μ` in
the dominance order on partitions of `G`, then `Ψ_λ ≤ Ψ_μ` pointwise; hence
`λ ≺ μ`, and the pair does not cross.

*Proof.* By the Muirhead/Brylawski transfer lemma `μ` is obtained from `λ` by
finitely many moves `(a, b) ↦ (a−1, b+1)` with `a > b` (parts `0` allowed, which
is how a new part appears). Concavity means the increments
`Δ_g = κ_g − κ_{g−1}` are non-increasing, so `a ≥ b+1` gives `Δ_a ≤ Δ_{b+1}`,
i.e. `κ_a + κ_b ≤ κ_{a−1} + κ_{b+1}`: each move does not decrease `Σ_i κ_{λ_i}`.
Hence `K_λ ≤ K_μ` and `Ψ_λ ≤ Ψ_μ` for all `τ`. ∎

**Concavity, restated.** `M_g = H_g(τ)/H_g(0)` with `H_g` the `g × g` Hankel
determinant of the tilted Chebyshev weight `w_τ(x) = e^{2τx}√(1−x²)`, and
`H_{g+1}H_{g−1}/H_g² = b_g²`, the squared off-diagonal Jacobi coefficient. Since
`b_g(0) = 1/2` exactly (Chebyshev `U`),

```
κ_{g+1} + κ_{g−1} − 2κ_g  =  log( 4 b_g(τ)² ) ,
```

so **`κ` is concave in `g` iff `b_g(τ) ≤ 1/2`**: the exponential tilt only lowers
the Jacobi coefficients below their free value. The tilt is the Toda flow, under
which `d/dτ log b_g² = 2(a_g − a_{g−1})` in the standard indexing, so — since
`b_g(0) = 1/2` — the inequality is equivalent to the *diagonal* Jacobi
coefficients being non-increasing in `g` all along the flow. That is very
plausible: the tilt pulls mass to `x = 1`, so the low-degree orthogonal
polynomials are centred furthest to the right, and `a_g → 0`, `b_g → 1/2` as
`g → ∞` by Rakhmanov's theorem. Plausible is not proved.

*Computed:* `max_τ 4 b_g(τ)² = 1`, attained only as `τ → 0`, for every `g ≤ 14`
over `τ ∈ [10⁻⁴, 10⁵]`; at 40 digits,

```
4 b_1(1)²   − 1 = −0.3491620314502367171856839336293172536365
4 b_3(1)²   − 1 = −0.001180306696243948850390380547425789174727
4 b_6(1)²   − 1 = −1.921965684292005536713148133107900063328e−9
4 b_9(10)²  − 1 = −0.237865576971722590624129258225866822525
```

and the margin `1 − 4b_g(τ)²`, which must be positive for every `τ > 0`:

| `g` | `τ = 0.1` | `τ = 1` | `τ = 5` | `τ = 30` |
|---:|---:|---:|---:|---:|
| 1 | `4.979·10⁻³` | `0.34916` | `0.94355` | `0.99835` |
| 3 | `1.387·10⁻⁹` | `1.180·10⁻³` | `0.62098` | `0.98844` |
| 6 | `2.086·10⁻²¹` | `1.922·10⁻⁹` | `0.05305` | `0.95709` |
| 8 | `4.777·10⁻³⁰` | `4.491·10⁻¹⁴` | `1.415·10⁻³` | `0.92524` |
| 11 | `8.893·10⁻⁴⁴` | `8.502·10⁻²²` | `6.622·10⁻⁷` | `0.86108` |

Both ends are proved: `b_g(0) = 1/2` exactly, and as `τ → ∞` the tilted weight
becomes the Laguerre weight `e^{−u}√u` at scale `1/(2τ)`, giving
`b_g² → g(g+½)/(2τ)² → 0`. The margin shrinks like `τ^{2g}` as `τ → 0` — which
is why it must be computed at high precision and not read off a plot. **This is
the one analytic gap of the session.**

*Also computed:* the increments `Δ_g` are convex in `g`, and the third
differences have the completely-monotone sign as well, but the pattern stops
there: `(−1)^j Δ^j Δ_g ≥ 0` fails from `j = 4`. `Δ` is not completely monotone,
so no Hausdorff representation is available.

### 4.2 The theorem, and where it stops

Dominance is a **total** order on partitions of `G` iff `G ≤ 5`. At `G = 6` there
are exactly two incomparable pairs, and they are exactly the two `t`-ties:

| genus 6 incomparable pair | `Δm₂` | `Δt` | what settles it |
|---|---:|---:|---|
| `(4,1,1)` vs `(3,3)` | `+1` | `0` | `Δ₂ + Δ₃ ≤ Δ₁ + Δ₄` |
| `(3,1,1,1)` vs `(2,2,2)` | `+1` | `0` | `2Δ₂ ≤ Δ₁ + Δ₃` |

Both follow from **convexity of the increment sequence** `Δ_g` (`{1,4} ⊵ {2,3}`
and `{1,3} ⊵ {2,2}` as multisets, and a convex sequence has the larger sum on
the more spread-out multiset) — the order-2 member of the same computed family
as concavity. Hence:

> **Theorem 3.** For genus `G ≤ 6` the multiplicity-free symplectic measures of
> genus `G` are totally ordered by *pointwise domination of `Ψ`*, hence by
> `mid`, and no pair among them crosses. The order is the dominance order on
> partitions, completed at `G = 6` by `(4,1,1) ≻ (3,3)` and
> `(3,1,1,1) ≻ (2,2,2)`.

*Computed confirmation:* 0 crossing pairs at genus `≤ 6`, and the dominance
prediction holds on every dominance-comparable pair up to genus 15 without
exception (12817 of the 15400 pairs at genus 15 alone).

**The proof breaks at genus 7, and it is easy to say where.** The minimum of
`Σλ_i²` over `ℓ`-part partitions of `G` first drops below the maximum over
`(ℓ+1)`-part partitions at `G = 7`: `(4,3)` has `Σλ² = 25` with `ℓ = 2`, while
`(5,1,1)` has `27` with `ℓ = 3`. That makes `(Δm₂, Δt) = (−1, −2)` comonotone,
so by §3.2 the pair **must** cross — the first same-genus crossing of the
multiplicity-free cone. Polished at 40 digits,

```
mid(Ψ_{511} − Ψ_{43}) = +0.2618589419637525451766927437570167868348
    sup +0.58958244911875033218 at τ = 2.152007156
    inf −0.065864565191245241825 at τ = 30.19342567
```

### 4.3 Genus 7 to 15: the cycles appear at 13 *(computed, exhaustive)*

`transitivity_dominance.py`, every partition of every genus, `Ψ` from exact
`κ_g` computed at working precision `60 + 1.5·g(g−1)·log₁₀(2τ)` digits:

| genus | partitions | pairs | dominance-comparable | crossing | **3-cycles** | closest approach |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 11 | 55 | 53 | 0 | **0** | `+0.2332` |
| 7 | 15 | 105 | 101 | 2 | **0** | `+0.1696` |
| 8 | 22 | 231 | 216 | 7 | **0** | `+0.0790` |
| 9 | 30 | 435 | 400 | 17 | **0** | `+0.0790` |
| 10 | 42 | 861 | 776 | 46 | **0** | `+0.0505` |
| 11 | 56 | 1540 | 1370 | 97 | **0** | `+0.0022` |
| 12 | 77 | 2926 | 2541 | 215 | **0** | `+0.0022` |
| **13** | **101** | **5050** | **4338** | **415** | **2** | **`−0.0080`** |
| 14 | 135 | 9045 | 7640 | 823 | 3 | `−0.0080` |
| 15 | 176 | 15400 | 12817 | 1526 | 6 | `−0.0080` |

The genus-11 near miss, `(6,2,1,1,1) → (4,3,3,1) → (5,2,2,2)`, is the first
triple whose three vertices are pairwise dominance-incomparable, and it is three
orders tighter than anything below it: polished, its three midranges are
`+0.0022007643065372019, −0.0217134944503198474, −0.0079979242503499643`, so it
misses by `2.2007643065372019·10⁻³`. **The search space:**
every partition of every genus `≤ 15`, i.e. every multiplicity-free symplectic
product up to genus 15 — exhaustive, not sampled, and the two cycles at genus 13
are the first that exist.

---

## 5. The genus-13 cycle *(computed, verified independently)*

```
(7,2,2,1,1) = USp14 × USp4 × USp4 × SU2 × SU2
(5,4,2,2)   = USp10 × USp8 × USp4 × USp4
(6,4,1,1,1) = USp12 × USp8 × SU2 × SU2 × SU2
```

all of genus 13, `α_max = 26`, `m₂ = 5, 4, 5` and `t = Σλ² + 13/2 =
65.5, 55.5, 61.5`. Two of the three edges are `(m₂, t)`-comonotone and cross for
that reason; the third, `(6,4,1,1,1)` against `(7,2,2,1,1)`, is an **`m₂`-tie**
(`Δm₂ = 0`, `Δt = −4`), one of the 185 crossing ties of §3.2. That is exactly
the configuration of the certified `F_101` finite-`q` witness, where two of the
three pencils share largest fibre, `m₂` and `ν(P)` and the edge is decided by
the multiplicity alone — the limit reproduces the finite-`q` mechanism.

| edge | `mid`, search grid | `mid`, independent grid + polished extrema | `sup` at `τ` | `inf` at `τ` |
|---|---:|---:|---|---|
| `(7,2,2,1,1) → (5,4,2,2)` | `−0.0134051` | `−0.013390001003294656` | `+0.58648` at `2.108` | `−0.61326` at `16.078` |
| `(5,4,2,2) → (6,4,1,1,1)` | `−0.0080000` | `−0.007997924250349964` | `+0.43931` at `13.427` | `−0.45531` at `1.5854` |
| `(6,4,1,1,1) → (7,2,2,1,1)` | `−0.0344364` | `−0.034436066427288228` | `+0.18812` at `21.048` | `−0.25699` at `3.5271` |

margin `7.9979242503499642791·10⁻³`. All six contacts are **interior**, at
`τ` between `1.6` and `21` — no endpoint contact anywhere, so nothing here
depends on how `Ψ(∞)` is supplied.

`transitivity_cycle13.py` re-derives it from scratch on a **wider and finer
grid** (`10⁻⁶…10⁸`, 3001 points, against `10⁻⁴…10⁵`, 1201) using only the ranks
the cycle actually contains (`g ≤ 8`, so the determinants are eight-fold rather
than fifteen-fold and the precision margin is far larger), polishes every
interior extremum by golden section on `log τ`, checks that `|D|` decays
monotonically at `τ = 10⁵, 10⁶, 10⁷, 10⁸` (it falls by a factor of about `8.3`
per decade, the `log τ / τ` law) so that nothing hides in the tail, and
re-validates every block against its cumulants:

```
g       1        2        3        4        5        6        7        8
κ₂   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000   1.0000     (exact value 1)
κ₄  -1.0000  -1.7e-7  -1.2e-14 -5.6e-22 -1.7e-29 -3.8e-37 -5.3e-45  9.3e-46    (exact -1, then 0)
```

A second genus-13 cycle,
`(8,1,1,1,1,1) → (4,4,4,1) → (6,3,2,2)`, has polished midranges
`−0.0070192921, −0.0004379082, −0.0335325399`, margin
`4.3790817996092885·10⁻⁴`, and is verified the same way.

**How it grows out of the genus-11 near miss.** Cancelling the common factor
`USp8` from the second edge gives the reduced pair `(5,2,2) / (6,1,1,1)`, which
is exactly the reduced form of the genus-11 near-miss edge
`(5,2,2,2) / (6,2,1,1,1)` — the two share a midrange to sixteen digits,
`−0.007997924250349964`. The genus-13 cycle is that near miss with its third
vertex replaced by one that only genus 13 makes available.

**Why this matters.** `FINDINGS.md` establishes that *every* cycle it found
needs a vertex with a repeated isogeny factor `k·G`, `k ≥ 2`, and concludes that
the `q → ∞` question "is now a single arithmetic question": whether such
families exist. **These cycles have no repeated factor at any vertex.** Repeated
*groups* are not repeated *factors*: `USp4 × USp4` is two **independent**
abelian surfaces (`α = tr g₁ + tr g₂`, `m₂ = 2`), whereas the forbidden `2·USp4`
is one surface counted twice (`α = 2 tr g`, `m₂ = 4`). Every vertex here is a
product of independent blocks, so the obstruction that `FINDINGS.md`'s reduction
ran into is simply absent; what they need instead is a curve family of genus 13
whose Jacobian splits into five (resp. four, five) independent factors of the
prescribed dimensions with full symplectic monodromy on each. That is a
construction problem of a different and much milder kind than `Jac ∼ A^k` —
cyclic covers and quotients by automorphism groups produce multi-factor
splittings routinely — but it is not carried out here, and until it is these
cycles are statements about the measure cone, exactly like the nine cycles of
`FINDINGS.md`.

**And it removes the genus constraint.** Brief F concluded that the limiting
comparison is "non-transitive **exactly** across genera", and brief J's witness
search is confined to mixed-genus triples on that basis. Both need amending:
same-genus cycles exist, they simply start at genus 13.

---

## 6. Measures with the same `α_max` are not totally ordered *(computed, verified two ways)*

A convex combination `p μ + (1−p) ν` of two measures with the same `α_max` is
again a mean-zero probability measure with that `α_max`. These are the trace
measures of two-component compact groups whenever such a group exists; they are
honest measures regardless, and they are the correct test of the statement as
posed.

*Computed* (`transitivity_mixtures.py`), inside each genus over the
multiplicity-free products and all uniform 2-mixtures of them, with cycles
counted exactly by `#cyclic = C(n,3) − Σ_i C(outdeg_i, 2)` so that nothing is
missed:

| genus | measures | triples | 3-cycles |
|---:|---:|---:|---:|
| 3 | 6 | 20 | 0 |
| 4 | 15 | 455 | 0 |
| 5 | 28 | 3276 | 0 |
| 6 | 66 | 45760 | 0 |
| 7 | 120 | 280840 | 0 |
| 8 | 253 | 2667126 | **38** |
| 9 | 355 | 7393585 | **236** |
| 10 | 367 | 8171255 | **201** |
| 11 | 381 | 9145270 | **790** |
| 12 | 402 | 10746800 | **414** |

The first uniform-mixture cycle is at genus 8 (`α_max = 16`):

```
½·USp16 + ½·(USp8×SU2⁴)   ≺   ½·USp16 + ½·USp4⁴   ≺   ½·(USp10×SU2³) + ½·(USp6×USp4²×SU2)   ≺  …
```

with midranges, at 40 digits and golden-section polished,

```
−0.02633875643733460272830901456544823309216      sup +0.35895 at τ = 1.406, inf −0.41163 at τ = 9.490
−0.04526323648456318677501288512742763277942      sup +0.33910 at τ = 5.655, inf −0.42963 at τ = 0.898
−0.03863526524068381663173215347728271905236      sup +0.20033 at τ = 0.624, inf −0.27760 at τ = 2.511
```

margin `2.633875643733460·10⁻²`. At genus 11 one has a **pure-product vertex**,

```
USp16×SU2³   ≺   ½·(USp18×SU2²) + ½·(USp12×USp8×SU2)   ≺   ½·(USp16×USp6) + ½·(USp14×USp4²)
```

midranges `−0.0433372226, −0.0265305251, −0.0609750245`, margin
`2.653052506327591·10⁻²`. Both reproduce to five digits on the independent grid
`10⁻⁵…10⁶`.

With the mixing weights free (`transitivity_freeweight.py`, differential
evolution plus Nelder–Mead on the three weights, seeded on the forty closest
triples of each genus) **a cycle appears already at genus 4**:

```
USp4 × SU2 × SU2
0.836189 · (USp4 × USp4)  +  0.163811 · SU2⁴
0.721649 · USp8           +  0.278351 · SU2⁴
```

`α_max = 8`, margin `8.674·10⁻³` on the search grid and `8.673·10⁻³` on the
independent grid. Genus 5, 6, 7 also cycle with free weights
(`4.5·10⁻³, 7.0·10⁻³, 2.5·10⁻⁴`).

**So the statement the brief asked to be proved is false**, and false at
`α_max = 8`, inside the range every library in this project already covers.
Same-genus transitivity was never a consequence of the shared endpoint; it is a
property of the rigidity of the *product* family, and it holds only up to
genus 12 even there. The mixtures are honest measures but they are not known to
be vertical Sato–Tate measures: that would need a compact group with two
components carrying exactly those coset trace measures and weights, and none is
exhibited — the same caveat `FINDINGS.md` attaches to its 3091 cross-genus hull
cycles. The genus-13 cycle of §5 does not carry that caveat.

---

## 7. The level lemma, proved — and why it does not close the gap *(proved + computed)*

Write `D₁, D₂, D₃` for the three differences around a triangle: `Σ_i D_i ≡ 0`
pointwise, and inside a genus `D_i(0) = D_i(∞) = 0`.

**Theorem 4 (two-level lemma).** Suppose there are two points `p, q` with

```
sup D_i = max(0, D_i(p), D_i(q))   and   inf D_i = min(0, D_i(p), D_i(q))
```

for each `i`. Then the triple is not a strict 3-cycle.

*Proof.* Put `a_i = D_i(p)`, `b_i = D_i(q)`, so `Σa_i = Σb_i = 0`; write
`M_i = max(0, a_i, b_i)`, `m_i = min(0, a_i, b_i)`, so `mid(D_i) = ½(M_i + m_i)`.
Suppose every `mid(D_i) < 0`. Then `M_i + m_i < 0`, which rules out
`a_i, b_i ≥ 0` (there `M_i + m_i = max(a_i,b_i) ≥ 0`). In the surviving cases
`M_i + m_i ≥ a_i + b_i`: equality when the signs are mixed, and
`M_i + m_i = min(a_i,b_i) ≥ a_i + b_i` when both are `≤ 0`. Summing,
`Σ(M_i + m_i) ≥ Σ(a_i + b_i) = 0`, contradicting `Σ(M_i + m_i) < 0`. ∎

**Theorem 5 (three-level cap).** If `|D_i| ≤ 1` and `Σ_i D_i ≡ 0` then
`min_i (−mid(D_i)) ≤ 1/4`, and `1/4` is attained.

*Proof.* With `P_i = sup D_i ≥ 0`, `Q_i = −inf D_i ∈ [0,1]`,
`−mid(D_i) = (Q_i − P_i)/2`. Suppose `Q_i − P_i > 1/2` for all `i` and let `i*`
maximise `Q`. Where `D_{i*}` attains its infimum the other two sum to `Q_{i*}`,
so `P_j ≥ Q_{i*}/2` for some `j ≠ i*`; then
`Q_{i*} ≥ Q_j > P_j + 1/2 ≥ Q_{i*}/2 + 1/2`, so `Q_{i*} > 1` — a contradiction.
Attained by `D₁ = (−1, ½, ½)`, `D₂ = (½, −1, ½)`, `D₃ = (½, ½, −1)`. ∎

Together they reproduce the whole `level_lemma.py` table without an LP and say
what it means:

| independent levels | optimum of `min_i(−mid)` |
|---|---|
| `≤ 2` (both ends pinned with `n ≤ 2` interior; or one end free with `n = 1`) | `0` — **no cycle** |
| `≥ 3` | exactly `1/4` |

A free endpoint *is* a level: that is the precise content of "the endpoint gap is
worth a whole level".

**But the hypothesis is badly false.** *Computed* (`transitivity_levels.py`), on
all 119 surviving same-genus triples of the symplectic library:

* the three interior **minima** sit at three distinct `τ` on every triple —
  never two — spread by a median factor of `57` and up to `1.5·10⁵`;
* the three interior **maxima** likewise, spread by a median factor of
  `4.5·10⁴`.

So the library uses six levels, not two, and Theorem 4 never applies to it. Route
2's plan — bound the number of independent levels — cannot work here: the level
count is already maximal and the triples are transitive anyway. Whatever holds
genus `≤ 12` together is a cancellation between six genuinely different levels,
not a shortage of them. The same conclusion transfers to brief G §4.1: counting
extrema is the wrong invariant; by Theorem 1 the extrema of `Ψ_μ − Ψ_ν` are the
parallel tangencies of two concave inverse-rate functions, and there are as many
of them as the shapes allow.

---

## Corrections

* **To `FINDINGS.md`'s headline reduction.** "The `q → ∞` question is now a
  single arithmetic question. Cycles persist for fixed genus profile as `q → ∞`
  **iff** there are one-parameter families of curves over `F_q` … whose
  Jacobians have an isogeny factor of multiplicity `k ≥ 2`." The "iff" is false.
  The genus-13 cycle of §5 is **multiplicity-free at every vertex**: no repeated
  isogeny factor anywhere. What its vertices need instead is a Jacobian splitting
  into four or five *distinct* independent factors — a much weaker demand. The
  correct statement is that a repeated factor is needed only if one insists on
  genus `≤ 12`.
* **To `FINDINGS.md`'s summary, and to brief J's scope.** "The limiting
  comparison is transitive inside a genus and non-transitive across genera …
  non-transitive **exactly** across genera." True for genus `≤ 12` in the
  multiplicity-free cone and for the whole `α_max ≤ 12` library, false in
  general: same-genus 3-cycles exist from genus 13. Brief J's confinement of the
  witness search to mixed-genus triples is therefore a convenience, not a
  theorem.
* **To `FINDINGS.md`'s crossing reduction.** "A pair whose `Ψ`s do not cross is
  ordered by pointwise domination, and a tournament all of whose edges are
  pointwise dominations is transitive. So **a 3-cycle needs all three pairs to
  cross**." The premise is right, the conclusion one edge too strong: two
  *non-crossing* edges of a triangle are consecutive and chain, so a 3-cycle
  needs **at least two** crossing edges. On the symplectic library 119 triples
  survive the correct test, 13 the incorrect one.
* **To brief K and to `lex_exceptions.py`.** "The fitted `0.668` is suspicious —
  it is close to `2/3`, and `2/3` would suggest a clean derivation." It is
  exactly `2/3`, and there is no derivation: `2/3` is `Δm₂/Δt` for one family of
  library pairs and hence the boundary of the interval of thresholds that
  classifies them correctly. The full constraint system is **infeasible**, and
  two `(Δm₂, Δt)` values carry both signs of `mid`, so no rule in those two
  statistics exists at all.
* **A sharpening, not a correction, of `level_lemma.py`.** Its docstring calls
  the `n = 2` value "an explanation and not a proof". It is now a three-line
  theorem (Theorem 4), and the `n ≥ 3` value is exactly `1/4` (Theorem 5). What
  remains conjectural is not the lemma but its hypothesis — and §7 shows the
  hypothesis is false on the library.

---

## Open

* **The one inequality.** `b_g(τ) ≤ 1/2` for the Jacobi coefficients of
  `e^{2τx}√(1−x²)`, equivalently `M_{g+1}M_{g−1} ≤ M_g²` for the `USp(2g)` trace
  MGFs, equivalently concavity of `g ↦ K_{USp(2g)}(τ)`. Verified to 40 digits
  over eleven decades and every rank to 14, with both endpoint regimes proved. A
  Toda-flow proof needs `a_g(τ)` non-increasing in `g`; a combinatorial proof
  would need the right total positivity of the `USp(2g)` trace-moment array.
  **With it, Theorem 3 is unconditional** — and Theorem 2 makes the whole
  dominance-comparable part of every genus unconditional as well, which is 83 %
  of the pairs even at genus 15.
* **Realising the genus-13 cycle.** Curve families of genus 13 with
  `Jac ∼ A₇ × A₂ × A₂' × E × E'`, `A₅ × A₄ × A₂ × A₂'` and
  `A₆ × A₄ × E × E' × E''`, all factors having independent full symplectic
  monodromy. Cyclic covers and quotients by automorphism groups give multi-factor
  splittings; the work is in getting the dimension vector and the independence
  simultaneously. This is now **the** arithmetic input, and it replaces
  `FINDINGS.md`'s repeated-factor question as the thing to settle.
* **Is genus 13 sharp?** Genus 12 has no cycle among 73150 triples and genus 13
  has two; the closest approach at genus 11–12 is `+0.0022`. The proof that
  genus `≤ 6` cannot cycle is Theorem 3; genus 7–12 is exhaustive computation
  with no structural reason attached. A structural explanation of the interval
  7–12 would be worth more than the computation.
* **A criterion for the crossing pairs.** §3 reduces everything to the crossing
  pairs, but §7 shows the level count is the wrong handle on them. The Legendre
  picture of Theorem 1 — parallel tangencies of concave inverse-rate functions
  with common endpoints — is the natural place to look for the missing
  inequality, and it is untouched.
* **Are the mixture cycles realisable?** The genus-4 counterexample uses
  weighted mixtures of products. A compact group whose two cosets carry those
  trace measures with those weights would push the arithmetic counterexample down
  from genus 13 to genus 4. `FINDINGS.md`'s swap construction
  (`½(μ∗μ) + ½δ₀`) is the only coset measure known to be realisable here, and it
  is not of this shape.

---

## Files

| file | what |
|---|---|
| `kappa_lib.py` | `K_{USp(2g)}`, the Hankel determinants and the Jacobi coefficients on a grid, one Bessel call per order per `τ`; agrees with `st_lib.mgf_classical` exactly |
| `transitivity_pairs.py` | every same-genus pair with `(Δm₂, Δt)`, `mid`, extrema locations and sign changes; the infeasibility certificate for the linear rule; writes `transitivity_pairs.csv` |
| `transitivity_dominance.py` | concavity of `κ_g` in `g`, the higher-order differences, and the multiplicity-free order genus by genus to 15 including the genus-13 cycles; writes `transitivity_dominance.csv` |
| `transitivity_certificate.py` | the crossing census and the reduction of 6840 triples to 119; writes `transitivity_certificate.csv` |
| `transitivity_cycle13.py` | the genus-13 cycles re-derived on a wider grid with polished extrema, tail bounds and cumulant re-validation; writes `transitivity_cycle13.csv` |
| `transitivity_levels.py` | the two proved level statements checked at random, and the measured spread of the extrema on the surviving triples; writes `transitivity_levels.csv` |
| `transitivity_duality.py` | Legendre invariance of `sup`, `inf`, `mid`, with a refinement study; writes `transitivity_duality.csv` |
| `transitivity_mixtures.py` | the same-`α_max` cycle search over products and uniform mixtures, cycles counted exactly; writes `transitivity_mixtures.csv` |
| `transitivity_freeweight.py` | free mixing weights, and the genus-4 counterexample; writes `transitivity_freeweight.csv` |
| `transitivity_verify.py` | the mixture cycles and the tightest certificates at 40 digits with golden-section polished extrema; writes `transitivity_verify.csv` |
