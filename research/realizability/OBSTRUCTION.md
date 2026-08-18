# Findings — there is no obstruction: `C_4` **is** realisable

Answer to session brief I. Every claim below is marked *proved* or *computed*;
the two headline items are certified at 40 digits (`i_certify.py`,
`i_witness.py`).

## Summary

The brief asked for a proof that no four signatures realise `s·C_4`. **The
opposite is true.** `C_4` is realisable, brief G §4.3's `1.255692` is a stalled
search and not an obstruction, and the reduction the brief proposed cannot work
because the relaxation it uses already realises `C_4`.

* **`C_4` is realisable** *(proved; certified to `6.6·10⁻⁴¹` at 40 digits)*.
  Explicitly, `d = s·C_4` with
  `s = 0.201980198313297395324406847972…`, realised by four points of the cone
  that is the projective closure of `{F_a = log Z_a}`. Every *strict*
  inequality of the construction holds with margin `≥ 6.757·10⁻³` (one
  monotonicity constraint is met with equality — two lines of `Φ_1` merge,
  which is legal); `i_certify.py`, `i_certify.json`.
* **Consequently the infimum of the `C_4` distortion over signature 4-tuples is
  exactly `1`** *(proved)*, approached as `1 + 4.54/log r`. An explicit integer
  family with `log₁₀ r ≈ 10.9` already gives certified distortion **`1.2254`**,
  below brief G's best; the ladder reaches `1.00019` (`i_witness.py`).
* **The brief's co-peaked reduction is a dead end** *(proved by explicit
  counterexample)*. Co-peaked unimodal 1-Lipschitz bumps under the envelope
  realise `s·C_4` exactly for every `s ≤ log 2 / 3`, so no relation is "forced
  by co-peakedness" (`i_copeaked.py`).
* **The right structure theorem is a Hilbert projective metric** *(proved)*.
  `d(a,b)` is exactly the Hilbert metric between `F_a` and `F_b` on the cone
  `C = {Φ convex, nondecreasing, Φ ≥ Λ_Φ·β}`, and `C` is exactly the projective
  closure of the achievable set. Cartesian powers are the projective rescaling
  the Hilbert metric quotients out — which reproves FINDINGS §1.6.
* **The defect constant is halved and is sharp** *(proved)*:
  `d ≤ |Δσ| + log(1 + e^{−|Δσ|})`, i.e. `ε ≤ log(1+e^{−Δ})`, not
  `2 log(1+e^{−Δ})`. Hence `d ≤ ℓ + log 2`. This settles the first companion
  question — but by a different mechanism than brief G guessed (see
  Corrections 3).
* **Convexity of `U` in `s` is proved**, in two lines:
  `U'' = U'·(F − βF')/F + β²F''/F > 0`. This settles the third companion
  question.
* **The maximum triangle curl is `(log 2)/2`** *(computed to 9 digits, in the
  cone, where the search is exact and unconstrained)*, one third of the proved
  bound `3(log 2)/2`. Not proved.

---

## Notation

That of `FINDINGS.md`. `Z_a(β) = Σ a_i^β`, `F_a = log Z_a`, `U_a(s) = log F_a(e^s)`,
`d(a,b) = osc_s(U_b − U_a)`, `A = mid`, `R = log r`, `Λ = log max a_i`,
`σ = log(R/Λ)`, `ψ = ½log(RΛ)`, `ε = P + Q`, `D = (P−Q)/2`,
`sp(t) = log(1+e^t)`, `sigm = sp'`. `φ = (1+√5)/2`. `(1,)` and all-ones
signatures are excluded throughout, so `r ≥ 2` and `Λ > 0`.

---

## 1. The exchange metric is a Hilbert projective metric *(proved)*

### 1.1 The identity

`d(a,b) = osc_s log(F_b/F_a) = log sup_β(F_b/F_a) + log sup_β(F_a/F_b)`, taken
over `β ∈ [0,∞]` with the two endpoints read as `R` and `Λ`. That is exactly the
**Hilbert projective metric** `d_H(F_a, F_b)` on the cone of positive functions
on `[0,∞]`. Nothing here is new except the name — but the name is the whole
point, because a Hilbert metric is insensitive to positive rescaling of either
argument, and `F_{a^{⊗k}} = k·F_a`. So:

> **The Cartesian power acts on `F` by the scalar `k`, and the Hilbert metric
> quotients scalars out.** This is FINDINGS §1.6's "`d` and `D` are unchanged,
> `ψ → ψ + log k`" with no computation.

### 1.2 The cone *(proved)*

**Theorem 1.** Let
`C = { Φ : [0,∞) → (0,∞) convex, nondecreasing, with Φ(β) ≥ Λ_Φ·β for all β }`,
`Λ_Φ = lim Φ'`. Then

1. `F_a ∈ C` for every signature `a`;
2. every `Φ ∈ C` is a limit of `(1/K)·F_{a^{(K)}}`, uniformly on `[0,∞)`, at rate
   `O(log k / K)`; hence `d(a^{(K)}, b^{(K)}) → d_H(Φ_a, Φ_b)` at the same rate.

*Proof.* (1) `F = log Σ m_i e^{βx_i}` is convex and nondecreasing (`x_i ≥ 0`);
`F ≥ Λβ` is FINDINGS §1.1. (2) A piecewise-linear `Φ ∈ C` is
`max_j (c_j + x_j β)` with every `c_j ≥ 0`: along the envelope the slopes
increase and the intercepts decrease, so `c_j ≥ c_last`, and `c_last ≥ 0` is
exactly `Φ ≥ Λβ` for large `β`; conversely `c_j ≥ 0` gives `Φ ≥ Λβ`. Taking `m_j = ⌈e^{Kc_j}⌉` copies of the atom `e^{Kx_j}` gives
`(1/K) log Σ m_j e^{Kβx_j} = max_j(c_j + βx_j) + O(log k / K)` uniformly, by the
standard `max ≤ LSE ≤ max + log k`. Since `Φ ≥ Φ(0) > 0`, the same bound
transfers to `log Φ` and hence to the oscillation. General `Φ ∈ C` are uniform
limits of piecewise-linear ones. ∎

*Verified* (`i_validate_cone.py`): 40 random cone pairs, `K·|d − d_H| = 1.35,
1.46, 1.49, 1.43, 1.40` at `K = 4, 8, 16, 32, 64` — the predicted `O(1/K)` with
constant `≈ 1.4`. Also `d` computed directly from `F` agrees with `common.py` to
`1.7·10⁻⁸` (grid resolution) on 200 random integer pairs.

**Corollary 1.1.** For any target metric `δ`, the infimum of the scale-free
distortion over signature `n`-tuples equals the infimum over `C`. So a
realisation *in the cone* settles realisability up to distortion `1 + o(1)`.

Inside `C` the metric is an **exact finite maximum**: `Φ_b/Φ_a` is a ratio of
affine functions between consecutive breakpoints, hence monotone there, so

```
d(a,b) = osc over {0, ∞} ∪ {breakpoints of Φ_a, Φ_b} of log(Φ_b/Φ_a).
```

No grid, no Lipschitz bracket. This is why every search in this session is
sharper than brief G's.

---

## 2. Two exact structure theorems *(proved)*

### 2.1 `U` is strictly convex in `s` — the open item of FINDINGS §1.2

**Theorem 2.** For every signature, `U(s) = log F(e^s)` satisfies

```
U'(s) = ν := βF'/F ∈ (0,1),      U''(s) = ν·(F − βF')/F + β²F''/F  >  0 .
```

*Proof.* `U' = βF'/F = ν`. Then `U'' = β dν/dβ = ν + β²F''/F − ν²
= ν(1−ν) + β²F''/F`, and `1 − ν = (F − βF')/F ≥ 0` because `F' ≤ Λ` and
`F ≥ βΛ` (FINDINGS §1.1); `F'' = Var_{Gibbs}(log a) ≥ 0`. Both terms vanish only
if `ν ∈ {0,1}` and all atoms are equal; a flat signature has
`ν = βΛ/(R+βΛ) ∈ (0,1)` strictly, so `U'' > 0` everywhere. ∎

*Verified* (`i_verify.py`): 214 signatures × 4001 points, no negative value; the
formula is cancellation-free once `F − βF'` is evaluated as
`log S + β·Σm_i y_i e^{−βy_i}/S` with `y_i = Λ − x_i ≥ 0`.

### 2.2 The scale function, and the sign of `(U_b − U_a)'`

**Theorem 3 (scale function).** Put `T(β) = (F(β) − βF'(β))/F'(β) > 0` — the
ratio *intercept / slope* of the tangent to `F` at `β` — and `S(s) = log T(e^s)`.
Then

```
T is strictly decreasing,     U'(s) = β/(β + T) = sigm( s − S(s) ),
```

and therefore, for any two signatures,

```
sign (U_b − U_a)'(s)  =  sign ( T_a(β) − T_b(β) )  =  sign ( S_a(s) − S_b(s) ).
```

Equivalently `U` is the upper envelope of translates of the softplus:
`U = sup_j [ γ_j + sp(· − σ_j) ]`.

*Proof.* `F − βF'` has derivative `−βF'' ≤ 0` and `F'` has derivative `F'' ≥ 0`,
so `T` = (nonincreasing)/(nondecreasing) is nonincreasing, strictly unless
`F'' ≡ 0`. It is positive by `F ≥ βΛ ≥ βF'`. The identity
`β/(β+T) = βF'/F = U'` is immediate. For the envelope form: on the `j`-th piece
of a piecewise-linear `Φ = max_j(c_j + x_j β)`, `log Φ = log x_j + sp(s − σ_j)`
with `σ_j = log(c_j/x_j)`, and `σ_j` is strictly decreasing along the envelope. ∎

*Verified* (`i_verify.py`): on the same 214 signatures the maximum *relative*
increase of `T` along the grid is `4.5·10⁻¹⁶` and
`max |U' − β/(β+T)| = 7.8·10⁻¹⁶`.

**Remark (the conjugate picture).** `βF' − F = F*(F'(β))` is the Legendre
transform of `F` evaluated at the slope, so `T = −F*(p)/p` with `p = F'(β)`.
Geometrically, the tangent to `F` at `β` meets the `β`-axis at `−T(β)`, so:

> `U_b − U_a` is critical exactly at those `β` where the tangents to `F_a` and
> to `F_b` **hit the `β`-axis at the same point**,

and it increases where `F_a`'s tangent-root is further left. That is the
"tangency" formulation brief K found for the rate functions, in the form it
takes here, and it is what makes the four-node analysis of §5 finite: `T` is
monotone, so a step-function `T` is the general case.

**Corollary 3.1 (free parametrisation).** Up to the projective scaling, an
element of `C` with `k` pieces is *any* pair of unconstrained sequences

```
σ_1 > σ_2 > … > σ_k        (the tangent scales)
s_1 < s_2 < … < s_{k−1}    (the breakpoints)
```

because `x_{j+1}/x_j = (e^{σ_j} + β_j)/(e^{σ_{j+1}} + β_j) > 1` for any
`β_j > 0`. Equivalently: **every nonincreasing step function `S` occurs**, and
`U` is recovered from `S` by `U' = sigm(s − S(s))`. This is the parametrisation
all of §5 uses.

**Corollary 3.2 (oscillation count).** In the cone, `U_b − U_a` has at most
`k_a + k_b − 2` interior extrema, `k` being the number of linear pieces
(= distinct tangent scales), because `S_a − S_b` is a step function with that
many jumps and each jump changes the sign at most once. This is FINDINGS §4.1's
conjecture, in the cone. *It is not the mechanism behind anything* — brief K's
level-count verdict transfers, and §5 shows there is no obstruction to be
mechanised anyway. Recorded as colour only.

---

## 3. The sharp defect constant *(proved)* — brief G's first open question

**Theorem 4.** For any two signatures, with `Δ = |σ_a − σ_b|`,

```
d(a,b) ≤ Δ + log(1 + e^{−Δ}),      i.e.   ε = P + Q ≤ log(1 + e^{−Δ}),
```

and the constant is attained in the limit. In particular `d ≤ ℓ + log 2`, half
of FINDINGS Corollary A2.

*Proof.* Normalise `β` and `Φ_a` so that `R_a = Λ_a = 1` (so `σ_a = 0`) and
`Λ_b = 1`, `R_b = ρ = e^{σ_b}`; assume `σ_b ≥ σ_a`, i.e. `ρ ≥ 1` (otherwise swap
`a` and `b`). Let `β_1` maximise and `β_2` minimise `Φ_b/Φ_a`, so
`e^d = (Φ_b(β_1)/Φ_b(β_2))·(Φ_a(β_2)/Φ_a(β_1))`. Only the sandwich
`max(R,Λβ) ≤ Φ ≤ R + Λβ` and `Φ' ≤ Λ` are used.

*Case `β_1 = t < T = β_2`.* `Φ_b(β_2) ≥ max(Φ_b(β_1), T)` and
`Φ_b(β_1) ≤ ρ + t`, so `Φ_b(β_1)/Φ_b(β_2) ≤ min(1, (ρ+t)/T)`; and
`Φ_a(β_2) ≤ Φ_a(β_1) + (T − t)` with `Φ_a(β_1) ≥ max(1,t)`, so
`Φ_a(β_2)/Φ_a(β_1) ≤ 1 + (T−t)/max(1,t)`. If `t ≤ 1` the product is at most
`1 + ρ`: for `T ≤ ρ + t` it is `≤ 1 + T − t ≤ 1 + ρ`, and for `T > ρ + t` the
function `T ↦ (ρ+t)(1+T−t)/T` has derivative of the sign of `t − 1 ≤ 0`, so it is
maximal at `T = ρ + t`, where it equals `1 + ρ`. If `t > 1` the product is at
most `max( 1 + ρ/t , (ρ+t)/t ) = 1 + ρ/t < 1 + ρ`.

*Case `β_2 = t < T = β_1`.* Now `Φ_b(β_1)/Φ_b(β_2) ≤ 1 + (T−t)/max(ρ,t)` and
`Φ_a(β_2)/Φ_a(β_1) ≤ min(1, (1+t)/T)`. If `T ≤ 1 + t` the product is at most
`1 + 1/max(ρ,t) ≤ 1 + 1/ρ`. If `T > 1 + t` and `t ≤ ρ` the product is
`(ρ+T−t)(1+t)/(ρT)`, whose `T`-derivative has the sign of `t − ρ ≤ 0`, so it is
maximal at `T = 1 + t`, where it equals `1 + 1/ρ`. If `T > 1+t` and `t > ρ` the
product is at most `(T/t)·((1+t)/T) = 1 + 1/t < 1 + 1/ρ`. Since `ρ ≥ 1`,
`1 + 1/ρ ≤ 1 + ρ`.

Hence `e^d ≤ 1 + ρ = 1 + e^{Δ}`, i.e. `d ≤ log(1+e^{Δ}) = Δ + log(1+e^{−Δ})`. ∎

*Sharpness.* `Φ_a = max(1, β)` against the flat `Φ_b = ρ + β` gives exactly
`d = Δ + log(1+e^{−Δ})`. In signatures those are `a_r = (r,1,…,1)` with
`r → ∞` and a flat `b` — FINDINGS §1.4's own ladder. `i_verify.py` runs it:
`d/bound = 0.9978, 0.9985, 0.9994` at `r = 10²⁰⁰` for `Δ = 0, 0.5, 1.5`, rising
monotonically to 1.

*Verified:* over 22 791 pairs from three random integer pools plus adversarial
thin/flat ones the bound is never violated (worst slack `−1.2·10⁻³`).

*Computed independently:* a hill-climb over the cone at fixed `Δ`
(`i_constants.py` and the `Δ`-ladder) returns

| `Δ` | 0 | 0.10 | 0.25 | 0.50 | 0.75 | 1.00 | 2.00 |
|---|---|---|---|---|---|---|---|
| `sup ε` found | 0.69315 | 0.63111 | 0.57539 | 0.47319 | 0.38670 | 0.31324 | 0.12692 |
| `log(1+e^{−Δ})` | 0.69315 | 0.64440 | 0.57594 | 0.47408 | 0.38687 | 0.31326 | 0.12693 |

**Corollary 4.1.** `|D| ≤ ε/2 ≤ ½log(1+e^{−Δ}) ≤ (log 2)/2`, unchanged and still
sharp (`i_constants.py` attains `0.346573590` in the cone).

---

## 4. The large-scale bound for `C_4` *(proved)* — the brief's LP argument

**Theorem 5.** If four signatures realise `s·C_4` then `s ≤ log φ = 0.4812118…`,
`φ` the golden ratio. With only FINDINGS' unsharpened `2log(1+e^{−Δ})` the bound
is `s ≤ 0.7644901717…`, the root of `s = 2log(1+e^{−s})`.

*Proof.* Write `ℓ_ab = |σ_a − σ_b|`, a **line metric** on the four points. By
Theorem 4, `ℓ_ab ≤ d_ab ≤ g(ℓ_ab)` with `g(ℓ) = ℓ + log(1+e^{−ℓ})` increasing,
so `g^{-1}(d_ab) ≤ ℓ_ab ≤ d_ab`. For four reals `x_1 ≤ x_2 ≤ x_3 ≤ x_4` the three
perfect-matching sums are `x_3+x_4−x_1−x_2` (twice) and `x_2+x_4−x_1−x_3` (once),
so **the largest matching sum of a line metric is attained at least twice**.
Hence whichever matching is the diagonal one `{13, 24}`, some side matching
`{12,34}` or `{14,23}` has `ℓ`-sum at least as large. Therefore

```
2s ≥ ℓ_12 + ℓ_34   or   ℓ_14 + ℓ_23   ≥   ℓ_13 + ℓ_24   ≥   2 g^{-1}(2s),
```

i.e. `g(s) ≥ 2s`, i.e. `log(1+e^{−s}) ≥ s`, i.e. `e^{2s} − e^{s} − 1 ≤ 0`,
i.e. `e^s ≤ φ`. ∎

The realised scale of §5 is `0.2020`, comfortably inside. The bound is about
scale, not about realisability, and — as brief I anticipated — it says nothing
about small `s`. What it does say is that **any** realisation is small-scale.

---

## 5. `C_4` is realisable *(proved; certified)*

### 5.1 The pattern at the four contact points is forced

Suppose `d = s·C_4` on the 4-cycle `1–2–3–4`. Let `β⁺, β⁻` maximise and minimise
`φ_13` and `γ⁺, γ⁻` do the same for `φ_24`. Since
`φ_13 = φ_12 + φ_23 = φ_14 + φ_43` with each summand of oscillation exactly `s`
and `φ_13` of oscillation `2s`, `β⁺` maximises and `β⁻` minimises **each** of
`φ_12, φ_23, φ_14, φ_43`; likewise `γ±` for `φ_21, φ_14, φ_23, φ_34`. Reading off
the six values, and using that `U_a` is free up to an additive constant (the
projective scaling) and that a common profile may be subtracted:

```
U_a(P_k) = κ_k + c_a + s·T[a][k],
T = [[0,0,0,0],[1,0,0,1],[2,1,0,1],[1,1,0,0]],
(P_1,…,P_4) = (β⁺, γ⁺, β⁻, γ⁻).
```

The four points are pairwise distinct (each pair is separated by some `φ` taking
its max at one and its min at the other), so at most two of them are the
endpoints `0, ∞` and at least two are interior.

*(Aside: the pattern is an isometric copy of `s·C_4` in `ℓ_∞⁴/ℝ1` with the
oscillation norm, so the obstruction — if there were one — could only come from
the shape of the admissible `U`, never from the `osc` structure.)*

### 5.2 Feasibility, and hence realisability

By Corollary 3.1 one may take `S_a` constant on each cell cut by
`P_1 < … < P_4`; then every `U_b − U_a` is monotone on each cell, the
oscillation is exactly the oscillation over the four nodes, and the increments
are

```
U_a(θ_{k+1}) − U_a(θ_k) = δ_{L_k}( S_{a,k} − θ_k ),
δ_L(v) = sp(L − v) − sp(−v) : R → (0, L) decreasing.
```

Setting `S` to a common value on the two outer cells makes `φ` constant there, so
`φ(0)` and `φ(∞)` add nothing. Feasibility is therefore the finite system

```
u_{a,k} = ρ_k + s(T[a][k+1] − T[a][k]) ∈ (0, L_k),
S_{a,k} = θ_k + δ_{L_k}^{-1}(u_{a,k})   nonincreasing in k.
```

**It is feasible.** `i_pattern.py` maximises `s` over the seven unknowns
`(L_1,L_2,L_3, ρ_1,ρ_2,ρ_3, s)` for each of the 24 orderings of
`(β⁺, γ⁺, β⁻, γ⁻)`; the run is in `i_pattern_output.txt` (a partial run: 15 of the 24
orderings, each recorded row complete). **Most orderings are feasible**: of the
15 recorded, 14 are feasible and only `γ⁺β⁻β⁺γ⁻` is not. Separate
higher-effort runs on individual orderings push the scale up:

| `θ`-ordering | largest scale `s` found |
|---|---:|
| `β⁺ γ⁺ β⁻ γ⁻` | 0.2619 |
| `β⁺ γ⁺ γ⁻ β⁻` | 0.2135 |
| `β⁺ β⁻ γ⁺ γ⁻` | 0.1117 |
| `β⁺ β⁻ γ⁻ γ⁺` | 0.1117 |

Feasibility is the whole point: **any** `s > 0` gives an exact realisation, and
the numbers only bound how large the scale can be made (see Open). The witness
of §5.3 uses the ordering `β⁺ γ⁺ γ⁻ β⁻` at `s = 0.20198`.

### 5.3 The certified witness

`i_certify.py` rebuilds the `β⁺ γ⁺ γ⁻ β⁻` solution in `mpmath` at 60 digits
(400-step bisection for `δ_L^{-1}`), forms the four `Φ_a = max_j(c_j + x_jβ)`
exactly, and evaluates `d` at 40 digits. Result:

```
s     = 0.201980198313297395324406847972
θ     = 0,  3.838690666329304568,  5.2289587638930064,  7.1710085468037781
S_0   = 8.92527974788611096, 5.25556557243888581, 4.61856421579323452, 1.13867743754035860, 0.63867743754035860
S_1   = 8.92527974788611096, 8.42527974788611096, 4.00338417509673069, 4.00338417509672104, 0.63867743754035860
S_2   = 8.92527974788611096, 8.42527974788611096, 4.61856421579323452, 4.00338417509672104, 0.63867743754035860
S_3   = 8.92527974788611096, 5.25556557243888581, 5.24880825574554980, 1.13867743754035860, 0.63867743754035860

d     = s · C_4   with   max_ij | d_ij − s·(C_4)_ij |  =  6.60·10⁻⁴¹
```

An **independent cross-check** inside the same script recomputes the whole
`d` matrix in double precision on a `4·10⁶`-point grid in `log β`, sharing no
code with the exact breakpoint evaluation: it returns distortion `1.0000055`,
the residual `2.2·10⁻⁶` being the grid's own resolution at the one contact that
sits exactly on a breakpoint.

Margins of every strict inequality the construction uses: `min(u, L−u) ≥
9.917·10⁻³`, node gaps `≥ 1.390`, `min_ij d_ij = 0.2020`, and the monotonicity
gaps `S_{a,k} − S_{a,k+1} ≥ 6.757·10⁻³` **except** for the single pair
`S_{1,2} = S_{1,3}` which is exactly `0` — two lines of `Φ_1` merging, which is
legal (`Φ_1` then has four pieces, not five) and is not a strict inequality of
the construction. All margins are far above the brief's `10⁻⁶` threshold and
above the `10⁻¹⁰` tie threshold. Full data in `i_certify.json`.

### 5.4 Signatures: the infimum is `1`

Theorem 1(2) turns the cone witness into signatures. With `B = λ·max_j c_j` the
multiplicity budget (`B ≈ log r`), `i_witness.py` reports

| `B` | `log₁₀ r` | distortion | realised scale | `B·(dist − 1)` |
|---:|---:|---:|---:|---:|
| 25 | 10.9 | 1.2316762697 | 0.180018 | 5.79 |
| 50 | 21.7 | 1.1177880571 | 0.191476 | 5.89 |
| 100 | 43.4 | 1.0547370779 | 0.197259 | 5.47 |
| 200 | 86.9 | 1.0239965103 | 0.200022 | 4.80 |
| 400 | 174.0 | 1.0110722192 | 0.201220 | 4.43 |
| 800 | 347.7 | 1.0057294381 | 0.201675 | 4.58 |
| 1600 | 695.1 | 1.0028383082 | 0.201838 | 4.54 |
| 3200 | 1390.0 | 1.0014191829 | 0.201909 | 4.54 |
| 6400 | 2779.7 | 1.0007101055 | 0.201945 | 4.54 |
| 12800 | 5559.1 | 1.0003578482 | 0.201962 | 4.58 |
| 25600 | 11118.0 | 1.0001908432 | 0.201969 | 4.89 |

> **Theorem 6.** `inf { distortion(d, s·C_4) : four signatures, s > 0 } = 1`,
> and the infimum is approached like `1 + 4.54/log r`.

An **exact integer family** at `B = 25` (five distinct atoms each, atoms of 7 to
574 digits, multiplicities up to `7.2·10¹⁰`, listed in `i_witness.json`) has
`d`-matrix evaluated at 40 digits

```
0            0.205200414  0.407167843  0.221721372
0.205200414  0            0.210231227  0.361868955
0.407167843  0.210231227  0            0.194032905
0.221721372  0.361868955  0.194032905  0
```

certified distortion **`1.225424`**, already better than brief G's `1.255692`
with a fraction of the search effort — and it is the *first rung* of the ladder.
(It differs slightly from the ladder's `1.2317` because flooring the
multiplicities and rounding the atoms moves the configuration; here it happens
to move it downhill.)

*Not settled:* whether distortion exactly `1` is attained by an actual signature
family rather than only in the closure. A damped Gauss–Newton on the five
scale-free equations, over 20 real log-atoms (and over all 40 parameters) from
the `B = 100` and `B = 300` rungs, stalls at residual `1.8·10⁻³` and `6.2·10⁻⁴`
(distortion `1.01586`, `1.00594`) — a genuine stationary point of the residual,
not a failure of the solver. The five-line ansatz is presumably too rigid; more
atoms give more freedom. See Open.

---

## 6. Why the brief's small-`s` route could not have worked *(proved)*

Brief I proposed: with all `σ_a` equal, `d(a,b) = osc(w_b − w_a)` for co-peaked
unimodal bumps, and asked for "a linear relation forced by co-peakedness". There
is none, because the class

```
W = { w ≥ 0 : w(±∞) = 0, unimodal with peak at 0, 1-Lipschitz, w ≤ E(t) = log(1+e^{−|t|}) }
```

— exactly what FINDINGS Theorem A gives — **realises `s·C_4` exactly**.
`i_copeaked.py` builds four explicit members of `W`, piecewise linear on
`(−1.5, −0.9, −0.6, 0)` with value pattern

```
w_1 : 0, 2ς, 2ς, 3ς      w_2 : 0, ς, 2ς, 2ς
w_3 : ς, ς, 2ς, 3ς       w_4 : 0, ς, ς, 3ς
```

(clipped by `E`, tapered by `E`-shaped tails, right half `h_a·E(t)/log 2`), and
checks unimodality, the Lipschitz constant and `w ≤ E` on a grid of `8·10⁵`
points:

| `ς` | 0.05 | 0.10 | 0.15 | 0.20 |
|---|---:|---:|---:|---:|
| distortion | 1.0000000000 | 1.0000000000 | 1.0000000000 | 1.1724914 |

so `s·C_4` is realised in `W` for every `s ≤ log 2/3 = 0.2310` (the envelope
begins to bind at `ς = 0.20`, where `max w = 3ς = 0.6 → log 2`).

The same relaxation is strictly weaker than the truth in a second, quantitative
way: inside `W`, `sup osc(w_b − w_a) = 2h* = 0.9624236…` where `h*` solves
`h = log(1+e^{−h})` (a plateau against a tent, and their mirror images), whereas
the truth for signatures is `log 2 = 0.6931` (Theorem 4). **The bump
relaxation loses too much; the cone of §1 is the right object.**

---

## 7. The remaining companion question — the triangle curl *(computed)*

`curl A = Σ_cyc D_ij`, and `|D| ≤ (log 2)/2` gives `3(log 2)/2 = 1.0397`.
Hill-climbing in the cone (`i_constants.py`), where the objective is an exact
finite maximum and the parameters are unconstrained,

| pieces `k` | 2 | 3 | 4 | 6 |
|---|---:|---:|---:|---:|
| `sup \|D\|` | 0.346573590 | 0.346573590 | 0.346573590 | 0.346573590 |
| `sup \|curl\|` | 0.346573590 | 0.346573590 | 0.346573590 | 0.346573590 |

`0.346573590 = (log 2)/2` to nine digits, from four independent parametrisations,
and never `0.9889·(log 2)/2` as in brief G. So:

> *Computed, not proved:* the maximum triangle curl is exactly `(log 2)/2`,
> one third of the proved bound. Equivalently, around a 3-cycle the **total**
> rate asymmetry — not just the geometric mean — is at most a factor 2.

A proof does not follow from Theorem 4 alone: `Σ_cyc log(1+e^{−Δ_ij})/2` over a
line-metric triple is still `3(log 2)/2` at `Δ ≡ 0`. It needs a joint
constraint on the three `(P, Q)` pairs.

---

## 8. The other targets of §4.3 — a first pass *(computed)*

The `C_4` argument is a **sufficient test for any target**, and it is cheap:

1. find an isometric copy of `δ` in the oscillation norm on `m+1` nodes, i.e.
   `Y ∈ ℝ^{n×(m+1)}` with `osc_k (Y[b,k] − Y[a,k]) = δ_ab`;
2. test the node system of §5.2 for that `Y`.

A feasible point realises `δ` **exactly** in the cone, hence with distortion
`1 + O(1/log r)` by signatures (Theorem 1). `i_targets.py` runs it; it
reproduces `C_4` as a self-check. For the uniform metric step 1 is free: any
**antichain of subsets** works, since `osc(1_{A_b} − 1_{A_a}) = 2` whenever
neither set contains the other — in particular `Y = I_n` (the singletons) and
the `k`-subsets of an `m`-set.

Results so far (each "REALISABLE" is a *proof*; each "infeasible" is only a
statement about that ansatz):

| target | ansatz | verdict | brief G §4.3 |
|---|---|---|---|
| `C_4` | forced 4-node pattern | **REALISABLE** | 1.2557 |
| `K_4` | `Y = I_4` | **REALISABLE** | 1.0002 |
| `K_5` | `Y = I_5` | this ansatz infeasible | not tested |
| `K_5` | 2-subsets of 4 nodes | **REALISABLE**, `s = 0.0372` | not tested |
| `K_6` | `Y = I_6` | this ansatz infeasible | 1.4500 |

(the `K_6` singleton row is from a separate earlier run; `i_targets_output.txt`
holds a partial re-run of the same sweep.)

**`K_5` is realisable**, and it is realisable on *four* nodes although the
singleton ansatz on five nodes is not — so the node count, not the number of
points, is what binds, and FINDINGS §4.3's "an `n`-point equilateral set in
`ℓ∞^k` needs `n ≤ 2^k`, and the effective dimension here is 2–3" is not the
right heuristic. `K_6` also has an antichain embedding on four nodes (the six
2-subsets of a 4-set); whether it is feasible was still running when this file
was written — see `i_targets_output.txt` for the completed table.

---

## Corrections

1. **FINDINGS §4.3 is wrong: `C_4` is realisable.** "Only the 4-point
   equilateral metric is realisable; nothing else tried is, even up to scale" and
   "this is the one place where realisability looks genuinely obstructed" must be
   withdrawn. `1.255692` was a stalled 24-parameter search, and this session's
   very first *integer* rung already beats it (`1.2254`). **The obstruction was
   inadequate search, not geometry.**
2. **The other §4.3 verdicts are now unsafe too, and one of them is already
   contradicted.** `C_5`, `K_{2,3}`, `K_{3,3}`, `K_6`, `K_8`, Petersen and
   Wagner `V_8` were declared unrealisable by the same search on the same
   objective. §8 shows `K_5` (untested by brief G, and the natural companion of
   its `K_6 = 1.4500`) **is** realisable, on four nodes. The rest should be
   treated as *open*, not as negative results, and §4.3's `ℓ∞`-dimension
   heuristic should be dropped.
3. **"`ε ≤ log(1+e^{−Δ})` because one of `P, Q` is always zero" — the conclusion
   is right, the reason is wrong.** `min(P,Q)` is *not* always zero: a hill-climb
   in the cone reaches `max min(P,Q) = 0.2231435513 = log(5/4)`, stably at
   `k = 2,…,6` pieces (`i_constants.py`). Brief G's "always with one of `P, Q`
   exactly zero" is an artefact of hill-climbing over atom multiplicities. The
   bound itself is now proved (Theorem 4) by a route that never needs the claim.
4. **`sup ε = 0.9889·log 2` is a search artefact; the true supremum is exactly
   `log 2`.** Likewise `max|D| = 0.9889·(log 2)/2` → exactly `(log 2)/2` and
   `max|curl| = 0.9889·(log 2)/2` → exactly `(log 2)/2`. FINDINGS §1.4/§3.5's
   `0.9889` is what a bounded-atom climb reaches, not a supremum.
5. **FINDINGS §4.3's distortions are biased low by ~`10⁻⁵`.**
   `g2_metrics.distortion` takes the *raw* grid max/min of `φ` on a step-`0.01`
   grid with no parabolic refinement (unlike `common.matrices`), which
   underestimates `d` by a different amount for each pair. Re-running the same
   search and re-scoring with `common.certified_extrema` gives `1.264058`,
   `1.271621`, `1.255702` at `r = 3,4,6` against the reported `1.264051`,
   `1.271617`, `1.255692` (`i_check_g2.py`). Small, but the reported numbers are
   not upper bounds for the configurations they name.
6. **FINDINGS §1.2's "computed, not proved: `U` is convex in `s`" is now
   proved** (Theorem 2), and the constant `3.6·10⁻⁹` was second-difference noise.
7. **FINDINGS §1.6 ("the Cartesian power acts on the potential and on nothing
   else") is a one-line corollary of the Hilbert-metric identity**, not a
   separate verification.

---

## What this decides for the papers

The rigidity of FINDINGS is real and Theorem A stands. What must change is the
sentence after it. The framework is **not** metrically obstructed at `n = 4`;
what it is, is

> the **Hilbert projective metric on the cone of convex nondecreasing `Φ`
> with `Φ ≥ Λ_Φ β`**, whose elements are exactly the `log Z` of signatures up to
> Cartesian power and tropical limit; equivalently, `U = log F` runs over the
> upper envelopes of translates of `log(1+e^t)`, and `d` is the oscillation of
> differences.

That is a sharper and more useful statement than "one-dimensional plus a bounded
defect", it contains FINDINGS §1.1–§1.6 as corollaries, and it says what the
framework *can* do: it can realise `C_4`, and the price is `log r`. The honest
headline for the metric side is therefore **not** an obstruction theorem but a
rate:

> `C_4` is realisable with distortion `1 + O(1/log r)` — measured constant
> `4.54` — and cannot be realised at any scale above `log φ`. (Whether the rate
> is also a lower bound is open.)

The `A`-side conclusions of FINDINGS (§2, §3) are untouched, and §3.4's
"geometric-mean rate asymmetry at most 2 around any cycle" is unchanged and
still the recommended one-sentence claim.

---

## Open

* **Is distortion exactly `1` attained by a signature family**, or only in the
  projective closure? Gauss–Newton on a five-line ansatz stalls at residual
  `6·10⁻⁴`; the natural next step is 7–9 distinct atoms per signature and a
  least-squares solve seeded from a `B ≈ 10³` rung.
* **The largest realisable scale for `C_4`.** The four-node model gives
  `s ≥ 0.2619`; Theorem 5 gives `s ≤ log φ = 0.4812`. More nodes should raise the
  lower bound. Closing the gap would give a genuinely sharp constant.
* **Which metrics *are* realisable?** With Corollary 1.1 and Corollary 3.1 this
  is a finite feasibility question for each target (§8), and `i_targets.py`
  runs it. `C_5`, `K_{2,3}`, `K_{3,3}`, `K_6`, `K_8`, Petersen and Wagner `V_8`
  should all be pushed through it, with several oscillation embeddings and
  several node counts each. The evidence so far (`C_4`, `K_4`, `K_5` all
  realisable; `K_5` and `K_6` infeasible only for the *singleton* ansatz) points
  towards the exchange metric being **universal in the closure**, but that is a
  conjecture, not a result, and the right way to settle it is a proof that every
  oscillation embedding can be re-timed into a monotone `S` — or a target for
  which no `Y` at all is feasible.
* **Prove `max |curl A| = (log 2)/2`** (§7). The missing ingredient is a joint
  constraint on the three `(P,Q)` pairs of a triangle.
* **Is the `4.54` in `1 + 4.54/log r` intrinsic?** It should be
  `c·log(#atoms)/(s·log r)` for a `c` of order 1; a matching lower bound
  `distortion ≥ 1 + c'/log r` would be the correct replacement for the withdrawn
  obstruction — a *rate* theorem in place of an impossibility theorem.
* **A lower bound for the `2 log 2` metric budget** (still open from FINDINGS).

---

## Files

All under `research/realizability/`.

| file | what |
|---|---|
| `i_cone.py` | the tropical cone `C`: `Trop` (max of lines), exact Hilbert metric, `from_centers` (the free parametrisation of Cor. 3.1), tropical → signature translation |
| `i_validate_cone.py`, `_output` | `d = osc log(F_b/F_a)` on genuine signatures; the `O(1/K)` tropical limit |
| `i_verify.py`, `_output` | numerical verification of Theorems 2, 3, 4 and of the Theorem 5 constants; the sharpness ladder for Theorem 4 |
| `i_constants.py`, `_output` | hill-climbs in the cone: `sup ε`, `max min(P,Q)`, `sup \|D\|`, `sup \|curl\|` |
| `i_check_g2.py`, `_output` | audit of FINDINGS §4.3's `1.255692` against certified extrema |
| `i_nodes.py` | the exact finite (common-node) model of the cone; `dmat` from a nonincreasing step scale function |
| `i_deep.py` | basin-hopping search in the node model, seeded from brief G's own best `C_4` signatures; it stalls at `1.2618` at four nodes, which is why the *forced pattern* of `i_pattern.py` was needed |
| `i_c4.py`, `i_search.py` | earlier searches in the `(centres, breakpoints)` and `(intercept, slope)` parametrisations; kept because they show where the naive searches stall (`1.2558–1.2624`) |
| `i_pattern.py`, `_output` | the forced four-node pattern of §5.1 and its feasibility over all 24 orderings — **the proof that `C_4` is realisable** |
| `i_certify.py`, `_output`, `i_certify.json` | the witness rebuilt at 60 digits and certified at 40: `‖d − s·C_4‖_∞ = 6.6·10⁻⁴¹`, all margins, plus an independent double-precision cross-check on a `4·10⁶`-point grid |
| `i_witness.py`, `_output`, `i_witness.json` | cone witness → signatures: the `1 + 4.54/log r` ladder and an exact integer family with certified distortion `1.2254` |
| `i_logsig.py` | signatures held as `(log multiplicity, log atom)`, so the ladder can run to `r = 10^{11118}` |
| `i_copeaked.py`, `_output` | the co-peaked relaxation realises `C_4` exactly — §6 |
| `i_solve_sig.py`, `_output` | damped Gauss–Newton on the five scale-free equations: the "is distortion `1` attained by an actual signature family?" test |
| `i_targets.py`, `_output` | §8: the general sufficient test — oscillation embedding, then the monotone-`S` node system — applied to the uniform metrics and the graph metrics of §4.3 |

Reproduce in order: `i_validate_cone.py`, `i_verify.py`, `i_constants.py`,
`i_check_g2.py`, `i_pattern.py`, `i_certify.py`, `i_witness.py`,
`i_copeaked.py`, `i_solve_sig.py`, `i_targets.py --quick`. The last two and
`i_pattern.py` are the slow ones (tens of minutes each); `i_targets.py` without
`--quick` adds `K_8` and the graph metrics and runs for hours.
