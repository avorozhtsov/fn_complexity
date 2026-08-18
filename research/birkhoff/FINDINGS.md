# Findings — a processor is **not** a Birkhoff contraction, and `log 2` **is** a projective diameter

Answer to session brief N. Transcribed from the reporting agent's verified
output (the harness blocked it from writing this file). Every claim is marked
*proved* or *computed*; the headline constants are certified at 40 digits
(`n2_diameter.py`, `n3_curl.py`).

## Summary

The brief's stated failure mode is the right one, but the truth is sharper — in
both directions.

* **No resource operation of this framework is a Birkhoff contraction**
  *(proved)*, for two independent reasons that should not be conflated.
  Tensoring is a **translation** `Φ ↦ Φ + F_c` on the cone that carries the
  metric — not linear, and not even homogeneous of degree 1, so it falls outside
  Birkhoff **and** outside Nussbaum's nonlinear Perron–Frobenius theory. The
  operations that *are* linear — Cartesian power, fibre power, their tensor
  mixtures, and the flattening retraction `L` — all shift `σ` by a constant,
  hence have infinite projective diameter and Birkhoff ratio `tanh(∞/4) = 1`;
  and the value `1` is *attained*, because each is an isometry on the flat
  subcone.
* **Tensoring does not merely fail to contract: it expands `d` without bound**
  *(proved; exact family, 40 digits)*. `d(a, a^{⊗k}) = 0` for every `k`, yet with
  `Φ_a = max(1,β)` and a flat `Φ_c = N + εβ`,

  ```
  d(a⊗c, a^{⊗k}⊗c) = log [ (k+ε)(1+N) / ((1+ε)(k+N)) ]  ⟶  log k .
  ```

  The disjoint union does the same. So conversion is not a map on the metric
  space at all: `d` is a pseudometric on resources and `⊗` does not descend to
  its quotient. Any statement of the form "conversion loses distinguishability at
  a definite rate" is **false**.
* **The correct positive statement is in Thompson's metric, not Hilbert's**
  *(proved)*. All four resource operations — `⊗`, `⊔`, Cartesian power, fibre
  power — are **nonexpansive in Thompson's part metric**
  `d_T = max(log sup Φ_b/Φ_a, log sup Φ_a/Φ_b)`. Hilbert and Thompson differ
  exactly along the projective direction, and the projective direction is what
  `⊗` trades in.
* **A clean sufficient condition** *(proved)*: if `a` and `b` are **mutually
  non-dominating** — `C(a→b) ≤ 1` and `C(b→a) ≤ 1`, i.e. neither is the better
  resource at *every* temperature — then `d(a⊗c, b⊗c) ≤ d(a,b)` for every `c`.
  On the `F_3` quadratic, `F_3` cubic and `F_5` quadratic map classes **every**
  pair is non-dominating and tensoring strictly contracts on every triple
  (largest contraction `0.0649`).
* **`log 2` is a projective diameter** *(proved)*. With
  `C_0 = {Φ ∈ C : Φ(0) = Λ_Φ}` the `σ = 0` fibre,

  ```
  diam_d C_0 = log 2 exactly, attained between max(1,β) and 1+β,
  ```

  the two extreme rays of the fibre. Consequently `(C,d)` fibres over the
  `σ`-line: `σ` is 1-Lipschitz, the `β`-dilation group acts by isometries with
  **geodesic** orbits, and every fibre has diameter `log 2`. The attached
  Birkhoff ratio is `tanh((log 2)/4) = 3 − 2√2 = 0.1715728752538099024…`, and it
  is achieved by explicit positive linear maps into `C_0` (§3.4).
* **The `(log 2)/2` curl lower bound is now proved, not computed**
  *(proved; 40 digits)*. `Φ_1 = max(1,β)`, `Φ_2 = 1+β`, `Φ_3 = 1+e^Tβ` give
  `curl A = ½(log 2 − log(1+e^{−T})) ↗ (log 2)/2`. So `sup |curl A| ≥ (log 2)/2`,
  the supremum is **not attained**, and only the upper bound remains open.
* **The upper bound is reformulated.** `curl A = ½ log Ω` with `Ω` the cycle
  asymmetry ratio *(proved)*, so "`max curl = (log 2)/2`" is exactly "`Ω ≤ 2`
  around every triangle" — the forward arbitrage exceeds the backward one by at
  most **one bit**. The candidate `|curl| ≤ max_e |D_e|` is **refuted** (ratio `3`
  attained); `|curl| ≤ max_e ε_e` survives and would give `log 2`.
* **Where the curl extremum lives** *(computed)*. Inside a single `σ`-fibre the
  supremum is `½ log(4/3) = 0.14384104` (9 digits), and `½ log(5/4) = 0.11157178`
  (11 digits) with two lines — against `(log 2)/2` globally. **No argument
  confined to one fibre can produce the global constant.**
* **Refuted: the exchange framework and quantum-channel contraction are not the
  same theorem on two cones** *(proved)*. A channel is a linear map on the object
  carrying the metric; here conversion is an *inequality between objects*, and
  the only maps are the monoidal operations, which are translations.
* **Correction to the brief's own premise** *(proved)*. `d` is the Hilbert metric
  of the ambient cone `K` of **positive functions** (pointwise order), restricted
  to `C` — *not* the intrinsic Hilbert metric of `C`. Birkhoff's `Δ` must be
  measured in the cone whose order defines the metric, and the domain-restricted
  form of the theorem is **false** (§2.5, one-line counterexample). This is the
  step the brief warned was easy to wave at, and waving at it is exactly what
  would produce a spurious contraction theorem.

---

## Notation

That of `research/realizability/OBSTRUCTION.md`. `Z_a(β) = Σ a_i^β`,
`F_a = log Z_a`, `U_a(s) = log F_a(e^s)`, `L(a,b) = sup_s (U_b − U_a) = −log C(a→b)`,
`d = osc`, `A = mid`, `R = F(0) = log r`, `Λ = lim F' = log M`,
`σ = log(R/Λ)`, `ψ = ½log(RΛ)`, `ε = P + Q`, `D = (P−Q)/2`,
`sp(t) = log(1+e^t)`, `E(t) = log(1+e^{−|t|})`.

Two cones and two metrics, which must be kept apart:

```
K   = { Φ : [0,∞] → (0,∞) },  pointwise order
      d(Φ,Ψ)   = log sup(Ψ/Φ) + log sup(Φ/Ψ)          Hilbert
      d_T(Φ,Ψ) = max( log sup(Ψ/Φ), log sup(Φ/Ψ) )    Thompson
C   = { Φ ∈ K : convex, nondecreasing, Φ ≥ Λ_Φ β } ⊂ K,   the achievable set
C_0 = { Φ ∈ C : Φ(0) = Λ_Φ }                             the σ = 0 fibre
```

`C` is a subcone of `K`, but the metric is `K`'s, not `C`'s.

---

## 1. The dictionary *(proved)*

| resource operation | signature | action on `Φ = F` | action on `Z` |
|---|---|---|---|
| parallel composition `a ⊗ c` | `(a_i c_j)` | `Φ ↦ Φ + F_c` (translation) | `Z ↦ Z·Z_c` |
| Cartesian power `a^{⊗k}` | — | `Φ ↦ kΦ` (scalar) | `Z ↦ Z^k` |
| fibre power `a^{[m]}` | `(a_i^m)` | `Φ ↦ Φ(mβ)` (dilation) | `Z ↦ Z(mβ)` |
| labelled alternative `a ⊔ c` | `sort(a,c)` | `Φ ↦ log(e^Φ + e^{F_c})` | `Z ↦ Z + Z_c` |
| flattening `L` | `r` copies of `M` | `Φ ↦ Φ(0) + Λ_Φ β` | — |

All five verified against genuine integer signatures on a `β`-grid `[0,10³]`:
errors `1.4·10⁻¹²`, `2.7·10⁻¹²`, `0`, `4.5·10⁻¹³`.

**The two linear structures are exchanged by `log`.** `⊔` is addition on `Z`;
`⊗` is addition on `Φ = log Z`. The exchange rate is a ratio of `F`'s, so the
metric lives one logarithm above the cone on which the resource semiring is
linear. Everything in §2 follows from that sentence.

### 1.1 The linearity audit *(computed)*

Relative additivity / homogeneity defects over 400 random cone pairs:

| operation | `T(Φ+Ψ) − TΦ − TΨ` | `T(λΦ) − λTΦ` | linear? |
|---|---:|---:|---|
| tensor with a fixed `c` | `4.6e−1` | `1.8e+0` | **no** |
| Cartesian power `k=3` | `4.6e−16` | `4.4e−16` | yes |
| fibre power `m=2` | `4.4e−16` | `4.3e−16` | yes |
| dilation mixture | `4.7e−16` | `6.1e−16` | yes |
| flattening `L` | `4.3e−16` | `4.3e−16` | yes |
| corner `N : Φ ↦ max(Φ(0), Λ_Φβ)` | `1.8e−3` | `3.3e−16` | **no** |
| disjoint union with a fixed `c` | `5.0e−1` | `2.3e+0` | **no** |

Note the pair `L` / `N`: brief G's structure theorem is written with the tropical
corner `N`; only the flat retraction `L` is available to Birkhoff.

---

## 2. Birkhoff's verdict *(proved)*

**Proposition 2.1.** A map `T` on `C` is additive and homogeneous iff
`T(a ⊗ b) = T(a) ⊗ T(b)` and `T` commutes with Cartesian powers — i.e. a linear
induced map is exactly a `⊗`-homomorphism commuting with powers. Tensoring with a
fixed `c` is not one. ∎

**Proposition 2.2.** The positive linear operators on `C` commuting with the
`β`-dilation group are exactly

```
T_μ Φ (β) = μ({0})·Φ(0) + ∫_{(0,∞)} Φ(cβ) dμ(c) + μ({∞})·Λ_Φ β ,
```

`μ` positive on `[0,∞]`. Every finitely supported `μ` on `(0,∞)` is a resource
operation, namely `⊗_j (a^{[c_j]})^{⊗ w_j}`; `μ = δ_0 + δ_∞` is flattening. ∎

**Theorem 2.3.** For every `T_μ`, `σ(T_μΦ) = σ(Φ) + const`. Since
`d(a,b) ≥ |σ_a − σ_b|` and `σ` is unbounded on `C`, `Δ(T_μ) = ∞` and the Birkhoff
ratio is `tanh(∞/4) = 1`. Moreover `T_μ` maps flats to flats and two flats have
`d = |Δσ|` exactly, so `T_μ` is an **isometry on the flat subcone**: its Lipschitz
constant is exactly `1`. ∎

*Verified:* over `2·10⁵` random cone triples every linear operation has
`max (d(Ta,Tb) − d(a,b)) = 0.000000` and `sup d(Ta,Tb)/d(a,b) = 1.000000000`;
`d(T flat_0, T flat_W) = W` exactly for `W = 1, 10, 100, 600`.

> **Corollary 2.3.1.** No processor of this framework is a Birkhoff contraction.
> The Birkhoff coefficient of every linear resource operation is exactly `1`, and
> the Hilbert-metric identification of brief I is **descriptive, not dynamical**.

### 2.3 Tensoring: nonlinear, and expanding

**Proposition 2.4.** With `Ξ = F_c`:

1. If `Φ_a/Φ_b` takes the value `1` somewhere — equivalently `C(a→b) ≤ 1` **and**
   `C(b→a) ≤ 1` — then `d(a⊗c, b⊗c) ≤ d(a,b)` for every `c`.
   *Proof.* `(Φ_a+Ξ)/(Φ_b+Ξ)` is a mediant of `Φ_a/Φ_b` and `1`, hence lies
   pointwise between them; if `1` is in the range, the mediant's range is
   contained in it. ∎
2. The condition cannot be dropped: `d(a, a^{⊗k}) = 0` while
   `d(a⊗c, a^{⊗k}⊗c) = log[(k+ε)(1+N)/((1+ε)(k+N))] ↗ log k`.

*Certified:* at `N = 10⁴⁰`, `ε = 10⁻⁴⁰` the closed form agrees with `log k` to
`1.5·10⁻⁴⁰` (`k=2`), `2.7·10⁻⁴⁰` (`k=3`), `7.9·10⁻⁴⁰` (`k=8`). A hill-climb of the
expansion restricted to non-dominating pairs returns `+0.000000000000`, so
2.4(1) is sharp. The disjoint union does the same on genuine signatures. Over
`2·10⁵` triples the largest expansion under `⊗` is `+1.7308`, the largest ratio
`2.256·10³`.

### 2.4 Thompson nonexpansiveness *(proved)*

**Theorem 2.5.** For every `c`, `d_T(a⊗c, b⊗c) ≤ d_T(a,b)` and
`d_T(a⊔c, b⊔c) ≤ d_T(a,b)`; Cartesian and fibre powers are `d_T`-isometries.

*Proof.* For `⊗`: with `t = sup(Φ_b/Φ_a)`,
`Φ_b + Ξ ≤ tΦ_a + Ξ ≤ max(t,1)(Φ_a + Ξ)`. For `⊔`: `g(x) = log(e^x + e^{ξ})` with
`ξ ≥ log 2 > 0` is increasing and `g(tx) ≤ t g(x)` for `t ≥ 1` by
superadditivity of `y ↦ y^t`. ∎

*Verified:* `max (d_T(Ta,Tb) − d_T(a,b)) = −7.2·10⁻⁴` for `⊗` over `10⁵` triples,
`+0.000000` for `⊔`.

> **Conversion is nonexpansive for the absolute (Thompson) geometry and expansive
> for the projective (Hilbert) one**, because Cartesian powers — the projective
> direction — are exactly what conversion trades in. Nussbaum's theory needs
> degree-1 homogeneity, and 2.4(2) is a quantitative witness that it cannot be
> dropped.

### 2.5 Which cone carries `Δ` *(proved)*

**The domain-restricted form of Birkhoff's theorem is false.** If `Δ(T)` could be
the `d`-diameter of `T(K')` for a subcone `K' ⊆ K` on which `T` is linear, take
`K' = C_0` and `T = Id`: the diameter is `log 2` (§3.1), so the "theorem" would
give Lipschitz constant `tanh((log 2)/4) = 0.1716 < 1` for the identity. Hence
`Δ(T)` must be the diameter of the image of the **whole** cone whose order
defines the metric.

---

## 3. `log 2` as a projective diameter *(proved)*

**Theorem 3.1.** `C_0` is a subcone of `C` and `diam_d C_0 = log 2`, attained
between `Φ_min(β) = max(1,β)` and `Φ_max(β) = 1 + β`.

*Proof.* `R` and `Λ` are additive, so `R = Λ` is preserved by `+`. Normalise
`R = Λ = 1`; the sandwich `max(R,Λβ) ≤ Φ ≤ R + Λβ` reads `Φ_min ≤ Φ ≤ Φ_max`. The
upper bound is brief I's Theorem 4 at `Δ = 0`. For attainment,
`Φ_max/Φ_min = (1+β)/max(1,β)` has supremum `2` at `β = 1` and infimum `1` at
both ends. ∎

*Verified:* `d(Φ_min,Φ_max) = 0.693147180559945`; over `2·10⁵` random pairs the
largest distance is `0.990772 log 2`; hill-climbs at `k = 2,3,5` lines return
`0.693147180449`.

**Remark.** The generic order-interval bound is `2 log 2` — exactly FINDINGS
Corollary A2. Theorem 3.1 says the fibre is *half* as wide as its order interval,
which is what brief I's Theorem 4 adds.

**Theorem 3.2.** The `β`-dilation group acts by `d`-isometries,
`σ(Φ(c·)) = σ(Φ) − log c`, each orbit is an isometric copy of `ℝ`, and

```
|σ_a − σ_b| ≤ d(a,b) ≤ |σ_a − σ_b| + log(1 + e^{−|σ_a−σ_b|}) ≤ |σ_a−σ_b| + log 2 .
```

So `σ : (C,d) → (ℝ,|·|)` is surjective 1-Lipschitz, all fibres isometric to `C_0`
of diameter `log 2`, with a global geodesic section.

> **The exchange geometry is the real line thickened by exactly one bit.**

*Verified:* `max |d(Φ,Φ(c·)) − |log c|| = 8.9·10⁻¹⁶`; the hill-climbed
`sup{d : Δσ = Δ}` matches `Δ + log(1+e^{−Δ})` to `1.1·10⁻¹⁰`.

### 3.3 Constants *(computed, 40 digits)*

```
log 2           = 0.6931471805599453094172321214581765680755
(log 2)/2       = 0.3465735902799726547086160607290882840378
tanh((log 2)/4) = 0.1715728752538099023966225515806038428607  =  3 − 2√2
```

exactly, since `e^{2·(log2)/4} = √2` and `(√2−1)/(√2+1) = 3−2√2`; agreement
`−1.8·10⁻⁴⁶` at 45 digits.

### 3.4 Where a genuine contraction lives *(proved; computed)*

`R` and `Λ` are positive linear functionals, so for positive `μ_1, μ_2`

```
T Φ = ( ∫ Φ dμ_1 )·max(1,β) + ( ∫ Φ dμ_2 )·(1+β)
```

is positive, linear on all of `K`, with image in `C_0`. Hence `Δ(T) ≤ log 2` and
Birkhoff applies: `d(TΦ,TΨ) ≤ (3 − 2√2)·d(Φ,Ψ)`. *Verified* for
`μ_i = δ_{0.7}, δ_{4.0}`: computed `Δ(T) = log 2` to 12 digits, empirical
Lipschitz ratio `0.140033` over `10⁵` pairs, inside `0.171573`.

**But `T` is not a processor.** It reads the free energy at two fixed
temperatures and rebuilds a flat resource — a measurement-and-rebuild map, not a
conversion. A Birkhoff contraction here must *forget* `σ`, and no resource
operation forgets `σ`.

---

## 4. The triangle curl

### 4.1 The lower bound, in closed form *(proved; 40 digits)*

**Theorem 4.1.** With `Φ_1 = max(1,β)`, `Φ_2 = 1+β`, `Φ_3 = 1+e^Tβ`:
`D_12 = (log 2)/2` exactly, `D_23 = 0`, `D_31 = −½ log(1+e^{−T})`, so

```
curl A = ½ ( log 2 − log(1 + e^{−T}) )  ↗  (log 2)/2 .
```

*Proof.* `U_1 = max(0,s)`, `U_2 = sp(s)`, `U_3 = sp(s+T)`. For `(1,2)`,
`φ = E(s)` gives `P = log 2`, `Q = 0`. For `(3,1)`,
`φ' = 1_{s>0} − sigm(s+T)` is negative on `s<0` and positive on `s>0`, so the
minimum is at `s = 0`, giving `P = 0`, `Q = log(1+e^{−T})`. ∎

*Certified:* closed form and exact tropical computation agree to `0` /
`2.8·10⁻¹⁶` at `T = 1, 5, 20, 100, 300`; at `T = 100, 300` the value equals
`(log 2)/2` to all 40 digits.

> This upgrades OBSTRUCTION §7's nine-digit hill-climb to
> `sup |curl A| ≥ (log 2)/2`, **proved**, with the supremum **not attained**.

### 4.2 Curl is a cycle asymmetry ratio *(proved)*

`curl A = ½ log Ω` with

```
Ω = ( sup(Φ_2/Φ_1)·sup(Φ_3/Φ_2)·sup(Φ_1/Φ_3) ) / ( sup(Φ_1/Φ_2)·sup(Φ_2/Φ_3)·sup(Φ_3/Φ_1) ) ,
```

both cycle sums being `≥ 0` by the quasi-metric triangle inequality — the
**forward and backward arbitrage** — with sum the perimeter.

> `max |curl A| = (log 2)/2` ⟺ `Ω ≤ 2`: **around a 3-cycle the forward arbitrage
> exceeds the backward one by at most one bit.**

*Verified* to `8.9·10⁻¹⁶` over `2·10⁵` triples; largest sampled `Ω` is `1.888944`.

### 4.3 Two exact reductions *(proved)*

**(a) Circulation identity.** With `p_e, q_e` the argmax/argmin of `U_j − U_i`,
`2 curl A = J(p) + J(q)` where
`J(x) = ∫_{x_3}^{x_2} u_1 + ∫_{x_1}^{x_3} u_2 + ∫_{x_2}^{x_1} u_3`, `u_i = U_i'`
nondecreasing from `0` to `1`. *Verified:* residual `1.9·10⁻¹⁵` over `4·10⁴`
triples.

**(b) Envelope/cocycle form.** `U_i = n_i − g_i` with `n_i` flat and
`0 ≤ g_i ≤ E(s − σ_i)`; then `U_j − U_i = m_{ij} + h_{ij}` with `m` monotone and
`h` a pointwise cocycle, `Σ_cyc h = 0`, and
`D_ij = ½[ sup(h_ij − ν_ij) − sup(−h_ij − ν̃_ij) ]` with `ν + ν̃ ≡ |Δσ|`.
Dropping `ν, ν̃` gives the per-edge `(log 2)/2` in two lines — but summing around
a cycle gives only `3(log 2)/2`, and the relaxation is genuinely lossy. **Any
proof of the sharp bound must retain `ν`, i.e. couple the per-edge defect to the
`σ`-gaps of the other two edges.**

### 4.4 The upper bound *(computed)*

Hill-climbs, `k` lines per point, 5–8 restarts of differential evolution plus a
ten-scale pattern search:

| `k` lines | 2 | 3 | 4 | 6 | 8 |
|---|---:|---:|---:|---:|---:|
| `sup \|curl A\|` | 0.346573590 | 0.346573590 | 0.346573590 | 0.346573590 | 0.346573590 |

Five parametrisations agreeing with `(log 2)/2` to nine digits, never exceeding.

Two candidate joint constraints:

* **`|curl A| ≤ max_e |D_e|` is refuted** — the ratio reaches exactly `3`. This is
  the route that would have closed §7 in one line.
* **`|curl A| ≤ max_e ε_e` survives** every search (saturating at `1.0000000000`);
  it would give `log 2`, still a factor `2` above the truth. Not proved.

**The extremum is not inside one fibre** *(computed)*:

| `k` lines | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| `sup \|curl A\|`, equal `σ` | 0.111571776 | 0.143841029 | 0.143831369 | 0.143835850 | 0.143773758 |

```
k = 2 :  0.111571775650   vs   ½ log(5/4) = 0.111571775657   (7·10⁻¹²)
k ≥ 3 :  0.143841028898   vs   ½ log(4/3) = 0.143841036226   (7·10⁻⁹)
```

so intra-fibre `Ω = 4/3` (and `5/4` with two lines) against `Ω = 2` globally; the
extremal intra-fibre triangle is equilateral in the defect.

> Any proof of `|curl A| ≤ (log 2)/2` must use the `σ`-spread. The fibre diameter
> `log 2` is the right constant but the wrong mechanism on its own: it bounds a
> single edge, and the extremal triangle is one edge at the fibre diameter plus
> two edges asymptotically defect-free because their `σ`-gap is infinite.

---

## 5. Placement against Birkhoff, Bushell and Reeb–Kastoryano

Birkhoff (1957): a linear cone-preserving map contracts the Hilbert metric with
ratio `tanh(Δ/4)`. Bushell (1973): Banach-space form and sharpness.
Reeb–Kastoryano–Wolf (2011): the PSD cone, relating `tanh(Δ/4)` to trace-norm and
relative-entropy contraction coefficients. Nussbaum: nonexpansiveness for
order-preserving **degree-1 homogeneous** maps.

**The structural difference.** In all of that the channel is a *linear map on the
object carrying the metric*. Here the object is `F = log Z`, the resource
semiring is linear one logarithm below, conversion is an *inequality between
objects* rather than a map, and the only genuine maps — the monoidal operations —
are **translations**.

**Verdict: refuted.** The two are not the same theorem on two cones. What
survives:

| quantum channels | exchange framework |
|---|---|
| cone of PSD operators | cone `K` of positive functions of `β` |
| state `ρ` | `Z_a`, or `F_a = log Z_a` |
| channel `T` (linear, positive) | *no analogue*; conversion is an order relation |
| Hilbert contracted, ratio `tanh(Δ/4)` | contracted only by measure-and-rebuild maps (§3.4) |
| `Δ(T) < ∞` for full-rank output | `Δ = ∞` for every linear resource operation |
| Dobrushin / relative-entropy coefficients | Thompson nonexpansiveness (Thm 2.5) |

The analogy does bite for `⊔` on the `Z`-cone: `Z ↦ Z + Z_c` is the mediant map
and contracts `δ(a,b) = osc_β(F_b − F_a)`, the `Z`-cone Hilbert metric, under the
same non-domination hypothesis. But `δ` is `⊗`-invariant, hence blind to
Cartesian powers, and is not the exchange metric.

---

## Corrections

1. **Brief N's premise is misstated, and the misstatement is the trap.** `d` is
   the Hilbert metric of the ambient cone of **positive functions**, restricted to
   `C`; not the intrinsic Hilbert metric of `C`. Birkhoff's `Δ` is the diameter of
   the image of the cone whose *order* defines the metric, and the
   domain-restricted version is false (§2.5).
2. **OBSTRUCTION §7's `max |curl A| = (log 2)/2` is no longer purely computed.**
   The lower bound is proved in closed form (Theorem 4.1) and the supremum is
   shown not to be attained. The upper bound remains open between `(log 2)/2` and
   `3(log 2)/2`.
3. **The natural guess `|curl A| ≤ max_e |D_e|` is false** — the ratio attains `3`.
4. **No correction to brief I or brief G is required.** Every statement leant on
   here reproduced exactly, and Theorem 4's `log 2` is confirmed independently as
   a projective diameter.

---

## Open

* **Prove `Ω ≤ 2`**, equivalently `|curl A| ≤ (log 2)/2`. §4.3 gives two exact
  reductions; §4.4 says the argument must use the `σ`-spread; §4.1 pins the
  extremal configuration.
* **Prove or refute `|curl A| ≤ max_e ε_e`** — would give `log 2`, a factor `3/2`
  better than the proved bound.
* **Prove the intra-fibre constants `½ log(4/3)` and `½ log(5/4)`.** Bounded
  problem, possibly within reach; would be the project's first exactly-solved
  cycle constant. Is this `log(5/4)` the same as `max min(P,Q)` in OBSTRUCTION
  Correction 3?
* **Cycles of length `n > 3`.** `curl A = ½ log Ω` holds verbatim; does the
  supremum grow with `n`?
* **A Birkhoff-type theorem for the Thompson metric** — a nontrivial contraction
  coefficient for `⊗` with fixed `c` on a bounded region.
* **The `Z`-cone metric `δ`.** Does it have an operational reading — a one-shot
  rather than asymptotic conversion cost? If so, Reeb–Kastoryano transfers.

---

## Files

All under `research/birkhoff/`; they import `i_cone.py` and `optimizers.py` from
`research/realizability/`.

| file | what |
|---|---|
| `n1_maps.py`, `_output` | the dictionary; the linearity audit; the contraction audit in Hilbert and Thompson; exact unbounded-expansion families for `⊗` and `⊔`; the non-dominating hill-climb; the `F_3`/`F_5` map classes |
| `n2_diameter.py`, `_output` | `diam C_0 = log 2`; the `2 log 2` comparison; the two-fibre law at 40 digits; `3 − 2√2`; the counterexample to domain-restricted Birkhoff; an explicit contracting linear map; the dilation action |
| `n3_curl.py`, `_output` | the exact `curl → (log 2)/2` family at 40 digits; the `Ω` reformulation; the circulation identity; the upper-bound search; the equal-`σ` slice; the two candidate constraints |
| `n4_fibre_curl.py`, `_output` | intra-fibre curl constants `½log(4/3)` and `½log(5/4)`, with the extremal configuration |
