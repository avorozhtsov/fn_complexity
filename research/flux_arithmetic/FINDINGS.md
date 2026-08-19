# Findings — is the flux arithmetic?

Answer to session brief E. **No, and the negative is sharp enough to be
useful.** The curl of the flux `A` — the part of the exchange comparison that no
scalar can reproduce — is *not* an arithmetic phenomenon. A random pool of
signatures matched to the arithmetic pool in the only things that matter (`q`
fibers, `Σ N_c = q²`, and the trace distribution) has the same curl to within a
few percent, and the tightest control — one whose multiset of largest fibers is
*identical* to the arithmetic pool's — has slightly **more** curl and **more**
cycles. Over 108 pools, arithmetic and control alike, `‖curl‖/‖A‖` is a function
of the trace spread to `R² = 0.93`, and every arithmetic pool is an ordinary
residual of that regression (`−0.53σ` to `+1.24σ`, against a control spread of
`±1.07σ`). The brief's seed measurement, `0.091` against `0.041–0.051`, compared
`q`-entry curve signatures with 2-to-7-entry random integers; the factor of two
is the shape of the signature, not the arithmetic.

Two things do come out positively, and they are the paper if there is one.

* **`ψ_opt` is what brief D's flat-connection theorem predicts between
  `φ`-classes, and it is the moment ladder inside one.** The endpoint potential
  `½ log φ` is within `1.7·10⁻⁴` of the best possible function of the largest
  fiber at every `q` from 11 to 23. Inside a class of equal largest fiber —
  where `φ` is exactly tied and where every cycle lives — the potential is
  `m₂, m₃, m₄, log μ` with **`m₂` entering with the opposite sign to brief B's
  addendum**, reaching `R² = 0.978` at `q = 11` and `0.985` at `q = 13`.
* **The certified cycle descends to the curves, and it descends *too well*.**
  The three signatures of the `F_11` cycle are realised by 2, 2 and 3
  non-isomorphic pencils, so the cycle is a cycle among `2·2·3 = 12` explicit
  triples of genus-two fibrations. But the two pencils realising the signature
  of `B` have 1 and 4 fibers whose Jacobian splits, and the comparison cannot
  tell them apart. Over the complete `F_11` enumeration, **every one of the 209
  signatures realised by more than one pencil is realised by pencils with
  different fiberwise isogeny data** — 419 of 420 at `q = 13`. The flux is an
  invariant of the pair of signatures and of nothing finer.

Everything below is either proved or computed; each statement says which.
Five earlier claims are corrected in [Corrections](#corrections).

---

## Notation

Conventions of `research/m_and_e_and_a_c/PLAN.md` and brief D Part 0, with signs
fixed once here because the brief's seed script carried the opposite ones.

`Z_a(β) = Σ_i a_i^β`, `u_a(β) = log log Z_a(β)`, `C(g→f) = inf_β log Z_g/log Z_f`
with the implementer first, `L(a,b) = −log C(a→b)`. Then with `g = u_b − u_a`

```
L = S + A,     S(a,b) = ½ osc g = d(a,b)/2,     A(a,b) = mid g,
a ≺ b  ⟺  C(a→b) < C(b→a)  ⟺  A(a,b) > 0.
```

`A` is antisymmetric; `S` is the metric and is comparison-blind.

**HodgeRank.** For antisymmetric `A` on the complete graph the least-squares
gradient part is `grad ψ_opt` with

```
ψ_opt = −rowmean(A),      (grad ψ)(a,b) = ψ(b) − ψ(a),
```

and `curl = A − grad ψ_opt`. Norms are Frobenius over ordered pairs; "energy" is
the square of the norm ratio.

Signatures: `σ(f) = {N_c}` with `N_c = #f⁻¹(c)` the **affine** plane count,
`Σ_c N_c = q²`, `a_c = q − N_c`, `m_k = (1/q^k) Σ_c a_c^k`,
`M = max_c(−a_c) = max_c N_c − q`, `μ` the multiplicity of the largest fiber,
`ν(P) = #{(x,x′) : P(x) = P(x′)}`, `φ(a) = log(#fibers)·log(max fiber)`, and the
endpoint potential `ψ_end = ½ log φ`.

**A caution used throughout.** At a fixed `q` every pool member has exactly `q`
fibers, so `M`, `log max_c N_c` and `½ log φ` are three *link functions of the
same integer* `max_c N_c`. They induce **the same order** and differ only as
approximations to the potential. Any table that ranks them ranks shapes, not
statistics.

---

## Reproduction of the seed, and the engine *(computed)*

`build_f11_pool.py` and `hodge_split.py` reproduce exactly: 296 signatures,
`‖grad‖/‖A‖ = 0.9959`, `‖curl‖/‖A‖ = 0.0908`, 132 strict 3-cycles, and brief E's
E2 `R²` column to six decimals.

`pools.py` here is a third, independent enumeration. The signature depends on
`P` only through its value-multiplicity vector `n_u`, because
`N_c = q + Σ_u n_u χ(u+c)`, so the whole enumeration is one histogram and one
`q × q` matrix multiply. It gives **296** positive signatures at `q = 11` from
175461 pencils and **698** at `q = 13` from 399763 — brief B's census numbers,
now confirmed by three unrelated codes. `flux.py` then gives **132** and
**1475** strict 3-cycles, again brief B's numbers exactly.

Accuracy of the engine:

| check | result |
|---|---|
| grid `A` against `fn_complexity.exchange_rate`, 200 random `q=11` pairs | max deviation `4.2·10⁻⁸` |
| grid self-convergence, density `×2` then `×2` again | `1.66·10⁻⁸` then `4.0·10⁻⁹` (second order in the step) |
| `cycle_count` against brute force over all triangles, `n = 40, 80, 120` | identical (0, 4, 7) |
| smallest `|A(a,b)|` in the `q = 11` pool | `1.30·10⁻⁷`, three orders above the grid error |
| ties at the `10⁻¹⁰` floor | **0** in both pools; every pair is strict |

Runs below use the `×4` grid (`|ΔA| < 5·10⁻⁹`) with horizon `360 q`, well past
the `9.3√q` at which the `F_101` witness has its deepest contact.

The certified `F_11` cycle, re-verified against the package to twelve digits:

```
A->B  0.990213498322   B->A  0.994908823027   margin 4.695e-03
B->C  0.981637513410   C->B  0.983352120195   margin 1.715e-03
C->A  0.979537663762   A->C  0.981637513410   margin 2.100e-03
```

and the two endpoint contacts check out — `C(B→C) = C(A→C) = log 18/log 19`
exactly — so the spurious-`β = 258.9` trap does not bite here. The `F_101` cycle
of brief B's addendum 2 is re-verified in `verify_f101.py`: signatures
recomputed by point count (each summing to `10201`), flux
`A = (−1.619720, −1.301503, −1.164720)·10⁻⁴` on the three edges, reproducing the
addendum's table to all quoted digits, `r = |curl A|/Σ|A| = 1.000000000000`.

---

## E1 — the Hodge decomposition on arithmetic pools *(computed; negative)*

### Why the split is the right measurement, and it is a theorem

**Lemma (proved).** *Let `A` be antisymmetric on the complete graph on a set
`V`. Then `A(a,b) + A(b,c) + A(c,a) = 0` for every triple iff `A = grad ψ` for
some `ψ : V → ℝ`.*

*Proof.* If `A = grad ψ` the sum telescopes. Conversely fix `o ∈ V` and put
`ψ(x) = A(o,x)`; the triangle `(o,x,y)` reads `A(o,x) + A(x,y) + A(y,o) = 0`,
i.e. `A(x,y) = ψ(y) − ψ(x)`. ∎

So **every scalar invariant, at every order of any expansion, contributes
exactly zero to the curl**: `φ`, the addendum's `φ̃`, and whatever comes next are
annihilated without being computed. That is the formulation of
`research/m_and_e_and_a_c/curl_on_curve_families.py` (brief D's session), and it
is a better justification for E1 than the brief's own; the energy split
`‖curl‖/‖A‖` is the `L²` version of the same statistic.

### The calibration run, and a discrepancy inside the brief

`research/realizability/tournament_seed.py`, re-run unmodified: `‖curl‖/‖A‖` =
`0.0511, 0.0510, 0.0407` at `n = 8, 16, 24` with **0, 0, 0** strict 3-cycles.
That reproduces the `0.041–0.051` of brief E's *E1-seed* paragraph. It does
**not** reproduce the `0.065–0.088` with `0, 3, 2` cycles of brief E's *E1*
paragraph; those numbers are not what the committed script produces. Over ten
seeds the statistic is very noisy at these sizes anyway:

| `n` | `‖curl‖/‖A‖`, 10 seeds | 3-cycles |
|---:|---|---|
| 8 | `0.0153 … 0.0491` (median `0.0263`) | 0 … 0 |
| 16 | `0.0250 … 0.0641` (median `0.0389`) | 0 … 1 |
| 24 | `0.0283 … 0.0491` (median `0.0363`) | 0 … 2 |

### The matched controls

Every control has the two structural properties a curve signature has — `q`
positive entries summing to `q²`, equivalently `q` traces summing to zero — so
the `β = 0` endpoint is a tie across the pool exactly as it is for curves. They
differ in how much of the trace *distribution* they are told to match
(`pools.control_pool`):

| control | what it matches |
|---|---|
| `loose` | traces uniform on the genus-two Weil box `[−4√q, 4√q]`; deliberately over-spread (`m₂ ≈ 4–5` against the pool's `1.0–1.3`) |
| `m2matched` | each signature's `m₂` drawn from the pool's empirical `m₂` distribution |
| `marginal` | traces iid from the pooled empirical trace law of the arithmetic pool — the Sato–Tate law of that pool; destroys only *which* traces co-occur inside one pencil |
| `sigshuffle` | each order statistic drawn from that order statistic's empirical law |
| `maxmatched` | **the multiset of largest fibers is identical to the arithmetic pool's**, so `ψ_end = ½ log φ` agrees signature by signature and the whole gradient part is matched by construction |

A construction trap worth recording: the obvious way to force `Σ a_c = 0` after
sampling is to move units off the extreme coordinate. The typical excess is
`√(q·var)`, so that shaves a dozen units off the largest entries and destroys
`m₂` — the first version of these controls came out at `m₂ = 0.73` against the
pool's `1.33`. Spread the correction evenly instead.

Full pools, five replicates each:

| `q` | `n` | arithmetic | `loose` | `m2matched` | `marginal` | `sigshuffle` | `maxmatched` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 11 | 296 | **0.09084** | 0.11499 | 0.07896 | 0.08961 | 0.07967 | 0.09701 |
| 13 | 698 | **0.08517** | 0.13068 | 0.08140 | 0.08566 | 0.07457 | 0.09411 |

(control standard deviations `0.0014`–`0.0043`.) Strict 3-cycles, same pools:

| `q` | arithmetic | `loose` | `m2matched` | `marginal` | `sigshuffle` | `maxmatched` |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | **132** | 1688 | 90 | 75 | 84 | 145 |
| 13 | **1475** | 31278 | 1223 | 1112 | 628 | 1820 |

**The arithmetic pool is in the middle of its own controls, on both statistics.**
The `marginal` control — which knows only the one-point trace law — reproduces
the curl fraction to 1.4% at `q = 11` and 0.6% at `q = 13`. The `maxmatched`
control, the tightest, *exceeds* it on both. The one control far away is
`loose`, and it is far in the direction of *more* curl, because it is more
spread out: the curl is set by how far the interior excursion beats the
endpoint, which is exactly brief B's mechanism.

### It is not an `n` effect

The seed compared `n = 296` with `n = 8, 16, 24`. That could have been a size
effect; it is not. `‖curl‖/‖A‖` is flat in `n` for every pool (8 replicates per
cell, `q = 11`):

| `n` | arithmetic | `loose` | `m2matched` | `marginal` | `sigshuffle` | `maxmatched` |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 0.1116 | 0.1141 | 0.0848 | 0.1129 | 0.0719 | 0.0985 |
| 24 | 0.0935 | 0.1168 | 0.0809 | 0.0801 | 0.0898 | 0.0980 |
| 100 | 0.0921 | 0.1226 | 0.0811 | 0.0912 | 0.0805 | 0.0976 |
| 296 | 0.0908 | 0.1212 | 0.0784 | 0.0931 | 0.0816 | 0.0952 |

(`e1_scaling.csv` runs the `q = 13` version to `n = 698`, same story.) That is
what the noise model predicts: if `A = grad ψ + ε` with edge noise of scale `σ`,
then `‖curl‖/‖A‖ ≈ σ/(2 sd ψ)`, independent of `n` once `n ≫ 1`.

### Larger `q`, and what does set the curl

Sampled pools, `n = 280`, three replicates (`e1_supplement.py`):

| `q` | arithmetic | `loose` | `m2matched` | `marginal` | `sigshuffle` | `maxmatched` |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | 0.0908 | 0.1258 | 0.0811 | 0.0930 | 0.0795 | 0.0962 |
| 13 | 0.0873 | 0.1332 | 0.0809 | 0.0833 | 0.0718 | 0.0937 |
| 17 | 0.0837 | 0.1461 | 0.0792 | 0.0831 | 0.0698 | 0.0924 |
| 19 | 0.0883 | 0.1425 | 0.0730 | 0.0895 | 0.0654 | 0.0979 |
| 23 | 0.0843 | 0.1346 | 0.0737 | 0.0790 | 0.0605 | 0.0871 |
| 29 | 0.0836 | 0.1316 | 0.0683 | 0.0806 | 0.0610 | 0.0933 |

Regressing `‖curl‖/‖A‖` on the spread statistics of the trace vectors
(`mean m₂`, `sd m₂`, `sd α_max`, `sd log log max N`) over **all 108 pools,
arithmetic and control alike**:

```
R² = 0.9309,   residual sd = 0.00592
```

and the arithmetic pools land where an unremarkable pool lands:

| | residuals of the arithmetic pools, in units of the residual sd |
|---|---|
| `q = 11` | `+0.42  +0.30  +0.62` |
| `q = 13` | `+0.31  +0.17  +1.24` |
| `q = 17` | `−0.53  −0.08  −0.28` |
| `q = 19` | `+0.66  +0.41  −0.06` |
| `q = 23` | `−0.24  +0.63  +0.07` |
| `q = 29` | `+0.13  +0.39  +0.78` |
| controls | mean `−0.06`, sd `1.07`, range `−2.45 … +2.95` |

Not one arithmetic pool is beyond `1.3σ`.

### Reconciliation with `curl_on_curve_families.py`

That script measures `r = |curl A|/Σ|A|` and `max|curl A|` on 90 sampled
signatures per field. Its numbers and mine are compatible; the comparable ratio
is `max|curl A|` against the typical edge:

| source | `q` | `n` | `max r` | `max curl` | `median A` | ratio |
|---|---:|---:|---:|---:|---:|---:|
| that script | 101 | 90 | 1.000000 | `3.051·10⁻³` | `2.835·10⁻³` | 1.08 |
| that script | 211 | 90 | 1.000000 | `1.811·10⁻³` | `1.527·10⁻³` | 1.19 |
| here | 11 | 296 | 1.000000 | `2.893·10⁻²` | `1.725·10⁻²` | 1.68 |
| here | 13 | 698 | 1.000000 | `2.058·10⁻²` | `1.499·10⁻²` | 1.37 |

`max r = 1` in every case, which is only the statement that strict 3-cycles
exist — known exhaustively here, and by that script at `q = 101, 211`. Its
controls table (`e1_curl.csv`) reproduces the same ordering as the energy split:
`loose` above arithmetic, `sigshuffle` below.

### E1's answer

> **The curl of the flux is a generic feature of nearly flat signatures.** Given
> a pool of `q`-entry signatures with `Σ N_c = q²` and the trace spread of a
> genus-two pool, `‖curl‖/‖A‖ ≈ 0.06–0.10` whether or not the signatures come
> from curves, at every `q` from 11 to 29; regressing on the spread puts every
> arithmetic pool inside `1.3σ` of the control cloud. The arithmetic is in
> *which* signatures occur, not in how much the comparison curls.
>
> That is the negative half of brief E's dichotomy, and the brief's seed
> evidence for the positive half does not survive a matched control.

---

## E2 — what `ψ_opt` is *(computed)*

### Between classes: the endpoint potential, and it is optimal in its family

Read the table with the caution in the Notation: at fixed `q`, `M`,
`log max_c N_c` and `½ log φ` are three link functions of one integer, so the
honest ceiling for that family is the *categorical* fit — the best function of
`max_c N_c`, whatever it is.

| model | `q=11` | `q=13` | `q=17` | `q=19` | `q=23` | `q=13`, genus 2–4 |
|---|---:|---:|---:|---:|---:|---:|
| `M` | 0.975382 | 0.981648 | 0.984819 | 0.982081 | 0.987751 | 0.983966 |
| `log max N` | 0.988641 | 0.991632 | 0.990231 | 0.987897 | 0.990548 | 0.991483 |
| **`½ log φ`** | **0.990106** | **0.992521** | **0.990561** | **0.988195** | **0.990804** | **0.992193** |
| *ceiling: any function of `max N`* | *0.990225* | *0.992692* | *0.990680* | *0.988357* | *0.990857* | *0.992316* |
| `m₂` alone | 0.423200 | 0.410227 | 0.334310 | 0.344008 | 0.381698 | 0.422778 |
| `M, m₂` (addendum, free fit) | 0.982352 | 0.986143 | 0.990409 | 0.989947 | 0.993137 | 0.989192 |
| `φ̃ = M − 0.0858 m₂` (its own coefficient) | 0.971797 | 0.979254 | 0.982720 | 0.979936 | 0.986300 | 0.981364 |
| `½ log φ, m₂` | 0.995447 | 0.996149 | 0.995839 | 0.995174 | 0.995976 | 0.996217 |
| `½ log φ, m₂, m₃, m₄, log μ` | 0.998062 | 0.998196 | 0.997806 | 0.997918 | 0.998278 | 0.998701 |
| *ceiling: any function of `(max N, μ)`* | *0.995522* | *0.996088* | *0.993501* | *0.992893* | — | *0.995643* |

> **(E2.a)** `ψ_opt` is a function of the largest fiber alone to `R² ≈ 0.99`,
> and **`½ log φ` is within `1.7·10⁻⁴` of the best possible such function** at
> every `q` tested. The potential of the flat-connection regime survives as the
> between-class potential where the connection is not flat.

> **(E2.b)** The remaining `1%` is carried by the moment ladder and by `μ`;
> `½ log φ` plus `m₂` reaches `0.9954`–`0.9960`, the full ladder
> `0.9978`–`0.9987`.

The slope of `ψ_opt` on `½ log φ` is `0.967, 0.973, 0.979, 0.984` at
`q = 11, 13, 17, 19` — approaching 1, as it must if the flat-connection theorem
is the `q → ∞` picture.

### `R²` on the potential badly overstates what a scalar explains

`½ log φ` is a function of an integer taking about `2g√q` values, so it is
*exactly tied* on a sixth of the pairs. The fraction of ordered pairs given the
right sign, ties counted as failures:

| scalar | `q=11` | `q=13` | `q=17` | `q=19` | `q=23` |
|---|---:|---:|---:|---:|---:|
| `M`, `log max N`, `½ log φ` (identical orders) | 0.82494 | 0.83316 | 0.83554 | 0.83187 | 0.85602 |
| `φ̃ = M − 0.0858 m₂` | 0.83889 | 0.84846 | 0.85255 | 0.85089 | 0.87220 |
| `½ log φ` + fitted `m₂` | 0.95692 | 0.95839 | 0.95666 | 0.95467 | 0.96107 |
| `ψ_opt` | **0.98314** | **0.98651** | **0.99044** | **0.99056** | **0.99046** |

and the misses split exactly as brief B's census says. At `q = 11`, `φ` is wrong
on `15286` ordered pairs of which `14798 = 2 × 7399` are exact `φ`-ties and
`488 = 2 × 244` are genuine violations; at `q = 13`, `79296 = 2 × 39648` ties and
`1872 = 2 × 936` violations. **Brief B's census columns are reproduced here by a
third pipeline, to the unit.** So is its `φ̃` audit: of the 43660 pairs at
`q = 11`, `φ̃` is wrong on `7034` of which `1028` are exact `φ̃`-ties, leaving
**6006** decided-and-wrong; it decides **6371** of the 7399 `φ`-blind pairs and
is wrong on **90.4%** of them — brief B's three numbers exactly.

### Inside a `φ`-class, where all the content is

Between classes `φ` decides and is right 99.33% of the time; inside a class it
is exactly tied and, by brief B, not one rate is attained at an endpoint. So the
within-class flux is the whole of what a scalar has to explain
(`e2_within_class.py`, classes with three or more members):

| | `q = 11` | `q = 13` | `q = 101` |
|---|---:|---:|---:|
| classes used | 9 of 11 | 11 of 14 | 19 of 20 |
| within-class `‖curl‖/‖A‖` | **0.09447** | **0.10368** | **0.08534** |
| whole-pool `‖curl‖/‖A‖` | 0.09084 | 0.08517 | 0.07922 |
| `R²` on `m₂` | 0.765253 | 0.732081 | 0.562681 |
| … on `m₂, m₃` | 0.958741 | 0.951275 | 0.829007 |
| … on `m₂, m₃, m₄` | 0.973071 | 0.978202 | 0.932169 |
| … on `m₂, m₃, m₄, log μ` | 0.975369 | 0.983979 | 0.947023 |
| … plus `log min N` | 0.978334 | 0.985301 | — |
| pooled slope of `ψ` on `m₂` | **`+0.00885`** | **`+0.00823`** | **`+0.00419`** |
| "the larger `m₂` precedes" is right on | **9.6%** | **10.9%** | — |

> **(E2.c)** Inside a `φ`-class the potential is the moment ladder
> `m₂, m₃, m₄, log μ`, reaching `R² = 0.95`–`0.98`; the slope on `m₂` is
> **positive**, i.e. *the larger `m₂` follows*, and the addendum's rule that the
> larger `m₂` precedes is right 9.6% of the time at `q = 11` and 10.9% at
> `q = 13`.
>
> The within-class curl fraction is **not** smaller than the whole-pool one — at
> `q = 11` and `13` it is larger. Removing `φ` does not remove the obstruction;
> it exposes it.

The 9.6% and 10.9% agree to the digit with the curve-family session's
independently measured 90.4% and 89.1% wrong. The positive `m₂` slope is the
same finding seen from the regression side.

---

## E3 — does the cycle descend to the curves? *(computed, plus two proofs)*

### E3.1 The arithmetic of the eleven fibers of each certified `F_11` pencil

For a genus-two curve over `F_p`, counting the smooth projective model over
`F_p` and `F_{p²}` determines the L-polynomial and hence, by Tate, the isogeny
class of the Jacobian:

```
L(T) = 1 − s₁T + e₂T² − p s₁T³ + p²T⁴,
s₁ = p+1−#C(F_p),   s₂ = p²+1−#C(F_{p²}),   e₂ = (s₁²−s₂)/2,
```

the Jacobian is isogenous over `F_p` to a product of elliptic curves iff the real
Weil polynomial `x² − s₁x + (e₂−2p)` has integer roots, and the `p`-rank is
`2, 1, 0` according as `p ∤ s₁`, `p | s₁ ∤ e₂`, `p | s₁` and `p | e₂`.

Points at infinity: one for `deg P = 5`, two for monic `deg P = 6`. So **the
trace `a_c = q − N_c` of `m_and_e_and_a_c/FINDINGS.md` equals the genuine
Frobenius trace `s₁` for degree 5 and `s₁ + 1` for degree 6** — a one-unit shift
that matters in E3.3.

Validation of the `F_{121}` arithmetic, independent of the L-polynomial theory:
for 60 random degree-5 pencils and the twisting substitution `x ↦ 2x` (`2⁵` is a
non-residue mod 11), 60/60 have their `F_{121}` fiber counts preserved and 60/60
satisfy `#C′(F_11) = 2(q+1) − #C(F_11)`. Every one of the 33 certified fibers
has `(s₁, e₂)` inside the Honda–Tate/Rück admissible region.

| | `A` (deg 5) | `B` (deg 6) | `C` (deg 6) |
|---|---|---|---|
| signature | `{18,16,15,15,14,12,9,6,6,5,5}` | `{18,18,14,13,12,9,9,9,8,7,4}` | `{19,14,12,11,11,10,10,10,9,9,6}` |
| all fibers smooth | yes | yes | yes |
| `ν(P)` | 33 | 29 | 21 |
| `m₂ = ν(P)/q − 1` | 2.000000 | 1.636364 | 0.909091 |
| split Jacobians | 1 / 11 | 1 / 11 | **3 / 11** |
| `p`-rank | all 2 | all 2 | **three fibers of rank 1** |
| distinct isogeny classes | 11 | 10 | 10 |
| `P′(x)` factors as | one irreducible quartic | `2 + 3` | one irreducible quintic |
| critical values in `F_11` | 0 | 0 | 0 |

The last row *proves* every fiber is smooth: `P(x)+c` is squarefree for all
`c ∈ F_11` exactly when no critical value of `P` is rational. It also confirms
T2.3's identity `m₂ = ν(P)/q − 1` on the nose.

**The branch maps have full monodromy, all three (proved).** `P(x) − t` is
linear in `t`, hence irreducible over `F_q(t)`, so the arithmetic monodromy
group `G` is transitive. The observed Frobenius cycle types force `G = S_d`:

* `A`, `d = 5`: the type `(2,3)` occurs, an element of order 6. The transitive
  subgroups of `S_5` are `C_5, D_5, F_20, A_5, S_5`, with element orders
  `{1,5}, {1,2,5}, {1,2,4,5}, {1,2,3,5}, {1,…,6}`. Only `S_5` has order 6.
* `B`, `d = 6`: the type `(1,5)` makes the point stabiliser transitive on the
  other five points, so `G` is 2-transitive, hence primitive; the type
  `(1,1,1,1,2)` is a transposition; a primitive group containing a transposition
  is symmetric (Jordan). So `G = S_6`.
* `C`, `d = 6`: `(1,5)` again gives primitivity, so `G` is one of `A_5`
  (`= PSL(2,5)`), `S_5` (`= PGL(2,5)`), `A_6`, `S_6`. The type `(1,1,4)` is odd,
  killing `A_5` and `A_6`. In `PGL(2,5)` on `P¹(F_5)` an element of order 3 has
  no fixed point (`3 ∤ 5−1`, so it is non-split and acts as `(3,3)`), so the
  observed type `(1,1,1,3)` rules `PGL(2,5)` out. So `G = S_6`. ∎

**So the three pencils differ by a great deal classically** — degree, `ν(P)`,
the number of split Jacobians, the `p`-rank distribution — and share almost
nothing: `A` and `B` have exactly one isogeny class in common among their 22
fibers, `A` and `C` and `B` and `C` none. What they do *not* differ by is
monodromy: all three branch maps are full symmetric. The cycle is not a cycle
among arithmetically indistinguishable objects; the point is the opposite, that
no ordering by any of those invariants can be consistent, since the comparison
closes a loop.

### E3.2 The map `pencil ↦ signature` and its fibers

The affine group acts on the enumeration: `x ↦ ax + b` and renormalising gives
`Q(x) = a^{−d}(P(ax+b) − P(b))`, and the fiber of the `Q`-pencil over `c` is the
quadratic twist by `a^d` of a fiber of the `P`-pencil. So the subgroup
`G₀ = {(a,b) : a^d is a square}` acts by isomorphisms of the fibration, and the
complement sends the signature to its twist `N_c ↦ 2q − N_c`. Both verified:
the signature is constant on `G₀`-orbits, and the twisting elements do send it
to `2q − N` (14641/14641 degree-5 polynomials at `q = 11`). Counting pencils in
`G₀`-orbits is the honest count.

| `q` | order of `G₀`, deg 5 / deg 6 | pencils (orbits) | signatures | mean fiber | max fiber | injective on |
|---:|---|---:|---:|---:|---:|---|
| 11 | 55 / 110 | 1744 | 296 | 5.89 | 82 | 87 / 296 |
| 13 | 78 / 156 | 2796 | 698 | 4.01 | 117 | 278 / 698 |

**The collapse is not mild**, contrary to what brief E expected from T2.2. The
tension resolves once the two measurements are separated: T2.2 counted distinct
signatures among 400 *randomly sampled* genus-`≥2` fibrations at
`q = 101, 211, 503`, which is a birthday statistic and depends on `q`. Redoing
exactly that measurement for genus-two pencils:

| `q` | 11 | 13 | 17 | 19 | 23 | 31 | 53 | 101 | 211 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| distinct signatures among 400 random pencils | 143 | 217 | 313 | 358 | 381 | 387 | 399 | 399 | 400 |

T2.2's `398–400 of 400` is right at `q ≥ 53` and badly wrong at `q = 11, 13` —
which is where every certified cycle in the programme lives.

The three certified signatures are realised by **2, 2 and 3** pencils.

### E3.3 The signature collapses arithmetic that the curves make

For every signature realised by two or more pencils, compare the multisets of
`(s₁, e₂)` over the eleven (thirteen) fibers:

| `q` | signatures with ≥2 pencils | all pencils isogeny-identical | **pencils with different isogeny multisets** | differing already within one degree |
|---:|---:|---:|---:|---:|
| 11 | 209 | 0 | **209** | 190 |
| 13 | 420 | 1 | **419** | 419 |

The certified triple exhibits it directly:

* `σ(A)` is realised by brief B's degree-5 pencil **and by a degree-6 pencil**,
  `y² = x⁶ + 7x⁴ + x³ + 2x² + c`. Their fibers have *different point counts over
  `F_11`* — the trace multisets differ by the one-unit shift above,
  `(−7,−5,−4,−4,−3,−1,2,5,5,6,6)` against `(−8,−6,−5,−5,−4,−2,1,4,4,5,5)`.
* `σ(B)` is realised by brief B's pencil and by
  `y² = x⁶ + 9x⁵ + 3x⁴ + x³ + x² + x + c`. Same trace multiset (both degree 6),
  but **1 split Jacobian against 4**.
* `σ(C)` is realised by three pencils with 2, 1 and 3 split Jacobians.

Within one degree the signature pins the multiset of `F_q` point counts exactly
and leaves the `F_{q²}` counts free. The extreme case is the flat signature
`(q,…,q)`, which by T2.3 is exactly the permutation-polynomial locus: at
`q = 13` it is realised by `y² = x⁵ + c`, all thirteen of whose fibers are
supersingular with `L(T) = 1 + q²T⁴`, and by at least six other pencils whose
fibers have `s₁ = 0` but `e₂` ranging over `{−22, …, 20}` — every one a
different isogeny class, all with the same flux to every other signature.

> **(E3)** The flux `A(f,g)` is an invariant of the pair of *signatures* and of
> nothing finer. The certified cycle descends to curve families in the strongest
> possible sense — it holds for **all** `2·2·3 = 12` triples of pencils realising
> the three signatures — and in the weakest possible sense, because those pencils
> are not isogenous and the comparison never knew.

### E3.4 Brief C's claim, with the numbers in it

| | `q = 11` | `q = 13` |
|---|---:|---:|
| signatures / unordered pairs / ties | 296 / 43660 / **0** | 698 / 243253 / **0** |
| gradient energy of `A` | **99.175%** | **99.275%** |
| curl energy of `A` | **0.825%** | **0.725%** |
| ordered pairs `ψ_opt` gets right | **98.314%** | **98.651%** |
| ordered pairs `½ log φ` gets right | 82.494% | 83.316% |
| strict 3-cycles | 132 of 4278680 | 1475 of 56434696 |
| edge-disjoint 3-cycles (greedy) | 42 | 214 |
| ⇒ pairs **any** scalar must misorder | ≥ 42 (0.0962%) | ≥ 214 (0.0880%) |

**Lemma (proved).** *A set of `k` pairwise edge-disjoint strict 3-cycles forces
every scalar `ψ` to misorder at least `k` pairs.* Each cycle contains at least
one edge whose direction disagrees with `grad ψ`, and edge-disjointness makes
those `k` edges distinct. ∎

> **The statement to make, and no more.** There is no `φ` with
> `a ≺ b ⟺ φ(a) < φ(b)` on the genus-two pencils over `F_11`, because `A` has
> non-zero curl: any scalar misorders at least 42 of the 43660 pairs. The best
> scalar found, the least-squares potential `ψ_opt`, gets **98.31%** of the
> ordered pairs right and captures **99.18%** of the flux energy; **1.69%** of
> the comparisons and **0.82%** of the energy are irreducibly pairwise. At
> `q = 13` the numbers are 98.65%, 99.28%, 1.35% and 0.72%. The hedge the
> addendum's §7 anticipated — "no scalar beyond the `φ̃` truncation" — is not
> needed: `φ̃` is a *worse* order than `φ` (§E2), so the unqualified statement
> stands.

### E3.5 The `F_101` witness: what separates `f_1` from `f_2`

Brief B's addendum 2 exhibits three pencils over `F_101` in which `f_1` and
`f_2` agree on the largest fiber (123), on `m₂` (`0.851485`) and on `ν(P)`
(187), hence on `φ` and on `φ̃`, so their edge is decided by the interior alone.
The question posed was whether *anything* classical separates them. It does, and
`e3_f101.py` names it:

| | `f_1` | `f_2` | `f_3` |
|---|---:|---:|---:|
| max fiber | 123 | 123 | 122 |
| **multiplicity `μ` of the largest fiber** | **2** | **1** | 3 |
| smallest fiber | 72 | 77 | 71 |
| `m₂` | 0.851485 | 0.851485 | 0.990099 |
| `m₃` | **+0.017645** | **−0.002941** | — |
| `m₄` | +0.022060 | +0.019936 | — |
| smooth fibers | **100 / 101** | 101 / 101 | 101 / 101 |
| split Jacobians | 14 / 101 | 9 / 101 | 6 / 101 |
| `p`-rank 1 fibers | 6 | 4 | 3 |
| `P′(x)` factors as | `1 + 3` | `2 + 2` | irreducible quartic |
| isogeny classes shared by `f_1` and `f_2` | 4 of 101 | | |

So `f_1` and `f_2` are tied on exactly three statistics and differ on everything
below them in the ladder: `μ`, `m₃` (opposite signs), `m₄`, the smallest fiber,
the ramification type of the branch map, the split count, the `p`-rank
distribution, and — a fact the addendum does not mention — **`f_1` has a
singular fiber and `f_2` does not.** `P′_1` has a rational root, so one critical
value of `P_1` lies in `F_101` and the fiber over it is a singular curve. The
`F_101` cycle is therefore *not* an all-smooth witness, unlike the `F_11` one.

Out-of-sample test of (E2.c). Fitting the within-class potential on a
702-signature `q = 101` pool and then predicting the `f_1`/`f_2` edge, which the
fit does not see as a special pair:

| within-class model | predicted `ψ(f_2) − ψ(f_1)` | edge |
|---|---:|---|
| `m₂` | `0` (exactly tied) | wrong |
| **`log μ`** | `−4.452·10⁻⁴` | **right** |
| `m₃` | `+9.768·10⁻⁵` | wrong |
| `m₂, m₃` | `+4.601·10⁻⁴` | wrong |
| `m₂, m₃, m₄` | `+4.392·10⁻⁴` | wrong |
| `m₂, m₃, m₄, log μ` | `+2.116·10⁻⁴` | wrong |

(the truth is `A(f_1,f_2) = −1.6197·10⁻⁴`, i.e. `f_2 ≺ f_1`.)

> **(E3.5)** Something classical *does* separate `f_1` from `f_2`, and it is the
> multiplicity of the largest fiber — addendum 1's §4(b) term, the one its own
> §3 scalar `φ̃` omits. On this pair the multiplicity alone predicts the edge
> correctly and the moment ladder does not. So the `F_101` witness is not
> evidence that the flux reads something beyond the classical statistics; it is
> evidence that the ladder needs one more rung than `φ̃` has, exactly where
> addendum 1 said it would. **One pair is one bit; this is a spot check, not a
> measurement.**

---

## Corrections

1. **Brief E's E1 conclusion is withdrawn.** Its seed reads "the arithmetic pool
   has about twice the curl of a random one — `0.091` against `0.041–0.051`".
   That compared `q`-entry curve signatures against 2-to-7-entry random
   integers. All five matched controls put the arithmetic pool in the middle of
   its own range, the tightest control exceeds it, and over 108 pools the curl
   fraction is a function of the trace spread to `R² = 0.93` with every
   arithmetic pool inside `1.3σ`. There is no evidence that the flux is
   arithmetic in this sense.

2. **Brief E's E2 conclusion is half right, and its refutation of the addendum
   is stated wrongly.** `½ log φ` does survive as the between-class potential,
   and it is essentially the *optimal* function of the largest fiber — which the
   brief did not check and which is the sharper statement. But:
   * `R² = 0.9901` for `½ log φ` translates into only **82.5%** of ordered pairs
     given the right sign, against `ψ_opt`'s 98.3%. The `1%` of variance the
     brief dismisses is the entire content of the comparison inside a `φ`-class.
   * The brief refutes the addendum by comparing free fits of `(M, m₂)` at
     `R² = 0.982` against `0.990`. That comparison **fails** at `q = 17, 19, 23`,
     where the free fit of `(M, m₂)` scores `0.9904, 0.9899, 0.9931` against
     `½ log φ`'s `0.9906, 0.9882, 0.9908`. The refutation that holds at every `q`
     is about the coefficient, not the fit: `φ̃` with its own stated coefficient
     scores *below `M` alone* at every `q` tested, and the freely fitted ratio of
     the `m₂` to the `M` coefficient is `+0.40, +0.38, +0.52, +0.71, +0.72` — the
     **opposite sign** to the addendum's `−0.0858`. The within-class slope on
     `m₂` is positive too: `+0.0089`, `+0.0082`, `+0.0042` at `q = 11, 13, 101`.

3. **T2.2's "genus ≥ 2 gives 398–400 of 400 distinct signatures" does not
   transfer to `q = 11, 13`.** Redoing that measurement gives 143 and 217 of
   400. The statement is `q`-dependent and holds only for `q ≳ 53`. Brief E's
   expectation that "the collapse should be mild" at `q = 11, 13` is wrong: the
   exhaustive count is 1744 pencils onto 296 signatures, mean fiber 5.89.

4. **The `F_101` witness of brief B's addendum 2 is not an all-smooth cycle.**
   `P′_1` has a rational root over `F_101`, so one critical value of `P_1` is
   rational and one of the 101 fibers of `f_1` is singular (100/101 smooth;
   `f_2` and `f_3` are 101/101). The `F_11` cycle of
   `curve_family_cycles/FINDINGS.md` remains the all-smooth witness. Also, the
   addendum's phrasing that `f_1` and `f_2` are "exactly tied" should be read as
   tied *on `M`, `m₂` and `ν(P)`*: they differ on `μ` (2 vs 1), on `m₃` (opposite
   signs), on `m₄`, on the smallest fiber, and on every isogeny statistic.

5. **Sign convention in the seed script.** `hodge_split.py` computes
   `A(a,b) = mid(u_a − u_b)`, the negative of brief D Part 0's convention, so its
   `psi` is `−ψ_opt`. Every statistic it prints is invariant under the flip, so
   its numbers stand unchanged; `flux.py` here uses brief D's signs.

---

## Open

* Whether the within-class potential's `m₂, m₃, m₄, log μ` ladder has a closed
  form. Brief B's `β ≍ √q` scaling says the right variable is `Ψ(τ) = Λ(τ)/τ`
  and the within-class flux should be a midrange of `Ψ_a − Ψ_b`; nobody has
  written that expansion, and (E2.c) says it would explain 95–98% of the
  within-class potential. The `F_101` spot check (E3.5) says the `μ` term is not
  yet correctly weighted relative to the moments.
* Whether the curl fraction has a limit as `q → ∞`. It sits in `0.083–0.091` for
  the arithmetic pools at `q = 11 … 29` with no visible trend, but the pools at
  `q ≥ 17` are samples and the `m2matched` and `sigshuffle` controls drift
  downwards over the same range.
* Whether any *pairwise* classical invariant correlates with the curl. E1 rules
  out the one-family statistics; brief E's own argument was that `A` is a
  function of a pair and no classical invariant of pairs is on offer. None was
  found here, and none was systematically searched for.
* **The programme's arithmetic claim.** E1 says the flux correlates with nothing
  beyond the trace spread; E3 says it *cannot* correlate with anything finer,
  because it is a function of the signature, which at `q = 11` merges on average
  5.89 pencils with different isogeny data. What survives is (E2.c): the
  within-class comparison is a function of `m₂, m₃, m₄, μ` to `R² ≈ 0.95–0.98`,
  i.e. of the trace moments, which are classical. **On the present evidence the
  honest position is that the exchange comparison of curve families is a
  functional of the trace moments plus an irreducible `0.7–0.8%` of pairwise
  structure, and there is no arithmetic in it beyond the moments — but that
  fraction is provably not a functional of anything scalar, and it is what the
  cycles are made of.**

---

## Files

| file | what |
|---|---|
| `pools.py` | third independent enumeration of the genus-two pools; the five matched random controls |
| `flux.py` | `A = mid(u_b − u_a)`, HodgeRank split, cycle counting, verification against the package |
| `curves.py` | point counts over `F_p` and `F_{p²}`, L-polynomials, splitting, `p`-rank, factorisation and ramification types |
| `orbits.py` | the affine action on pencils and its `G₀`-orbits |
| `e1_hodge.py` | E1: smoke test, full pools, controls, `n`-scaling → `e1_smoke.csv`, `e1_pools.csv`, `e1_scaling.csv` |
| `e1_supplement.py` | E1 at `q = 11 … 29`, and the regression of curl on spread → `e1_supplement.csv` |
| `e1_curl.py` | the triangle statistic `r`, for reconciliation → `e1_curl.csv` |
| `e2_potential.py` | E2 regressions, ceilings, order agreements → `e2_regressions.csv`, `e2_potential.csv` |
| `e2_within_class.py` | E2 inside a `φ`-class → `e2_within_class.csv` |
| `e2_f101_class.py` | the out-of-sample test of E2.c on the `F_101` witness → `e2_f101_class.csv` |
| `e3_curves.py` | E3.1–E3.3 → `e3_certified_fibers.csv`, `e3_signature_fibers.csv`, `e3_collisions.csv` |
| `e3_brief_c.py` | E3.4, with the certified lower bound → `e3_brief_c.csv` |
| `e3_f101.py` | E3.5, the arithmetic of the `F_101` witness → `e3_f101.csv` |
| `verify_f101.py` | independent re-verification of the `F_101` cycle |
| `build_f11_pool.py`, `hodge_split.py` | brief E's seed scripts, unmodified |
| `*_output.txt` | the console transcript of each run |
