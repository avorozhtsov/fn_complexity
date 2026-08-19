# Findings — what can an exchange geometry be?

Answer to session brief G. Authored from the reporting agent's verified output;
`structure.py` and `g3_bound.py` were re-run independently before committing.

## Summary

Neither branch of the brief's dichotomy is right. **The framework is rigid, and
the rigidity is an exact structure theorem — but it does not forbid tournaments;
it quantifies them, and it does appear to forbid metrics.**

Writing `R = log r`, `Λ = log(max fiber)`, `σ = log(R/Λ)`, `s = log β`:

```
u_a(e^s) = log Λ_a + max(σ_a, s) + w_a(s),   0 ≤ w_a(s) ≤ log(1 + e^{−|s−σ_a|})
```

with `w_a` **unimodal, peaked exactly at `s = σ_a`, and 1-Lipschitz**. Every `u`
is a kink plus a single bump of height ≤ `log 2`. Hence

```
d(a,b) = |σ_b − σ_a| + ε,      0 ≤ ε ≤ 2 log(1 + e^{−|Δσ|})
A(a,b) = ψ(b) − ψ(a) + D,     |D| ≤ ½ log(1 + e^{−|Δσ|}) ≤ (log 2)/2 = 0.34657
```

with `ψ = ½ log(log r · log M)` the brief-D(c) potential. Both bounds are
independent of `n` and of everything about the signatures except `σ`. The
constant `(log 2)/2` is **sharp**.

* **G1 — every tournament on `n ≤ 6` is realisable** (2 + 4 + 12 + 56 = 74
  isomorphism classes, all certified). Two fibers already suffice for a 3-cycle.
* **G2 — the metrics do *not* fill MET.** Only the uniform 4-point metric was
  realised (distortion 1.0002); `C_4`, `C_5`, `K_{2,3}`, `K_{3,3}`, `K_6`, `K_8`,
  Petersen and Wagner `V_8` all resist, best distortions **1.256 to 3.58**.
  Proved: every exchange metric is within additive `2 log 2 = 1.386` of a
  **one-dimensional** metric, uniformly in `n`.
* **G3 — `‖curl‖/‖A‖` is *not* bounded away from 1** (attained exactly at
  `n = 3,4,5,6`), but there is a sharp *inequality*, which is the theoretical
  half of brief E's `R² = 0.93` measurement and of brief F's `q → ∞` constant.

---

## Notation

Conventions of `PLAN.md` and brief D Part 0. `Z_a(β) = Σ a_i^β`, `F_a = log Z_a`,
`u_a = log F_a`, `L(a,b) = −log C(a→b) = sup_β(u_b−u_a)`, `d = osc`, `A = mid`,
`a ≺ b ⟺ A(a,b) > 0`. Work in `s = log β`, `U_a(s) = u_a(e^s)`.
`R_a = log r_a = F_a(0)`, `Λ_a = log M_a = lim F_a'`, `τ = R/Λ`, `σ = log τ`,
`ψ = ½ log(RΛ) = ½ log φ`. For `φ = U_b − U_a`: `P = max φ − max(φ(±∞))`,
`Q = min(φ(±∞)) − min φ`, `ε = P+Q`, `D = (P−Q)/2`. Signatures have `r ≥ 2`
fibers and `M ≥ 2`; `(1,)` and all-ones signatures are excluded throughout.

**`σ` and `ψ` depend on a signature only through `(r, M)`.** Fixing `(r,M)` makes
`A = D` exactly. This is used constantly.

---

## 1. The structure theorem *(proved; verified to 3.6e−15)*

### 1.1 The sandwich

**Lemma.** `max(R, βΛ) ≤ F(β) ≤ R + βΛ` for all `β ≥ 0`.

*Proof.* `Z ≥ M^β` gives `F ≥ βΛ`; `F` non-decreasing gives `F ≥ F(0) = R`;
`Z ≤ r M^β` gives `F ≤ R + βΛ`. ∎

The two bounds differ by exactly the factor `1 + e^{−|s−σ|}`.

### 1.2 Kink plus bump

**Theorem A.** With `w_a(s) = U_a(s) − log Λ_a − max(σ_a, s)`:

1. `0 ≤ w_a ≤ log(1 + e^{−|s−σ_a|})`, `w_a(±∞) = 0`;
2. `U_a'(s) ∈ (0,1)`, so `w_a` is 1-Lipschitz in `s`;
3. `w_a` strictly increases on `s < σ_a`, strictly decreases on `s > σ_a` —
   **unimodal, peak exactly at `σ_a`** — of height `log(F_a(τ_a)/R_a) ∈ (0, log 2]`,
   `= log 2` **iff `a` is flat**.

*Proof.* (1) is §1.1 divided by `max(R,βΛ)`. (2) `U'(s) = βF'/F`; `F,F' > 0` gives
`U' > 0`; `F'(β) = Σa_i^β log a_i / Σa_i^β ≤ Λ` with `F ≥ βΛ` gives `U' ≤ 1`;
`w' = U' − 1_{s>σ}`. (3) For `s<σ`, `max(R,βΛ)=R` so `w' = U' > 0`; for `s>σ`,
`w' = U'−1 < 0` since `F > βΛ` strictly when `r ≥ 2`. Peak height `log(F(τ)/R)`,
and `F(τ) ≤ R + τΛ = 2R` with equality iff `Σa_i^τ = rM^τ`, i.e. iff every
`a_i = M`. ∎

*Verified* (`structure.py`) on 340 signatures — random integer pools to 20 fibers
and values `10⁵`, plus an adversarial thin/flat pool with up to `10⁸` fibers:
worst violation `3.6e−15`, unimodality `1.8e−12`, `|w(σ) − max w|` at grid
resolution `5e−4`.

*Computed, not proved:* `U` is **convex** in `s` on every pool tested,
`max(−U'') = +3.6e−9` against second-difference noise `~1e−6` (`g2_metrics.py`).

### 1.3 The two halves

**Corollary A1.** `φ(−∞) = log R_b − log R_a`, `φ(+∞) = log Λ_b − log Λ_a`, so
`φ(−∞) − φ(+∞) = σ_b − σ_a` and `½(φ(−∞)+φ(+∞)) = ψ_b − ψ_a`. Hence
`d = |σ_b−σ_a| + P + Q` and `A = (ψ_b−ψ_a) + (P−Q)/2`.

**Corollary A2.** With `Δ = |σ_b − σ_a|`: `P, Q ≤ log(1+e^{−Δ})`, so
`ε ≤ 2log(1+e^{−Δ}) ≤ 2log 2` and `|D| ≤ ½log(1+e^{−Δ}) ≤ (log 2)/2`.

*Proof.* By §1.1, `φ ≤ m(s) + E_b(s)` with `m = log max(R_b,βΛ_b) −
log max(R_a,βΛ_a)` monotone with endpoint values `φ(±∞)`, and
`E_b(s) = log(1+e^{−|s−σ_b|})`. Take `σ_b ≥ σ_a`. On `s ≤ σ_a`: `m = φ(−∞)`,
`E_b ≤ log(1+e^{−Δ})`. On `σ_a ≤ s ≤ σ_b`, with `t = σ_b − s`:
`m + E_b = φ(+∞) + log(1+e^t)`, increasing, sup `φ(−∞) + log(1+e^{−Δ})` at
`t = Δ`. On `s ≥ σ_b`: `≤ φ(+∞) + log 2`, smaller since
`2e^{−Δ} ≤ 1+e^{−Δ}`. So `P ≤ log(1+e^{−Δ})`; apply to `−φ` for `Q`. ∎

*Verified:* over 12 000 pairs `d = |Δσ| + P + Q` to `8.9e−16` with `P,Q ≥ 0`
exactly; `|D| ≤ ½log(1+e^{−Δ})` violated by at most `2.2e−16` on every pool in
`g3_bound.csv`, including brief E's arithmetic pools.

### 1.4 `(log 2)/2` is sharp *(proved and computed)*

Take `a_r = (r, 1, …, 1)` with `r` fibers (`R = Λ = log r`, `τ = 1`, `σ = 0`) and
`b` any flat signature with `σ_b = 0`. Then `w_b = log(1+e^{−|s|})` exactly (peak
`log 2`), while `0 ≤ w_{a_r} ≤ log(1 + log 2/log r) → 0` because
`F_{a_r}(β) = log(r^β + r − 1) ≤ max(R,βΛ) + log 2`. Since `σ_a = σ_b`, `φ` is a
constant plus `w_b − w_a`, so `P → log 2`, `Q → 0`, `D → (log 2)/2`.

| `r` | `P` | `Q` | `D` |
|---:|---:|---:|---:|
| `10` | 0.45535103 | 0 | 0.22767551 |
| `10³` | 0.59823756 | 0 | 0.29911878 |
| `10⁶` | 0.64428005 | 0 | 0.32214002 |
| `10¹²` | 0.66838193 | 0 | 0.33419097 |
| `10⁴⁰` | 0.68564991 | 0 | 0.34282496 |
| `10⁸⁰` | 0.68939141 | 0 | 0.34469570 |
| `→ ∞` | `log 2 = 0.69314718` | 0 | `(log2)/2 = 0.34657359` |

An independent hill-climb over free multiplicities (`extremes.py`) reaches
`max|D| = 0.34273769` = **98.89 %** of the bound, at exactly this shape. The same
climb gives `sup ε = 0.68547538 = 0.9889·log 2`, always with `Q = 0` — so the
true supremum of `ε` looks like `log(1+e^{−Δ})`, half the proved bound (see Open).

### 1.5 Two exactly-exact families *(proved)*

* **One fiber.** `u_{(a)}(e^s) = s + log log a`, so `u_b − u_a` is constant:
  `d ≡ 0`, `A = dψ` exactly, total order.
* **Flat signatures.** For two flat signatures
  `φ = log((R_b+βΛ_b)/(R_a+βΛ_a))` is a ratio of affine functions, hence
  monotone: `P = Q = D = 0`. **A pool of flat signatures has zero curl and can
  never contain a cycle.** Confirmed: the `k = 1` rows of `extremes.py` return
  `max|D| = 0` and `max|curl| = 0` after a full hill-climb. This generalises the
  endpoint-regime theorem of brief B / D(c): it is not only the endpoint regime
  that is curl-free, it is the whole flat locus.

### 1.6 What the Cartesian power does *(proved; verified to 1.1e−13)*

`a^{⊗k}` has `r^k` fibers and largest fiber `M^k`, so `R → kR`, `Λ → kΛ`:

```
σ unchanged,  w unchanged,  ψ → ψ + log k,  d unchanged,  D unchanged
A(a^{⊗j}, b^{⊗k}) = A(a,b) + log(k/j)   exactly.
```

**The Cartesian power acts on the potential and on nothing else.** This refines
FINDINGS T1.5 Theorem B and is the exact statement of why `d` is the invariant of
the three homogeneity degrees recorded in brief H. Verified: `σ` fixed to `0`,
`ψ` shift equal to `log k` to `1e−15`, `d` and `A − log(k/j)` fixed to `1.1e−13`.

---

## 2. G1 — which tournaments are realisable

### 2.1 The seed, reproduced exactly

`tournament_seed.py`, `random.seed(11)`: pool 298; `(0,1,2,3) → 1499`,
`(0,2,2,2) → 1`, other two classes `0`; Hodge `|grad|/|A| = 0.9987, 0.9987,
0.9992`, residual `0.0511, 0.0510, 0.0407`, **zero** 3-cycles at `n = 8, 16, 24`.

Random pools are almost transitive for a now-provable reason: `|D| ≤ (log 2)/2`
while the `ψ`-differences in such a pool have rms `0.31`, and every
intransitivity must sit inside a `ψ`-window of width `< log 2` (§2.6).

### 2.2 Where the cycles live

Fixing `(r, M)` kills `dψ`, so `A = D` and the comparison is pure defect. These
pools are cycle-rich (`g1_mine.py`):

| pool | signatures | 3-cycles | triples |
|---|---:|---:|---:|
| `r=6, M=12` | 700 | 6 802 | 56 921 900 |
| `r=5, M=12` | 700 | 5 012 | 56 921 900 |
| `r=7, M=12` | 700 | 6 262 | 56 921 900 |
| `r=6, M=20` | 700 | 6 262 | 56 921 900 |
| random 2–7 fibers, values ≤ 40 | 24 | 0 | 2 024 |

**Recipe: to build a prescribed comparison, hold `(r, M)` fixed and vary the
interior.**

### 2.3 The exact realisability criterion *(proved)*

**Theorem D.** Fix signatures `a_1,…,a_n` with defect matrix `D`. Replacing `a_i`
by `a_i^{⊗k_i}` changes `A_ij` to `D_ij + (ψ_j + log k_j) − (ψ_i + log k_i)` and
leaves `d` and `D` untouched (§1.6). Writing `c_i` for the free potential
(realisable to arbitrary precision since `{log(k/j)}` is dense in `ℝ`), the system
`sign(c_j − c_i + D_ij) = T_ij` is a system of difference constraints along the
edges of `T`. By Bellman–Ford feasibility it is solvable with margin `m` iff every
**directed cycle** of `T` has `Σ_{e∈cycle} D_e ≥ m·(length)`. Hence

> **A tournament `T` is realisable from a defect matrix `D` iff every directed
> cycle of `T` has positive `D`-sum, and the largest achievable margin is exactly
> the minimum mean cycle weight of `D` over `T`'s directed cycles.**

Since `Σ_{cycle} D_e = curl A` around that cycle, **realisability of tournaments
= realisability of curl patterns**, made exact. It also gives an `O(n³)` test in
place of a search.

### 2.4 Realisation, certified

Using two-fiber bases (whose Cartesian powers have only `k+1` distinct values, so
they stay certifiable) plus the powers of Theorem D (`g1_potential.py`):

| `n` | classes | realised | certified min abs A, range | powers `k` used |
|---:|---:|---:|---|---|
| 3 | 2 | **2/2** | ≥ 1.4e−3 | ≤ 250 |
| 4 | 4 | **4/4** | ≥ 1.75e−5 | ≤ 250 |
| 5 | 12 | **12/12** | 3.7e−3 … 1.8e−2 | ≤ 250 |
| 6 | 56 | **56/56** | 1.7e−3 … 1.4e−2 | ≤ 250 |

Every witness re-verified at 40 digits (`|double − mpmath| ≤ 3e−16`), and the
isomorphism class of the certified matrix re-canonicalised and confirmed. The
`n = 5` doubly-regular tournament (scores `2-2-2-2-2`) is realised with certified
margin `5.277e−3`. Independently, focused `(r,M)` pools alone (no powers) already
give 2/2, 4/4 and 9/12 at `n = 3,4,5` with margins `1.1e−4` to `3.8e−3`.

**Answer to G1: no obstruction was found, and none is expected.**

### 2.5 A certified 3-cycle among **two-fiber** signatures

A dense sweep of the whole two-fiber family (3600 signatures, `x ∈ [0.25,4]`,
`7.77·10⁹` triples) contains **2493** 3-cycles. Refined and integerised by the
exact power symmetry:

```
a = (5180584820303422554112, 2234356748607932858368)
b = (31539938692971208223227904, 1)
c = (405578209186060584353792, 1000000000000)

A(a,b) = +0.001382452670    A(b,c) = +0.001382452701    A(c,a) = +0.001382452672
directed 3-cycle;  margin 1.3825e−3;  |curl| = Σ|A| = 0.004147358
σ = (−4.278520, −4.439179, −4.362134),  spread 0.1607
ψ = ( 1.772747,  1.853077,  1.814554),  spread 0.0803   (bound log 2 = 0.6931)
mean |A| = 0.0013825                                    (bound (log2)/2 = 0.3466)
```

mpmath agreement `7.3e−17`.

### 2.6 What is proved impossible

* One-fiber and flat families: total orders (§1.5).
* **Every 3-cycle has `ψ`-spread `< log 2`.** On a directed cycle each edge has
  `ψ_{i+1} − ψ_i + D_i > 0` with `|D_i| ≤ (log 2)/2`, and the `ψ`-increments sum
  to zero, so the total descent is under `(#descending edges)·(log 2)/2`; a
  3-cycle has at most two descending edges. The known cycle
  `{(6,3,3),(7,2,1),(6,5,1)}` has `ψ`-spread `0.0413`. Since `ψ` is a function of
  `(r,M)` alone, this is the exact form of brief B/E's "cycles live inside
  `φ`-classes", and the analogue of brief F's endpoint-level LP argument.

---

## 3. G3 — how much curl is possible

### 3.1 The inequality *(proved)*

`curl A = curl D` around every cycle, because `dψ` is exact; and the HodgeRank
gradient projection fixes gradients, so `‖A − grad(A)‖_F = ‖D − grad(D)‖_F ≤ ‖D‖_F`.
Entrywise,

```
(B1)  |D(a,b)| ≤ ½ (d(a,b) − |σ_a − σ_b|)          exact, computable from d and σ
(B2)  |D(a,b)| ≤ ½ log(1 + e^{−|σ_a−σ_b|})         a priori, from σ alone
(B3)  |D(a,b)| ≤ (log 2)/2                          universal
```

so `‖curl A‖/‖A‖ ≤ ‖bound‖_F/‖A‖_F` in each case. **There is no `n` anywhere —
the bounds are per-edge, which is exactly why brief E measured the statistic to be
flat in `n` from 8 to 698.** In the form of (B3):
`‖curl A‖/‖A‖ ≤ (log 2)/2 / rms|A|`, and `rms|A|` is dominated by the spread of
`ψ` — the "function of the spread" brief E measured to `R² = 0.93`.

### 3.2 How tight *(computed, `g3_bound.py`)*

| pool | `n` | `‖curl‖/‖A‖` | `B1` | `B1`/actual | rms abs A | `sd(σ)` |
|---|---:|---:|---:|---:|---:|---:|
| random short | 8 | 0.0240 | 0.0338 | 1.41 | 0.334 | 0.342 |
| random short | 16 | 0.0383 | 0.0655 | 1.71 | 0.160 | 0.216 |
| random short | 24 | 0.0309 | 0.0389 | 1.26 | 0.269 | 0.303 |
| random short | 48 | 0.0287 | 0.0413 | 1.44 | 0.323 | 0.327 |
| random short | 96 | 0.0387 | 0.0575 | 1.49 | 0.261 | 0.293 |
| arithmetic `q=11` | 60 | 0.0780 | 0.1567 | 2.01 | 0.0256 | 0.0366 |
| control `q=11` | 60 | 0.1031 | 0.2270 | 2.20 | 0.0217 | 0.0313 |
| arithmetic `q=13` | 60 | 0.0913 | 0.1752 | 1.92 | 0.0199 | 0.0284 |
| control `q=13` | 60 | 0.0881 | 0.1599 | 1.81 | 0.0178 | 0.0257 |
| arithmetic `q=17` | 60 | 0.0876 | 0.1680 | 1.92 | 0.0174 | 0.0248 |
| control `q=17` | 60 | 0.0898 | 0.1608 | 1.79 | 0.0154 | 0.0220 |
| bucket `r=6,M=12` | 60 | 0.0636 | 1.0400 | 16.4 | 0.0428 | 0 |

`B1` is tight to a factor **1.3–2.2** on every pool with any `σ`-spread,
reproduces the arithmetic/control ordering, and is itself **flat in `q`** (0.157,
0.175, 0.168 at `q = 11,13,17`) exactly as brief F's `q → ∞` constant requires. It
is vacuous only when `σ` is constant across the pool — the `(r,M)`-buckets — which
is the regime where the curl is maximal. `B2`/`B3` are 8–22 there and useless;
**`B1` is the bound to quote.** It also explains the two regimes: `q`-entry
near-flat signatures have `σ` almost constant, so `rms|A|` is ~15× smaller than
for short random signatures while the defect budget is unchanged — hence `≈ 0.09`
versus `≈ 0.04`.

### 3.3 The ratio is **not** bounded away from 1 *(computed, `g3_hodge.py`)*

With a hard floor `min|A| > 1e−6` (four orders above the tie threshold), so no
degenerate `A ≈ 0` families are accepted:

| `n` | `r` | max ratio | `‖A‖` | min abs A | 3-cycles | `ψ` spread |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 3 | **1.000000** | 6.23e−3 | 2.54e−3 | 1/1 | 8.9e−3 |
| 4 | 4 | **1.000000** | 1.81e−2 | 1.81e−4 | 2/4 | 3.0e−2 |
| 5 | 5 | **1.000000** | 1.73e−2 | 5.63e−4 | 4/10 | 4.9e−2 |
| 6 | 6 | **1.000000** | 2.15e−2 | 2.67e−5 | 6/20 | 3.2e−2 |

Purely rotational flows — with no potential component at all — are attained.
**The rigidity is the inequality of §3.1, not a universal constant.**

### 3.4 The strength of cyclic preference *(proved)*

Around a directed cycle all `A`'s share a sign, so `|curl A| = Σ|A|`, and
Corollary A2 gives
`(1/k) Σ|A(a_i,a_{i+1})| ≤ (1/k) Σ ½log(1+e^{−|Δσ_i|}) ≤ (log 2)/2`.
Since `2A(a,b) = log(C(b→a)/C(a→b))`:

> **The geometric mean of the rate asymmetry `C(b→a)/C(a→b)` around any
> preference cycle is at most 2.**

Cycles exist in abundance, but they are necessarily weak: the strong preferences
in the matrix are always the transitive ones. This strictly improves brief D(d)'s
`Σ S ≥ |Σ A|` whenever the cycle's metric perimeter exceeds `k log 2 / 2`. On the
known cycle: `|curl A| = 0.015814`, `Σ S = 0.201398` (slack ×12.7), sharpened
`Σ S − spread(σ) = 0.118896` (slack ×7.5), universal `3(log2)/2 = 1.039721`.

### 3.5 The largest triangle curl *(computed)*

Hill-climbing over free multiplicities gives
`max |curl A| = 0.34273769 = 0.9889·(log 2)/2`, attained where one edge carries
the extremal defect and the other two essentially none. The proved bound is
`3(log2)/2 = 1.0397`; the numerics say the truth is `(log 2)/2`. Bounded-atom
climbs give `0.1288, 0.1653, 0.1848, 0.1976, 0.2067` for `r = 2,…,6`.

---

## 4. G2 — which metrics are realisable

### 4.1 The oscillation count

Over 5310 pairs from three random integer pools (to 20 fibers, values `10⁵`),
interior local extrema of `φ = u_b − u_a`:

| extrema | 0 | 1 | 2 | 3 | 4 | ≥5 |
|---|---:|---:|---:|---:|---:|---:|
| 2–7 fibers, ≤40 | 756 | 794 | 202 | 16 | 1 | 0 |
| 2–12 fibers, ≤400 | 820 | 630 | 291 | 25 | 4 | 0 |
| 2–20 fibers, ≤10⁵ | 1027 | 589 | 139 | 13 | 2 | 0 |

**Never more than four, usually zero or one, even with forty atoms in play** — so
at most six values of `β` (four extrema plus two endpoints) at which the sup-norm
can be active per pair. Mechanism: `φ' = U_b' − U_a'` is a difference of two
increasing maps of `(0,1)` onto itself (numerically `U` is convex in `s`), of total
variation ≤ 2, changing sign only where the two Gibbs measures swap which
atom-scale dominates.

### 4.2 The metric rigidity statement *(proved)*

With `ℓ_ab = |σ_a − σ_b|`, a **line metric** (isometric in `ℝ`, hence hypermetric
and of negative type):

```
ℓ(a,b) ≤ d(a,b) ≤ ℓ(a,b) + 2 log(1 + e^{−ℓ(a,b)}) ≤ ℓ(a,b) + 2 log 2 = ℓ + 1.3863 .
```

**Every exchange metric on any number of points is a one-dimensional metric plus
an additive correction bounded by `1.3863`.** Corollaries: `a ↦ σ_a` embeds
`(X,d)` in `ℝ` with distortion ≤ `1 + 2log(1+e^{−t})/t`, `t = min ℓ`, so
`ℓ₂`-distortion `→ 1` exponentially as the `σ`-separation grows; and every
negative-type failure, every pentagonal violation and all of the distortion above
1 in T1.1–T1.3 live in a window of diameter `2 log 2` and **cannot be amplified by
adding points**. That is why T1.3 saw `c₂ ≈ 1.1` with `c₂/log n` *falling*:
Bourgain's `O(log n)` is unreachable because the non-line budget does not grow
with `n`.

### 4.3 Realisation up to scale — the answer is **no** *(computed)*

Minimising the scale-free distortion `max(d/δ)/min(d/δ)` over `n` signatures with
`r` fibers each (`g2_metrics.py`; per target: differential evolution with
population `min(20·nr, 120)` for 250 generations plus three compass-search
polishes, two restarts, at `r = 3, 4, 6`; the parameter count `nr − 1` is 11–29
against 5–44 essential distance ratios, so the search is not dimension-starved):

| target | `n` | best distortion | realised scale `d/δ` |
|---|---:|---:|---:|
| uniform `K_4` | 4 | **1.000207** | 0.066 |
| uniform `K_6` | 6 | 1.450031 | 0.054 |
| uniform `K_8` | 8 | 1.977406 | 0.089 |
| `C_4` | 4 | 1.255692 | 0.148 |
| `C_5` | 5 | 1.272589 | 0.050 |
| `K_{2,3}` | 5 | 1.546122 | 0.062 |
| `K_{3,3}` | 6 | 2.138484 | 0.073 |
| Petersen | 10 | 3.578379 | 0.084 |
| Wagner `V_8` | 8 | 2.737475 | 0.081 |

**Only the 4-point equilateral metric is realisable; nothing else tried is, even
up to scale.** The pattern is coherent with §4.1: an `n`-point equilateral set in
`ℓ∞^k` needs `n ≤ 2^k`, and the effective `ℓ∞`-dimension available here is the
number of active `β`-points, `≈ 2–3` — so `K_4` fits and `K_6`, `K_8` do not, and
the graph metrics, which need genuinely two-dimensional structure, fail at `1.26`
upward. **This is the one place where realisability looks genuinely obstructed,
and it is where the framework's content is.** *Computed, not proved:* there is no
proof that `C_4` is unrealisable, only a 24-parameter search that stalls at
`1.256` at three independent `r`.

---

## Corrections

1. **The brief's dimension count for G1 is wrong.** Two fibers already suffice for
   a 3-cycle (§2.5), below the predicted `r ≳ (n−1)/2` threshold; sign conditions
   are open, so a dimension count neither guarantees nor forbids them.
2. **The brief's dichotomy is false.** G1 comes out universal (for `n ≤ 6`) *and*
   the framework has content, because the content is quantitative (§1.3, §3.1,
   §4.2) — and separately G2 comes out **not** universal.
3. **The seed script's curl residual is pool-dependent and the dependence is now
   explained**: it equals `‖½(d−ℓ)‖/‖A‖` to within a factor of 2. The brief's
   numbers (`0.0511, 0.0510, 0.0407`) reproduce exactly from `random.seed(11)`; a
   different subset draw from the same 298-signature pool gives `0.0240, 0.0383,
   0.0309`. Quote the script, not a number.
4. **The withdrawn seed comparison ("arithmetic ≈2× random curl") is not used**;
   §3.2 gives the mechanism for the `0.09` vs `0.04` contrast in terms of
   signature *shape*.
5. **"Exchange metrics leave `HYP` at the first opportunity" (T1.1) is correct but
   must now be read with §4.2**: they leave it and cannot go far, because the
   whole non-line part is bounded by `2 log 2` uniformly in `n`.

---

## What this decides for the papers

**The framework is rigid, and the rigidity is Theorem A.** An exchange geometry is

> a **one-dimensional** geometry (the line metric `|σ_a − σ_b|` and the potential
> `ψ`), plus a correction universally bounded by `log 2` in the metric and
> `(log 2)/2` in the comparison, whose entire `β`-dependence is a **single
> unimodal bump of height ≤ `log 2` sitting at `β = τ_a = log r_a/log M_a`**.

That is a theorem about the framework with an explicit sharp constant, and the
observed `c₂ ≈ 1.1`, the mild pentagonal failure, the `4–9 %` curl, the rarity of
cycles and the `q → ∞` constant of brief F are all instances of it. It should be
the opening theorem of whatever paper carries the exchange geometry.

**But it does not answer brief A's objection in the negative, and
`A_paper_or_section.md`'s verdict stands.** The objection is about the map
`family ↦ signature`; the rigidity is about the framework's shape. Two things
change in how that section is written:

1. **The escape route stays open, but now has a cap.** Cycles are realisable in
   abundance (§2.2), so `Question~\ref{q:arith-cycle}` cannot be closed
   negatively; but §3.4 caps what any cycle can ever say — geometric-mean rate
   asymmetry at most 2, and `ψ`-spread under `log 2`. State the cap next to the
   question.
2. **`A` is not merely "unexamined".** The Szegedy discriminant cannot see `A` at
   all; flat families have `A` exactly a potential (§1.5); and the symmetric part
   *is* the line metric to within `2 log 2` (§4.2). The honest framing: *`A`
   differs from a potential by a quantity with an explicit universal bound, and
   that difference is where every non-scalar phenomenon in the programme lives.*

**Recommended one-sentence claim.** *The exchange comparison is an exact 1-form
plus a defect bounded by `(log 2)/2`, sharply; equivalently, the exchange rates
around any preference cycle differ by a geometric-mean factor of at most two.*

---

## Open

* **Is `ε ≤ log(1+e^{−Δ})` rather than `2log(1+e^{−Δ})`?** Every search returns
  `0.9889·log 2` at `Δ = 0`, always with one of `P, Q` exactly zero. A proof would
  also give `d ≤ ℓ + log 2`.
* **Is the maximum triangle curl `(log 2)/2` rather than `3(log 2)/2`?** Three
  independent parametrisations reach `0.9889·(log2)/2` and no further.
* **Prove `C_4` is unrealisable.** The single most valuable open item: it would
  convert §4.3 from a stalled search into the obstruction theorem. The peak
  heights `h_a = w_a(σ_a)` satisfy `|h_a − h_b| ≤ d(a,b) ≤ h_a + h_b` (proved,
  when the `σ` agree), which is the natural starting point.
* **Convexity of `U` in `s`.** Verified to `3.6e−9`; a proof would turn §4.1 into
  "at most `(#distinct atom-scales of a) + (of b) − 2` interior extrema".
* **Every tournament on every `n`?** Settled for `n ≤ 6` (74 classes). Theorem D
  reduces `n = 7` to exhibiting a defect matrix with the right cycle-sum signs.
* **A matching lower bound for the `2 log 2` metric budget.**

---

## Files

All under `research/realizability/`.

| file | what |
|---|---|
| `common.py` | `Sig` (multiplicity-compressed, exact big integers), `U`, `w`, `d`, `A`, `P`, `Q`, `D`; certified extrema by the Lipschitz bracket; 40-digit mpmath re-evaluation |
| `validate_against_package.py` | agreement with `fn_complexity.exchange_rate` on 600 integer pairs: `1.8e−13` |
| `validate_certification.py` | mpmath certification against the package on the known cycle |
| `gpools.py`, `optimizers.py` | pools; dependency-free DE and compass search |
| `structure.py`, `.csv` | Theorem A verification and the sharpness ladder |
| `extremes.py`, `.csv` | hill-climbs of the defect, of `ε`, and of triangle curl |
| `g1_mine.py`, `g1_classes.csv` | focused `(r,M)` pools, cycle counts, class census |
| `g1_target.py`, `g1_target_n5.*` | annealed subset search per isomorphism class |
| `g1_potential.py`, `g1_potential_n{5,6}.csv` | Theorem D: min-mean-cycle criterion + Cartesian-power realisation |
| `g1_atoms.py`, `.csv` | dense sweep of the two-fiber family |
| `g1_r2_witness.py`, `.json`, `.csv` | the certified two-fiber 3-cycle |
| `g2_metrics.py`, `.csv` | oscillation count; realisation of `K_n`, `C_4`, `C_5`, `K_{2,3}`, `K_{3,3}`, Petersen, Wagner `V_8` |
| `g3_curl.py`, `g3_hodge.py`, `g3_bound.py` and their CSVs | seed reproduction; hill-climbed curl fraction; the rigidity inequality |
| `*_output.txt` | captured stdout of each run |
| `tournament_seed.py` | the brief's seed measurement, unmodified |

Reproduce in order: `structure.py`, `extremes.py`, `g1_mine.py 4`, `g1_atoms.py`,
`g1_r2_witness.py`, `g1_potential.py 5`, `g1_potential.py 6`, `g2_metrics.py`,
`g3_bound.py`, `g3_hodge.py`.
