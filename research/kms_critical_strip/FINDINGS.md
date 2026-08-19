# Findings — does the KMS condition give a comparison of resources on the critical strip?

Answer to session brief H, the first (and, on the evidence below, the last)
session of brief D's stage 5.

**No. Stop the programme.** The reason is a theorem, not a numerical failure:
for `0 < β ≤ 1` the Bost–Connes system has **exactly one** KMS`_β` state, so on
the critical strip there is no pair of states to compare and the `Ẑ*` symmetry
that was supposed to do the work acts trivially. Above `β = 1`, where a
`Ẑ*`-family of extremal states does exist, the Galois action is a **relabelling**
of the energy levels, and every monotone of the exchange framework is a
symmetric function of the level multiset — so the family is invisible to the
framework at every temperature where the framework is defined.

Everything below is either proved, or computed and independently re-verified;
the distinction is marked in each case. Two of the brief's own statements are
corrected: the direction of the Bost–Connes phase transition, and the claim that
the KMS route is "the only route in the whole programme that reaches `Re β = ½`".

**Placement, in the sentence the brief asks for.** Sections 1, 2 and 4 are
Julia's Riemann gas and Bost–Connes restated in exchange language; nothing in
them is new to a referee who knows that literature, and the exchange framework
neither shortens nor is shortened by it. The one item here that is not a
restatement is §3's decomposition `ξ = A − B` of the completed resource, which
is Riemann's own 1859 computation read as a statement about the semiring, and it
is a negative: `ξ` is not reachable from `⨂_p P_{p,∞}` by any operation the
framework has.

**Standing obstructions, unchanged by anything in this note.** Both are recorded
in `research/m_and_e_and_a_c/FINDINGS.md` and are repeated here because the KMS
route removes neither:

* **atomic measures `Σ_i δ_{a_i}` are not admissible Weil test functions**, so
  the matrix `E` is a finite-rank truncation and `|Z_a(½+iγ)|` does not decay;
* **the exchange monotone diverges in the critical strip.**

The KMS route addresses the second and not the first. What §§3–5 below show is
that it does not really address the second either: it *names* the divergence
(Hagedorn), *classifies* what replaces the Gibbs state (Bost–Connes), and then
supplies nothing that a resource comparison can be built from.

---

## Correction to the brief

Session brief H states the Bost–Connes phase transition backwards. It says

> for `β > 1` there is a unique KMS`_β` state, and it is the Gibbs state […] at
> `β ≤ 1` the Gibbs state does not exist; KMS`_β` states still do, they form a
> simplex, and for `β ≤ 1` the symmetry group `Ẑ* ≅ Gal(Q^{ab}/Q)` acts on them —
> spontaneous symmetry breaking

The theorem (Bost–Connes 1995, Thm 5; Connes–Marcolli, *Noncommutative Geometry,
Quantum Fields and Motives*, Thm 3.1) is the other way round, in the direction
spontaneous symmetry breaking always goes — symmetry breaks at **low**
temperature:

| range | KMS`_β` states | type | symmetry |
|---|---|---|---|
| `0 < β ≤ 1` | **unique** | III`_1` factor state | `Ẑ*` fixes it — unbroken |
| `β > 1` | simplex whose extreme points are a free transitive `Ẑ*`-torsor | I`_∞`, partition function `ζ(β)` | broken |

So the Galois-permuted *family* lives exactly where the exchange framework
already lives, and at `β = ½` there is one state and nothing to compare it with.
The brief's success criterion — "if the extremal KMS states below `β = 1` give a
*family* of comparisons permuted by `Gal(Q^{ab}/Q)`" — is excluded by the
uniqueness half of the theorem before any computation starts.

The brief's second claim, that this is "the only route in the whole programme
that reaches `Re β = ½` at all", is also wrong, and in the project's favour. §2
below locates the critical line at `β = 0` in the completed-zeta normalisation
the companion paper already uses. The framework reaches the critical *line*
trivially; what it does not reach is a *nondegenerate comparison* there, and the
reason is normalisation, not divergence.

---

## Notation — three inverse temperatures, and they are not the same

The brief warns about two unrelated `t`'s. There are three unrelated `β`'s, and
conflating them is what makes "the framework lives at `β > 1`, RH lives at
`β = ½`" sound like a near miss when it is not one.

```
s        the Dirichlet exponent / primon-gas inverse temperature.
         Z(s) = Σ_n n^{-s}.  Abscissa s = 1, strip 0 < s < 1, line s = ½.

β_x      the exchange framework's temperature: Z_a(β) = Σ_i a_i^β with entries
         a_i ≥ 1 and β ∈ [0,∞].  For a resource in the cost convention β_x = −s.

β_ξ      the completed-zeta temperature of the two-positivities note:
         Z_ξ(β) = ξ(½+β)/ξ(½).  Here β_ξ = s − ½, so the critical LINE is β_ξ = 0.
```

Under `β_x = −s` the primon gas `{1,2,3,…}` converges only for `β_x < −1`. The
framework's admissible range `β_x ∈ [0,∞]` therefore lies **entirely inside the
divergence half-line**, and the critical strip sits at `β_x ∈ (−1,0)` — at
*negative* exchange temperature. The gap is not "the framework stops at 1 and RH
starts at ½"; in this normalisation the framework and the zeta resource never
overlap at all.

Under `β_ξ = s − ½` the picture is the opposite: `log Z_ξ` is finite for every
real `β_ξ`, the critical line is the framework's own `β = 0` endpoint, the
functional equation is the reflection `β_ξ ↦ −β_ξ` about it, and the zeros sit at
`β_ξ = iγ` — purely imaginary, i.e. real time. That last row is brief D's Wick
rotation, now with the axis fixed.

Both conventions appear in §3 of `exchange_positivity_and_weil.md`, in adjacent
paragraphs, without being distinguished.

---

## 1. The primon gas, and what breaks *(computed; independent verification)*

`H|n⟩ = log n |n⟩` on `ℓ²(ℕ)`, `Tr e^{−sH} = Σ_n n^{−s}`.

* `Tr e^{−sH} = ζ(s)`: partial trace to `n = 10⁵` plus an Euler–Maclaurin tail,
  computed without calling `zeta`, agrees with `ζ(s)` to `1.6·10⁻³³` at
  `s = 1.1` and better beyond.
* **The Euler product is the tensor factorisation of the resource.**
  `Π_{p ≤ 5·10⁵}(1−p^{−s})^{−1}` against `ζ(s)`: relative gap `1.4·10⁻⁷` at
  `s = 2`, `3.0·10⁻²⁵` at `s = 5`, and `0.124` at `s = 1.1` — the convergence of
  the product degrades exactly as the transition is approached.
* `U·(s−1) → 1` and `S/U → 1` at 40 digits, reproducing
  `research/m_and_e_and_a_c/primon_gas_hagedorn.py` (commit `c6848fa`) from an
  independent code path. Added here: `log ζ(s)/(−log(s−1)) → 1`, so the **free
  energy diverges only logarithmically** while energy and entropy diverge like a
  simple pole. In the isometry coordinate `u = log log Z` of Theorem 1 the
  divergence is `u ~ log log 1/(s−1)`: doubly logarithmic. The resource does run
  away, but at the slowest rate the framework's coordinates can express.
* **The Hagedorn statement in counting-function form.** `N(E) = #{n : log n ≤ E}
  = ⌊e^E⌋`, so `log N(E)/E → 1` (measured `0.999999999958` at `E = 20`). The
  abscissa of convergence of an exchange monotone equals the exponential growth
  rate of the level counting function, which is `1` here and `0` for every finite
  signature. **No object of the exchange semiring has a phase transition at all**
  — `Z_a` is entire for every finite `a`. The transition belongs to a limit point
  that is not in the semiring.

Reproduce: `primon_gas.py`.

## 2. The framework *is* the `s > 1` Gibbs restriction *(computed)*

* The exchange monotone of a truncated Euler factor is its Gibbs free energy:
  the repo's `exchange_rate` on `(1,p,…,p^K)` against the closed form
  `min(inf_β log Z_P/log Z_Q, log(K+1)/log(L+1), K log p/(L log q))` agrees to
  **`4.4·10⁻¹⁶`** over 16 pairs — at the solver's stated `~10⁻¹³` accuracy, an
  exact match.
* `analysis/xi_versus_euler_factors.py` reproduced at `mp.dps = 40` with true
  extrema (golden-section refinement, not a 240-point grid). **The published
  table is good to about `10⁻⁵`, and here is the error budget it needs stated
  alongside it:**

  | window `β_ξ` | published `d_W` (`K→∞`) | 40-digit `d_W` | gap |
  |---|---:|---:|---:|
  | `[0.5, 5]` | `2.283688` | `2.283692005907073682` | `4.0·10⁻⁶` |
  | `[1, 10]` | `2.236350` | `2.236350227673908866` | `2.3·10⁻⁷` |
  | `[2, 20]` | `2.114300` | `2.114256576533592508` | `−4.3·10⁻⁵` |

  Per-`(p,K)` entries agree to `3.5–4.3·10⁻⁶` on `[0.5,5]` and to
  `7·10⁻⁸–4.3·10⁻⁷` on `[1,10]`. The larger error on `[0.5,5]` has an
  identifiable cause: that window's left endpoint is `s = 1`, the pole of `ζ`,
  which the published script sidesteps by shifting the grid point by `10⁻⁶`.
  Here `(s−1)ζ(s)` is evaluated as a regular function instead. **The agreement
  claim is therefore: the published numbers are correct to five decimals and no
  further, and the discrepancy is grid plus pole-shift, not method.**
* **Where the critical line actually is.** `log Z_ξ(0) = 0` exactly: the
  completed resource is a probability measure, one unit of resource, while any
  finite signature has `log Z_a(0) = log r > 0`. So `u_ξ = log log Z_ξ → −∞` at
  `β_ξ = 0` and both rates against any signature vanish there. **At the critical
  line the axiom that fails is the normalisation, not the finiteness of the
  monotone** — the opposite of what happens in the primon-gas normalisation,
  where finiteness fails and normalisation is fine. Which axiom breaks first is
  therefore a property of the chosen embedding, not of the resource.
* The only nondegenerate quantity surviving the `β_ξ → 0` limit is a ratio of
  *curvatures*: for centred resources
  `lim_{β→0} log Z_a(β)/log Z_b(β) = (log Z_a)''(0)/(log Z_b)''(0)`, and for `ξ`
  the numerator is `(log Z_ξ)''(0) = 0.04620998623083794158`, i.e.
  `2 Σ_{γ>0} γ^{−2}`. **This is an endpoint limit of `C`, defined where `C`
  already is, so by the brief's own criterion it is a renaming, not a new
  functional.**

Reproduce: `gibbs_restriction.py`.

## 3. What the completion does, and why it is not a semiring operation *(proved + computed)*

The two-positivities note says the completion "repairs `ζ` as a partition
function" and calls the archimedean factor "the real place's own processor". The
first half is true; the second is imprecise in a way that decides this brief.

**Proposition.** Riemann's positive measure `Φ`, with
`ξ(½+β) = ∫Φ(u)e^{βu}du`, is a *difference* of two positive measures,

```
Φ = A − B,   A(u) = 4π² e^{9u/2} Σ_n n⁴ e^{−πn²e^{2u}},
             B(u) = 6π  e^{5u/2} Σ_n n² e^{−πn²e^{2u}},
```

whose two-sided Laplace transforms are, with `s = β + ½`,

```
Â(β) = 2 π^{−¼−β/2} Γ(9/4 + β/2) ζ(s),
B̂(β) = 3 π^{−¼−β/2} Γ(5/4 + β/2) ζ(s),
Â/B̂  = (2/3)(5/4 + β/2) = (s+2)/3,
Φ̂(β) = B̂(β)·(s−1)/3 = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s) = ξ(s).
```

Both `A` and `B` are honest resources — positive measures — and both have
**abscissa of convergence exactly `s = 1`**: as `u → −∞` each behaves like
`(3/2)e^{−u/2}`, so `∫A e^{βu}du` converges iff `β > ½`. Their price ratio
`Â/B̂ = (s+2)/3` equals `1` **exactly at `s = 1`**. The completion is the
cancellation of two divergent resources that become equal in price precisely at
the Hagedorn temperature; the surviving factor `(s−1)` is, on the measure side,
the first-order differential operator `f ↦ −f′ − f/2`.

**Consequence.** `ξ` is a legitimate resource, but multiplication by `s(s−1)` is
neither `⊗` nor `⊕`: it needs subtraction and differentiation, and the signature
semiring has neither. **There is no chain of framework operations carrying
`⨂_p P_{p,∞}` to the `ξ`-resource.** Tensoring with the archimedean measure
`π^{−s/2}Γ(s/2)` *is* a `⊗`, and it leaves the abscissa at `s = 1`; the step
that moves it to `−∞` is the one that is not a resource operation. The
completion does not move the phase transition — it cancels it between two copies
of itself.

Verified at 40 digits: `Φ > 0` and even (to `3·10⁻⁴⁰`); `∫Φ e^{βu}du = ξ(½+β)`
to `2.6·10⁻⁴¹` relative at `β = 0,½,1,2,5`; the closed form for `Â` against
quadrature-plus-analytic-tail to `1.4·10⁻⁴¹`; `A(u)e^{u/2}` and `B(u)e^{u/2}`
both `→ 3/2` with difference `6.1·10⁻²⁴` already at `u = −1.5`;
`Â/B̂ − (s+2)/3 = 0` exactly.

Also verified, because the note quotes it without a source:
`½(log Z_ξ)''(0) = Σ_{γ>0} γ^{−2}`. From 1200 zeros plus a Riemann–von Mangoldt
tail, `0.023105132519`, against `0.023104993115` from direct differentiation —
agreement `1.4·10⁻⁷`, the size of the tail estimate. Note that this is **not**
the classical `Σ_ρ 1/ρ = 1 + γ_E/2 − ½ log 4π = 0.023095708966`; the two differ
in the fifth digit and the note's `0.0231` is ambiguous between them.

Reproduce: `xi_as_positive_measure.py`.

## 4. The KMS replacement, and the Galois action *(theorem cited + computed)*

Restricted to the group algebra `C*(Z/q) ⊂ C*(Q/Z)`, the extremal KMS`_s` state
labelled by `g ∈ (Z/q)*` is the residue distribution of the zeta measure:

```
P_j^{(g)}(s) = ζ(s)^{−1} Σ_{n : gn ≡ j (q)} n^{−s} = q^{−s} ζ(s, {g^{−1}j}/q)/ζ(s).
```

Hurwitz form against direct summation with an Euler–Maclaurin tail: `5·10⁻¹⁰` at
`s = 1.5`, `4.3·10⁻¹⁸` at `s = 3`; total mass `1` to `2·10⁻⁴¹`.

### (a) The Galois action is a relabelling, hence invisible *(proved)*

**Proposition.** No exchange monotone can detect the `Ẑ*` symmetry breaking, at
any temperature.

*Proof.* The `Ẑ*` action commutes with `σ_t` and fixes every isometry `μ_n`,
hence fixes `H`, hence fixes `Tr e^{−sH} = ζ(s)` and every functional of the
level multiset. Equivalently `σ_t(e(a)) = e(a)`, so the whole Galois orbit lies
in the **zero-energy sector**, the fixed-point algebra of the dynamics. By the
first paper's representation theorem the `⊗`-multiplicative, `⊕`-additive,
`≼`-monotone functionals are exactly `a ↦ Z_a(β) = Σ_i a_i^β`, symmetric
functions of the level multiset. A symmetric function cannot see a permutation
of labels. ∎

Measured at `q = 12`, `s = 2`: `P^{(g)}` sorted agrees with `P^{(1)}` sorted to
`0`, while unsorted they differ by `0.586`, `0.598`, `0.606` for `g = 5,7,11`.
The order parameter of Bost–Connes is `φ_β(e(a))` for `a ∈ Q/Z` — a zero-energy
observable, outside the framework's entire vocabulary, which consists of
partition functions.

This is the brief's own contemplated honest negative: *"If the Galois action
turns out to be invisible to any resource-theoretic quantity, that is the honest
negative and is also worth writing."* It is invisible, and the proof is one line.

### (b) The `Ẑ*`-average, and the state that actually exists at `β = ½` *(computed)*

Averaging the extremal states over `Ẑ*` turns the polylogarithm into a Ramanujan
sum and `ζ` cancels:

```
F_s(q) = (1/φ(q)) Σ_{d|q} d^{1−s} μ(q/d) = Π_{p|q} p^{−(k_p−1)s} (p^{1−s} − 1)/(p−1).
```

Divisor sum against the average of extremal polylogarithms: `7·10⁻⁴²`;
Euler-product form against divisor sum: `< 10⁻³⁰` (asserted in code). Limits:
`F_∞(q) = μ(q)/φ(q)` (the Möbius function), `F_1(q) = 0` for `q > 1` (Haar
measure on `Ẑ`), `F_0(q) = 1`. This function, continued below `s = 1`, **is**
Bost–Connes' unique KMS`_s` state — and it is a genuine state: its Fourier
coefficients on `Z/q` at `s = ½` are nonnegative for every modulus tested, from
`0.2929` at `q = 2` down to `0.004277` at `q = 60`.

**So there is an object at `β = ½`.** There is exactly one, it is `Ẑ*`-invariant,
and it is the *same* object whatever resource one started from, because there is
only one Bost–Connes system. It is a state of one algebra, not a comparison of
two resources.

### (c) The naive continuation of the Gibbs weights is not a state *(computed)*

`q^{−s}ζ(s,j/q)/ζ(s)` continues below `s = 1`, keeps total mass `1` exactly, and
acquires negative entries:

| `q` | `min_j P_j` at `s = ½` | `s_c` (first sign change) | `q(1−s_c)` |
|---:|---:|---:|---:|
| 2 | `+0.2928932188` | none in `(0,1)` | — |
| 3 | `+0.0466841721` | `0.4312925683` | `1.706` |
| 4 | `−0.0821593398` | `0.6108185579` | `1.557` |
| 8 | `−0.2934350250` | `0.8351269682` | `1.319` |
| 20 | `−0.4513444790` | `0.9422584857` | `1.155` |
| 100 | `−0.5838742250` | `0.9895692295` | `1.043` |

`q(1−s_c) → 1`: the first negative weight appears at `1 − s_c ≍ 1/q`, so every
modulus `q ≥ 4` already fails at `s = ½`. The "analytically continued Gibbs
state" on the critical strip is a signed measure of total mass one, and the
resource-theoretic reading of a negative weight is a negative multiplicity. This
is the quantitative form of "the Gibbs state ceases to exist": not that the sum
diverges — the continuation is perfectly finite — but that it stops being
positive, at rate `1/q` in the modulus.

Reproduce: `kms_states.py`.

## 5. The Connes cocycle: the reduction that holds, and the one that fails *(computed + proved)*

The brief nominates `[Dφ : Dψ]_t` and says the reduction check "is the entire
scientific content of the checkpoint". Here it is.

### (a) The type-I reduction holds *(computed)*

For the truncated gas at `s₁, s₂ > 1`,
`[Dω₁:Dω₂]_t = ρ₁^{it}ρ₂^{−it} = diag(n^{−i(s₁−s₂)t})·(ζ(s₂)/ζ(s₁))^{it}`,
matched to `7.2·10⁻⁴¹`, with the cocycle identity
`u_{t+t'} = u_t σ^{ω₂}_t(u_{t'})` to `7.3·10⁻⁴¹`. Relative entropy is the
Bregman divergence of the log partition function,
`S(ω_{s₁}‖ω_{s₂}) = log Z(s₂) − log Z(s₁) − (s₂−s₁)(log Z)'(s₁)`, exact to
`5.7·10⁻⁴²` on a common finite truncation. **So the cocycle does reduce to a
ratio of Gibbs weights, as the brief expected.**

### (b) The reduction to `C(g→f)` fails, structurally *(proved, with numbers)*

A ratio of Gibbs weights is not an exchange rate. Four independent obstructions,
each measured, at `β = 2` on `a = (12,10,8,8,2,1)`, `b = (11,9,7,7,4,1)`.

1. **Difference, not ratio.** On a common algebra with two Hamiltonians,
   `S(ω_a^β‖ω_b^β) = β⟨H_b − H_a⟩_a + log Z_b(β) − log Z_a(β)` — verified exact
   to `1.6·10⁻⁴⁰`. The log-partition functions enter through their *difference*.
   `C(a→b) = inf_β log Z_a/log Z_b` needs their *ratio*, which is what makes `C`
   dimensionless, i.e. a rate.
2. **Three different homogeneity degrees under one operation.** With `a^{⊗k}` the
   Cartesian power:

   | quantity | behaviour | measured |
   |---|---|---|
   | `S(ω_{a^{⊗k}}‖ω_{b^{⊗k}})` | extensive, `= k S₁` | ratio `1.0` for `k = 1…5` |
   | `−log C(a^{⊗k}→b)` | log-extensive, `= −log k − log C` | residual `5·10⁻⁴²` |
   | `d(a^{⊗k}, b^{⊗k})` | intensive, invariant | drift `1.1·10⁻¹⁶` |

   No function of one can be the other.
3. **Relative entropy is not defined on the objects the framework compares.** It
   is an invariant of the *algebra*; the exchange framework retains only the
   *spectrum*, the multiset of entries. Merging equal products is a
   coarse-graining and by data processing strictly decreases relative entropy:
   on the spectral algebra `S_k/(kS₁)` falls `1, 0.9784, 0.7652, 0.6691, 0.5930`
   for `k = 1…5`, while `C` is unchanged by construction.
4. **The zero sets differ.** `d(a, a^{⊗2}) = 0` (measured `6.7·10⁻¹⁶`, with
   `C(a→a^{⊗2}) = 0.5` and `C(a^{⊗2}→a) = 2` exactly), while the two Gibbs
   states are different measures and every relative entropy between them is
   positive.

**What `C(g→f)` would have to mean in KMS language, stated as the brief asks.**
It cannot be a functional of a *pair of states*. `C` is a functional of a *pair
of Hamiltonians* — of two C\*-dynamical systems, not of two states of one. In
the quantum-statistical-mechanical register the exchange framework is therefore
comparing *QSM systems*, and the developed theory of that comparison is
Cornelissen–Marcolli: morphisms and isomorphisms of Bost–Connes systems, where
the answer is a **rigidity theorem** — two number fields are isomorphic iff their
BC systems are, and an isomorphism of abelianised-Galois character groups
matching all `L`-series forces isomorphism of the fields — not a rate. There is
no quantitative "exchange rate between QSM systems" in that literature, and this
session found no reason to expect one.

### (c) The type-III obstruction, named *(cited)*

Not computable: these are theorems about a type III`_1` factor with nothing
finite-dimensional to evaluate. They are the decisive part.

* **Uniqueness.** For `0 < β ≤ 1` the KMS`_β` state is unique, so for any two
  KMS`_β` states `φ = ψ`, hence `[Dφ:Dψ]_t = 1` and `S(φ‖ψ) = 0` identically.
  *There is no pair to compare.*
* **Across temperatures.** For `β₁ ≠ β₂` the two KMS states have modular groups
  `σ_{−β₁t}` and `σ_{−β₂t}`, differing by `σ_{(β₂−β₁)t}`, which is not inner for
  the BC dynamics; the states are not quasi-equivalent and Araki's relative
  entropy is `+∞`. The comparison takes only the values `{0, +∞}`.
* **The invariant vanishes for III`_1`.** The flow of weights of a type III`_1`
  factor is the trivial flow on a point; Connes' invariants are `S(M) = [0,∞)`
  and `T(M) = {0}`, both content-free. And by the Connes–Størmer transitivity
  theorem the unitary orbit of any faithful normal state is norm-dense in the
  normal state space of the hyperfinite III`_1` factor — so *any unitarily
  invariant comparison of states there is constant*.

That is the obstruction `X` the brief asked to have named, and it is sharper than
"there is no trace, so the arithmetic is hard". It is: **there is exactly one
state; even if there were two, the Radon–Nikodym invariant that would compare
them is precisely the one type III`_1` annihilates; and the comparison the
framework actually needs was never a comparison of states in the first place.**

Reproduce: `connes_cocycle.py`.

---

## Go / no-go, with its reason

**NO-GO. Do not continue the KMS programme (brief D stage 5).**

The brief set two criteria and neither is met, and both fail for reasons that are
theorems rather than difficulties:

* *"A comparison that is defined at `β = ½` and reduces to `C` at `β > 1`."*
  There is none. Below `β = 1` there is one state and one algebra: no pair of
  resources, no partition function, and the only Radon–Nikodym invariant is
  killed by type III`_1`. The one thing genuinely defined at the critical line
  (§2, the curvature ratio `(log Z_a)''(0)/(log Z_b)''(0)`) is an endpoint limit
  of `C` in the `ξ`-normalisation — defined where `C` already is, hence by the
  brief's own rule a renaming.
* *"The `Ẑ*` symmetry doing work."* It cannot. The symmetry breaks at `β > 1`,
  not `β ≤ 1`; where it does break it acts by relabelling energy levels; and
  every monotone the framework possesses is a symmetric function of that
  multiset. Proved in §4(a) in one line.

Add that no framework operation reaches the completed resource at all (§3), and
the honest summary is that the KMS route is a **precise renaming of the
obstruction, not a way past it** — the same verdict brief D reached for the
quantum language generally: *"The quantum language renames both obstacles
precisely and removes neither."*

**What is worth keeping**, and it is one paragraph, not a programme:

1. The `ξ = A − B` proposition of §3, as the exact statement of why the
   completion is outside the semiring. This belongs in §3 of
   `exchange_positivity_and_weil.md`, replacing the looser "the completed `ξ`
   adds the archimedean factor, which is the real place's own processor".
2. The three-`β` notation table, because that section currently uses two of the
   three conventions in adjacent paragraphs without saying so.
3. The `Ẑ*`-invisibility proposition of §4(a), as a one-line closure of the
   Bost–Connes item in brief D Part 3 — which called it "the most promising
   single item in the whole brief".
4. The error budget of §2 for `analysis/xi_versus_euler_factors.py`: the
   published table is good to `10⁻⁵`, and the larger error on `[0.5,5]` is the
   `ζ`-pole shift at that window's left endpoint.

Brief D commits the project to claiming only two new things on the quantum side —
the gauge decomposition and the `t = ½` Szegedy point. **Nothing here should be
added as a third.**

## Open

* Whether the `ξ`-resource `Φ` admits a *positive* factorisation into local
  pieces at all. §3 shows Riemann's is a difference; it does not prove that no
  positive one exists. A negative there would be the sharp statement, and it was
  not attempted.
* Whether `Ẑ*` acts on anything the framework could be *enlarged* to see. It
  would have to be a monotone that is not a symmetric function of the level
  multiset, i.e. outside the first paper's representation theorem — so this is a
  question about changing the framework, not about the KMS route.
* Cornelissen–Marcolli rigidity as a source of *qualitative* comparison of
  arithmetic resources (isomorphism, not rate). Not pursued; it is a different
  question from the exchange rate and would need a different brief.

## Files

| file | what |
|---|---|
| `primon_gas.py` | step 1: `Tr e^{−sH} = ζ`, Euler product as `⊗`, Hagedorn rates, counting-function form; independent check of `U(s−1)→1` and `S/U→1`. Writes `primon_gas.csv` |
| `gibbs_restriction.py` | step 2: `exchange_rate` against the closed form on truncated Euler factors; `xi_versus_euler_factors.py` at 40 digits with the error budget; the three-`β` table. Writes `gibbs_restriction.csv` |
| `xi_as_positive_measure.py` | §3: `Φ > 0`, `∫Φe^{βu}du = ξ(½+β)`, the `A − B` decomposition and its closed-form transforms, `½(log Z_ξ)''(0) = Σ γ^{−2}`. Writes `xi_as_positive_measure.csv` |
| `kms_states.py` | step 4: extremal KMS`_s` states as residue distributions, the Galois-relabelling measurement, the `Ẑ*`-average and its positivity at `s = ½`, the sign-change threshold of the continuation. Writes `kms_states.csv` |
| `connes_cocycle.py` | step 5: the type-I reduction of `[Dφ:Dψ]_t`, the Bregman identity, and the four obstructions to reading `C` as a relative entropy. Writes `connes_cocycle.csv` |

Prior art assumed and used: Bost–Connes (1995); Julia's Riemann gas (1990) and
its Hagedorn reading; Connes' trace formula and the adele class space;
Connes–Marcolli on quantum statistical mechanics and number theory;
Cornelissen–Marcolli on QSM, `L`-series and anabelian geometry; Connes' type
classification, the flow of weights, and the Connes–Størmer transitivity theorem;
Berry–Keating `xp` and Bender–Brody–Müller (2017) for the Hilbert–Pólya register
of brief D Part 1. None of it is shortened by anything in this repository.
