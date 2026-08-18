# Exchange rates as quantitative structure on a hierarchy of maps

*Working draft. Every numbered claim is either proved or computed; the
distinction is marked. Sources are the `research/` findings files.*

---

## 1. The framework is general, and it adds structure rather than objects

Fix any class of maps `P` and any notion of **implementation** — a relation
`f ⪯ g` meaning "`g` can be wired up to compute `f`". Classification of `P` under
`⪯` is the classical activity: it produces a preorder, and its Hasse diagram is
the hierarchy. For polynomial maps over a field the standard choice is *affine
implementation*, `f = b ∘ g ∘ a` with `a, b` affine, and the resulting posets of
quadratic and cubic maps over `F_q` are the objects this project started from.

A preorder answers one question — *can `g` do `f` at all?* — and is silent on
the quantitative one. The exchange rate supplies the missing structure. Writing
`f^{×k}` for the `k`-fold Cartesian power and

```
k_{g→f}(r) = max{ k : f^{×k} ⪯_aff g^{×r} },
```

block-diagonal composition gives `k(r+s) ≥ k(r) + k(s)`, so Fekete's lemma
yields

```
C_aff(g→f) = lim_{r→∞} k_{g→f}(r) / r .                              (1)
```

**This is the object.** It refines the hierarchy without changing it: `f ⪯ g`
already implies `C_aff(g→f) ≥ 1`, and `C_aff` grades every comparable pair by
*how many* copies, and every incomparable pair by an exchange ratio the preorder
cannot see. The construction is generic — it needs only a monoidal product and
an implementation relation — so it applies verbatim to any hierarchy of maps a
reader already cares about.

The computable proxy is the **signature rate**. The signature `σ(f)` is the
multiset of fiber cardinalities; for Cartesian powers fiber sizes multiply, and
the entropy method gives the closed form

```
C(g→f) = inf_{0 ≤ β ≤ ∞}  log Z_g(β) / log Z_f(β),   Z_a(β) = Σ_i a_i^β .   (2)
```

Here `β` is an inverse temperature and `Z` a partition function with energy
levels `−log n_i`; (2) is a theorem about the operational limit, not a second
definition. **Throughout, `C` means (2) and `C_aff` means (1).**

**The two are not the same, and the gap is the framework's main unaudited
liability.** Affine implementation forces signature implementation, so

```
C_aff ≤ C          always,                                            (3)
```

and no computation in this programme has measured the slack. It is not a
formality: over `F_11`, 1744 genus-two pencils `y² = P(x)+c` fall onto only 296
signatures, and **every one of the 209 signatures realised by more than one
pencil is realised by pencils with different fiberwise isogeny data**
(419 of 420 at `F_13`) [flux_arithmetic]. The signature is a genuinely lossy
shadow. Every arithmetic statement below is a statement about `C`, and the word
"signature" belongs in each of them.

---

## 2. What the rate reads off an arithmetic family

Let `f : A² → A¹` over `F_q` have geometrically irreducible fibers,
`N_c = #f^{-1}(c)`, `a_c = q − N_c` the Frobenius traces, `m₂ = q^{-2} Σ_c a_c²`,
and let `L = (q,…,q)` be the flat signature of a linear map. Two propositions
[m_and_e_and_a_c]:

**Proposition 1** *(proved).* `C(L→f) = log q / log max_c N_c`, attained at
`β = ∞`.

**Proposition 2** *(proved).* `1 − C(f→L) = (3−2√2)·m₂/(2q log q) + O(q^{−3/2})`,
with the infimum at

```
β* = √2 − 1 = 0.414213562…
```

independent of family, genus and `q`.

So the two ends of the temperature axis read two classical statistics: the
**extreme Frobenius trace** at `β = ∞` and the **second moment** at `β* = √2−1`.
The interior carries the whole moment ladder, each order damped by `≈ 0.6/√q`.
Flatness is equivalent to a permutation-polynomial condition, which subsumes the
known `q ≡ 2 mod 3` result and reaches non-congruence conditions such as
supersingular primes.

**Two things are provably invisible**, and they are invisible for structural
reasons rather than by accident: the *smallest* fiber, because isolating it
needs `β < 0` where `Z` is no longer order-preserving; and the Katz–Sarnak
symmetry type, because `Σ_c a_c = 0` identically and that first moment is
exactly the separating statistic.

**The objection, stated plainly.** `C(L→f)` is a strictly monotone function of
`max_c N_c`, and `1 − C(f→L)` is to leading order affine in `m₂`. The rate is a
function of the signature and the signature is the list `{N_c}`, so nothing the
rate reports about arithmetic could fail to be computable from `{N_c}` directly.
Information flows into the matrix, not out of it. What survives is not a new
theorem about curves but a statement about the *comparison*: of all ways to
compress `{N_c}`, resource-theoretic comparison against a linear map returns
exactly the extreme value and the second moment, and it returns them at two
specific temperatures. That is a consistency check on the framework, and it is
how it should be written.

---

## 3. The main result: the geometry is one-dimensional plus a bounded defect

This is, in my view, the most important thing the programme has produced,
because it is the only statement about the *framework* rather than about
particular resources.

Write `F_a = log Z_a`, `u_a = log F_a`, `R = log r` (`r` = number of fibers),
`Λ = log(max fiber)`, `σ = log(R/Λ)`, and `s = log β`.

**Theorem A** *(proved; verified to `3.6·10⁻¹⁵`)* [realizability].

```
u_a(e^s) = log Λ_a + max(σ_a, s) + w_a(s),      0 ≤ w_a ≤ log(1 + e^{−|s−σ_a|}),
```

with `w_a` **unimodal, peaked exactly at `s = σ_a`, and 1-Lipschitz**, of height
`≤ log 2` with equality iff `a` is flat. Every profile is a kink plus one bump.

Splitting `−log C` into symmetric and antisymmetric parts — `S = d/2` the
exchange metric, `A` the comparison, with `a ≺ b ⟺ A(a,b) > 0` — Theorem A gives

```
d(a,b) = |σ_a − σ_b| + ε,     0 ≤ ε ≤ log(1 + e^{−|Δσ|}) ≤ log 2,       (4)
A(a,b) = ψ(b) − ψ(a) + D,     |D| ≤ ½ log(1+e^{−|Δσ|}) ≤ (log 2)/2,     (5)
```

with `ψ = ½ log(log r · log M)`. Both constants are **sharp** and both bounds are
**independent of the number of points**. [The constant in (4) is half the first
published bound; the sharpening is due to the obstruction session.]

**Theorem B** *(proved)* [realizability/OBSTRUCTION]. `d(a,b)` is exactly the
**Hilbert projective metric** between `F_a` and `F_b` on the cone

```
C = { Φ : Φ convex, nondecreasing, Φ ≥ Λ_Φ·β },
```

and `C` is exactly the projective closure of the achievable set. Cartesian powers
are the projective rescaling the Hilbert metric quotients out.

**Everything previously observed is a corollary.** The `ℓ₂`-distortion `≈ 1.1`
with `c₂/log n` *falling* — Bourgain's `O(log n)` is unreachable because the
non-line budget does not grow with `n`. The mild pentagonal failure, and every
negative-type violation, confined to a window of diameter `2 log 2`. The measured
`4–9 %` curl of the comparison and its flatness in `n` from 8 to 698 points — the
bounds in (5) are per-edge. The same constant reappearing in the `q → ∞` limit.

One consequence deserves its own line, because it is what §4 turns on. Around a
directed preference cycle every `A` shares a sign, so `|curl A| = Σ|A|`, and (5)
gives:

**Corollary C** *(proved).* The geometric mean of the rate asymmetry
`C(b→a)/C(a→b)` around **any** preference cycle is at most **2**.

---

## 4. Are the cycles important?

The exchange comparison is not transitive. Over `F_11` there are 132 strict
3-cycles among genus-two pencils `y² = P(x)+c` in a complete enumeration, 1475
over `F_13`, and certified witnesses exist over `F_11` and `F_101` with margins
`10⁻³`–`10⁻⁴` verified three ways [curve_family_cycles, session D]. Substantial
effort in this programme went into finding and certifying them. **Are they
worth it?**

The case against has to be stated first, because three separate results weaken
the obvious reading.

* **Cycles are not arithmetic.** A random pool of signatures matched to a curve
  pool in the only things that matter (`q` fibers, `Σ N_c = q²`, matched trace
  distribution) has the *same* curl, and the tightest control has slightly more.
  Over 108 pools, `‖curl‖/‖A‖` is a function of the trace spread to `R² = 0.93`
  [flux_arithmetic]. Cycles are a property of near-flat signature *shape*, not of
  the arithmetic that produced the signature.
* **Cycles are weak.** Corollary C caps the asymmetry around any cycle at a
  geometric-mean factor of 2, and every 3-cycle has `ψ`-spread below `log 2`.
* **Cycles are a zero-temperature phenomenon.** Soften both extrema over a
  density `ρ` on `s = log β`:
  `A_λ = ½(softmax_λ + softmin_λ)` of `u_a − u_b`. Then `A_λ → A` as `λ → ∞`,
  and as `λ → 0` **both** soft-extrema tend to `∫(u_a−u_b)ρ`, so

  ```
  A_0(a,b) = Ψ(a) − Ψ(b),      Ψ(a) = ∫ u_a ρ,
  ```

  an **exact potential difference for every `ρ`** — hence a total order with no
  cycles. Measured on the standard 3-cycle: no cycle below `λ ≈ 200`, cycle above
  `λ ≈ 300` [maslov]. **The framework supplies a continuum of perfectly good
  scalar complexities, one per temperature prior, and the tropical limit is the
  unique member that fails to be one.**

That third point looks fatal and is not, for one reason: **the infimum in (2) is
not a modelling choice.** It is what the asymptotic conversion rate *is* — the
constraint is that *every* temperature's monotone be satisfied simultaneously, and
the entropy method's theorem is that the operational limit equals the infimum.
The soft family `A_λ` for `λ < ∞` currently has no operational reading at all.

So the defensible claim, and the only one I would put in a paper, is:

> **Asymptotic convertibility is not governed by any scalar complexity — and the
> failure is bounded: rates around any preference cycle differ by a geometric-mean
> factor of at most two.** The non-scalarity is a property of *worst-case*
> conversion, which is what conversion means; averaged readings of the same data
> are scalar and cycle-free, but do not describe conversion.

That is a complete narrative — promise, obstruction, sharp bound on the
obstruction — and it is the narrative I would continue. What it does **not**
support is the stronger story the cycle hunt was implicitly aimed at: that cycles
reveal something arithmetic. They do not; that reading is refuted.

**What would raise their status.** Two concrete things, both open. (i) A cycle in
`C_aff` rather than `C`: every certified cycle so far is a cycle of the signature
shadow, and by (3) it need not survive. (ii) An operational meaning for some
`λ < ∞` — average-case conversion, or conversion with error — which would turn
the `λ`-family from a deformation into a physical phase diagram, and `λ_c` into
a real transition rather than an artefact of the prior.

---

## 5. Status

Proved and independently re-verified: Propositions 1–2, Theorems A–B,
Corollary C, the flatness/permutation-polynomial theorem, the invisibility of the
smallest fiber and of symmetry type, superadditivity of `C` and `C_aff` under
Cartesian products. Computed and certified: the cycle censuses and their
witnesses. **Open and load-bearing:** the size of the gap in (3); whether
`C_aff` cycles; whether any `λ < ∞` is operational.

Earlier drafts of this material claimed an obstruction — that the exchange
metrics do not fill `MET` — on the strength of a stalled search. `C_4` is
realisable exactly, and that claim is withdrawn.
