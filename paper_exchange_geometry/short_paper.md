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
over `F_13`, with certified witnesses over `F_11` and `F_101` verified three ways
[curve_family_cycles, session D]. Substantial effort went into finding and
certifying them. **Are they worth it?** Three results weaken the obvious reading,
one strengthens it, and one bounds how much can be claimed.

**Against: cycles are not arithmetic.** A random pool of signatures matched to a
curve pool in the only things that matter (`q` fibers, `Σ N_c = q²`, matched
trace distribution) has the *same* curl, and the tightest control has slightly
more. Over 108 pools, `‖curl‖/‖A‖` is a function of the trace spread to
`R² = 0.93` [flux_arithmetic]. Cycles are a property of near-flat signature
*shape*, not of the arithmetic that produced the signature.

**Against: cycles mostly do not descend to the maps.** The `F_3` tensor
three-cycles are cycles of *signatures*, and since `C_aff` is not a signature
invariant (§1) the question is only well posed after choosing affine orbits. The
seven signature cycles lift to twelve orbit triples, and **eight of the twelve
are refuted under `C_aff`** [affine_rate] — in each case because one vertex is
affinely implementable from another, forcing the opposite comparison on that
edge. The remaining four are undecided. Where the descent can be tested, it
mostly fails.

**For, and this replaces an argument I had wrong.** An earlier draft claimed that
cycles are a zero-temperature artefact — that averaging over temperature destroys
them. That is false [maslov]. Write `c_β(a→b) = F_a(β)/F_b(β)` for the
fixed-`β` exchange rate and let `M_λ` be the `λ`-power mean over a temperature
prior `ρ`. Then the whole comparison is

```
A_λ(a,b) = ½ log( C_λ(a→b) / C_λ(b→a) ),   C_λ(a→b) = M_{−λ}(c_·(a→b)),    (6)
```

with `λ = ∞` the framework's infimum, **`λ = 1` the plain arithmetic average of
the transition rates**, and `λ = 0` the geometric mean. Two facts:

* At a *single* temperature the comparison is the total order of the scalar
  `F_·(β)`, since `c_β(a→b)·c_β(b→a) = 1`. **No pool cycles at any one
  temperature.** All non-scalar content comes from the `β`-dependence of the
  argmin.
* Cycles survive averaging. There are certified 3-cycles of **flat** signatures —
  which brief G *proves* cannot cycle at `λ = ∞` — living in bounded bands
  `λ ∈ (1.5211, 2.5376)` and `(0.8032, 1.3139)`, the second **a cycle at `λ = 1`
  exactly**. Cycles need `λ ≠ 0`; they do not need `λ` large, and some exist only
  at finite `λ`.

The comparison is a potential **iff** `λ = 0`, for a one-line reason: the
geometric mean is the unique power mean commuting with reciprocals, so
`C_0(a→b)C_0(b→a) = 1` — the trade becomes reversible. But the same identity
gives `S_λ = λ·Var_ρ/2 + O(λ³)`, so at `λ = 0` **every distance is zero and every
pair freely interconvertible**. The scalar complexity available there is the
potential of an empty theory.

> So non-scalarity is not an artefact of the infimum. It is a consequence of
> **irreversibility**: the comparison fails to be a potential exactly because a
> power mean of order `λ ≠ 0` does not commute with reciprocals. The infimum is
> the extreme case, not the source.

**Where the cycles come from, mechanically.** Brandão–Horodecki–Ng–Oppenheim–
Wehner govern conversion under thermal operations by a continuum of Rényi
`α`-free energies — the same `β`-indexed family. Their scalar is an infimum of a
**difference**, hence a partial order, acyclic by construction; ours is an
infimum of a **ratio**, hence a tournament, which can cycle [quantum]. And the
ratio is *forced*: these resources are unnormalised, so the difference diverges
like `β(Λ_a − Λ_b)` and the ratio is the only finite scalar. **Cycles are the
price of comparing resources with no common normalisation** — the situation for
maps, where nothing plays the role of a fixed particle number.

**The hedge, and it is the number to quote.** A single *unfitted* temperature
prior — uniform in `s = log β` on `[0,6]` — reproduces **99.39 %** of the tropical
order on the `F_11` arithmetic pool and **99.51 %** on `F_13`; a fitted prior
reaches **99.62 %** and transfers across `q` without loss, against an
information-theoretic ceiling of `99.90 %` (at least 42 edge-disjoint 3-cycles
must be misordered). Both beat the HodgeRank potential `ψ_opt` at 98.31 %
[maslov]. So the irreducibly non-scalar part of the comparison is **about four
parts in a thousand**, not the 1.7 % an earlier estimate suggested.

**The defensible claim, and the only one I would put in a paper:**

> Asymptotic convertibility is governed by no scalar complexity, because the
> comparison is a non-reversible aggregation over temperature. The failure is
> real but small: a single temperature prior recovers `99.4 %` of the order, and
> at `λ = ∞` rates around any preference cycle differ by a geometric-mean factor
> of at most two.

What this does **not** support is the story the cycle hunt was implicitly aimed
at — that cycles reveal something arithmetic. They do not; that reading is
refuted, and so is the "signature cycles descend to maps" reading, in eight of
the twelve cases where it can be tested.

**One genuinely quantum instance exists.** Replacing `Tr A^β` by a sandwiched
Rényi divergence against a background operator gives a non-spectral rate under
which three qubit states form a certified 3-cycle (`|curl A|/Σ|A| = 1` exactly,
40 digits) **whose decohered shadow is transitive** [quantum]. Two caveats keep
it honest: the plain quantisation `Tr A^β` is *empty* — it factors through the
spectrum identically, so a quantum cycle without a classical shadow cannot exist
— and every sandwiched profile still lies in the classical cone. What is
certified is the narrower physical statement that decoherence destroys that
cycle, not that the geometry is new.

**What would raise their status.** (i) A cycle in `C_aff` itself, on affine
orbits rather than signatures — the eight refutations above make this look harder
than it did. (ii) An operational reading of some `λ < ∞`, which would turn (6)
from a deformation into a physical phase diagram.

---

## 5. Status

Proved and independently re-verified: Propositions 1–2, Theorems A–B,
Corollary C, the power-mean identity (6) and its `λ = 0` reversibility, the
flatness/permutation-polynomial theorem, the invisibility of the smallest fiber
and of symmetry type, superadditivity of `C` and `C_aff` under Cartesian
products — and, for `C_aff`, **strict** superadditivity, witnessed by
`x² ⊗ x → x²+y` over `F_3` with a gap of at least `1/2` [affine_rate]. That
settles the one form of the synergy question no variational formula could have
decided: one affine processor draws a degree-2 part from `x²` and a linear part
from `x`, which neither factor supplies alone.

The first exact affine rates are in hand — `1/2` for `x² → x`, `x² → xy`,
`x² → x²+y²`, and `2/3` for `x²+y² → x` — and in the last the Fekete supremum is
attained at no `r ≤ 2`, so **no single-shot computation determines `C_aff`**.

Computed and certified: the cycle censuses and their witnesses; the flat-signature
band-cycles; the qubit cycle.

**Open and load-bearing:** whether `C_aff` cycles on affine orbits; whether any
`λ < ∞` is operational; the sharp finite-`λ` cycle constant off the flat locus
(there is no finite universal one — the flat-locus value grows like
`0.149·log(width of ρ)`); six of the twenty-five affine brackets over `F_3`.

**Withdrawn.** Three claims of earlier drafts. That the exchange metrics do not
fill `MET` — `C_4` is realisable exactly, and the obstruction rested on a stalled
search. That `C_aff ≤ C` — false, with the singular-processor counterexample now
in §1. And that cycles are a zero-temperature artefact — false, with certified
cycles at `λ = 1` now in §4. Each was refuted by the session sent to verify it.
