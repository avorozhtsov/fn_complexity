# Session brief K — is same-genus transitivity a theorem?

**Repo:** `fn_complexity`, work on a branch.

**Read first:** `research/sato_tate_limit/FINDINGS.md` §§ on the crossing
criterion, the lexicographic order and `level_lemma.py`. Then
`research/realizability/FINDINGS.md` §1 (the structure theorem) and §2.6.

## The question

Brief F established by **exhaustive computation over its libraries** that the
limiting comparison is transitive inside a genus and cycles only across genera.
Exhaustive computation over a library is not a theorem. Make it one, or find the
counterexample.

> **Prove: for measures with the same `α_max`, the midrange comparison
> `mid_τ(Ψ_μ − Ψ_ν)` is a total order.**

Everything downstream leans on this. It is why the `q → ∞` cycles must cross
genus classes, why the endpoint gap is "worth a whole level", and why the search
for witnesses (brief J) is confined to mixed-genus triples.

## What is already established, and where the gap is

**The crossing criterion** *(computed, 98.63%)*. From the two asymptotics

```
Ψ_μ − Ψ_ν  ≈  (m₂(μ) − m₂(ν))·τ/2        (τ → 0)
Ψ_μ − Ψ_ν  ≈  (t(ν) − t(μ))·log τ/τ      (τ → ∞)
```

the small-`τ` order is by `m₂` and the large-`τ` order is by `−t` (`t` = edge
exponent), so `Ψ_μ − Ψ_ν` changes sign **iff `m₂` and `t` are ordered the same
way**. Agrees with the computed sign-change count on 3954 of 4009 same-`α_max`
pairs; all 55 disagreements are exact `m₂` ties, where the criterion is vacuous.

**Why that nearly settles it.** A pair whose `Ψ`s do not cross is ordered by
pointwise domination, and a tournament all of whose edges are pointwise
dominations is transitive. So **a 3-cycle needs all three pairs to cross**, hence
needs `(m₂, t)` comonotone on the triple.

**The gap, precisely.** Inside a genus the comparison is the lexicographic order
`m₂ ascending, ties broken by t descending` on **745 of 765** pairs. The 20
exceptions all have `m₂` differing by 1 and `t` differing by a large amount — the
edge exponent overturning a small `m₂` gap. A one-parameter trade-off
`sign(Δm₂ − 0.668·Δt)` raises agreement to 731 of 735. **So the order is not
lexicographic, and the true rule is some trade-off between `m₂` and `t` that has
not been identified.** Find it.

## Three routes, in order of promise

1. **Identify the exact trade-off.** The fitted `0.668` is suspicious — it is
   close to `2/3`, and `2/3` would suggest a clean derivation. Compute
   `mid_τ(Ψ_μ − Ψ_ν)` as a functional of the pair and expand it in `(Δm₂, Δt)`
   near the tie locus. If the leading behaviour is `Δm₂ − c·Δt` with an explicit
   `c`, the order is a genuine linear functional and transitivity follows
   immediately, because any order by a linear functional of a fixed pair of
   statistics is total. **This is the most likely route to a theorem and should
   be tried first.**
2. **The level lemma, made exact.** `level_lemma.py` samples the three
   differences at `n` interior levels and solves the resulting LP over all
   argmax/argmin patterns with `|D| ≤ 1`, getting `0` at `n = 2` with both ends
   pinned (one genus) and `+1/4` with one end free (two genera). That is already
   the shape of the theorem: **with both endpoints pinned, two interior levels
   admit no cycle.** What is missing is the step from "two sampled levels" to
   "the whole `τ` axis" — i.e. a bound on how many independent levels a genuine
   `Ψ`-difference can have. Brief G's §4.1 is the finite-`q` analogue of exactly
   this ("never more than four interior extrema") and is also unproved. **The two
   questions may be the same question; try to prove one bound that serves both.**
3. **Look harder for a counterexample.** The 20 exceptions are the natural place.
   If same-genus cycling is possible at all, it lives among measures with
   `Δm₂` small and `Δt` large — construct such triples deliberately rather than
   sampling a library, and push to genus 8 where the `USp(2g)` edge exponents
   spread further. A same-genus cycle would be a *bigger* result than the
   theorem, because it would remove the genus constraint from brief J entirely.

## What a proof buys

* Brief F's headline becomes unconditional: the limiting comparison is
  non-transitive **exactly** across genera.
* Brief J's search space is provably confined to mixed-genus triples.
* It would be the limiting counterpart of brief G §2.6 ("every 3-cycle has
  `ψ`-spread `< log 2`") — the finite-`q` statement that cycles live inside a
  `φ`-class. Stating the two side by side is worth a section.

## Traps

* The `Ψ` endpoints must be supplied analytically (`Ψ(0) = 0` because
  `Σ_c a_c = 0` identically; `Ψ(∞) = α_max`), not read off a grid — grids
  underestimate `Ψ(∞)` badly, as an earlier low-resolution run showed.
* `τ` grids must reach at least 10: the certified `F_101` cycle uses contacts up
  to `τ = 9.3`.
* Weyl integration for `USp(2g)` at `g ≥ 6` needs care; verify every measure
  against its known moments (`E[α²] = 1` for `USp(2g)` and `SU(2)`) before using
  it. Brief F's `validate_library.py` already does this — reuse it.
* Do not reuse the coordinator's original seed library: 7 of its 10 measures are
  **unrealisable** as vertical Sato–Tate measures (Deligne semisimplicity plus a
  unipotent argument for elliptic pencils). Use brief F's validated library.

## Working standards

Proved vs computed marked on every statement. "Searched and did not find" is
acceptable only with the search space stated and an argument that it was
adequate. Verify headline numbers at 40 digits.

## Deliverable

`research/sato_tate_limit/TRANSITIVITY.md` in the house style (summary first,
sections marked proved/computed, corrections, open, files table), scripts and
CSVs alongside. **If the harness blocks you from writing that file, put its
complete body in your final message instead and say so.**
Build on `level_lemma.py`, `validate_library.py`, `symplectic_search.py`,
`lex_exceptions.py`.
