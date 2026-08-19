# Session brief F — do the cycles survive q → ∞? (Sato–Tate measures)

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (branch; commit or
stash first).

**Read first:** `research/curve_family_cycles/FINDINGS.md` — the whole file, and
especially the section "Why the regime permits it, quantitatively";
`research/m_and_e_and_a_c/FINDINGS.md` Notation and T2.1;
`research/session_briefs/D_quantum_mechanics_of_the_exchange_matrix.md` Part 0.

## The one thing brief B left open

> Not settled: whether cycles persist for a *fixed* genus as `q → ∞`. Everything
> here scales as `1/(√q log q)`, so the margins shrink; the counts stay in the
> hundreds in every sample taken, but at `q ≥ 31` the search is a sample and not
> a census, so that is evidence and not a theorem.

**This brief says the limit is computable exactly, with no curves in it at all**,
and that settling it is a finite analytic problem about a short list of
probability measures.

## The reduction

Brief B established, exactly:

```
a_c = q − N_c,   α_c = −a_c/√q,
Λ_f(β) = log( (1/q) Σ_c (N_c/q)^β ),      log Z_f(β) = (1+β) log q + Λ_f(β)
```

Substituting `β = τ√q` and using `(N_c/q)^β = exp(β log(1 + α_c/√q))
→ exp(τ α_c)`:

```
Λ_f(τ√q) → K_μ(τ) := log E_μ[e^{τα}],      the cumulant generating function
Ψ_f(τ) := Λ_f(τ√q)/τ → K_μ(τ)/τ =: Ψ_μ(τ)
```

where `μ` is the limiting empirical distribution of the normalised traces
`α_c = (N_c − q)/√q`. **Katz–Sarnak says what `μ` is:** for a family with big
monodromy group `G`, the `α_c` equidistribute according to the trace measure of
`G` in its standard representation. So `μ` is not an unknown — it is the pushforward
of Haar measure on `USp(2g)`, or `SU(2)`, or a CM torus, or an orthogonal group,
under the trace.

Two boundary values, both forced and both already recorded by brief B:

```
Ψ_μ(0) = K'_μ(0) = E_μ[α] = 0        (because Σ_c a_c = 0 identically)
Ψ_μ(∞) = ess sup supp(μ) = α_max     (the only value φ reads)
```

And the comparison, by brief D's Part 0, is the midrange:

```
f ≺ g  ⟺  mid_τ (Ψ_μ − Ψ_ν) < 0,      mid = ½(sup + inf)
```

> **Therefore: cycles persist as `q → ∞` if and only if the midrange comparison
> on the scaled cumulant generating functions of Sato–Tate measures is
> non-transitive.**

That is a self-contained question about three probability measures on a compact
interval. No fields, no polynomials, no exchange-rate solver. It is the theorem
brief B could not reach by searching.

## What to compute

1. **Build the library of `Ψ_μ`.** At minimum:
   * `SU(2)` — the semicircle on `[−2,2]`, `dμ = (1/2π)√(4−α²) dα`. (Non-CM
     elliptic; Sato–Tate.)
   * `U(1)` — the arcsine on `[−2,2]`. (CM elliptic, split case.)
   * The CM measure proper: `½δ_0 + ½·arcsine` at the inert primes.
   * `USp(4)`, `USp(6)` — Weyl integration formula, `g = 2, 3`.
   * `SU(2) × SU(2)` and `SU(2)`-with-multiplicity — split Jacobians, genus 2.
   * `O(2g)` / `SO(odd)` orthogonal families, for contrast with the symplectic ones.
   * The Weyl integration formulas are standard; `K_μ(τ) = log ∫ e^{τ tr} dHaar`
     is a one-dimensional integral for `SU(2)`/`U(1)` and a `g`-dimensional one
     for `USp(2g)` — do them numerically to 12 digits and check against the
     known moments (`E[α²] = 1` for `USp(2g)` and `SU(2)`, `E[α²] = 2` for the
     split CM case; T2.2's `m₂` table is the cross-check).

2. **Compute `Ψ_μ(τ) = K_μ(τ)/τ` on a τ-grid** with the two endpoints handled
   analytically (`Ψ(0) = 0`; `Ψ(∞) = 2g` for `USp(2g)`, `2` for `SU(2)`).

3. **Search for a 3-cycle in the midrange comparison** over the library, and
   over convex combinations and products of library measures (a product of
   families multiplies the MGF, hence *adds* the `K`s — so the reachable set is
   a convex cone and is easy to search).

4. **If a cycle is found:** it is a theorem that cycles occur for all large `q`,
   provided each of the three measures is realised by an actual family of curves
   of the given genus over `F_q` for infinitely many `q`. Check realisability
   family by family — this is where the arithmetic re-enters, and it is the only
   place it does. Give explicit families.

5. **If no cycle is found:** that is equally sharp, and it says the `F_11`/`F_13`
   cycles measure the *fluctuation* around Katz–Sarnak rather than the limit.
   Then quantify: the limiting comparison is transitive with gaps of size `G`,
   the finite-`q` fluctuation of `Ψ_f` around `Ψ_μ` is `O(q^{−1/2})` by
   Katz–Sarnak (or whatever the correct exponent is — derive it), so cycles die
   once `G ≫ q^{−1/2}`, i.e. beyond `q ≈ …`. Predict the `q` at which the census
   counts should collapse, and then **test the prediction** by running brief B's
   search at that `q`.

## Seed results — already computed, reproduce them first

I built the library at low resolution and ran the search. **No cycle was found**,
and the reason is structural and worth having before you start.

The `mid` matrix over the standard `α_max = 2` measures (`U(1)`, `SU(2)`, `CM`)
and the `α_max = 2g` measures (`USp(4)`, `USp(6)`, `SU(2)×SU(2)`) is a **total
order**, dominated by `α_max = 2g`: the endpoint gap between genus classes is
`O(1)` and so is the interior, so the endpoint wins and `φ` decides. Cycles
therefore cannot come from mixing genera. **Fix `α_max` and vary the measure.**

At fixed `α_max = 4` I then took ten genus-two monodromy measures — `USp(4)`,
`SU(2)²`, `U(1)²`, `SU(2)×U(1)`, `CM²`, `SU(2)×CM`, `U(1)×CM`, and the
multiplicity-two measures `2·SU(2)`, `2·U(1)`, `2·CM` (isogenous factors, i.e.
`α = 4cos θ`). Results:

```
total order:  USp4 < SU2xSU2 < SU2xCM < CMxCM < SU2xU1
                   < U1xCM < SU2_mult2 < U1xU1 < CM_mult2 < U1_mult2
3-cycles: 0 of 120 triangles
Hodge residual ||mid - grad psi|| / ||mid|| = 6.8e-2      (NOT zero)
crossing pairs: 3 of 45
```

Three things to take from this.

* **`Ψ` is additive under independent products**, trivially:
  `Ψ_{μ*ν} = (K_μ + K_ν)/τ = Ψ_μ + Ψ_ν`. So the reachable set is a **convex
  cone** and the whole library is generated by `{U(1), SU(2), CM, USp(2g)}` under
  addition. Searching convex combinations and products is therefore cheap and
  should be exhaustive over the cone, not sampled.
* **The decisive criterion is crossing.** If `Ψ_μ − Ψ_ν` has a fixed sign on
  `(0,∞)` then `mid` has that sign and the pair is ordered by pointwise
  domination — and a tournament in which every edge is a pointwise domination is
  automatically transitive. **Only 3 of 45 pairs cross:** `SU2×U1` vs `CM²`,
  `CM²` vs `2·SU(2)`, and `U(1)×CM` vs `2·SU(2)`. A cycle needs a triangle whose
  edges cross. **Start the search there and build outward.**
* **The curl is not zero — it is 6.8%**, the same order as the `≈ 7%` measured on
  random integer signatures (brief G). So the limiting comparison is *not* a
  potential; it is merely too close to one for the residual to flip a sign in
  this library. That is a quantitative near-miss, not a structural impossibility,
  and it says a cycle is plausible with a richer library.

**So the likely outcome is the negative branch, and the negative branch is the
valuable one.** Do not spend the session trying to force a positive. Spend it on:
(i) widening the library until either a crossing triangle cycles or the cone is
exhausted; (ii) then executing step 5 below properly — predicting the crossover
`q` and testing the prediction against brief B's census.

**Caveat on my numbers.** They were computed on binned measures (4000 bins) with
`τ` grids to 400 and the `τ = ∞` endpoint supplied analytically as `α_max`. The
`USp(6)` grid was coarse (`60³`). Margins of `10⁻²` are safe; anything below
`10⁻³` in my tables is grid noise. **Recompute at proper resolution before
building on any individual entry.** The qualitative facts — total order,
6.8% curl, 3 crossing pairs — were stable across two resolutions.

## Why a positive answer is still possible

`mid` is odd and positively homogeneous but not additive, so three functions
whose pairwise differences sum to zero pointwise can each have negative
midrange — brief B's model is `(−2,1,1), (1,−2,1), (1,1,−2)`. What that needs is
"deep dip at one temperature, mild rise at the others", i.e. crossing `Ψ`s. The
standard library is too rigid: its measures are nested by tail weight. Richer
candidates, in order of promise:

* **Measures with equal `α_max` and deliberately opposite tail shapes** — a
  measure concentrated near `±α_max` versus one concentrated near `0` with the
  same support. Both are reachable: `2·U(1)` is the former, `USp(4)` the latter,
  and they are the two extremes of the order above.
* **Non-connected monodromy** — the CM/inert mixtures already introduce an atom
  at `0`, and atoms are what bend `Ψ` at small `τ`. Try heavier atoms.
* **Higher genus at fixed `α_max`** — `USp(4)` versus a genus-3 family whose
  `α_max` is forced down to 4 by a constraint. If none exists, say so.

## The trap that will bite

Brief B's leading-order midrange law **gets the sign of one certified edge wrong
at `q = 11`**, because it drops `Λ_v/((1+β)log q) → log(max_v/q)/log q`, nominally
`1/(√q log q) = 0.13` but actually `0.23` at that field, and because
`log(1 + α/√q) ≈ α/√q` fails when `α/√q = 0.72`. **So the limiting law is not a
validation of the `F_11` cycle and the `F_11` cycle is not evidence for the
limiting law.** They are different regimes and must be kept apart in writing.
Verify the limit numerically at increasing `q` (`q = 101, 401, 1601, 6401` with
sampled families) and show `Ψ_f → Ψ_μ` before drawing any conclusion from `Ψ_μ`.

Other traps: the `β` horizon scales with `q` (grids to `360q`, brief B);
`exchange_rate` is good to `1e−13`, tie threshold `1e−10`; genus one gives only
`gcd(4,q−1)+1` distinct signatures so it is useless for a *census* but is fine
here, because here the measure is the object, not the signature.

## Success criterion

Either three named Sato–Tate measures with a computed midrange cycle, their
margins, the curves realising them and the resulting theorem statement; or a
proof/strong numerical demonstration that the limiting comparison is transitive,
with the crossover `q` predicted and then tested.

## Reproduce / build on

`research/curve_family_cycles/regime.py` (the `Ψ` construction and the two-scale
table), `common.py`; `research/m_and_e_and_a_c/t2_1_genus_scaling.py` (the
`USp(2g)` extreme-value work, which already has the Weyl integration machinery).
