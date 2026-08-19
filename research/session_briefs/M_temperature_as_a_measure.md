# Session brief M — temperature as a measure, and the death of cycles

**Read first:** `research/maslov/lambda_seed.py` and run it;
`research/realizability/FINDINGS.md` §1 and §2.6;
`research/session_briefs/D_quantum_mechanics_of_the_exchange_matrix.md` Part 3
item 1 (Maslov dequantisation, where `ħ = 1/t` — this brief is that item made
concrete).

## The question

The exchange rate is an **infimum** over inverse temperature. An infimum is a
zero-temperature operation. What happens if one *averages* over `β` instead —
i.e. puts a distribution on temperature and reads `M` as an average rather than
an extremum?

## The seed result, already computed — reproduce it first

Put `s = log β`, let `ρ` be a probability density in `s`, and soften both
extrema:

```
softmax_λ(f) = (1/λ) log ∫ e^{λf} ρ ,   softmin_λ(f) = −(1/λ) log ∫ e^{−λf} ρ
A_λ(a,b)     = ½(softmax_λ + softmin_λ) of f = u_a − u_b
```

`A_λ → A = mid` as `λ → ∞` (the current framework). And as `λ → 0` **both**
soft-extrema tend to `∫ f ρ`, so

```
A_0(a,b) = ∫ (u_a − u_b) ρ = Ψ(a) − Ψ(b),   Ψ(a) = ∫ u_a ρ
```

— an **exact potential difference, for every `ρ`**. Hence:

> **At `λ = 0` the comparison is a total order with no cycles; at `λ = ∞` it has
> them. Every cycle in this framework lives above a critical `λ_c`.**

Measured on the known 3-cycle `{(6,3,3),(7,2,1),(6,5,1)}`: no cycle at
`λ ≤ 100`, cycle from `λ ≈ 300`, and `λ_c = 230, 287, 323, 350, 185` for uniform
priors on `s ∈ [−8,8], [−12,12], [−16,16], [−6,20], [−20,6]`.

**So the answer to "is `M` an average of transition probabilities?" is: it is
the `λ → ∞` member of a family whose `λ → 0` member is exactly an average — and
averaging destroys precisely the non-scalar content.**

## What to establish

1. **Is `λ_c` intrinsic?** It moves with `ρ` by a factor ~2 above. Find a
   formulation in which it is not arbitrary — e.g. normalise `λ` by the width of
   `ρ`, or use a canonical `ρ`. **A prior-free critical exponent is the result to
   aim for.** If none exists, say so and report the scaling law `λ_c(ρ)`.
2. **Is there a canonical `ρ`?** `Ψ(a) = ∫ u_a ρ` is a *bona fide* scalar
   complexity for every `ρ`, so this construction manufactures scalar
   complexities at will. Which `ρ` gives the best one — e.g. the one whose `Ψ`
   best reproduces the tropical order (brief E measured `½ log φ` at 98.3% of
   ordered pairs; beat it, or fail to)?
3. **The phase diagram.** As `λ` grows, cycles appear. Do they appear one at a
   time, or in a burst? Is `λ_c` the same for all cycles in a pool? Is there an
   order parameter (the curl fraction `‖curl A_λ‖/‖A_λ‖` is the obvious
   candidate — it is exactly `0` at `λ = 0` by the potential identity).
4. **Does the structure theorem deform?** Brief G bounds the tropical defect by
   `(log 2)/2`, sharp. What is the bound on `|A_λ − dΨ|`? It must vanish as
   `λ → 0` and reach `(log 2)/2` as `λ → ∞`.
5. **Transition probabilities, literally.** The user's phrasing was "`M` as an
   average of transition probabilities". At fixed `β` there is a Gibbs
   transition kernel on the fibers. Write down that kernel, define the fixed-`β`
   comparison it induces, and check whether `∫ (fixed-β object) ρ(β) dβ`
   reproduces `A_λ` for some `λ`, or is a genuinely third construction.

## Why it matters

It reframes every cycle result in this project. If cycles exist only in the
tropical limit, then "complexity is not a scalar" is a statement about the
`inf`, not about the resources — and the framework has a one-parameter family of
scalar complexities available at any `λ < λ_c`. That is either a serious
limitation of the cycle narrative or the right way to state it; decide which,
with numbers.

## Traps

The `λ → 0` transitivity is an identity, not a measurement — do not "verify" it
numerically and report a tolerance. `soft` extrema overflow badly; subtract the
max/min before exponentiating (the seed script does). β horizons to `10³`;
tie threshold `1e-10`.

## Deliverable

`research/maslov/FINDINGS.md`, house style. **If the harness blocks the write,
paste the full body in your final message.**
