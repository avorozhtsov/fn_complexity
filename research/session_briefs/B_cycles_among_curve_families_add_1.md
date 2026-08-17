# Addendum 1 to brief B — where on the β axis a cycle could possibly live

Read together with `B_cycles_among_curve_families.md`. Written in a later
session; **everything here is an asymptotic derivation, not a computation.** It
reproduces both recorded propositions of T2.1 and T2.2 as special cases, which is
the only evidence it has. Verify §2–§3 numerically before letting it steer the
search — the check is cheap and is spelled out in §6.

The point of the addendum: the brief says "find pairs where the computed
comparison disagrees with `φ`". That instruction is right but under-specified,
because **the first correction to `φ` is still a scalar**, so it produces no
cycles either. The search has to aim one order deeper, and this note says where.

---

## 1. Two exact blind points, and what they force

```
Z_f(0) = #{fibers}  = q       for every f : A² → A¹ over F_q
Z_f(1) = Σ_c N_c    = q²      for every f
```

Both are identities, not approximations. Consequences worth having explicit:

* `R(β) := u_f(β)/u_g(β)` equals **1 at β = 0 for every pair**, where
  `u = log Z`. Hence `C(f→g) = inf R ≤ 1` and `C(g→f) = 1/sup R ≤ 1`
  automatically, and

  ```
  d(f,g) = log(sup_β R / inf_β R)
  ```

  is literally the oscillation of `R`. A cycle is a statement about the *shapes*
  of three such oscillations, not about their sizes.
* The leading correction carries the factor `β(β−1)`, which vanishes exactly at
  the two blind points. That is not a coincidence — `β = 0` is `Σ_c 1 = q` and
  `β = 1` is `Σ_c N_c = q²`, the two constraints that kill the zeroth and first
  moments. The interval `(0,1)` is therefore **the only place where a fibration
  can beat the flat reference**, and `β* = √2 − 1` is the minimum inside it.

## 2. The deviation function and its two scales

Write `M_f = max_c(−a_c)`, `T_k = Σ_c a_c^k`, `m₂ = T_2/q²`. For `β ≪ √q`,

```
Z_f(β) = q^{β+1}·(1 + Σ_{k≥2} C(β,k)·T_k/q^{k+1})
u_f(β) = (β+1)log q + Σ_{k≥2} C(β,k)·T_k/q^{k+1} + …
```

Define the **deviation function** `D(β) = (R(β) − 1)·q·log q`. Then

```
D(β) ≈ ψ₂(β)·Δm₂ + O(q^{−1/2})      for β = O(1),    ψ₂(β) = β(β−1)/(2(β+1))
D(∞)  =  ΔM                          exactly, to this order
```

with `Δm₂ = m₂(f) − m₂(g)` and `ΔM = M_f − M_g`. `ψ₂` is negative on `(0,1)`,
attains `ψ₂(β*) = −(3−2√2)/2 = −0.08578644` at `β* = √2 − 1`, and grows without
bound afterwards.

**The two scales.** In the common denominator `q log q`: the endpoint
contributes `ΔM`, an **integer** of size up to `≈ 4√q`; the interior contributes
`|ψ₂|·|Δm₂| ≤ 0.0858·|Δm₂|` with `Δm₂ = O(1)`. So the endpoint is larger by a
factor `√q` unless the maxima coincide. This — not any property of `L` — is why
everything degenerates to endpoints. `L` is only the extreme case (`M_L = 0`, so
`ΔM > 0` always).

**Crossover.** The moment expansion is valid while `β|a_c|/q ≪ 1`, i.e.
`β ≪ √q`; the endpoint expansion starts when `β log(max N)` dominates. Neither
description covers `β ~ √q`, where every term of the moment sum is the same
order. For `q = 10³–10⁶` that is `β ≈ 30–10³` — inside the range the brief
already warns must be gridded densely, and worth watching for its own sake.

## 3. The first correction to `φ` is still a scalar

Take the two candidates for the infimum (`β*` and `β = ∞`) and the two for the
supremum (`0` and `∞`). For `Δm₂ > 0`:

```
C(f→g) = 1 + min(−0.0858·Δm₂, ΔM)/(q log q)      inf at β* when ΔM > 0
C(g→f) = 1 − max(0, ΔM)/(q log q)                sup at ∞
```

Note the two directions of one pair are attained at **different temperatures** —
one interior, one at infinity. That asymmetry is the real content of a non-flat
pair, and it does not occur against `L`.

Substituting into `f ≺ g ⟺ C(f→g) < C(g→f)` and checking the sign cases gives

> **f ≺ g ⟺ φ̃(f) < φ̃(g),  where  φ̃ = M_f − ((3−2√2)/2)·m₂(f)**
>
> equivalently, the comparison flips exactly when `Δm₂ = (6 + 4√2)·ΔM`,
> `6 + 4√2 = 11.6568542`.

**This is still a total order, so it still forbids cycles.** `φ̃` is the
endpoint invariant `φ` (monotone in `M`) with the second moment added as a
tie-breaker — one order deeper, same obstruction. Searching for
`φ`-disagreements will surface pairs where `m₂` overturns `M`, and those pairs
are *not* cycle candidates; they are ordered by `φ̃`.

**Consistency check (why I believe the algebra).** Setting `g = L`
(`M_L = 0`, `m₂(L) = 0`) reproduces both recorded propositions with no extra
input:

```
C(L→f) = 1 − M/(q log q)          = log q/log(q+M) + O(q^{−2})     [T2.1]
C(f→L) = 1 − (3−2√2)·m₂/(2q log q)                                 [T2.2]
```

## 4. Where scalarity can actually break — three candidates, ranked

A cycle needs at least one edge on which the `φ̃` gap is smaller than the terms
`φ̃` drops. Three places, in decreasing order of how promising they look:

**(a) Exact `(M, m₂)` ties, decided by `m₃`/`m₄`.** This is the regime the
recorded `q = 211` pair already lives in: identical image size, identical largest
fiber, identical `m₂`, all endpoint probes agreeing to `0.000e+00`, interior
probes separating by `2.42e−7` and predicted in advance by
`c₃Δm₃ + c₄Δm₄ = 2.32e−7`. Exact `m₂` ties are *constructible*, not accidental:
for `f = y² − P(x)`,

```
m₂ = ν(P)/q − 1,      ν(P) = #{(x,x′) : P(x) = P(x′)}
```

so `m₂` lies in `(1/q)·Z` and equality is an integer condition on `ν(P)`. The
minimal nonzero gap `|Δm₂| = 1/q` moves `C` by `≈ 0.0858/(q² log q)`, which at
`q = 211` is `3.6e−7` — the same order as the observed `m₃`-driven `2.42e−7`.
**So `q ≈ 200` is where the two effects are comparable, and the window closes as
`q` grows.** Do not push to large `q` here; push to designed degeneracy at small
`q`.

**(b) Multiplicity of the largest fiber.** At large `β`,
`u_f = β·log max N_f + log μ_f + …` where `μ_f` is how many fibers attain the
maximum. So

```
R(β) = L_f/L_g + (log μ_f − (L_f/L_g)·log μ_g)/(β·L_g) + O(β^{−2}),   L = log max N
```

When the maxima tie (`L_f = L_g`) the multiplicity becomes the leading endpoint
datum and enters at order `1/(β log q)` — which for moderate `β` is much larger
than any moment term. `C(L→f)` cannot see `μ` at all (it is a pure `β = ∞`
limit), so this is genuinely new information available only to non-flat pairs,
and it is the term most likely to be doing the work in case (a) as well.
**Unverified; check it first, it is one line of code.**

**(c) The crossover region `β ~ √q`.** Neither expansion describes it, so no
scalar prediction exists there by construction. If an infimum is ever attained
at `β ≈ √q` rather than at `β*` or `∞`, that pair is the best cycle candidate in
the whole search. Log the argmin `β` of every computed rate and histogram it —
if the histogram has mass anywhere other than `{0} ∪ [0,1] ∪ {∞}`, that is the
finding.

## 5. Consequence for how to search

Random families are the wrong sample. In a random pair `ΔM ≠ 0` with probability
close to 1, so the endpoint wins by `√q` and the comparison is `φ`. The brief's
"look for `φ`-disagreements" will then find only `m₂`-driven flips, which are
ordered by `φ̃`.

**Search over designed degeneracies instead:**

1. Compute `(M, ν(P), μ)` for a large pool of families at a small `q`
   (`q ≈ 101, 211, 503`). These are cheap integer statistics — no exchange rates
   needed.
2. Bucket by `(M, ν(P))`. Only buckets with ≥ 3 members can contain a cycle.
   Genus ≥ 2 is required for variety: 400 random elliptic fibrations give only
   `gcd(4,q−1)+1` distinct signatures, so elliptic pencils collapse into the same
   bucket for a trivial reason and will not produce anything.
3. Within a bucket, compute the full exchange matrix on a dense grid to
   `β ~ 10³` and search for 3-cycles. Record the argmin `β` of every rate.
4. Only then widen to near-ties (`|Δm₂| = 1/q`, `|ΔM| = 1`) and see whether the
   `m₃` term ever overturns them.

Step 2 is the whole point: it turns a `O(n²)` blind sweep into an integer
bucketing that isolates the only pairs where the interior is not overwhelmed.

## 6. What to verify before trusting any of this

Cheap, and it either validates §2–§3 or kills them:

* Take 20 curve families at `q = 211` and `q = 1009`. Plot measured
  `(C(f→g) − 1)·q log q` against the predicted `min(−0.0858·Δm₂, ΔM)`. The
  claim is agreement to `O(q^{−1/2})` relative.
* Check that `φ̃ = M − 0.08578644·m₂` reproduces the observed `≺` on every pair
  in that pool. Every exception is a cycle candidate and should be recorded with
  its argmin `β`.
* Check the multiplicity term (b) directly: build two signatures with identical
  `{N_c}` except that the maximum is attained twice instead of once, and confirm
  the `1/(β log q)` behaviour.

Precision, from the parent brief: `exchange_rate` is good to `~1e−13`, so treat
differences below `1e−10` as ties. Nothing in §4(a) is safe below that floor —
`2.42e−7` is three orders above it, `0.0858/(q² log q)` at `q = 1009` is
`1.2e−8`, which is only two.

## 7. If no cycle is found

The addendum sharpens what the negative result would say. Not merely "the regime
forbids it", but:

> The comparison of curve families over `F_q` is governed to two orders by the
> scalar `φ̃ = max_c(−a_c) − ((3−2√2)/2)·m₂`, and the corrections that could
> overturn it are `O(q^{−1/2})` relative to the `m₂` term. Cycles are therefore
> confined to pairs whose `φ̃` gap is below that threshold, i.e. to exact
> `(M, m₂)` coincidences — of which there are `N` in a pool of `n` families at
> `q = …`, and none of them closes a triangle.

That is a quantitative theorem-shaped statement, and it is what brief C needs:
"no scalar invariant" is false at this order, so brief C's claim must be
formulated as *no scalar invariant beyond the `φ̃` truncation* — or abandoned if
`φ̃` turns out to be exact.
