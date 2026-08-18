# Session brief D — is the exchange matrix a quantum system?

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (branch; commit or
stash first).

**Read first:** `research/m_and_e_and_a_c/FINDINGS.md` (Notation, T1.2, T1.3,
T1.4, T1.5), `paper_finite_fields_maps/docs/exchange_positivity_and_weil.md`
(§1, §2, §3, §3b), and the output of
`research/m_and_e_and_a_c/gauge_decomposition.py`.

## The three questions this brief answers

1. Can the Weil matrix `E` be phrased with a Hamiltonian, commutators and a
   density matrix? Is there anywhere to go from there?
2. `M` is real and non-symmetric. Can it be mapped into the quantum world —
   e.g. is `MᵀDM` or `MᵀM` a density matrix, or related to one?
3. Can `M` be redefined and *quantised*? Are there methods for quantising
   stochastic matrices, so they become classical limits of quantum processes?

Short answers: (1) yes, and the formulation is the hundred-year-old
Hilbert–Pólya programme in which Weil positivity is already the statement that
a metric operator is positive definite — nothing here shortens that road, but
one theorem already in the repo turns out to be a statement in that language;
(2) yes but the naive version is empty, and the non-empty version is a gauge
decomposition that has immediate consequences for brief B; (3) yes, there are
four standard methods, and the informative fact is *which part of `M` each one
destroys*.

---

## Part 0 — what is already settled

Verified in `gauge_decomposition.py`; do not re-derive.

Write `L(a,b) = −log C(a→b)`, which by Theorem 1 is `sup_β g` for
`g = u_b − u_a`, `u_a = log log Z_a`. Split it:

```
L = S + A,     S(a,b) = (L_ab + L_ba)/2,     A(a,b) = (L_ab − L_ba)/2
```

Because `sup(−g) = −inf g`:

```
S = (max g − min g)/2 = d/2      the half-RANGE of g      (the metric)
A = (max g + min g)/2            the MIDRANGE of g        (a 1-form)
```

**The two halves of the exchange matrix are the range and the midrange of one
and the same function.** Four consequences.

* **(a) `a ≺ b ⟺ A(a,b) > 0`.** The comparison is nothing but the sign
  pattern of the antisymmetric part. `S` — the whole content of the metric,
  Theorem 2, the `l2`-distortion, the negative-type certificate — is
  *comparison-blind*.
* **(b) A triangle is a strict 3-cycle iff `|curl A| = Σ|A|` over its edges and
  `min|A| > 0`,** where `curl A = A(a,b) + A(b,c) + A(c,a)`. The ratio
  `|curl A|/Σ|A|` lies in `[0,1]` and hits `1` exactly on cycles. Measured:
  `1.0000000000` on the known cycle `{6,3,3}→{7,2,1}→{6,5,1}`, against
  `0.0176`, `0.0078`, `7.03·10⁻⁴` and `0.0000000000` on control triangles — the
  last being a triangle wholly inside the endpoint regime, where (c) forces the
  curl to vanish exactly.
* **(c) If both infima of a pair are attained at an endpoint, `A` is exact:**

  ```
  A(a,b) = ψ(b) − ψ(a),      ψ = ½ log[ log(#fibers) · log(max fiber) ] = ½ log φ
  ```

  verified to `1.1·10⁻¹⁶` on 347 endpoint pairs, against a median deviation of
  `1.215·10⁻²` on 1423 interior pairs. An exact 1-form has zero curl, so **no
  cycle can live in the endpoint regime** — which is the endpoint-regime theorem
  of brief B, re-derived as *zero curvature*, with `φ` emerging as the potential
  rather than being guessed.
* **(d) No arbitrage says `Σ S ≥ |Σ A|` around every loop:** the metric part
  dominates the flux part. On the known 3-cycle, `0.201398 ≥ 0.015814`.

**Correction to FINDINGS T1.3.** The minimal five-point certificate has PSD
threshold `t* = 1.0918390486`, not the `0.124` quoted there — that value belongs
to the *superseded* thirteen-point family (`cert13` in
`t1_2_part1_thresholds.csv`). The sentence "the certificate's PSD threshold is
only `t* = 0.124`, while artificially hill-climbed metrics reach `t* ≈ 12–17`"
must be re-stated: the gap is `1.09` against `12–17`, an order of magnitude
narrower than advertised. `xᵀDx = 9.8126948851·10⁻⁴` reproduces unchanged.

---

## Part 1 — `E` in quantum-mechanical language

The formulation exists and is standard. The Hilbert space is
`L²(ℝ_{>0}, d×x)`; dilations act unitarily by `(D_a ψ)(x) = a^{−1/2} ψ(x/a)`,
and **the `a^{−1/2}` is the `½` of the critical line**. The generator is
`H = −i(x d/dx + ½)`, self-adjoint, the Berry–Keating `xp` with
`[log x, H] = i`. Its spectral representation is the Mellin transform
restricted to `Re s = ½`. So:

| exchange theory | quantum mechanics |
|---|---|
| dilation portfolio `T_g = ∫ g(a) D_a d×a` | `ĝ(H)`, functional calculus of the Hamiltonian |
| monotone `Z_a(β)`, `β` real | imaginary time / inverse temperature |
| `Z_a(½ + iγ)` | real time; a transition amplitude `⟨γ\|a⟩` |
| signature `a` | state `\|a⟩ = Σ_i δ_{a_i}` |
| `E_ab = Σ_ρ Z_a(ρ) conj(Z_b(ρ))` | `⟨a\| Π_T \|b⟩`, `Π_T` the projector on the first `N` zero modes |
| `E/N` | the maximally mixed state on the zero subspace |
| `M ↔ E` | **Wick rotation `β ↦ ½ + iγ`** |

That last row is the honest summary of the whole `M`/`E` relation: `M` is the
Euclidean (thermodynamic, Gibbs, tropical) section of the monotone family and
`E` is the Lorentzian (unitary, Hilbert) section, and they are the same analytic
object read on two lines in the `s`-plane.

**Where RH sits.** If RH fails, `E_ab = Σ_ρ Z_a(ρ) conj(Z_b(1−ρ̄))` and the
involution `θ : s ↦ 1−s̄` is no longer the identity on zeros. `E` stops being
`⟨a|Π|b⟩` for a projector and becomes `⟨a|Θ|b⟩` for an indefinite `Θ`. This is
exactly Krein-space / pseudo-Hermitian quantum mechanics: `H` is
`η`-pseudo-Hermitian, `H† = ηHη^{−1}`, and its spectrum is real iff `η` can be
chosen positive definite. **Weil positivity is the positivity of the metric
operator**; nothing in this project shortens the distance to it. Relevant
literature to read before writing anything: Berry–Keating on `xp`, Connes'
trace-formula reformulation, Bender–Brody–Müller (2017) on a `PT`-symmetric
Hamiltonian for the zeros, de Branges spaces. Assume all of it is known to a
referee.

**The one thing this project looked like it contributed to that picture —
settled, and negatively.** T1.5 Theorem A says `Z_{λa}(ρ) = λ^ρ Z_a(ρ)`, hence
`E_{λa,λb} = λ^{2Re ρ} E_{ab}`, *and the factor is constant across zeros only
because every `Re ρ = ½`* — verified to `6·10⁻¹⁴`, broken to `1.3·10⁻¹` by
moving 60 of 1200 zeros to `Re = 0.7`. In the quantum language this reads:

> **RH ⟺ the dilation group acts on the span of the zero modes by a scalar,
> i.e. the Weil Gram matrix is covariant under dilations.**

**This is false, and the error is instructive.** There are two matrices:

```
W_ab = Σ_ρ Z_a(ρ) Z_b(1−ρ)          the Weil pairing
G_ab = Σ_ρ Z_a(ρ) conj(Z_b(ρ))      what it becomes under RH
```

`riemann_hypothesis_exchange_matrices.md` defines the pairing as `W`; the
research thread computed `G` throughout, which is harmless while the zeros used
are on the line and wrong the instant one is moved off it. Since
`Z_{λa}(s) = λ^s Z_a(s)` and `ρ + (1−ρ) = 1`, **`W` is homogeneous of degree
exactly 1 for any multiset of zeros whatever.** Scale covariance is the
*functional equation*, not RH. Verified to `10⁻¹²` with zeros in quadruples
`½ ± δ ± iγ` out to `Re ρ = 0.99`; `G` over the same zeros breaks by a factor of
1319 at `λ = 210`. The earlier counterfactual also moved zeros singly, which
breaks closure under `ρ ↦ 1−ρ` and so is not a legal zero set for the pairing.

The salvage is a measurement, not a criterion. What responds to `δ` is
*positivity* of `W` — Weil's actual criterion — and on atomic measures it has a
detection threshold: the 16-signature family of T1.4 stays PSD even at
`Re ρ = 0.99`, and a 40-signature family first goes negative between `Re = 0.7`
and `Re = 0.8`. **The `Re = 0.7` of the original counterfactual sits just below
detection.** Since atomic measures are inadmissible this bounds the truncation
and says nothing about `ζ` — but it closes the question of whether signature
families could be a numerical RH probe. They cannot.

Record: `research/m_and_e_and_a_c/t1_5_scale_covariance.py`, FINDINGS T1.6.

---

## Part 2 — `M` in quantum-mechanical language

**The trap, first.** `MᵀM` is positive semidefinite for *every* real matrix, so
`MᵀM/tr(MᵀM)` is a density matrix for free and carries no information about `M`
whatever. Likewise `MᵀDM ⪰ 0` for any `D ⪰ 0`. Anything built this way will
look like quantum mechanics and mean nothing. Do not build the programme on it.
The question worth asking is not "is some PSD matrix available" but "**which**
PSD matrix is canonical, and what does its failure to be PSD cost".

**The canonical one already exists in the repo.** It is the Gibbs family

```
ρ_t = e^{−t·d} / tr e^{−t·d},        e^{−t·d}_{ab} = (C(a→b)C(b→a))^t
```

a state whose Hamiltonian is the exchange metric and whose inverse temperature
is `t`. By Schoenberg this is a legitimate quantum state for every `t > 0` iff
`d` is of negative type. Theorem 2 says it is not. Therefore:

> **The exchange data admits a Hilbert-space description only above a critical
> inverse temperature `t*`.** Below it, `ρ_t` has a negative eigenvalue — a
> negative probability. T1.2 measured `t*` and found exactly one sign change of
> `λ_min(t)` in every one of 6886 cases: a single, sharp classical/quantum
> transition, with `λ_min(t)` as order parameter.

This reframes T1.2 from a technical positivity scan into a phase diagram, and
"positivity is worst in the middle" becomes a statement about where the
transition is sharpest. It costs nothing to state — the numbers exist.

**`t = ½` is a distinguished point, and it is the quantisation point.** For a
stochastic matrix `P`, Szegedy's quantum walk is governed by the *discriminant*
`√(P ∘ Pᵀ)` — the entrywise geometric mean of forward and backward rates. For
the exchange matrix that is

```
√(M ∘ Mᵀ) = e^{−d/2},
```

i.e. exactly the `t = ½` member of the Gibbs family. So `t* ⋛ ½` decides, family
by family, whether the exchange market has a Szegedy quantisation whose
discriminant is a Gram matrix. Measured:

| family | `n` | `λ_min(√(M∘Mᵀ))` | `t*` |
|---|---:|---:|---:|
| minimal certificate | 5 | `−2.62·10⁻⁴` | `1.0918` |
| cert13 (superseded) | 13 | `+1.53·10⁻³` | `0.1241` |
| greedy | 25 | `−2.13·10⁻²` | `14.742` |
| greedy | 30 | `−2.77·10⁻²` | `15.862` |
| greedy | 40 | `−3.21·10⁻²` | `16.713` |

**The minimal five-signature certificate already fails at `t = ½`.** That is a
crisp, falsifiable, previously unstated fact about the objects of the second
paper.

**The right decomposition of `M` is additive, not multiplicative.** Part 0: `L`
splits into a symmetric part `S = d/2` — a Hamiltonian — and an antisymmetric
part `A` — a lattice gauge field. In that reading `S` is the energy, `A` is the
vector potential, `curl A` is the magnetic flux through a triangle, and
**the 3-cycles of the first paper are Aharonov–Bohm phases**: loops through
which the flux is so large that it saturates `|curl A| = Σ|A|`. The
endpoint-regime theorem is the statement that the field is a pure gauge,
`A = dψ` with `ψ = ½ log φ`, so the flux vanishes and the tournament is
transitive.

In stochastic thermodynamics the same `A` is the *affinity*, its cycle sums are
the entropy production, and `A = dψ` is detailed balance. This makes the
obstruction to quantising `M` completely explicit and is the answer to Q3 below.

---

## Part 3 — quantising `M`, and what each method destroys

Four standard routes. None is exotic; the content is the accounting.

1. **Maslov dequantisation / idempotent analysis (Litvinov–Maslov).** The `Lᵗ`
   interpolation of §3 of the two-positivities note *is* this, with `ħ = 1/t`:
   `C(a→b)^{−1} = lim_{t→∞}(∫(φ_b/φ_a)^t dμ)^{1/t}`, and `(max,+)` at `ħ → 0`.
   So the framework already has a quantisation parameter with a principled name,
   and the exchange rate is a classical limit by construction rather than by
   analogy. **Keep this family strictly separate in writing from the Schoenberg
   family `e^{−td}`** — they share the letter `t` and nothing else, and the note
   currently invites the confusion even while disclaiming it.
2. **Rokhsar–Kivelson / stochastic-matrix form (Castelnovo–Chamon–Mudry–Pujol).**
   The textbook way to make a Hamiltonian out of a stochastic matrix: for a
   *reversible* `P` with stationary `π`, `H = 1 − π^{1/2} P π^{−1/2}` is
   Hermitian and PSD, with zero-energy ground state `Σ_x √π(x)|x⟩` and
   relaxation spectrum. It requires detailed balance, i.e. `A = dψ`. **So the
   standard quantisation of `M` exists only after discarding `A` — precisely the
   part that carries the comparison, the cycles and the entire arithmetic
   content.** This is a negative result of the useful kind: the absence of a
   Hermitian quantisation of `M` is not a gap to be filled but a theorem, and
   `d` is the quantitative obstruction. It also vindicates §3b(ii)'s claim that
   the correspondent of an off-line zero is irreversibility.
3. **Szegedy walk.** Works for any stochastic `P`, reversibility not required,
   which is why it is the right tool here. Needs a stochastic matrix built from
   `M` — the honest choice, and a decision the session must make and defend.
4. **Lindblad / classical embedding.** Any stochastic generator is the diagonal
   restriction of a GKSL generator; "quantising" is choosing the off-diagonal
   coherences, and the moduli of those choices is the object to study. Most
   open-ended, least likely to produce a theorem quickly.

**And the one genuinely arithmetic quantum system: Bost–Connes.** §3 of the
two-positivities note observes that `⨂_p P_{p,∞} = {1,2,3,…}` with partition
function `ζ` — this is the primon gas / Riemann gas (B. Julia, 1990), whose
non-commutative completion is the Bost–Connes C\*-dynamical system. Its
partition function is `ζ(β)` and **it has a KMS phase transition at `β = 1`**.
That is the same `β = 1` as the note's caveat "the Euler product converges only
for `Re s > 1`, so the critical strip is exactly where the monotone diverges".
The caveat is therefore not a defect of the framework but a known physical
phenomenon with a developed theory: below `β = 1` the Gibbs state does not
exist and is replaced by the KMS states, which Bost–Connes classify and on which
the idele class group acts. **This is the most promising single item in the whole
brief**, because it says what to do when `inf_β Z_a/Z_b` is meaningless: replace
the infimum of ratios by a KMS condition.

---

## The plan

Five stages, ordered by certainty of payoff. Stages 1–3 are safe and will
produce publishable propositions; stage 4 is a restatement whose value is
unknown; stage 5 is the only one aimed at RH and is a multi-session programme.

**Stage 1 — the gauge decomposition. DONE, and it answered brief B.**
Feeding `|curl A|/Σ|A|` into brief B's search found a certified strict 3-cycle
among genus-2 hyperelliptic pencils over `F_101` — 9 of them among 117480
triangles, margins to `1.16·10⁻⁴` against a `10⁻¹⁰` floor, verified three ways,
and cycles again at `q = 211`. Record:
`research/session_briefs/B_cycles_among_curve_families_add_2.md`,
`analysis/certify_curve_family_cycle.py`. The write-up below is still owed to
the second paper. Original text of the stage:

**Stage 1 — write up the gauge decomposition (certain, one session).**
Proposition: `L = S + A` with `S` the half-range and `A` the midrange of
`u_b − u_a`; `a ≺ b ⟺ A > 0`; strict 3-cycle `⟺ |curl A| = Σ|A|` with
`min|A| > 0`; `A = dψ`, `ψ = ½ log φ`, in the endpoint regime, hence no cycles
there; no arbitrage `⟺ Σ S ≥ |Σ A|`. All four are one-line proofs and all four
are verified. *Deliverable:* a short section for the second paper, and a smooth
`[0,1]`-valued search objective for brief B where sign-hunting currently gives
none. **Feed `|curl A|/Σ|A|` and `|A − dψ|` back into brief B's curve-family
search immediately** — they replace the `φ`-violation heuristic with a
continuous, differentiable one and they say precisely how far a family is from
the regime that forbids what brief B is looking for.

**Stage 2 — the phase diagram. DONE.** `t1_2_phase_diagram.py`: one sign change
of `λ_min(t)` per family on a 400-point log grid over `[10⁻³, 10³]`, `t*` marked
against the Szegedy point `t = ½`, CSV and figure written. The `t* = 0.124`
error in FINDINGS T1.3 is fixed. Original text:

**Stage 2 — the phase diagram (certain, one session).** Fix the `t* = 0.124`
error in FINDINGS T1.3. Plot `λ_min(t)` for the minimal certificate, the greedy
families and the random ensemble; state the single-sign-change observation as
the conjecture it is; identify `t*` as a critical inverse temperature and
`λ_min` as an order parameter. *Success criterion:* one figure and a paragraph
that survives a referee asking "so what is `t*`".

**Stage 3 — the Szegedy quantisation. DONE, and negatively.** With the lazy
chain `P_ab ∝ C(a→b)^θ` the discriminant is exactly a positive diagonal
congruence of the Gibbs kernel, `√(P_ab P_ba) = Δ_a Δ_b e^{−θd(a,b)/2}` to
`10⁻¹⁶`, so by Sylvester's law its inertia — in particular whether it is a Gram
matrix — is the Schoenberg question at `t = θ/2` and does not involve `A` at all.
Flipping `A ↦ −A` moves the discriminant spectrum by `3.5·10⁻⁵` and its inertia
by zero signs. The walk is unitary to `6.7·10⁻¹⁶` and the chain is *not*
reversible, so Szegedy is genuinely the right tool and it still cannot see the
comparison; `A` reaches it only through the row sums, i.e. as a potential.
Record: `szegedy_walk.py`. Original text:

**Stage 3 — the Szegedy quantisation (one to two sessions).** Choose a
stochastic matrix from `M` and defend the choice; build `W(P)`; compute the
spectrum and compare `arccos` of the discriminant eigenvalues against it; report
where `t* > ½` makes the discriminant non-Gram. *Success criterion:* an explicit
signature family with a computed quantum-walk spectrum, and a statement of what
the walk mixes to. *Risk:* the choice of `P` may be arbitrary enough that
nothing invariant comes out. Kill the stage rather than defending an arbitrary
normalisation.

**Stage 4 — dilation covariance as RH. DONE, and dropped.** Covariance is
unconditional; it is the functional equation. Two corrections propagated: the
"this step uses RH" callout in FINDINGS T1.5, and the bullet in §3b(ii) of the
two-positivities note claiming that "the scale-invariance of the Weil geometry
is a manifestation of the critical line". The impossibility argument of T1.5 is
unaffected — it needs only that the two invariance groups are transverse, and
both invariances are now unconditional, which if anything strengthens it. See
Part 1 above for what replaced it.

**Stage 5 — checkpoints 1 and 2 DONE; the shape of the answer is now known.**
The question was: is there a comparison of resources built from KMS states
rather than from partition functions, continuing below `β = 1`? The answer
splits three ways, and the split is the useful part.

*Known, no work needed.* Bost–Connes has a unique KMS state for `0 < β ≤ 1` and
extremal states parametrised by embeddings for `1 < β ≤ ∞`; the low-temperature
ones are type `I_∞`, the `β ≤ 1` one is type `III₁`. **So `β = 1` is a transition
in the *type* of the algebra.** In type III there is no trace — hence no density
matrix, no entropy, no partition function, and no thermomajorisation, so the
"second laws" of quantum thermodynamics have no formulation there at all. What
does survive is everything relative: Araki's relative entropy, the Connes
cocycle, and the Araki–Masuda `L^p` spaces, which *are* the sandwiched Rényi
divergences at `α = p/2` and satisfy data processing for `α ≥ ½` over arbitrary
von Neumann algebras (Berta–Scholz–Tomamichel; Jenčová). Separately, in
reversible resource theories the asymptotic conversion rate is known to be a
*ratio of regularised relative entropies* — the same shape as
`C(a→b) = inf_β log Z_a/log Z_b`.

> **The comparison does continue below `β = 1`, but only in relative form: the
> value of a resource is gone, the rate between two of them is not.** That is
> the project's own methodological choice — exchange rates primitive, no scalar
> value — arriving as a structure theorem instead of a taste.

*Derivable in a session — done, in `kms_comparison.py`.*

1. **The first paper's monotones are divergences against the trace.** Put `H_a`
   with eigenvalues `−log a_i`, so `Tr e^{−βH_a} = Z_a(β)`; then
   `D_α(ρ_a(β)‖tr) = (log Z_a(αβ) − α log Z_a(β))/(α−1)`, verified to `10⁻³⁰`,
   and the whole family `{Z_a(β)}` is recoverable from it. The trace is explicit,
   which is exactly why the construction cannot cross into type III.
2. **Between two primon-gas KMS states the Rényi divergence is exactly the
   convexity defect of `log ζ` along a chord:** with `L = log ζ` and
   `β_α = αβ₁ + (1−α)β₂`,
   `D_α(ω_{β₁}‖ω_{β₂}) = (L(β_α) − αL(β₁) − (1−α)L(β₂))/(α−1)`. Confirmed
   against direct summation to the truncation tail.
3. **The Hagedorn wall truncates the `α`-family.** `D_α` is finite exactly while
   `β_α > 1`, i.e. for `α < α* = (β₂−1)/(β₂−β₁)`; verified sharp to `±0.002` in
   `α` at three parameter pairs. So above the transition **only finitely many of
   the second laws are available**. `β = 1` is not merely a boundary in
   temperature — it cuts the family of monotones that governs conversion.

*Requires a theorem — the hard part, and it is one obstacle, not several.*
There is no asymptotic-interconversion (Strassen / Stein) theorem in type III,
so `D_α` there is a monotone without operational meaning. And the reason is
sharper than "nobody has done it": **the operational definition of the exchange
rate, `C(g→f) = lim max{k : f^k is implemented by g^n}/n`, is a counting
statement, and counting is a trace.** Type III has no trace, so the definition
does not merely fail to converge — it fails to parse. Any KMS comparison will be
a divergence first and a *rate* only if someone proves a coding theorem for it.
That is the same gap already recorded for `C_ζ` in §5B of the two-positivities
note, now located precisely.

*Literature caution.* Targeted searches found the ingredients (Araki–Masuda in
type III, ratio-of-relative-entropies conversion rates, the generalised quantum
Stein's lemma settled in finite dimensions by Hayashi–Yamasaki and Lami, 2024)
but nothing at the intersection of Bost–Connes and resource conversion. Absence
of search results is not absence of literature; check properly before claiming
novelty.

*Superseded planning text.*

**Stage 5 — first checkpoint DONE.** `primon_gas_hagedorn.py`: `β = 1` is the
**Hagedorn temperature** of the primon gas, not a coordinate artefact. The mean
energy is `U = −ζ'/ζ = Σ Λ(n)n^{−β}` (verified against direct summation to the
tail bound), and as `β → 1⁺`, `U·(β−1) → 1` and `S/U → 1`, so the entropy–energy
curve tends to the line `S = U` and the supporting line of slope 1 has no finite
contact point. That is exactly why `C(ξ→P) = C(P→ξ) = 0` unrestricted and why
the companion paper had to impose a temperature window. Every finite portfolio
is Hagedorn-free; only the infinite tensor product acquires the transition —
the same sentence as "the zeros of `ζ` are a phenomenon of the infinite tensor
product", now with a temperature attached. The KMS classification below `β = 1`
is the remaining programme. Original text:

**Stage 5 — the KMS route (a programme, not a session).** Below `β = 1` the
monotone diverges and `C(a→b)` is undefined. Ask what replaces it. Concretely:
does the Bost–Connes KMS classification give a comparison of resources on the
critical strip, and does the idele class group action have a resource-theoretic
reading? *First checkpoint, achievable in one session:* compute the KMS states
of the primon gas at `β > 1`, verify that the Gibbs state reproduces
`Z(β) = ζ(β)` and that the exchange rate against a truncated Euler factor
matches `analysis/xi_versus_euler_factors.py`, and identify what exactly breaks
at `β = 1` in the resource language. Only then decide whether to go further.

---

## What will not work, and should be said up front

* **`MᵀM`, `MᵀDM` and every cousin.** PSD for free, informative never.
* **A Hermitian quantisation of `M`.** Forbidden by irreversibility; `d` is the
  obstruction; this is a theorem, not a gap.
* **Any implication toward RH.** The two obstacles of FINDINGS stand unchanged:
  atomic measures are not admissible Weil test functions (so `E` is a
  finite-rank truncation), and the monotone diverges in the critical strip. The
  quantum language renames both obstacles precisely and removes neither.
* **Novelty of the QM picture itself.** Hilbert–Pólya, `xp`, Connes and
  Bender–Brody–Müller are all prior art. Claim only the two new things: the
  gauge decomposition of `M`, and the `t = ½` Szegedy point.

## Traps

* `−½JDJ` has the constant vector in its kernel; work in an orthonormal basis of
  `{Σx = 0}` (FINDINGS T1.1).
* Any windowed `β` computation on these objects needs `β ~ 10³`; grids truncated
  below `β ≈ 500` hide the phenomena (FINDINGS T1.2).
* `exchange_rate` is accurate to `~1e−13`; treat differences below `1e−10` as
  ties. The `S`/`A` identity is exact by algebra, so a `β`-grid check of it
  measures the grid, not the identity.
* Two unrelated one-parameter families are both called `t` in the existing
  notes. Rename one before writing.

## Reproduce

`research/m_and_e_and_a_c/gauge_decomposition.py` — all of Part 0, the `t = ½`
table, and the `t*` correction.
