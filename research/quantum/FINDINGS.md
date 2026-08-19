# Findings — the quantum exchange rate

Answer to session brief L. Every claim below is marked *proved* or *computed*;
the headline witness is certified at 40 digits (`q_certify.py`,
`q_certify.json`).

## Summary

The brief's core observation — `Z_a(β) = Tr A^β`, so the definition never used
commutativity — is correct, and the quantisation it produces is **empty**. That
is settled in one line and is the first section. The repair the brief proposes,
the sandwiched Rényi divergence against a background operator, is **not** empty:
it is a genuinely non-spectral functional, the whole structure theorem of briefs
G and I survives it with the same constants, and it produces a **certified
3-cycle among three qubit states that decoherence destroys**. But it produces no
new *geometry*: every sandwiched profile lies in the classical tropical cone, so
every quantum configuration is realised by ordinary integer signatures to
distortion `1 + O(1/log r)`. And against Brandão–Horodecki–Ng–Oppenheim–Wehner
the family is theirs; only the functional taken of it is not.

* **The spectral quantisation is empty** *(proved, one line)*. `Tr A^β` is a
  function of `spec A`, so `C(A→B) = C(spec A → spec B)` identically; `d`, the
  midrange `A`, the tournament and every curl are the pullback along
  `A ↦ spec A`, which is *onto* the positive `r`-tuples. **There is no quantum
  3-cycle without a classical shadow and there cannot be one** — the decisive
  negative the brief asked to be reported early. Verified to `2.0·10⁻¹⁵` on 300
  pairs and 2000 non-commuting triples (`q_spectral.py`).
* **The sandwiched exchange rate is the non-empty definition** *(proved)*. A
  *quantum signature* is a pair `(A,S)` of positive operators with `A ≽ S ≽ I`
  — the operator form of "multiplicities ≥ 1, fiber sizes ≥ 1" — with profile
  `F_{(A,S)}(β) = log Tr[(S^{(1−β)/2β} A S^{(1−β)/2β})^β]`. `S = I` gives back
  `log Tr A^β`; `[A,S] = 0` gives back a classical cone point exactly; and for
  `[A,S] ≠ 0` it does not factor through the two spectra (explicit isospectral
  witnesses with `d` up to `0.098`).
* **The whole structure theorem survives, with the same constants** *(sandwich
  proved; the rest computed)*. `max(R, βΛ) ≤ F ≤ R + βΛ` is **proved** for
  admissible pairs by three lines of operator monotonicity, with
  `R = log Tr S`, `Λ = log λ_max(S^{−1/2}AS^{−1/2})`. That is exactly the
  hypothesis briefs G and I say their proofs use, so FINDINGS Theorem A(1) and
  Corollary A2 transfer verbatim: `d = |Δσ| + ε`, `ε ≤ 2log(1+e^{−Δ})`,
  `|D| ≤ ½log(1+e^{−Δ}) ≤ (log 2)/2`. Computed on 600 random admissible pairs:
  `F` increasing and convex to the noise floor, `F − βΛ` nonincreasing to
  `2.7·10⁻¹²`, `U' ∈ (0,1)`, `w` unimodal with peak at `σ` and height
  `≤ 0.692101 < log 2`.
* **The mechanism is a shear about `β = 1`** *(proved)*. With `P_S` the pinching
  onto `S`'s eigenspaces (decoherence in the energy basis),
  `F_{(A,S)} − F_{(P_S A,S)}` vanishes at `β = 0` and `β = 1` and has the sign
  of `β − 1` on `[½,∞)`, by data processing for the sandwiched divergence.
  Coherence therefore tilts each profile about `β = 1`, and the tilt is what can
  flip an edge of the tournament. Verified with **zero** violations on 300 pairs.
* **A 3-cycle created by coherence alone** *(computed; certified at 40 digits)*.
  Three qubit states in one background `S = diag(2.9015…, 8.3891…)` have
  `A(1,2), A(2,3), A(3,1) = +0.0035114…, +0.0035135…, +0.0035093…` — a strict
  3-cycle with `|curl A|/Σ|A| = 1` exactly, brief D's criterion — while their
  decohered shadows give `+0.030054…, +0.025608…, −0.055205…`, a **transitive**
  tournament. Margins `3.5·10⁻³` and `3.0·10⁻²`, i.e. `3.5·10⁷` times the brief's
  `10⁻¹⁰` tie threshold. The separation survives taking the infimum over the
  full `β ∈ [0,∞]` instead of the data-processing range `[½,∞]`. A `3×3` witness
  with margin `6.6·10⁻³` is also recorded, and a longer search reaches
  `8.4·10⁻³` at `r = 2`.
* **But the quantisation adds no new geometry** *(computed)*. Every admissible
  quantum profile satisfies the cone conditions of OBSTRUCTION Theorem 1, so it
  is a classical cone point: its tangent-line envelope reproduces the quantum
  `d` to `1.2·10⁻⁸`, and pushing that envelope down to genuine **integer
  signatures** reproduces `d` at the predicted `O(1/K)`, measured constant
  `0.146`. **So the quantum cycle does have a classical shadow — just not the
  decohered one.** The honest statement of what quantisation buys is:
  non-spectral as a *functional on operators*, identical as a *realisable
  geometry*.
* **Placement against BHNOW: the family is theirs** *(proved by identity)*. With
  `γ = S/Tr S` and `Z = Tr S`, `F_{(A,S)}(β) = (β−1)·F̂_β(A,γ)/kT` where `F̂_α` is
  BHNOW's eq. (5), their **sandwiched** α-free energy. The `(β−1)` cancels in
  the ratio, so `C̃ = inf_β F̂_β(A)/F̂_β(B)` **is** an infimum of a ratio of their
  monotones. Their own scalar is the *work distance* `inf_α[F_α(ρ) − F_α(ρ')]`,
  an infimum of a **difference**. Ratio vs difference is the whole difference:
  the difference functional induces a partial order (acyclic by construction),
  the ratio functional a tournament (which cycles). Claim novelty only for the
  ratio-and-oscillation *geometry*, never for the family.

---

## Notation

That of `research/realizability/FINDINGS.md` and `OBSTRUCTION.md`.
`F_a = log Z_a`, `U_a(s) = log F_a(e^s)`, `φ = U_b − U_a`, `d = osc φ`,
`A = mid φ`, `a ≺ b ⟺ A(a,b) > 0`, `R = F(0)`, `Λ = lim F'`,
`σ = log(R/Λ)`, `w = U − log Λ − max(σ, s)`. `C(a→b) = inf_β F_a/F_b`
(`a` the implementer). Operators are real symmetric positive definite;
`P_S` is the pinching onto the eigenspaces of `S`. `Q̃_β(A‖S)` is the sandwiched
Rényi quantity, `D̃_β` the sandwiched divergence, so `log Q̃_β = (β−1)D̃_β`.

---

## 1. The naive quantisation is empty *(proved)*

**Proposition 1.** For positive operators `A, B` on a finite-dimensional Hilbert
space, `Tr A^β = Σ_k λ_k(A)^β`, hence

```
C(A→B) = C(spec A → spec B),   d(A,B) = d(spec A, spec B),   A(A,B) = A(spec A, spec B)
```

with the right-hand sides the classical quantities of the two positive
`r`-tuples. ∎

There is nothing more to it; the proof is the display. Three consequences, all
immediate, all worth stating because they close the brief's items 1 and 2.

**Corollary 1.1 (conjugation invariance).** `C(UAU*→B) = C(A→B)` for every
unitary `U`. So no arrangement of non-commuting operators is distinguishable
from a commuting one. *Verified* to `2.2·10⁻¹⁶` over 50 random rotations.

**Corollary 1.2 (no quantum cycle without a classical shadow).** The tournament,
the metric, the midrange and every curl of a family `A_1,…,A_n` of positive
operators equal those of the family `diag(spec A_1),…,diag(spec A_n)`. A quantum
cycle with no classical shadow does not exist. *Verified*: over 2000 random
non-commuting triples, `0` tournaments differ from their spectral shadow.

**Corollary 1.3 (nothing is even gained in range).** `A ↦ spec A` is onto the
positive `r`-tuples, so the achievable set of quantum profiles equals the
achievable set of classical *real-atom* signatures. That relaxation from integer
to real atoms is real, but it is not quantum: it is already inside OBSTRUCTION
Theorem 1's cone `C`, whose projective closure is the classical achievable set.

*Verified* (`q_spectral.py`, `q_spectral_output.txt`): on 300 random pairs with
`r = 2..7`, `max|d_quantum − d_classical| = 2.0·10⁻¹⁵` and
`max|A_quantum − A_classical| = 8.9·10⁻¹⁶`, over commutator norms up to `752`.
These are two ways of writing the same sum; the residual is LAPACK's eigenvalue
error.

**The only non-vacuous extension of Proposition 1 is not about
non-commutativity.** In infinite dimensions `Tr A^β` can diverge below some
`β₀` and converge above it — the primon-gas Hagedorn transition of brief D
Stage 5 and `research/kms_critical_strip/`. That is a genuine new regime, and it
is still a statement about a spectrum.

---

## 2. The sandwiched exchange rate *(definition; two proved reductions)*

**Definition.** A **quantum signature** is a pair `(A,S)` of positive definite
operators on the same space with

```
A ≽ S ≽ I                                              (admissibility)
```

— the operator form of FINDINGS §1.1's standing hypothesis. Its profile is

```
F_{(A,S)}(β) = log Q̃_β(A‖S) = log Tr[ (S^{(1−β)/2β} A S^{(1−β)/2β})^β ],
R = F(0) = log Tr S,        Λ = log λ_max(S^{−1/2} A S^{−1/2}),
C̃( (A,S) → (B,T) ) = inf_β  F_{(A,S)}(β) / F_{(B,T)}(β).
```

`S` is the background (a Hamiltonian `H = −log S`, an unnormalised thermal
state); it is part of the resource, exactly as `(ρ,H)` is in BHNOW.

**Proposition 2 (the two reductions).** *(proved)*

1. `S = I` gives `Q̃_β(A‖I) = Tr A^β`, so §1 is the case of a flat background.
2. If `[A,S] = 0` then, in a joint eigenbasis with `S = diag(s_i)`,
   `A = diag(a_i)`,

   ```
   Q̃_β(A‖S) = Σ_i a_i^β s_i^{1−β} = Σ_i s_i · (a_i/s_i)^β ,
   ```

   the classical profile of a cone point with **multiplicities** `m_i = s_i` and
   **atoms** `x_i = a_i/s_i`. Admissibility is exactly `m_i ≥ 1`, `x_i ≥ 1`,
   i.e. membership in the cone `C` of OBSTRUCTION Theorem 1. ∎

*Verified* (`q_sandwich.py` Q5): `max|F_pinched − classical log-sum-exp| =
2.0·10⁻¹⁰` over 200 draws on the widest available `β` grid.

**The domain of the infimum shrinks.** `D̃_α` is a monotone under quantum
channels only for `α ≥ ½` (BHNOW's quantum second laws are stated for `α ≥ ½`
for the sandwiched family, `0 ≤ α ≤ 2` for the Petz family). So the
operationally defensible object is

```
C̃_½ = inf_{β ∈ [½,∞]} F_a/F_b ,
```

and the classical `β = 0` endpoint `R = log r` — one of the two endpoints of the
Hilbert metric — is no longer available. **This is a real structural change that
the classical theory has no analogue of, and it is a consequence of data
processing, not of a numerical convenience.** Everything below is computed on
`[½,∞]`; §5 reports that the headline separation also survives on `[0,∞]`.

**Numerics.** The literal expression overflows as `β → 0`, since the exponent
`(1−β)/2β` diverges. The similarity `S^{−p}(S^p A S^p)S^p = A S^{2p}` turns the
`μ_k(β)` into the generalised eigenvalues of the pencil `(A, S^t)` with
`t = (β−1)/β ∈ [−1,1)` for `β ∈ [½,∞)` — perfectly conditioned on exactly the
data-processing range. Below `β = ½` we rescale `S` by its geometric mean
(`Q̃_β(A‖cS) = c^{1−β}Q̃_β(A‖S)`, an exact identity) and fall back to mpmath.

---

## 3. The structure theorem survives *(sandwich proved; the rest computed)*

### 3.1 The sandwich *(proved)*

**Theorem Q1.** For an admissible `(A,S)` and every `β > 0`,

```
max( R , βΛ )  ≤  F_{(A,S)}(β)  ≤  R + βΛ .
```

*Proof.* Three facts about `Q̃_β`, each one line.

*(i) Homogeneity.* `S^p(cA)S^p = c·S^pAS^p`, so `Q̃_β(cA‖S) = c^β Q̃_β(A‖S)`.

*(ii) Monotonicity in the first argument.* `A ≼ A'` gives
`S^pAS^p ≼ S^pA'S^p`, hence `λ_k(S^pAS^p) ≤ λ_k(S^pA'S^p)` for every `k` by
Weyl, hence `Tr X^β` is monotone for `β > 0`.

*(iii) The reference value.*
`Q̃_β(S‖S) = Tr[(S^{(1−β)/β+1})^β] = Tr[(S^{1/β})^β] = Tr S`.

*Upper.* `S^{−1/2}AS^{−1/2} ≼ e^Λ I` gives `A ≼ e^Λ S`, so by (ii), (i), (iii),
`Q̃_β(A‖S) ≤ e^{Λβ} Tr S`, i.e. `F ≤ βΛ + R`.

*Lower `R`.* `A ≽ S` gives `Q̃_β(A‖S) ≥ Q̃_β(S‖S) = Tr S`, i.e. `F ≥ R`.

*Lower `βΛ`.* The `μ_k(β)` are the eigenvalues of the pencil `(A, S^t)`,
`t = (β−1)/β`. Since `S ≽ I`, `t ↦ S^t` is nondecreasing in the PSD order, so
`μ_max(β) = max_x (x*Ax)/(x*S^t x)` is nonincreasing in `t`, hence in `β`, with
limit `λ_max(A,S) = e^Λ` as `β → ∞`. Therefore `μ_max(β) ≥ e^Λ` and
`F(β) = log Σ_k μ_k(β)^β ≥ β log μ_max(β) ≥ βΛ`. ∎

**Corollary Q1.1.** Because the proofs of FINDINGS Theorem A(1) and Corollary A2
use *only* this sandwich, they transfer verbatim to the sandwiched quantum
profile:

```
0 ≤ w ≤ log(1 + e^{−|s−σ|}),   d = |Δσ| + ε,   ε = P + Q ≤ 2 log(1 + e^{−Δ}),
|D| ≤ ½ log(1 + e^{−Δ}) ≤ (log 2)/2 = 0.34657… .
```

**The quantum defect obeys the same sharp constants as the classical one.**
This is the statement for the full domain `β ∈ [0,∞]`, which is where
`φ(−∞) = log R_b − log R_a` and `φ(+∞) = log Λ_b − log Λ_a` are the endpoints;
on the data-processing domain `[½,∞]` the upper bounds still hold, being
suprema over a subset, but the decomposition `d = |Δσ| + ε` does not, because
the `β = 0` endpoint is no longer available.

### 3.2 What is computed and not proved

The sharp halved constant of OBSTRUCTION Theorem 4 (`ε ≤ log(1+e^{−Δ})`) also
needs `F' ≤ Λ`, and the cone membership of §6 needs convexity, and its
tangent-line construction needs `F' ≤ Λ` again. Both are computed, not proved. On 600 random admissible pairs, `r = 2..6`, each on its
own widest double-precision `β` grid (3001 log-spaced points from `β_safe` to
`3000`):

| claim | worst violation |
|---|---:|
| `F ≥ R` | `−6.2·10⁻⁴` |
| `F ≥ βΛ` | `−2.8·10⁻²` |
| `F ≤ R + βΛ` | `−2.2·10⁻⁵` |
| `F` nondecreasing | `−3.1·10⁻⁶` |
| `F − βΛ` nonincreasing (i.e. `F' ≤ Λ`) | `−2.7·10⁻¹²` |
| `F` convex in `β` | `+1.4·10⁻⁶` |
| `U' ≥ 0` | `−2.4·10⁻⁴` |
| `U' ≤ 1` | `−3.0·10⁻⁶` |
| `U` convex in `s` | `−1.2·10⁻⁸` |
| `w ≥ 0` | `−3.0·10⁻⁶` |
| `w ≤ log(1+e^{−|s−σ|})` | `−8.2·10⁻⁶` |

(negative = the claim holds). Max bump height `0.692101` against
`log 2 = 0.693147`; `max|argmax w − σ| = 2.4·10⁻³`, below the grid step.

The one non-negative entry, `F` convex, is **noise**: the identical
second-difference statistic run on the *pinched* profiles, where convexity is a
theorem (Proposition 2.2 makes them log-sum-exps of affine functions), gives
`+1.2·10⁻⁶` — the same number. Convexity of the quantum `F` is therefore
*computed to the noise floor of the method*, and a proof is Open.

Coherence sampled over the ensemble: median `0.299`, max `0.670` (relative
Frobenius weight of the off-diagonal part of `A` in `S`'s eigenbasis).

### 3.3 The profile is genuinely non-spectral *(computed; explicit witnesses)*

Fix `S` and replace `A` by `OAOᵀ` for a random `O ∈ O(4)`. Both spectra are
unchanged; the profile is not:

| `max_s \|ΔU\|` | `d(base, rotated)` | `C̃(base→rot)` | `C̃(rot→base)` |
|---:|---:|---:|---:|
| `1.08·10⁻³` | `1.88·10⁻³` | `0.9989186213` | `0.9991993127` |
| `8.88·10⁻²` | `9.79·10⁻²` | `0.9149992133` | `0.9909628495` |
| `8.05·10⁻²` | `8.94·10⁻²` | `0.9226357021` | `0.9911482905` |
| `3.45·10⁻²` | `3.66·10⁻²` | `0.9660478450` | `0.9979061334` |
| `6.71·10⁻²` | `7.61·10⁻²` | `0.9350877960` | `0.9910655201` |
| `9.07·10⁻³` | `1.06·10⁻²` | `0.9985033461` | `0.9909618791` |

Every row is a pair of isospectral quantum signatures at nonzero exchange
distance. By Proposition 1 the spectral rate of §1 assigns `d = 0` to all of
them. **`C̃` does not factor through the spectra**, which is what the brief asked
for and is the whole reason §4–§6 exist.

---

## 4. The mechanism: coherence shears the profile about `β = 1` *(proved)*

**Theorem Q3.** Let `P_S` be the pinching onto the eigenspaces of `S`. Then for
an admissible `(A,S)`:

1. `F_{(A,S)}(0) = F_{(P_S A,S)}(0) = log Tr S`;
2. `F_{(A,S)}(1) = F_{(P_S A,S)}(1) = log Tr A`;
3. `sign( F_{(A,S)}(β) − F_{(P_S A,S)}(β) ) = sign(β − 1)` for `β ≥ ½`;
4. `Λ_{(A,S)} ≥ Λ_{(P_S A,S)}`.

*Proof.* (1) Both sides are `log Tr S` by Theorem Q1's (iii) argument as
`β → 0`. (2) At `β = 1` the sandwich exponent is `0`, so `Q̃_1(A‖S) = Tr A`, and
pinching preserves the trace. (3) `P_S` is a CPTP map with `P_S(S) = S`, so data
processing for the sandwiched divergence on `[½,∞]` gives
`D̃_β(P_S A‖S) = D̃_β(P_S A‖P_S S) ≤ D̃_β(A‖S)`; multiply by `(β−1)`, whose sign
flips at `β = 1`. (4) In `S`'s eigenbasis `λ_max(S^{−1/2}AS^{−1/2})` dominates
the diagonal entries `A_ii/s_i`, which are the eigenvalues of
`S^{−1/2}P_S(A)S^{−1/2}`. ∎

*Verified* (`q_sandwich.py` Q3) on 300 random admissible pairs: the worst
positive value of `F − F_pinched` on `[½,1)` is **exactly `0`**, the worst
negative value on `(1,∞)` is **exactly `0`**, and
`max|F(1) − F_pinched(1)| = 2.7·10⁻¹⁵`.

> Coherence acts on a quantum signature as a **shear of its profile about
> `β = 1`, pinned at `β = 0` and `β = 1`.** In a triple each node gets its own
> shear, and three independent shears are enough freedom to move a tournament.

---

## 5. A 3-cycle created by coherence alone *(computed; certified at 40 digits)*

### 5.1 The witness

Two qubits' worth of Hilbert space, one common background. All numbers below are
`q_certify.py` at 60 working digits, reported at 40.

```
S   = diag( 2.9015264522623596 , 8.38905609893065 )

A_1 = [[ 35.91614260201353 ,  28.889389303494585],
       [ 28.889389303494585,  61.075159852945504]]
A_2 = [[ 39.170607105541855, -16.214024560871955],
       [-16.214024560871955,  72.60949333751388 ]]
A_3 = [[  5.821732544969596,  12.856146279129781],
       [ 12.856146279129781, 113.25225733556378 ]]
```

Admissibility, at 40 digits: `min spec(S) − 1 = 1.9015…`, and the spectra of
`S^{−1/2}A_iS^{−1/2}` are `[3.4430…, 16.2157…]`, `[6.9949…, 15.1603…]`,
`[1.4433…, 14.0632…]` — every one bounded below by `1`, with the tightest margin
`0.4433`.

**Quantum midranges** (`β ∈ [½,10⁶] ∪ {∞}`), 40 digits:

```
A(1,2) = 0.003511424424370079214602686839085134010108
A(2,3) = 0.003513481085729650524813865820852736725468
A(3,1) = 0.003509281100295324562510779292328990370496
d(1,2) = 0.0607310831090470955895640647321708752441
d(2,3) = 0.06306753856513011382151901814125953526124
d(3,1) = 0.09792600378175595968222957945466435620793
```

All three midranges positive: `1 ≺ 2 ≺ 3 ≺ 1`, a **strict 3-cycle**, with
`|curl A|/Σ|A| = 1.0` exactly — brief D Part 0(b)'s criterion — and
`curl A = 0.010534186610…`. No-arbitrage holds:
`Σ S = 0.11086 ≥ |Σ A| = 0.01053`.

**Decohered shadow** `(P_S A_i, S)`, the same 40 digits:

```
A(1,2) = +0.03005373440485496379140519499985175926445
A(2,3) = +0.02560750584265668921820237521438404458845
A(3,1) = −0.05520470591844018514792241972630329451117
```

Signs `+ + −`: the tournament is `3 ≺ 1 ≺ 2`, **transitive**, with
`|curl A|/Σ|A| = 0.00412`.

```
quantum cycle margin  = +3.509281100295·10⁻³    (>0: strict 3-cycle)
shadow  cycle margin  = −3.005373440485·10⁻²    (<0: transitive)
separation            = +3.509281100295·10⁻³
                      = 3.5·10⁷ times the brief's 10⁻¹⁰ tie threshold
```

**Independent cross-check.** The double-precision path (`q_cycle.py`) computes
`μ_k` from the LAPACK Cholesky of the pencil `(A, S^t)`; the certification path
(`q_certify.py`) solves the `2×2` characteristic quadratic in mpmath. The two
share no code and agree on all three midranges to `10⁻⁹`, the double path's own
precision.

**Robustness to the domain.** Recomputing with the infimum over the *full*
classical domain `β ∈ [10⁻⁶,10⁶] ∪ {0,∞}` — i.e. dropping the `α ≥ ½`
data-processing restriction — moves the individual midranges but not the
verdicts: quantum `0.05225, 0.003513, 0.003509` (still a cycle, margin
`3.509·10⁻³`), shadow `+0.06567, −0.03727, +0.01762` (still transitive, margin
`−3.727·10⁻²`). **The cycle is not an artefact of the restricted range.**

A `3×3` witness is recorded as well (`q_cycle.json`, key `best_r3`), with
`S = diag(1.0706, 1.0498, 1.0624)`, quantum midranges
`+0.0065990, +0.0065974, +0.0065974` and shadow `+0.03915, −0.07299, +0.03384`
— margin `6.6·10⁻³` against shadow acyclicity `3.9·10⁻²`.

**The margin is not near a ceiling.** A longer search (10 restarts instead of 3,
`q_cycle_full_output.txt` — a partial run: the `r = 2` half is complete, the
`r = 3` half was still running when this file was written) reaches
**`8.408·10⁻³`** at `r = 2`, with `S = diag(1.1048, 2.9437)`, quantum midranges
`−0.0084078, −0.0084079, −0.0084076` (the cycle runs the other way) and shadow
`−0.010189, −0.018193, +0.009752`. Notably its coherences are *lower* —
`0.013, 0.149, 0.236` against the certified witness's `0.500, 0.268, 0.158` — so
the size of the effect is not monotone in the coherence, and 4 of the 10
restarts found no separation at all. The certified numbers above are those of
the `--quick` witness, which is the one `q_cycle.json` holds and `q_certify.py`
reads by default; `q_certify.py q_cycle_full.json best_r2` certifies the other.

### 5.2 What this is and is not

It **is** the statement that decoherence in the energy eigenbasis can turn a
strict preference cycle into a transitive order — coherence as the resource that
carries the intransitivity. It is *not* a violation of any second law: every one
of the six rates on the cycle is `< 1` (`0.9489…` to `0.9735…`), so no edge is a
feasible single-shot conversion, and BHNOW's relation "all free energies
decrease" is a partial order, acyclic by construction. The cycle lives in the
**rate tournament**, a different object (§7).

It is also **not** a cycle with no classical shadow. That is §6.

---

## 6. The sandwiched profile does not leave the classical cone *(computed)*

OBSTRUCTION Theorem 1: `C = {Φ convex, nondecreasing, Φ ≥ Λ_Φ β}` is exactly the
projective closure of the classical achievable set, and every `Φ ∈ C` is a
uniform limit of `(1/K)F_{a^{(K)}}` for genuine integer signatures.

Theorem Q1 proves `F ≥ βΛ`; §3.2 computes convexity and monotonicity to the
noise floor. So the quantum profile is in `C`. Made operational (`q_cone.py`),
with **every** metric an oscillation over the same domain `β ∈ [½,∞]`:

* **Cone membership.** The tangent-line envelope `Φ = max_j(c_j + x_jβ)` of a
  quantum profile has `min_j c_j = +7.07·10⁻²` over 200 random admissible pairs
  — `c_j ≥ 0` is exactly membership in `C`, and it follows from `F ≥ βΛ ≥ βF'`.
  With 60 tangents, `sup|Φ − F| ≤ 3.6·10⁻³`.
* **The metric agrees.** With 400 tangents, over 10 random quantum pairs,
  `|d_quantum − d_cone|` ranges from `9.4·10⁻¹⁰` to `2.9·10⁻⁶`.
* **Integer signatures.** Taking `⌈e^{Kc_j}⌉` copies of the atom `e^{Kx_j}`:

  | `K` | `log₁₀ r` | `d_signature` | `\|d_sig − d_cone\|` | `K·gap` |
  |---:|---:|---:|---:|---:|
  | 25 | 28.44 | 0.2689857771 | `3.79·10⁻³` | 0.095 |
  | 50 | 55.30 | 0.2673799549 | `2.18·10⁻³` | 0.109 |
  | 100 | 109.24 | 0.2664098586 | `1.21·10⁻³` | 0.121 |
  | 200 | 217.56 | 0.2658538460 | `6.53·10⁻⁴` | 0.131 |
  | 400 | 434.89 | 0.2655459838 | `3.45·10⁻⁴` | 0.138 |
  | 800 | 869.77 | 0.2653798422 | `1.79·10⁻⁴` | 0.143 |
  | 1600 | 1739.54 | 0.2652920392 | `9.13·10⁻⁵` | 0.146 |
  | 3200 | 3479.08 | 0.2652464801 | `4.57·10⁻⁵` | 0.146 |

  against `d_quantum = 0.2652007621`, `d_cone = 0.2652007504`
  (`|diff| = 1.2·10⁻⁸`). The `O(1/K)` rate of OBSTRUCTION Theorem 1(2) is
  reproduced with constant `0.146`.

> **Conclusion.** The sandwiched quantum exchange geometry embeds in the
> classical cone geometry. Every quantum configuration — the 3-cycle of §5
> included — is realised by ordinary integer signatures to distortion
> `1 + O(1/log r)`. The quantisation is **non-spectral as a functional on
> operators and empty as a source of new exchange geometry.**

This is the correct reading of the brief's success criterion. "A quantum cycle
with no classical shadow" fails at §1 for a trivial reason and at §6 for a
structural one; what survives, and is certified, is the strictly weaker and
physically meaningful "a quantum cycle with no *decohered* shadow".

---

## 7. Placement against BHNOW *(proved by identity; computed)*

Brandão–Horodecki–Ng–Oppenheim–Wehner, *The second laws of quantum
thermodynamics*, arXiv:1305.5278, PNAS **112** (2015) 3275. Their eq. (2)
defines `F_α(ρ,ρ_β) = kT D_α(ρ‖ρ_β) − kT log Z`; for states block-diagonal in
the energy basis, `ρ → ρ'` under catalytic thermal operations **iff**
`F_α(ρ) ≥ F_α(ρ')` for all `α ≥ 0`. For general states their eq. (4) uses the
Petz divergence (monotone `0 ≤ α ≤ 2`) and their eq. (5) the sandwiched one,

```
F̂_α(ρ,ρ_β) = kT/(α−1) · log Tr[ (ρ_β^{(1−α)/2α} ρ ρ_β^{(1−α)/2α})^α ] − kT log Z ,
```

monotone for `α ≥ ½`, and these quantum laws are **necessary but not
sufficient**.

**B1 — the classical exchange monotone is their family** *(identity; verified to
`4.8·10⁻¹⁶` relative)*. With `p_i = a_i/N`, `N = Σa_i`,
`log Z_a(β) = β log N + (1−β)H_β(p)` — the relation already flagged as unused in
`paper_exchange_rate/REVIEW.md` §5. So `log Z_a(β)` is an α-free energy of the
resource at trivial Hamiltonian, on an **unnormalised** state. Unnormalised is
what makes it strictly positive, and positivity is what makes the *ratio* of two
of them a rate. The index `β` is simultaneously the inverse temperature and the
Rényi order; for a flat background those coincide.

**B2 — the sandwiched exchange rate is literally theirs** *(identity; verified
to `5.5·10⁻¹⁶` relative)*. With `γ = S/Tr S` and `Z = Tr S`,

```
F_{(A,S)}(β) = (β−1) · F̂_β(A,γ) / kT ,      hence
C̃( (A,S)→(B,S) ) = inf_β  F̂_β(A,γ) / F̂_β(B,γ) ,
```

**the `(β−1)` cancelling identically in the ratio** — which is also why the
exchange rate is finite and continuous at `β = 1`, where the `(α−1)⁻¹` in the
definition of `F̂_α` has to be resolved by a limit. Because their quantum
laws are necessary only, `C̃` is an **upper bound on the achievable conversion
rate**, not the rate. (Classically, where the laws are necessary *and*
sufficient, `C` *is* the rate — which is the paper's own theorem, now with a
second proof.)

**B3 — ratio versus difference is the whole difference.** BHNOW's own scalar
functional of the family is the **work distance**
`D_work(ρ≻ρ') = kT inf_α [F_α(ρ) − F_α(ρ')]`, their eq. (7): an infimum of a
*difference*. The exchange rate is an infimum of a *ratio*. Differences are the
additive (single-shot, extractable-work) reading, ratios the multiplicative
(asymptotic, copies-per-copy) one. Two consequences:

* The difference functional induces a partial **order** — `ρ → ρ'` iff every
  `F_α` drops — which is transitive and cannot cycle. The ratio functional
  induces a **tournament** `a ≺ b ⟺ A(a,b) > 0`, which does cycle. **The
  3-cycles of this project are an artefact of forcing a total comparison on a
  partial order, and they contradict nothing.** On the §5 witness all six rates
  are `< 1`, so no edge is feasible in either direction.
* For unnormalised resources the difference functional is degenerate:
  `F_a(β) − F_b(β) ∼ β(Λ_a − Λ_b) → −∞` whenever `Λ_a < Λ_b`, so the work
  distance is `−∞` on most pairs (measured slopes `+0.0673, +0.0751, −0.1424`
  on the §5 witness). **The ratio is the only finite scalar available here**,
  which is a reason for the exchange rate rather than a preference.

**B4 — the exchange geometry is projective in the monotones** *(computed)*.
Replacing `F(β)` by `c(β)F(β)` for any positive `c` leaves `d` and `C` unchanged
(`0.042583151772` before and after) and changes the work distance
(`+0.2946 → +0.2034`). So the exchange geometry is a functional of the family of
monotones only **up to a positive ray at each `β`** — which is precisely
OBSTRUCTION §1's Hilbert projective metric, now with a resource-theoretic
reading: *the rate cannot see how you normalise your monotones.*

**What may be claimed as new, and what may not.**

* **May not:** the monotone family (BHNOW, and before them Horodecki–Oppenheim
  and the trumping literature); the sandwiched divergence and its `α ≥ ½` data
  processing; the observation that a family of additive monotones bounds an
  asymptotic conversion rate by an infimum of ratios — that is the standard
  resource-theory construction.
* **May:** the *geometry* of that infimum — the Hilbert projective metric on the
  cone `C`, the range/midrange (metric/1-form) split of brief D, the sharp
  constants `log 2` and `(log 2)/2`, the tournament and its cycles, and the
  arithmetic realisation of the family by fiber signatures of maps. None of
  those is in BHNOW, which never forms a ratio.

---

## Corrections

1. **Brief L's own framing of item 1 is confirmed and should be strengthened.**
   The brief asks to "say clearly whether the quantum object is therefore only a
   function of the two spectra". It is, and more: it is the *same function*, not
   merely determined by the spectra — `C(A→B)` is literally
   `C(spec A → spec B)` with no reparametrisation, so §1 is not a degeneracy to
   be measured but an identity (Proposition 1).
2. **The brief's success criterion is missing the case that actually holds.** It
   offers "either a quantum cycle with no classical shadow … or a proof that `C`
   factors through the spectra and the sandwiched version is the only
   non-trivial quantisation". The second horn is proved, but the sandwiched
   version is *also* free of cycles with no classical shadow (§6), for a
   completely different reason: it does not leave the cone. Non-spectral as a
   functional does not imply new geometry, and here it does not deliver it.
3. **`paper_exchange_rate/REVIEW.md` §5 understates its own suggestion.**
   "Rényi in one display … turns the Jensen analogy into a change of variables"
   is right but small: the same display identifies the exchange monotone with a
   named family of thermodynamic monotones whose completeness theorem is BHNOW's
   main result. The one-line change of variables is the citation hook for the
   whole second-laws literature, and belongs in the paper for that reason rather
   than for Jensen.
4. **No claim in this repo was found to be wrong by this session.** Brief D
   Part 2's warning about `MᵀM` and `MᵀDM` was respected and nothing here is
   built on them; brief D's Parts 1 and 3 are untouched. One addition to brief
   D's Traps: the `t = ½` Szegedy point and the `α = ½` data-processing
   threshold of §2 are unrelated coincidences of the numeral `½`, and are a
   fourth pair of one-parameter families in these notes that share a symbol and
   nothing else.

---

## Open

* **Prove convexity of `F_{(A,S)}` in `β`.** It is the one hypothesis of §6's
  cone-membership argument that is computed rather than proved, and everything
  in §6 rests on it. In the commuting case it is the log-convexity of
  `Σ m_i e^{βx_i}`; the quantum statement is that
  `β ↦ log Tr[(S^{(1−β)/2β}AS^{(1−β)/2β})^β]` is convex. A proof would upgrade
  §6 from computed to proved and would be worth having on its own as a fact
  about the sandwiched Rényi family.
* **Prove `F' ≤ Λ`**, equivalently that `β ↦ log Q̃_β(Ã‖S)` is nonincreasing when
  `Ã ≼ S`. Verified to `2.7·10⁻¹²`, the tightest of all the §3.2 checks, so it
  is very likely a clean lemma. It would give the *sharp* halved defect constant
  of OBSTRUCTION Theorem 4 in the quantum case.
* **How large can the coherence gap be?** §5 finds separation `3.5·10⁻³` at
  `r = 2` and `6.6·10⁻³` at `r = 3` from a modest search. The proved bound is
  `|D| ≤ (log 2)/2 = 0.3466`, about fifty times larger. Is there a bound on
  `|A_{(A,S)} − A_{(P_SA,S)}|` in terms of the coherence, and what is the
  supremum of the quantum-cycle margin over triples with transitive shadows?
* **Is the `α ≥ ½` restriction ever binding?** §5 shows the headline separation
  survives on `[0,∞]`. Is there a pair whose exchange rate genuinely differs
  between `[½,∞]` and `[0,∞]`, and is the resulting "data-processing exchange
  metric" a different metric on the cone? That is a purely classical question
  with a quantum origin, and it looks tractable in the cone model of
  `i_cone.py`.
* **The Petz alternative.** BHNOW's eq. (4) is monotone on the compact range
  `0 ≤ α ≤ 2`. The Petz exchange rate `inf_{β∈[0,2]}` is therefore a third
  quantisation, with a *bounded* `β` domain — the first object in this framework
  whose infimum does not reach `β = ∞`, so `Λ` stops being an endpoint. None of
  the geometry of briefs G and I applies to it unchanged.
* **Sufficiency.** BHNOW's quantum laws are necessary only, so `C̃` bounds the
  true rate from above. What *is* the achievable rate, and is the gap visible on
  the §5 witness?

---

## Files

All under `research/quantum/`.

| file | what |
|---|---|
| `q_core.py` | `QSig` — the quantum signature `(A,S)`, its profile `F`, the well-conditioned pencil form `μ_k(β) = eig(A, S^{(β−1)/β})`, the geometric-mean rescaling, the pinching, vectorised `F_grid`/`U_grid`, and the oscillation/midrange/rate machinery with golden refinement |
| `q_spectral.py`, `_output` | §1: `Tr A^β` sees only the spectrum. 300 pairs, 50 conjugations, 2000 non-commuting triples |
| `q_sandwich.py`, `_output` | §3: the structure theorem on 600 random admissible pairs, with the pinched-profile noise floor; §4's pinching check; the isospectral non-spectrality table; the commuting reduction |
| `q_cycle.py`, `q_cycle_quick.txt`, `q_cycle.json` | §5: the search for a 3-cycle with a transitive decohered shadow, `r = 2` and `r = 3`. `--quick` reproduces the committed certified witness; `--tag NAME` writes elsewhere |
| `q_cycle_full_output.txt` | the 10-restart search, **a partial run** — the `r = 2` half complete (best margin `8.408·10⁻³`), the `r = 3` half unfinished |
| `q_certify.py`, `_output`, `q_certify.json` | §5: the `r = 2` witness rebuilt in mpmath at 60 digits from the `2×2` characteristic quadratic and certified at 40, with admissibility margins, the curl criterion and the full-domain robustness check |
| `q_cone.py`, `_output`, `q_cone.json` | §6: tangent-line cone membership, the Hilbert metric of the approximant against the quantum `d`, and the integer-signature `O(1/K)` ladder |
| `q_bhnow.py`, `_output`, `q_bhnow.json` | §7: the two identities B1 and B2, the ratio-vs-difference table on the certified cycle, and the reweighting invariance |

Reproduce in order: `q_spectral.py`, `q_sandwich.py`, `q_cone.py`,
`q_cycle.py --quick`, `q_certify.py`, `q_bhnow.py`. The slow ones are
`q_spectral.py` (about ten minutes; the 2000-triple scan) and `q_cycle.py`
(about ten minutes with `--quick`, hours without). `q_certify.py` and
`q_bhnow.py` read `q_cycle.json`, so run the cycle search first.
