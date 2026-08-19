# Findings — how far is `C_aff` from `C_sig`?

Answer to session brief O. Transcribed from the reporting agent's verified
output; the harness blocked it from writing this file, and the omission was
caught later by the session for brief P. Scripts and data are committed
alongside.

**The brief's premise is false.** The brief, and `research/synergy/FINDINGS.md`
§4, assert that affine implementation forces signature implementation and hence
`C_aff ≤ C_sig`. It does not, and it is not. The two rates are
**incomparable**: on the five non-constant affine classes of quadratic maps
`F_3^2 → F_3` there are nine ordered pairs with `C_aff < C_sig` and four with
`C_aff > C_sig`, all proved.

The sharpest form:

> `x` and `x²+y` have the same fiber signature `(3,3,3)`, so
> `C_sig(x → x²) = C_sig(x²+y → x²) = log 3/log 6 = 0.6131471927654584…`
> exactly. But `C_aff(x → x²) = 0` and `C_aff(x²+y → x²) = 1`.
>
> **Two pairs, the same signature rate to all 40 digits, affine rates `0` and
> `1`.**

Along the way:

* the **first exact values of `C_aff`** in this project — `1/2` for `x² → x`,
  `x² → xy`, `x² → x²+y²`, and `2/3` for `x²+y² → x`, each with a matching
  proved upper bound and an explicit witness processor;
* a **proved strictly superadditive** triple for `C_aff`;
* **mixing copies strictly helps**: `k(2) = 1` but `k(3) = 2`, so the Fekete
  supremum is attained at no `r ≤ 2`;
* **eight of the twelve orbit realisations of the `F_3` tensor three-cycles are
  refuted** under `C_aff`.

---

## 0. The reduction lemma *(proved)*

**Lemma (reduction).** For `f : K^n → K^m` and `g : K^{n'} → K^{m'}`,
`f^{×k} ⪯_aff g^{×r}` iff there are affine `α_1,…,α_r : K^{nk} → K^{n'}` such
that every component `f_s(x_i)` lies in
`span_K{ g_t ∘ α_j } + K·1`.

*Proof.* Write `a = (α_1,…,α_r)` in blocks and `b(y) = By + u`. Component
`(i,s)` of `b ∘ g^{×r} ∘ a` is `Σ_{j,t} B_{(i,s),(j,t)} g_t(α_j(x)) + u_{is}`;
equating with `f_s(x_i)` is the stated span condition. ∎

Call `g ∘ α` a **`g`-atom** and set
`N_k(g→f) = min{ |T| : T a set of `g`-atoms, every `f_s(x_i) ∈ span(T) + K·1` }`.
Then `k_{g→f}(r) = max{k : N_k ≤ r}`, `N_{k+k'} ≤ N_k + N_{k'}`, and
```
C_aff(g→f) = sup_k k / N_k = 1 / inf_k (N_k / k).
```

*Independently checked:* `run_checks.py` compares brute force over all
`729 × 9` affine processor pairs against `N_1 ≤ 1` for all `196` ordered pairs
of `F_3` classes — **0 mismatches** — and reproduces `AFFINE_INPUT_COVERS` in
`src/fn_complexity/cubic_field_maps.py` exactly.

**Corollary A (dimension bound, proved).** If no non-zero combination of
`f_1,…,f_m` is constant then `mk + 1 ≤ m'r + 1`, so
`C_aff(g→f) ≤ m'/m`. In particular `C_aff ≤ 1` for maps `K^n → K`.

**Corollary B (zero rate, proved).** `C_aff = 0` iff `f ∉ span(S_g) + K·1`.

**Corollary C (attainment, proved).** `N_k = k` for some `k` ⟺ `N_1 = 1` ⟺
`f ⪯_aff g`. Hence `f ⪯_aff g ⟹ C_aff = 1`.

---

## 1. The pool *(computed, exhaustive)*

All `3^8 = 6561` degree-≤3 functions `F_3^2 → F_3` sorted into affine
left–right classes: **14 classes**.

| class | representative | degree | signature | orbit size | atoms on `F_3²` |
|---|---|---:|---|---:|---:|
| `constant` | `0` | 0 | `(9)` | 3 | 1 |
| `linear` | `x` | 1 | `(3,3,3)` | 24 | 27 |
| `rank1` | `x²` | 2 | `(6,3)` | 72 | 14 |
| `split` | `xy` | 2 | `(5,2,2)` | 324 | 183 |
| `anisotropic` | `x²+y²` | 2 | `(4,4,1)` | 162 | 105 |
| `parabolic` | `x²+y` | 2 | `(3,3,3)` | 144 | 135 |
| `cubic-711` | `xy+xy²` | 3 | `(7,1,1)` | 216 | 123 |
| `cubic-63` | `y+y²+x²y` | 3 | `(6,3)` | 432 | 98 |
| `cubic-522-a` | `y²+x²y` | 3 | `(5,2,2)` | 1296 | 291 |
| `cubic-522-b` | `xy²` | 3 | `(5,2,2)` | 648 | 291 |
| `cubic-441-a` | `xy+y²+x²y` | 3 | `(4,4,1)` | 1296 | 291 |
| `cubic-441-b` | `y+xy+xy²` | 3 | `(4,4,1)` | 432 | 195 |
| `cubic-333-a` | `y+x²y` | 3 | `(3,3,3)` | 648 | 291 |
| `cubic-333-b` | `x+y+y²+x²y` | 3 | `(3,3,3)` | 864 | 531 |

**This table is the whole point.** There are 4 quadratic and 5 cubic
*signatures*, but **5 and 13 affine classes**. `C_sig` sees the signature
column; `C_aff` sees the classes.

`N_1` for all `14 × 14` ordered pairs is in `n1_matrix.json`, computed by BFS
over the `3^9` function space; `inf` means the target lies outside
`span(S_g) + F_3·1`, so `C_aff = 0` by Corollary B.

---

## 2. Exact values *(proved, with witnesses)*

Lower bounds are Fekete; upper bounds come from one mechanism.

**The support monotone.** For `h` of degree ≤2 on `F_3^N` write `q(h)` for its
degree-2 part, `e(h)` for its degree-1 part, `supp(q) = (rad q)^⊥`. For a
subspace `Q` with basis `q_1,…,q_ρ`, every `q ∈ Q` has
`supp(q) ⊆ Σ_i supp(q_i)`.

**Proposition 1 *(proved)*.** For `g = x²` and `f ∈ {x, xy, x²+y², x²+y}`,
`N_k(g→f) ≥ 2k`, hence `C_aff(x² → f) ≤ 1/2`.

*Proof.* An atom is `h = L²`, `L = ℓ + c`, so `q(h) = ℓ²`, `e(h) = 2cℓ`,
`supp(q(h)) = span{ℓ} ∋ e(h)`. Write `f(x_i) = Σ_j λ_{ij} h_j + const` with `r`
atoms, `W = span{ℓ_j}`, `ρ = dim span{q_j} ≤ r`.
For `f = x`: `Σ λ_{ij} q_j = 0` and `Σ λ_{ij} e_j = x_i`; the `λ_i` are
independent so `ρ ≤ r − k`, and `x_i ∈ W` gives `k ≤ dim W ≤ ρ`.
For `f = xy`: `supp(x_iy_i) = span{x_i,y_i} ⊆ W`, so `2k ≤ ρ ≤ r`.
For `f = x²+y²`: identical. For `f = x²+y`: `x_i ∈ W` and `y_i ∈ W`. ∎

**Proposition 2 *(proved)*.** For `g = x²+y²`, `f = x` over `F_3`,
`N_k ≥ ⌈3k/2⌉`, hence `C_aff(x²+y² → x) ≤ 2/3`.

*Proof.* An atom is `h = L² + M²`, `q = ℓ²+m²`, `e = 2(cℓ + dm)`. **Because
`−1` is not a square in `F_3` the form `⟨1,1⟩` is anisotropic**, so `ℓ²+m² = 0`
forces `ℓ = m = 0`; hence `supp(q) = span{ℓ,m}` always and `e ∈ supp(q)`. Then
`ρ ≤ r − k` and `k ≤ dim W ≤ 2ρ`, so `r ≥ 3k/2`. ∎

That the form is anisotropic is load-bearing: for the isotropic `xy` the same
atom has `q = ℓ_1ℓ_2` but `e` need not lie in `supp(q)`, the argument collapses,
and indeed `x ⪯_aff xy` with `C_aff(xy → x) = 1`. *(This is the observation
brief P turned into a theorem over arbitrary fields — see `ISOTROPY.md`.)*

### Witnesses *(computed, re-verified on all 81 points of `F_3^4`)*

`N_2(x²+y² → x) = 3` with `h_1 = x_2²+x_1²`, `h_2 = x_2²+(x_1+1)²`,
`h_3 = (x_2+1)²+x_1²`, and
```
x_1 = h_1 + 2h_2 + 1 ,        x_2 = h_1 + 2h_3 + 1 .
```
`N_2(x²+y² → xy) = 3` with `x_1y_1 = h_1 − h_2`, `x_2y_2 = h_1 − h_3` (because
`4uv = uv` in `F_3`). `N_2(x² → x) = N_2(x² → xy) = N_2(x² → x²+y²) = 4`,
exhaustively.

### The table

| `g → f` | `N_1` | `N_2` | `C_aff` | `C_sig` (40 digits) | verdict |
|---|---:|---:|---|---|---|
| `x → x` | 1 | 2 | **1** | `1` | equal |
| `x → x²` | ∞ | ∞ | **0** | `0.6131471927654584131297538615321791235349` | `C_aff < C_sig` |
| `x → xy` | ∞ | ∞ | **0** | `0.682606…` | `C_aff < C_sig` |
| `x → x²+y²` | ∞ | ∞ | **0** | `0.792481…` | `C_aff < C_sig` |
| `x → x²+y` | ∞ | ∞ | **0** | `1` | `C_aff < C_sig` |
| `x² → x` | 2 | 4 | **1/2** | `0.6309297535714574370995271143427608542996` | `C_aff < C_sig` |
| `x² → x²` | 1 | 2 | **1** | `1` | equal |
| `x² → xy` | 2 | 4 | **1/2** | `0.6309297535714574370995271143427608542996` | `C_aff < C_sig` |
| `x² → x²+y²` | 2 | 4 | **1/2** | `0.6309297535714574370995271143427608542996` | `C_aff < C_sig` |
| `x² → x²+y` | 3 | 6 | `[1/3, 1/2]` | `0.630929…` | `C_aff < C_sig` |
| `xy → x` | 1 | 2 | **1** | `0.984018804151416701917751837935076442380` | **`C_aff > C_sig`** |
| `xy → x²` | 1 | 2 | **1** | `0.897377759644884381678216046246028718178` | **`C_aff > C_sig`** |
| `xy → xy` | 1 | 2 | **1** | `1` | equal |
| `xy → x²+y²` | 2 | 4 | `[1/2, 1]` | `0.996600…` | open |
| `xy → x²+y` | 2 | 4 | `[1/2, 1]` | `0.984018…` | open |
| `x²+y² → x` | 2 | 3 | **2/3** | `0.9753010616074255613128791871270129862381` | `C_aff < C_sig` |
| `x²+y² → x²` | 1 | 2 | **1** | `0.7737056144690831737404922769356417529303` | **`C_aff > C_sig`** |
| `x²+y² → xy` | 2 | 3 | `[2/3, 1]` | `0.861353…` | open |
| `x²+y² → x²+y²` | 1 | 2 | **1** | `1` | equal |
| `x²+y² → x²+y` | 2 | 4 | `[1/2, 1]` | `0.975301…` | open |
| `x²+y → x` | 1 | 2 | **1** | `1` | equal |
| `x²+y → x²` | 1 | 2 | **1** | `0.6131471927654584131297538615321791235349` | **`C_aff > C_sig`** |
| `x²+y → xy` | 2 | 4 | `[1/2, 1]` | `0.682606…` | open |
| `x²+y → x²+y²` | 2 | 4 | `[1/2, 1]` | `0.792481…` | open |
| `x²+y → x²+y` | 1 | 2 | **1** | `1` | equal |

**9 pairs strictly below, 4 strictly above, 6 undetermined, 6 equal.**
*(Brief P later closed or halved four of the six open brackets; see
`ISOTROPY.md`.)*

Every class here has fiber sizes summing to `9`, so `Z(1) = 9` and the ratio is
exactly `1` at `β = 1`: in this pool `C_sig ≤ 1` automatically.

---

## 3. The two rates are incomparable *(proved)*

Four pairs violate the brief's premise, the extreme one being
```
C_aff(x²+y → x²) = 1  >  C_sig((3,3,3) → (6,3)) = log 3/log 6 = 0.61314719…
```
with the implementation `x² = g(a(x,y))`, `g = x²+y`, `a(x,y) = (x,0)`,
`b = id`.

**Why the implication fails.** The affine model explicitly permits **singular**
input processors (`main.tex`: "The matrices `A` and `B` may be singular,
including zero"). Here `a(x,y) = (x,0)` is `3`-to-`1`: it collapses the
nine-point source onto a three-point line and thereby *manufactures* a fiber of
size `6` out of a map all of whose fibers have size `3`. The signature model
cannot do that — its processors inject each target fiber into an assigned source
fiber, so fiber sizes are inherited, never enlarged.

**And the failure is not a fringe effect.**

| `g` | `f` | `C_sig(g→f)` | `C_aff(g→f)` |
|---|---|---|---|
| `x` (linear) | `x²` | `0.6131471927654584131297538615321791235349` | `0` |
| `x²+y` (parabolic) | `x²` | `0.6131471927654584131297538615321791235349` | `1` |

Same signature, same `C_sig` to 40 digits, `C_aff` at the two extremes of its
range. **`C_aff` is not a function of the fiber signature**, so no functor
between the two processor categories, in either direction, can exist.

---

## 4. Strict superadditivity of `C_aff` *(proved)*

`research/synergy/FINDINGS.md` §4 proves
`C_aff(a⊗b→c) ≥ C_aff(a→c) + C_aff(b→c)` and leaves strictness open. It is
strict. Take `a = x²`, `b = x`, `c = x²+y` over `F_3`:

* `C_aff(a→c) ≤ 1/2` — Proposition 1.
* `C_aff(b→c) = 0` — the `x`-atoms are affine functions, whose span contains
  nothing of degree 2 (Corollary B).
* `C_aff(a⊗b→c) ≥ 1`, because `c ⪯_aff a⊗b` in one shot: with
  `a⊗b : F_3^4 → F_3^2`, `(u,v,x,y) ↦ (u², x)`, take `α(x,y) = (x,0,y,0)` and
  `B = (1,1)`; then `B·(a⊗b)(α(x,y)) = x² + y = c(x,y)`. Verified on all 9
  points.

Hence `C_aff(a⊗b→c) ≥ 1 > 1/2 ≥ C_aff(a→c) + C_aff(b→c)`, **gap ≥ 1/2**.

**The mechanism is the predicted one, not a temperature artefact.** `x²`
supplies a degree-2 part and cannot supply a linear part cheaply (Proposition 1:
two atoms per copy); `x` supplies a linear part and cannot supply a degree-2
part at all. The tensor product lets a single affine processor draw one
component from each. There is no variational formula in sight.

*Caveat.* `b = x` is degree-degenerate. Strictness with `a,b,c` all of degree
exactly 2 and full quadratic rank is **open**.

### 4b. Mixing copies of a single resource strictly helps *(computed + proved)*

```
N_1(x²+y² → x) = 2 ,   N_2(x²+y² → x) = 3 ,
so   k(1) = 0,  k(2) = 1,  k(3) = 2,  and   k(2)/2 = 1/2 < 2/3 = k(3)/3 .
```

`N_2 = 3 < 4 = 2·N_1`: three copies of `x²+y²` deliver two copies of `x`, while
two independent uses of the single-copy solution would need four. The witness
shares the atom `h_1` between the two output blocks — precisely "the block
processor is allowed to mix the `r` copies". With Proposition 2 this pins
`N_{2m} = 3m` exactly and `C_aff(x²+y² → x) = 2/3`.

> **In the smallest non-trivial example, `k(r)/r` is not constant and the Fekete
> supremum is attained at no `r ≤ 2`.** Any claim that a single-shot computation
> determines `C_aff` is wrong.

The same happens for `x²+y² → xy`. Two entries go the other way, both checked
exhaustively: `N_2(x² → x²+y) = 6 = 2N_1` and
`N_2(x² → x) = N_2(x² → xy) = N_2(x² → x²+y²) = 4 = 2N_1`. So subadditivity of
`N_k` is sometimes strict and sometimes not — the affine analogue of Strassen's
additivity question, with both answers present at this size.

---

## 5. Do the `F_3` tensor three-cycles survive? *(computed; eight of twelve refuted)*

`analysis/tensor_cycles_f3.py` reports seven distinct strict three-cycles in
case 3 (homogeneous quadratic maps `F_3^3 → F_3^3`, 50 orbits). **They are
cycles of signatures**, and §3 shows `C_aff` is not a signature invariant, so
the question is only well posed after choosing orbits: 50 orbits on 39
signatures give **12 orbit triples**.

Corollary A plus Corollary C forces the comparison on every `⪯`-comparable pair:
if `a` implements `b` then `C_aff(a→b) = 1 ≥ C_aff(b→a)`, so `a ≺_aff b` is
impossible.

| cycle (signatures, `C_sig` orientation) | orbit triples | verdict |
|---|---|---|
| `{4,2^11,1} < {4,4,3,2^8} < {12,6,6,3}` | `(c37,c35,c8)`, `(c37,c41,c8)` | **refuted** |
| `{4,2^11,1} < {4,4,4,2^7,1} < {12,6,6,3}` | `(c37,c28,c8)`, `(c37,c39,c8)` | **refuted** |
| `{4,2^11,1} < {4,4,4,3,2^6} < {12,6,6,3}` | `(c37,c27,c8)` | **refuted** |
| `{4,2^11,1} < {4^4,2^5,1} < {12,6,6,3}` | `(c37,c31,c8)`, `(c37,c43,c8)` | **refuted** |
| `{4,2^11,1} < {4^5,2^3,1} < {12,6,6,3}` | `(c37,c25,c8)` | **refuted** |
| `{5,4,2^9} < {6,6,4,4,2,2,2,1} < {7,6,4,4,2,2,2}` | `(c45,c12,c13)` | undetermined |
| `{6,3,2^9} < {6,5,2^8} < {9,4,4,4,2,2,2}` | `(c24,c22,c16)`, `(c26,c22,c16)`, `(c34,c22,c16)` | undetermined |

All five cycles through `{12,6,6,3}` are killed by the single fact that its
orbit `c8` is affinely implementable from both other orbits, forcing the
opposite comparison on one edge.

**Reading.** The certified `C_sig` cycles do not automatically descend to the
maps, and where they can be tested they mostly fail. This is brief C's crux one
level down, and it lands on the negative side: the arithmetic cycle results of
this project are results about the signature shadow.

---

## Corrections

1. **`research/synergy/FINDINGS.md` §4 is wrong** where it says "`C_aff ≤ C`:
   affine implementation implies signature implementation, so the signature rate
   is an upper bound". Session brief O repeats it as its premise. Affine
   implementation permits singular input processors, which enlarge fibers; the
   signature model's processors inject and cannot. The §2–§3 signature gaps in
   that document bound nothing about `C_aff`. *(Corrected in place.)*
2. **Brief O's framing of the pool is wrong.** The affine classes, not the
   signatures, are the objects: 5 non-constant affine classes on 4 quadratic
   signatures, 13 on 5 cubic ones. Any table indexed by signature cannot express
   `C_aff`.
3. *(Confirmation.)* The remark "No silent identification with the affine rate"
   in `paper_finite_fields_maps/main.tex` is exactly right and now has explicit
   counterexamples in both directions. It should be promoted from a remark to a
   proposition.

---

## Open

1. **Six brackets** — `xy → x²+y²`, `xy → x²+y`, `x²+y² → xy`,
   `x²+y² → x²+y`, `x²+y → xy`, `x²+y → x²+y²`. *(Four closed or halved by
   brief P; see `ISOTROPY.md`.)*
2. **Is `C_aff = 1` possible without `f ⪯_aff g`?** Corollary C says `N_k = k`
   forces one-step implementability, but `N_k = k+1` for all `k` would still give
   `C_aff = 1`. No example and no obstruction is known. This blocks the open
   brackets.
3. **Is `C_aff(x² → x²+y)` equal to `1/3` or `1/2`?** Needs `N_3 ∈ {7,8,9}`.
4. **Strict superadditivity with all three factors non-degenerate.**
5. **The four undetermined tensor triples**, needing `N_k` in dimension `3^{3k}`.
6. **Is `N_k` additive in general?** §4b exhibits both strict subadditivity and
   exact additivity at `k = 2`. The finite-field affine analogue of Strassen's
   additivity conjecture; what governs the two behaviours is not understood, but
   in both proved cases the deciding structure was `dim supp(q)` for the atoms —
   whether the resource's quadratic part has rank 1 or 2.

---

## Files

| file | what |
|---|---|
| `affine_atoms.py` | the reduction lemma made executable: affine classes, atom sets, `F_3` linear algebra, exhaustive `N_k` searches |
| `run_n1.py` | exact `N_1` for all `14 × 14` pairs by BFS over `3^9` → `n1_matrix.json` |
| `run_n2.py` | exhaustive `N_2` on the quadratic pool → `n2_matrix.json` |
| `run_n2_small.py` | the `N_2 ≤ 3` decision by hull enumeration → `n2_le3.json` |
| `run_rank1_parabolic.py` | depth-5 exhaustive search settling `N_2(x² → x²+y) = 6` |
| `run_witness.py` | explicit processors and coefficients for the computed `N_2`, re-verified on all 81 points |
| `run_checks.py` | independent checks using no reduction lemma: brute force over all `729 × 9` processor pairs vs `N_1 ≤ 1` (196 pairs, 0 mismatches) |
| `run_csig.py`, `run_digits.py` | `C_sig` at 60 digits, reported to 40, with contact temperatures |
| `run_summary.py` | assembles the `C_aff` brackets against `C_sig` |
| `run_tensor_cycles.py` | the seven tensor three-cycles against the affine order, orbit triple by orbit triple |
