# Findings — `C_aff` sees the square class, and it is not a `q`-dependent accident

Answer to session brief P. Every claim below is marked **proved** or
**computed**; "computed" always means an exhaustive search with an exact
termination argument, never a sample. `C_sig` values are certified at 40 digits
and their endpoints checked against closed forms.

## Summary

* **The congruence is settled, and the hypothesis is not merely true — it is a
  theorem about arbitrary fields** *(proved)*. For every binary quadratic form
  `Q` over any field `K` with `char K ≠ 2`,

  ```
  N_k(Q → x) = ceil(3k/2)  and  C_aff(Q → x) = 2/3     if Q is anisotropic,
  N_k(Q → x) = k           and  C_aff(Q → x) = 1       if Q is isotropic.
  ```

  Nothing in the proof is finite-field-specific. So over `F_q`,
  `C_aff(x²+y² → x)` is `2/3` for `q ≡ 3 mod 4` and `1` for `q ≡ 1 mod 4`; and
  the uniformly anisotropic family `x² − n y²` has `C_aff(→ x) = 2/3` at **every**
  odd `q`. The same statement holds over `Q_p` and over `R`.
* **`C_sig` does not jump, and the amplification factor grows without bound**
  *(computed, 40 digits)*. Along the uniform family, `C_sig(x² − n y² → x)`
  climbs smoothly `0.97530, 0.98435, 0.98858, …` to `1`, while
  `C_sig(xy → x)` climbs to `1` faster. `C_aff` separates the two classes by
  exactly `1/3` at every `q`; `C_sig` separates them by `0.0087` at `q = 3` and
  by `0.0014` at `q = 47`. The brief's "factor of 38" is the `q = 3` value of a
  quantity that is `245` at `q = 47` and appears to grow linearly in `q`.
* **Four of the six open brackets of `affine_rate` §2 are closed or cut in
  half** *(proved)*:

  ```
  C_aff(x²+y² → xy)   = 2/3        was [2/3, 1]
  C_aff(x²+y  → xy)   = 1/2        was [1/2, 1]
  C_aff(x²+y  → x²+y²)= 1/2        was [1/2, 1]
  C_aff(x²+y² → x²+y) in [1/2,2/3] was [1/2, 1]
  ```

  All four are **strictly below** their `C_sig` (`0.8614, 0.6826, 0.7925,
  0.9753`). The two brackets with the *split* resource `xy` remain open at
  `[1/2, 1]`. For one of them a third exact value is now known
  *(computed, exhaustive)*: `N_k(xy → x²+y²) = 2, 4, 6` at `k = 1, 2, 3`, so the
  lower bound `1/2` is not an artefact of stopping at `k = 2`.
* **The crux — "can `C_aff = 1` without `f ⪯_aff g`?" — is answered "no" on
  every pair of this pool where the answer is now known** *(proved)*, and the
  two surviving cases are isolated and explained. The tool that closes the four
  is a **specialisation lemma**: if `f` restricted to some affine line is a
  bijective affine function, then `N_k(g → f) ≥ N_k(g → x)` computed in `k`
  variables, so the `2/3` bound for an **anisotropic resource** transports from
  the linear target to `f`. Over `F_q` the targets admitting such a line are
  exactly `x`, `xy`, `x²+y`; the targets that do **not** are `x²` and the
  anisotropic forms. The two open brackets are exactly the two whose **resource
  is isotropic**, where the route is vacuous because an isotropic resource
  implements `x` in one step.
* **Rationality: every value proved here comes from an `N_k` that is exactly
  linear or exactly `ceil(3k/2)`, i.e. eventually exactly linear along
  arithmetic progressions** *(proved, case by case)*. Whether that is general
  is **open** and is the honest analogue of the rationality of the
  matrix-multiplication exponent; nothing here bears on it.
* **The anisotropic-pair question is empty over finite fields** *(proved,
  classical)*. Over `Q_p` it is non-empty, and it reduces to **the same
  combinatorial question** as the two brackets left open here — with the
  specialisation lemma failing for the opposite reason (there, the *target* is
  anisotropic, so it carries no affine line on which it is affine). What the
  `p`-adic version needs is stated in §6, together with one lever that exists
  over `Q_p` and not over `F_q`. As a by-product, Theorem 1 is new over `Q_p`:
  `C_aff(x² + d y² → x) = 2/3` for every anisotropic class.
* **Correction.** `research/affine_rate/FINDINGS.md` **does not exist in the
  repository.** Commit `d4ffdf1` says the harness refused the write and the body
  went into a session report only. Everything the brief asks to be read there
  was reconstructed from the committed code and re-derived from scratch; where a
  claim of that session could be checked, it checked out (§7).

---

## 0. Setup, and the reduction lemma restated

Throughout `K` is a field with `char K ≠ 2` and `|K| ≥ 3`; `q` denotes an odd
prime power. Maps are `f, g : K² → K` of total degree at most `2`, processors
are affine and may be singular
(`paper_finite_fields_maps/main.tex`, §"Processors, equivalence, and order"),
and

```
k_{g→f}(r) = max{ k : f^{×k} ⪯_aff g^{×r} },
C_aff(g→f) = lim_r k(r)/r = sup_r k(r)/r        (Fekete).
```

**Lemma 0 (reduction; restated).** *An affine `a : K^{2k} → K^{2r}` is a list of
`r` affine maps `α_j : K^{2k} → K²`, and an affine `b : K^r → K^k` is a list of
`k` affine functionals of the `r` outputs. Hence*

```
f^{×k} ⪯_aff g^{×r}  ⟺  ∃ α_1..α_r  with  f(x_i) ∈ span_K{g∘α_1,…,g∘α_r} + K·1
                                            for i = 1..k.
```

Write `N_k(g→f)` for the least such `r`. `N` is subadditive (block-diagonal
composition), `k(r) = max{k : N_k ≤ r}`, and

```
C_aff(g→f) = 1 / inf_k (N_k/k) = sup_k k/N_k ,
```

so **every computed `k/N_k` is a rigorous lower bound and every proved lower
bound on `N_k` is a rigorous upper bound on `C_aff`.** This is the lemma the
brief-O session proved; it is implemented in `affine_atoms.py` and restated in
that file's docstring.

**Lemma 1 (jet form).** *For `|K| ≥ 3` a polynomial of total degree `≤ 2` has
all individual exponents `< |K|`, so it is determined as a function by its
coefficients. Let*

```
J(h) = (quadratic part of h, linear part of h) ∈ Sym²(V*) ⊕ V*,   V = K^{2k},
```

*whose kernel is exactly the constants. Then `h ∈ span{h_j} + K·1` iff
`J(h) ∈ span{J(h_j)}`, and consequently*

```
N_k(g→f) = min { dim T : T ⊆ Sym²(V*) ⊕ V*  a subspace spanned by atom jets,
                         T ⊇ span{ J(f(x_1)), …, J(f(x_k)) } } .
```

*Explicitly, for `g(u,v) = a u² + b uv + c v² + d u + e v` and
`α(z) = (ℓ(z)+s, m(z)+t)`,*

```
J(g∘α) = ( a ℓ² + b ℓm + c m² ,  (2as + bt + d)·ℓ + (bs + 2ct + e)·m ).
```

*In particular, if the quadratic part `Q` of `g` is nondegenerate and `d = e = 0`,
the polarisation matrix `M = [[2a,b],[b,2c]]` is invertible, so as `(s,t)` runs
over `K²` the linear part runs over all of `span{ℓ, m}`:*

```
atoms(Q) = { ( Q(ℓ,m) , αℓ + βm ) : ℓ, m ∈ V*, α, β ∈ K }.        (*)
```

*And for `g = x² + y` the linear part `2sℓ + m` is unconstrained because `m` is
free, so*

```
atoms(x²+y) = { ( ℓ² , L ) : ℓ, L ∈ V* }.                          (**)
```

**Validation** *(computed)*. `run_isotropy_check.py` recomputes `N_1` and the
exhaustive "`N_2 ≤ 3`" decision for all `25` ordered pairs of `F_3` quadratic
classes from the jet representation and compares them with the committed
value-table results in `n1_matrix.json` and `n2_le3.json`. **All 50 entries
agree.** The jet form is a faithful, and much smaller, model: for `k = 3` the
value-table space has `3^729` functions, the jet space has `3^27` vectors.

**Lemma 2 (support).** *Let `Q` be anisotropic, `ℓ, m ∈ V*`, `W = span{ℓ,m}`,
`π = (ℓ,m) : V → K²`, `q = Q∘π`. Then `rad q = W^⊥` and `rank q = dim W`. In
particular `q = 0` iff `ℓ = m = 0`.*

*Proof.* `B_q(u,v) = B_Q(πu, πv)`, so `rad q = π^{-1}(Z ∩ Z^⊥)` with
`Z = im π`. If `0 ≠ z ∈ Z ∩ Z^⊥` then `2Q(z) = B_Q(z,z) = 0`, contradicting
anisotropy; so `Z ∩ Z^⊥ = 0` and `rad q = ker π = W^⊥`. ∎

This is the exact place where **anisotropy enters**, and it is the `F_3`
argument of the brief (`ℓ² + m² = 0 ⟹ ℓ = m = 0`) in invariant form. For
isotropic `Q` it fails: `xy` with `m = 0` gives `q = 0` but `W = ⟨ℓ⟩ ≠ 0`, and
those atoms are exactly the affine functions, which is why `xy` implements `x`
in one step.

**Lemma 3 (span controls support).** *If `q_1,…,q_N` satisfy `rad q_j = W_j^⊥`
and `S = span{q_j}` has dimension `s`, then `Σ_j W_j = Σ_{i=1}^{s} W_{j_i}` for
any independent spanning subset; hence `dim Σ_j W_j ≤ 2s` when every
`dim W_j ≤ 2`.*

*Proof.* If `q_j = Σ_i c_i q_{j_i}` then `rad q_j ⊇ ∩_i rad q_{j_i}`, i.e.
`W_j^⊥ ⊇ (Σ_i W_{j_i})^⊥`, i.e. `W_j ⊆ Σ_i W_{j_i}`. ∎

---

## 1. The congruence, settled

### 1.1 The lower bound (Proposition 2 of `affine_rate`, generalised)

**Proposition A** *(proved)*. *Let `K` have `char ≠ 2`, `|K| ≥ 3`, let `Q` be an
anisotropic binary quadratic form over `K`, let `V` be any finite-dimensional
`K`-space and let `z_1,…,z_k ∈ V*` be linearly independent. If `r` atoms
`h_1,…,h_r` of `Q` on `V` satisfy `z_i ∈ span{h_j} + K·1` for all `i`, then*

```
r ≥ ceil(3k/2).
```

*Proof.* Put `U = span{h_j} + K·1`, so `dim U ≤ r + 1`. Let `π_2` take the
degree-2 homogeneous part (well defined by Lemma 1). Then
`span{1, z_1,…,z_k} ⊆ ker(π_2|_U)` has dimension `k+1`, and
`π_2(U) = span{q_j}` where `q_j = Q(ℓ_j, m_j)`. Hence

```
dim span{q_j} = dim U − dim ker(π_2|_U) ≤ (r+1) − (k+1) = r − k.
```

Matching linear parts, `z_i = Σ_j λ_j L_j` with `L_j ∈ W_j = span{ℓ_j,m_j}`, so
`z_1,…,z_k ∈ W := Σ_j W_j` and `dim W ≥ k`. Lemmas 2 and 3 give
`dim W ≤ 2·dim span{q_j} ≤ 2(r−k)`. Therefore `k ≤ 2(r−k)`, i.e. `r ≥ 3k/2`. ∎

The brief asked whether Proposition 2 "generalises verbatim". It does more than
that: the proof never uses finiteness of `K`, never uses `n = 2` copies of
anything, and never uses a residue field. It uses exactly one arithmetic fact,
Lemma 2, and that fact is the definition of anisotropy.

### 1.2 The matching witness

**Proposition B** *(proved; verified computationally)*. *Let `Q` be any
nondegenerate binary quadratic form over `K`, `char K ≠ 2`, with polarisation
matrix `M`. Then `N_1(Q → x) ≤ 2` and `N_2(Q → x) ≤ 3`, hence
`N_k(Q → x) ≤ ceil(3k/2)` for all `k`.*

*Proof.* `Q(ℓ + s, m + t) = Q(ℓ,m) + (M(s,t))_1·ℓ + (M(s,t))_2·m + Q(s,t)`. For
`k = 1` take `α_1(x,y) = (x,y)` and `α_2(x,y) = (x+s, y+t)` with
`M(s,t) = (1,0)`; then `h_2 − h_1 − Q(s,t) = x`. For `k = 2` on
`K^4 = {(x_1,y_1,x_2,y_2)}` take

```
α_1 = (x_1, x_2),
α_2 = (x_1 + s_2, x_2 + t_2)   with  M(s_2,t_2) = (1,0),
α_3 = (x_1 + s_3, x_2 + t_3)   with  M(s_3,t_3) = (0,1),
```

giving `x_1 = h_2 − h_1 − Q(s_2,t_2)` and `x_2 = h_3 − h_1 − Q(s_3,t_3)`.
Block-diagonal composition then gives `N_{2m} ≤ 3m` and `N_{2m+1} ≤ 3m+2`. ∎

The `k = 2` witness **mixes the two copies** — all three atoms read `x_1` and
`x_2` and ignore `y_1, y_2` entirely. This is the concrete form of brief O's
"mixing copies strictly helps".

`run_isotropy_family.py` re-verifies both witnesses by evaluating the actual
functions on **all** `q²` resp. `q⁴` points, for `q = 3, 5, 7, 11, 13`, and
`run_isotropy_witness.py` extends the check to `q = 17, 19`. All pass.

### 1.3 The theorem and the congruence

**Theorem 1** *(proved)*. *For `K` with `char K ≠ 2`, `|K| ≥ 3`:*

* *`Q` anisotropic ⟹ `N_k(Q → x) = ceil(3k/2)` for every `k ≥ 1`, and*
  `C_aff(Q → x) = 2/3`.
* *`Q` isotropic and nondegenerate ⟹ `N_k(Q → x) = k` and `C_aff(Q → x) = 1`.*

*Proof.* First part: Propositions A and B. Second part: an isotropic
nondegenerate binary form is `λ·(uv)` after an invertible linear change, and
`α_i(z) = (x_i, 1)` in those coordinates gives `Q∘α_i = λ x_i`; so `N_k ≤ k`,
and `N_k ≥ k` because `1, x_1, …, x_k` are independent. ∎

**Corollary 1 (the congruence)** *(proved)*. *Over `F_q`, `q` odd:*

```
C_aff(x² + y² → x) = 2/3   if q ≡ 3 (mod 4),
C_aff(x² + y² → x) = 1     if q ≡ 1 (mod 4),
C_aff(x² − n y² → x) = 2/3 for every odd q  (n any non-residue).
```

*Proof.* `x² + y²` is anisotropic iff `−1` is a non-residue iff `q ≡ 3 (mod 4)`;
`x² − n y²` is anisotropic for every odd `q` because `x² = n y²` with `y ≠ 0`
would make `n` a square. Apply Theorem 1. ∎

**The hypothesis of the brief is therefore confirmed, and upgraded from a
`q`-by-`q` computation to a single proof.** The jump is not a feature of small
`q`; it is the definition of the square class of the discriminant.

### 1.4 Exhaustive confirmation

*(computed; `isotropy_family.csv`)*

| `q` | form | anisotropic | `N_1` exhaustive | `N_2 ≤ 3` exhaustive | witnesses re-verified | `C_aff` |
|---:|---|---|---:|---:|---|---|
| 3 | `x² − 2y²` | yes | 2 | 3 | on `3²` and `3⁴` points | `2/3` |
| 3 | `x² + y²`  | yes | 2 | 3 | on `3²` and `3⁴` points | `2/3` |
| 5 | `x² − 2y²` | yes | 2 | 3 | on `5²` and `5⁴` points | `2/3` |
| 5 | `x² + y²`  | **no** | **1** | **2** | one-step `α(x,y) = x·(1,2)+(0,1)` | `1` |
| 7 | `x² − 3y²` | yes | 2 | — | on `7²` and `7⁴` points | `2/3` |
| 7 | `x² + y²`  | yes | 2 | — | on `7²` and `7⁴` points | `2/3` |
| 11 | `x² − 2y²` | yes | 2 | — | on `11²` and `11⁴` points | `2/3` |
| 11 | `x² + y²` | yes | 2 | — | on `11²` and `11⁴` points | `2/3` |
| 13 | `x² − 2y²` | yes | 2 | — | on `13²` and `13⁴` points | `2/3` |
| 13 | `x² + y²` | **no** | **1** | — | one-step `α(x,y) = x·(1,5)+(0,1)` | `1` |

`N_1` is exhaustive over the full jet atom set at every `q` listed; the
`N_2 ≤ 3` decision is exhaustive at `q = 3, 5` (at larger `q` the `k = 2` atom
enumeration is `q^{10}` and was not run — but Theorem 1 makes it unnecessary).
The `N_1` column alone already exhibits the congruence at `q = 5` and `q = 13`.

---

## 2. `C_sig` crosses smoothly, and by less and less

*(computed, 40 digits; `isotropy_csig.csv`, `isotropy_csig_gap.csv`)*

Fibre signatures over `F_q` of a map `F_q² → F_q`: `x` has `q` fibres of size
`q`; an anisotropic `Q` is the norm form of `F_{q²}`, so it has one fibre of
size `1` and `q−1` of size `q+1`; an isotropic `Q` has one fibre of size `2q−1`
and `q−1` of size `q−1`. Every signature sums to `q²`, so
`log Z_g/log Z_f = 1` at `β = 1`, the `β = 0` endpoint is `log q / log q = 1`
exactly, and the `β = ∞` endpoints are `log(q+1)/log q` (anisotropic) and
`log(2q−1)/log q` (isotropic). The computed endpoints match those closed forms
to all printed digits; the infimum is attained at an interior contact near
`β ≈ 0.29–0.44`.

| `q` | `C_sig(xy → x)` | `C_sig(x²−ny² → x)` | `C_sig` gap | `C_aff` gap | amplification |
|---:|---|---|---:|---:|---:|
| 3 | `0.9840188042` | `0.9753010616` | `0.00871774` | `1/3` | `38.24` |
| 5 | `0.9931322258` | `0.9843464172` | `0.00878581` | `1/3` | `37.94` |
| 7 | `0.9958776174` | `0.9885812248` | `0.00729639` | `1/3` | `45.68` |
| 11 | `0.9978442410` | `0.9926324978` | `0.00521174` | `1/3` | `63.96` |
| 13 | `0.9982893949` | `0.9937552981` | `0.00453410` | `1/3` | `73.52` |
| 17 | `0.9988111705` | `0.9952264278` | `0.00358474` | `1/3` | `92.99` |
| 19 | `0.9989751802` | `0.9957337974` | `0.00324138` | `1/3` | `102.84` |
| 23 | `0.9992035183` | `0.9964868175` | `0.00271670` | `1/3` | `122.70` |
| 29 | `0.9994107434` | `0.9972291732` | `0.00218157` | `1/3` | `152.80` |
| 31 | `0.9994592306` | `0.9974129073` | `0.00204632` | `1/3` | `162.89` |
| 37 | `0.9995686868` | `0.9978447530` | `0.00172393` | `1/3` | `193.36` |
| 41 | `0.9996213344` | `0.9980621569` | `0.00155918` | `1/3` | `213.79` |
| 43 | `0.9996434417` | `0.9981555936` | `0.00148785` | `1/3` | `224.04` |
| 47 | `0.9996812035` | `0.9983184226` | `0.00136278` | `1/3` | `244.60` |

Forty-digit values are in the CSV; the `q = 3` entries are

```
C_sig(xy → x)      = 0.98401880415141670191775183793507644238
C_sig(x²+y² → x)   = 0.9753010616074255613128791871270129862381
```

reproducing the committed `csig_matrix.json` and the brief's numbers
`0.975301, 0.993132, 0.988581, 0.992632, 0.998289, 0.998811, 0.995734` at
`q = 3, 5, 7, 11, 13, 17, 19` exactly. (The brief's list mixes the two families:
at `q ≡ 1 mod 4` the entry is the *isotropic* value, because `x²+y²` is
isotropic there. The uniform anisotropic family is monotone.)

**The headline.** `C_aff` separates the isotropic from the anisotropic class by
`1/3` at every `q`; `C_sig` separates them by an amount that **tends to zero**.
The brief's factor of `38` is the value at `q = 3` of a quantity that grows
roughly linearly in `q`, reaching `245` at `q = 47`. In the large-`q` limit the
signature rate erases the square class completely, and the affine rate does not
move at all.

---

## 3. The six open brackets

`affine_rate` §2 left six pairs bracketed at `[·, 1]` with `N_1 = 2`. Four are
now settled or halved. The engine is:

**Proposition C (specialisation)** *(proved)*. *Suppose `γ : K → K²` is affine
with `f(γ(z)) = λz + μ`, `λ ≠ 0`. Then for every `g`,*

```
N_k(g → f) ≥ Ñ_k(g → x),
```

*where `Ñ_k` is the same minimum computed on `K^k` with targets `z_1,…,z_k`.*

*Proof.* `ι = γ^{×k} : K^k → K^{2k}` is affine, composition with `ι` is linear on
functions and fixes constants, `(g∘α_j)∘ι = g∘(α_j∘ι)` is again a `g`-atom, and
`f(x_i)∘ι = λ z_i + μ`. ∎

**Which `f` admit such a `γ`** *(proved)*. Writing `f = A + λ_1 + c` with `A` the
quadratic part: a direction `w ≠ 0` works iff `A(w) = 0` and
`B_A(w, w_0) + λ_1(w) ≠ 0` for some `w_0`. Over `F_q` this holds for `x`
(`A = 0`), `xy` (`A` isotropic of rank 2) and `x² + y` (`A = x²` of rank 1,
`λ_1 = y` nonzero on `ker x`), and **fails** for `x²` (`A = x²`, `λ_1 = 0`) and
for anisotropic `A`. Combining with Proposition A:

**Theorem 2** *(proved)*. *If `Q` is anisotropic over `K` (`char ≠ 2`) and `f`
admits such a line, then `N_k(Q → f) ≥ ceil(3k/2)` and `C_aff(Q → f) ≤ 2/3`. Over
`F_q` this covers `f ∈ {x, xy, x²+y}`.*

### 3.1 `x² + y² → xy` — **closed at `2/3`** *(proved)*

Upper bound `2/3` by Theorem 2. Lower bound `2/3` by an explicit witness with
`N_2 = 3`: with `u_i = x_i + y_i`, `v_i = x_i − y_i` and `Q = x² − n y²`,

```
h_1 = Q(u_1, u_2),   h_2 = Q(v_1, u_2),   h_3 = Q(u_1, v_2)
⟹  4·x_1y_1 = h_1 − h_2 ,   (−4n)·x_2y_2 = h_1 − h_3 ,
```

and `N_1 ≤ 2` from `4·xy = Q(u,0) − Q(v,0)`. Hence `N_k(Q → xy) = ceil(3k/2)`
and

```
C_aff(anisotropic → xy) = 2/3   over every F_q with q odd,
```

against `C_sig(x²+y² → xy) = 0.8613531…` at `q = 3`. Both witnesses are
re-verified on all `q²` and `q⁴` points for `q = 3,…,19`
(`run_isotropy_witness.py`). The `F_3` value `N_2 = 3` also reproduces the
committed exhaustive entry in `n2_le3.json`.

### 3.2 `x² + y → xy` and `x² + y → x² + y²` — **both closed at `1/2`** *(proved)*

By `(**)`, the atoms of `g = x²+y` are exactly `(ℓ², L)` with `ℓ` and `L`
ranging **independently** over `V*` — the resource's own linear term `y` buys an
unconstrained linear part. Consequently:

**Theorem 3** *(proved)*. *Let `g = x² + y` over `K`, `char K ≠ 2`, and let `f`
be nonconstant of degree `≤ 2` with quadratic part of rank `ρ`. Then*

```
N_k(g → f) = ρ k   (ρ ≥ 1),        N_k(g → f) = k   (ρ = 0),
C_aff(g → f) = 1/ρ (ρ ≥ 1),        C_aff(g → f) = 1 (ρ = 0).
```

*Proof.* Write the target jets as `(F_i, λ_i)`, `F_i` supported in block `i`.
*Lower bound:* projecting to `Sym²`, `F_i ∈ span{ℓ_j²}` for all `i`, so
`Σ_i F_i ∈ span{ℓ_j²}`; that form has rank `ρk`, and rank is subadditive over
rank-one forms, so `r ≥ ρk`. For `ρ = 0` project to `V*` instead: the `λ_i` are
independent, so `r ≥ k`. *Upper bound:* diagonalise `F_i = Σ_{t≤ρ} a_t ℓ_{it}²`
inside block `i`'s variables and take the `ρk` atoms `(ℓ_{it}², L_{it})`. The
coefficient vectors `λ^{(i)} ∈ K^{ρk}` are supported on disjoint blocks hence
independent, so the free `L_{it}` can be solved for the linear parts. ∎

With `ρ(xy) = ρ(x²+y²) = 2` this gives `C_aff = 1/2` for both, against
`C_sig = 0.6826062…` and `0.7924813…`. It also *explains* the rest of the
committed `x²+y` row: `ρ(x) = 0 ⟹ 1`, `ρ(x²) = ρ(x²+y) = 1 ⟹ 1`, and the
computed `N_2 = 4` for the two rank-two targets is `ρk` on the nose.

Witnesses (verified on all `q²` points for `q = 3,…,19`):
`4·xy = g(x+y,0) − g(x−y,0)`, `x²+y² = g(x,0) + g(y,0)`, `x²+y = g(x,y)`,
`x = g(0,x)`.

### 3.3 `x² + y² → x² + y` — **halved to `[1/2, 2/3]`** *(proved)*

`x² + y` admits the line `γ(z) = (0,z)`, so Theorem 2 gives `C_aff ≤ 2/3`. The
exhaustive `N_2 = 4` of `n2_le3.json` gives `C_aff ≥ 1/2`. Both endpoints of the
new bracket are below `C_sig = 0.9753011…`, so the pair is decided for the
purposes of the conjecture even though its exact value is not. Note that
`N_2 = 4 > 3 = ceil(3·2/2)`, so Proposition A is **not** tight here; deciding
`N_3 ∈ {5,6}` would separate `≥ 3/5` from the `N_k = 2k` pattern and was beyond
the exhaustive reach of this session (§8).

### 3.4 `xy → x² + y²` and `xy → x² + y` — **still `[1/2, 1]`** *(open)*

These are the two brackets whose **resource** is isotropic. For `g = xy` the
atoms with `m = 0` are exactly the affine functions, so Lemma 2 fails, and the
specialisation route of Proposition C is vacuous: it reduces `N_k(xy → f)` to
`Ñ_k(xy → x) = k`, which says nothing. What is proved:

* *(proved)* `N_k ≥ k + 1` for both. If `r = k` then the atom-spanned `T` of
  Lemma 1 equals the hull `H = span{target jets}`, so `H` would have to contain
  `k` independent atoms. For `xy → x²+y²`, a nonzero element of `H` is
  `Σ c_i F_i` of rank `2·#{i : c_i ≠ 0}`; rank `≤ 2` forces one term, and a
  single `c_i F_i` is anisotropic, hence not `ℓm`. For `xy → x²+y`, an element
  is `(Σ c_i x_i², Σ c_i y_i)` whose linear part never lies in the support
  `⟨x_i : c_i ≠ 0⟩` of its quadratic part, as `(*)` requires. So `H` contains no
  nonzero atom and `N_k ≥ k+1`. This matches, and reproves, the `Q_2` bound
  `k(n) ≤ n−1` of `p_adic_exchange_rate_attempts.md`.
* *(proved, restricted ansatz)* If every atom's quadratic part is **diagonal in
  the target's own basis** `{x_1,y_1,…,x_k,y_k}`, then `N_k(xy → x²+y²) = 2k`
  exactly over `F_3`. Diagonal products are `a·e_s²` or `a(e_s² − e_t²)`, i.e.
  vertices and edges of a graph on `2k` vertices; the edge-span on a component
  `C` is `{v : Σ_C v = 0}` of dimension `|C|−1`, and a vertex atom raises it to
  `|C|`. Every target `e_{2i−1}² + e_{2i}²` has nonzero coordinate sum on each
  component it meets, so **every** component must carry a vertex atom, and with
  `c` components `r ≥ (2k − c) + c = 2k`. This is a rigorous negative result for
  the ansatz only — exactly the status of "Attempt 3" in the `Q_2` document —
  and it says any construction beating `1/2` must be non-diagonal.
  **The ansatz is not obviously lossy**: run the same analysis with the roles of
  the two square classes exchanged (`Q`-form atoms `a(e_s² + e_t²)`, hyperbolic
  targets `e_{2i−1}² − e_{2i}²`) and the diagonal minimum is `ceil(3k/2)` — which
  is the *true* `N_k` of §3.1. `run_diagonal_ansatz.py` brute-forces both rows at
  `k = 1,2,3` and both match: `2, 4, 6` and `2, 3, 5`. The whole `1/2`-versus-
  `2/3` asymmetry between the two square classes is visible in this one
  two-line combinatorial model: "sum" atoms have nonzero coordinate sum and a
  spanning forest suffices, "difference" atoms do not and every component must
  buy a rank-one atom.
* *(computed, exhaustive)* **`N_3(xy → x²+y²) = 6`, exactly `2k`.** For a purely
  quadratic target the atoms may be taken with `s = t = 0`, so
  `N_k(xy → x²+y²)` is *exactly* the least dimension of a product-spanned
  subspace of `Sym²((F_3^{2k})*)` containing the `F_i` — a `21`-dimensional
  ambient space with `132 860` products at `k = 3`. `run_split_to_sum2.py`
  decides `N_k ≤ k+2` exhaustively: a solution spans `T ⊇ H` of dimension
  `≤ k+2`, so `T = H`, `H + ⟨p⟩` or `H + ⟨p,p'⟩` with `p, p'` products; the first
  product may be taken from a representative of each orbit of the stabiliser of
  `H` in `GL_{2k}(F_3)` (generated by the block permutations and the per-block
  similitudes of `x²+y²`), and the second is indexed by its residue modulo
  `H + ⟨p⟩` up to scalars. At `k = 2` the method returns the known `N_2 = 4`
  from `21` orbit representatives; at `k = 3` it returns `N_3 > 5` from `70`
  representatives in `163 s`, and `N_3 ≤ N_2 + N_1 = 6`, so `N_3 = 6`. Together
  with `N_1 = 2` and `N_2 = 4`, the sequence is `2, 4, 6` — the diagonal
  prediction, on the nose, three times.
* *(computed)* `run_n3.py` decides `N_3 ≤ 4` exhaustively for all `25` `F_3`
  quadratic pairs in the full jet model, by the argument that a four-atom
  solution has span of dimension `≤ 4` containing the three-dimensional hull,
  hence equal to `H` or to `H + ⟨a⟩` for a single atom `a`. Results in
  `n3_le4.json`; they agree with everything proved above.

The upper bound `1` is untouched: `C_aff ≤ 1` is all that is proved, so **the
crux question "can `C_aff = 1` without `f ⪯_aff g`?" survives exactly on these
two pairs**, and their `C_sig` values `0.9966000…` and `0.9840188…` sit close
enough to `1` that they do not decide it either.

### 3.5 Where that leaves the conjecture

| pair | `N_1` | `N_2` | `C_aff` | `C_sig` | verdict |
|---|---:|---:|---|---|---|
| `xy → x²+y²` | 2 | 4 | `[1/2, 1]` | `0.9966000` | open |
| `xy → x²+y` | 2 | 4 | `[1/2, 1]` | `0.9840188` | open |
| `x²+y² → xy` | 2 | 3 | **`2/3`** | `0.8613531` | `C_aff < C_sig` |
| `x²+y² → x²+y` | 2 | 4 | `[1/2, 2/3]` | `0.9753011` | `C_aff < C_sig` |
| `x²+y → xy` | 2 | 4 | **`1/2`** | `0.6826062` | `C_aff < C_sig` |
| `x²+y → x²+y²` | 2 | 4 | **`1/2`** | `0.7924813` | `C_aff < C_sig` |

Adding these to the nine pairs already decided in `affine_rate` §2, the tally on
the `F_3` quadratic pool is: **13 pairs with `C_aff < C_sig`, 4 with
`C_aff > C_sig` (all of them one-step implementable, where `C_aff = 1`), 6 with
both equal to 1, and 2 undetermined** (`13 + 4 + 6 + 2 = 25`). Every pair with
`C_aff > C_sig` has
`N_1 = 1`. So the statement

> off the one-step-implementable locus, `C_aff < C_sig`

now has **no counterexample and only two untested cases** on this pool, and it
is still exactly equivalent to "`C_aff = 1` forces `f ⪯_aff g`" there.

---

## 4. Rationality, stated honestly

`C_aff = 1/lim_k (N_k/k)` with `N_k` subadditive. A subadditive sequence's limit
is an infimum of rationals with unbounded denominators and **need not be
rational**. Nothing in this session bears on that question in general. What is
true is that every value this project can currently prove comes from an `N_k`
that is *exactly* linear, or exactly `ceil(3k/2)`, from `k = 1` on:

| pair | `N_k` | proved for | `C_aff` |
|---|---|---|---|
| `Q` anisotropic `→ x` | `ceil(3k/2)` | all `k`, any field `char ≠ 2` | `2/3` |
| `Q` anisotropic `→ xy` | `ceil(3k/2)` | all `k`, any `F_q`, `q` odd | `2/3` |
| `Q` isotropic `→ x` | `k` | all `k` | `1` |
| `x²+y → f`, `rank ρ ≥ 1` | `ρk` | all `k`, any `K`, `char ≠ 2` | `1/ρ` |
| `x²+y → f`, `rank 0` | `k` | all `k` | `1` |
| `x² → f` purely quadratic of rank `ρ` | `ρk` | all `k` | `1/ρ` |
| `x² → x` | `2k` | all `k` | `1/2` |

(The last two are re-proofs of brief O's `1/2` values by the same method: for
`g = x²` the atoms are `(ℓ², αℓ)`, the quadratic parts are rank one, and Lemma 3
degenerates to `dim Σ⟨ℓ_j⟩ ≤ dim span{ℓ_j²}`. For `x² → x`, `Σλ_jℓ_j² = 0` costs
`k` dimensions and `x_i ∈ Σ⟨ℓ_j⟩` costs `k` more, so `N_k = 2k`. The same
argument gives the new bound `N_k(x² → x²+y) ≥ 2k`, i.e.
`C_aff(x² → x²+y) ≤ 1/2`, consistent with the committed bracket `[1/3, 1/2]`.)

`ceil(3k/2)` is not linear, but it is linear on each residue class mod `2` and
`inf_k N_k/k = 3/2` is attained on the even ones, which is all Fekete needs. So:

* *(proved)* every `C_aff` known in this project is rational, with denominator
  `1, 2` or `3`;
* *(open)* whether `N_k` is eventually linear in general, equivalently whether
  `C_aff` is always rational. This is the affine analogue of the rationality of
  the matrix-multiplication exponent and should be labelled as such rather than
  guessed at. Note the analogy is not perfect: `N_k` here is a minimum over a
  *finite* set for each `k` (over `F_q`), so `C_aff` is at worst a limit of a
  computable sequence — but that gives no bound on the denominator.

---

## 5. The anisotropic-pair question is empty over finite fields

*(proved; classical)*

* Nondegenerate binary quadratic forms over `F_q` (`q` odd) are classified up to
  isometry by `disc ∈ F_q^×/(F_q^×)²`, which has two elements; scaling by
  `α ∈ F_q^×` multiplies `disc` by `α²` and so preserves the class. Hence there
  are exactly **two** similarity classes — hyperbolic and anisotropic — and
  exactly **one** anisotropic class, the norm form of `F_{q²}`. `C_aff(A → A) = 1`
  trivially.
* In three or more variables, Chevalley–Warning makes every form over a finite
  field isotropic, so there are **no** anisotropic classes at all.

So over `F_q` the question "what is `C_aff` between two anisotropic classes?"
has no instances. It should be recorded as empty rather than approximated.

---

## 6. What the `p`-adic version would need

*(exploratory; no new numbers claimed)*

Over `Q_p` the question is non-empty: `p_adic_quadratic_map_poset.md` records
eight square classes over `Q_2`, of which seven are anisotropic, and
`p_adic_exchange_rate_attempts.md` brackets both directions of the first pair
`(q_1, q_2) = (x²+y², x²+2y²)` in `[1/2, 1]`.

**What transports unchanged.** Lemmas 0–3 and Propositions A–C are stated and
proved above over an arbitrary field of characteristic `≠ 2`. Nothing in them
needs a residue field, so the brief's warning that "the reduction lemma of
`affine_rate` §0 does not transplant unchanged" is, for *these* statements, too
cautious — the reduction lemma is pure linear algebra on affine processors and
holds verbatim over `Q_p`. In particular, **new over `Q_p`**:

```
C_aff(x² + d y² → x) = 2/3   for every anisotropic class [d] over Q_p,
C_aff(x² − y²  → x) = 1.
```

That puts a number on the missing Hasse arrows of the `Q_2` poset: the seven
anisotropic classes do not implement a nonconstant linear map in one step, but
they implement it at rate exactly `2/3`, not `0`.

**What does not transport, and why.** Every bound proved in §3 that is strictly
below `1` goes through Propositions A and C together: A needs the **resource**
anisotropic, C needs the **target** affine along some line. For the `Q_p`
anisotropic pair `C_aff(q_d → q_e)` the resource is fine but the target is not:
**an anisotropic target has no such line, over any field.** For the two `F_3`
brackets of §3.4 it is the other way round: the target `x²+y` is fine but the
resource `xy` is isotropic, so Proposition A is unavailable. The two failures
land on the same combinatorial question, and the finite-field instance is small
enough to search exhaustively while the `p`-adic one is not. Concretely, both
ask:

> Given a class `D` of quadratic forms closed under pullback (the `q_d`-forms),
> what is the least dimension of a `D`-spanned subspace of `Sym²((K^{2k})*)`
> containing `k` copies of a form `q_e` on disjoint two-dimensional blocks?

with the trivial bounds `k+1 ≤ N_k ≤ 2k` from §3.4 and Attempt 1 respectively.
**Beating `k+1` asymptotically in either setting would settle the other's
shape.** That, and not a numerical experiment, is what the `p`-adic version
needs.

**One lever that exists over `Q_p` and not over `F_q`.** The rank-one atoms of a
form `Q` are `a·ℓ²` with `a` a value represented by `Q`. Over `F_q` the norm map
`F_{q²}^× → F_q^×` is **surjective**, so *every* rank-one form is an atom of
*every* anisotropic `Q` — no information. Over `Q_p` the represented values of
`x² + d y²` form the norm group of `Q_p(√−d)`, of **index 2** in `Q_p^×`. So
over `Q_p` the atom set is genuinely smaller and carries a `Z/2` invariant that
has no finite-field shadow. Any `p`-adic-specific improvement should start
there; the modulo-`4` search of Attempt 2 and the Schur-complement computation
of Attempt 3 are both special cases of using that invariant, and neither is
asymptotic. **This session claims no improvement on `[1/2, 1]` and does not
attempt one.**

---

## 7. Corrections

1. **`research/affine_rate/FINDINGS.md` is missing from the repository.** The
   brief instructs "read it in full"; it was never committed. Commit `d4ffdf1`
   ("Brief O: first data on the affine rate `C_aff` …") states in its message
   that "the harness refused the `research/affine_rate/FINDINGS.md` write" and
   that the findings live in a session report. This document therefore
   reconstructs the reduction lemma and the support monotone from
   `affine_atoms.py` (whose docstring states the lemma) and from
   `run_n2_small.py` / `run_summary.py`, and re-derives Propositions 1–2 rather
   than citing them. **Recommendation: commit the brief-O body.** Everything in
   it that could be checked here checked out — the `N_1`/`N_2` matrices, the
   `2/3` for `x²+y² → x`, the `1/2` values for the `x²` row, `N_2 = 3 < 4 = 2N_1`,
   and the `C_sig` matrix.
2. **`C_aff ≤ C_sig` remains refuted** (`research/synergy/FINDINGS.md` §4). It is
   not used anywhere above. The four pairs on this pool with `C_aff > C_sig` all
   have `N_1 = 1`; §3.5 gives the current, restricted, form of the conjecture.
3. **The brief's `C_sig` list is a mixed family.** `0.975301, 0.993132,
   0.988581, 0.992632, 0.998289, 0.998811, 0.995734` at `q = 3,…,19` is
   `C_sig(x²+y² → x)`, which switches signature at `q ≡ 1 mod 4`. It is not
   monotone for that reason, and the non-monotonicity is not evidence about a
   jump. The uniformly anisotropic family (§2) is monotone increasing to `1`.
4. **Proposition A is not tight for every target.** `N_2(x²+y² → x²+y) = 4`
   while `ceil(3·2/2) = 3`. Do not read `ceil(3k/2)` as a formula for `N_k`
   beyond the two targets where it is proved (§4).

---

## 8. Open

1. **`C_aff(xy → x²+y²)` and `C_aff(xy → x²+y)`**, both in `[1/2, 1]`. For the
   first, `N_1, N_2, N_3 = 2, 4, 6` is now exhaustively known, so the lower
   bound `1/2` is not going to improve at small `k`; an upper bound below `1`
   needs a lower bound on `N_k` growing faster than `k`, and no mechanism for
   that is known here — the two available mechanisms both fail (Proposition A
   needs an anisotropic resource, Proposition C needs a target that is affine on
   a line). This is the crux of `affine_rate`'s open item 2. `N_4 ≤ 6` for the
   first pair would force `C_aff ≥ 2/3` and is decidable by
   `run_split_to_sum2.py` if its `act` step is vectorised (the `k = 4` product
   set has about `10^7` elements).
2. **`C_aff(x²+y² → x²+y)` exactly**, in `[1/2, 2/3]`. `N_3 ∈ {5,6}` decides
   `≥ 3/5` versus consistency with `N_k = 2k`; the exhaustive `N_k ≤ k+2`
   decision needs a symmetry-reduced pair search that was not implemented.
3. **Rationality of `C_aff` in general** (§4).
4. **`C_aff` between anisotropic classes over `Q_p`**, `[1/2, 1]` (§6),
   equivalent in shape to item 1.
5. **Characteristic 2.** Lemma 1 needs `|K| ≥ 3` and Lemmas 2–3 need
   `char K ≠ 2`. `F_2`, `F_4`, `F_8` are untouched; over `F_2` the jet map is not
   injective (`x² = x` as functions) and the whole calculus needs restating.
6. **The `⌈3k/2⌉` phenomenon.** Two different targets (`x` and `xy`) give
   literally the same `N_k` for every anisotropic resource over every odd `F_q`.
   Whether `x` and `xy` are `C_aff`-equivalent as targets in general is not
   known.

---

## 9. Files

| file | what it is |
|---|---|
| `isotropy_atoms.py` | general-`q` jet calculus: atom generation, `F_q` linear algebra, the exhaustive `N_k` searches and the `N_k ≤ k+1` decision |
| `run_isotropy_check.py` | validates the jet model against the committed `F_3` value-table results (all 50 entries agree) |
| `run_isotropy_family.py` | `N_1`, `N_2 ≤ 3` and re-verified witnesses for `x²−ny²` and `x²+y²` at `q = 3,5,7,11,13` → `isotropy_family.csv` |
| `run_isotropy_witness.py` | the Theorem 2 and Theorem 3 witnesses re-verified on every point of `F_q^2` and `F_q^4`, `q = 3,…,19` |
| `run_isotropy_csig.py` | `C_sig(Q → x)` at 40 digits, endpoints against closed forms, and the amplification table → `isotropy_csig.csv`, `isotropy_csig_gap.csv` |
| `run_n3.py` | exhaustive `N_3 ≤ 4` decision for all 25 `F_3` quadratic pairs → `n3_le4.json` |
| `run_isotropy_summary.py` | the §3.5 bracket table, regenerated from `n1_matrix.json` and `csig_matrix.json` with the proved upper bounds folded in |
| `run_diagonal_ansatz.py` | brute-force check of the two diagonal-ansatz minima of §3.4 (`2k` and `ceil(3k/2)`) at `k = 1,2,3` |
| `run_split_to_sum2.py` | symmetry-reduced exhaustive `N_k(xy → x²+y²) ≤ k+2` decision in `Sym²`; gives `N_2 = 4` and `N_3 = 6` |
| `isotropy_family.csv` | the `q`-by-`q` congruence table of §1.4 |
| `isotropy_csig.csv` | 40-digit `C_sig(Q → x)` for the three families at `q ≤ 47` |
| `isotropy_csig_gap.csv` | the `C_sig` gap and the amplification factor of §2 |
| `n3_le4.json` | the `k = 3` decision table |

Pre-existing files this builds on: `affine_atoms.py`, `run_n1.py`,
`run_n2_small.py`, `run_witness.py`, `run_summary.py`, `n1_matrix.json`,
`n2_le3.json`, `csig_matrix.json`.

Reproduce with

```bash
python3 research/affine_rate/run_isotropy_check.py
python3 research/affine_rate/run_isotropy_family.py
python3 research/affine_rate/run_isotropy_witness.py
python3 research/affine_rate/run_isotropy_csig.py
python3 research/affine_rate/run_isotropy_summary.py
python3 research/affine_rate/run_diagonal_ansatz.py
python3 research/affine_rate/run_split_to_sum2.py
python3 research/affine_rate/run_n3.py
```

`run_n3.py` is the only slow one (about two hours; the `xy` atom set on
`F_3^6` has `1 192 101` elements).
