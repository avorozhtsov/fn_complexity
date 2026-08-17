# Findings — M ↔ E, and what M knows about {a_c}

Research record for two questions:

1. the connection between the directed exchange matrix `M_ab = C(a→b)` and the
   Weil matrix `E` whose positivity is equivalent to RH;
2. what information about Frobenius traces `{a_c}` is carried by `M`.

Six threads, run in parallel and revised as results arrived. Everything below is
either proved, or computed and independently re-verified; the distinction is
marked in each case. Conventions and the thread list are in `PLAN.md`.

---

## Summary

**Track 1 is answered, and the answer is negative but sharp.** `M` and `E` are
built from the same data — a signature *is* a Weil test measure — and the
pairing has a closed form. But no functional relation between the two induced
geometries can exist, because each is invariant under a group the other is not,
and one of those invariances is itself a manifestation of the critical line.
Separately, the geometry of `M` is now pinned exactly in the Deza–Laurent
hierarchy: it lies in `MET \ HYP`, with a minimal negative-type violation on
exactly five points — and that violation is quantitatively mild, the
`l2`-distortion being `1.3375` on the worst family found and about `1.1` on
random ones, far below the `O(log n)` a general metric may need.

**Track 2 is answered positively.** The two directions of the rate against a
linear map read different things: `C(L→f)` sees the largest fiber and nothing
else, exactly; `C(f→L)` has a universal bottleneck at `β* = √2 − 1` and its
interior carries the entire moment ladder of the traces. Flatness of a signature
is equivalent to a permutation-polynomial condition, which subsumes the known
congruence result and reaches non-congruence conditions such as supersingular
primes. Symmetry type is invisible, and provably so.

Five claims of mine were corrected in the process; they are listed in
[Corrections](#corrections-to-earlier-claims).


---

## Notation

Used throughout, and worth fixing because two of the results are statements
about which of these quantities the rate can see.

**`A² → A¹`** — affine spaces: `A^n` over `F_q` is just `F_q^n`, so
`f : A² → A¹` is a polynomial `f(x,y)` with values in `F_q`. The notation is
geometric rather than set-theoretic because what matters is the **fibers**
`f⁻¹(c)`, which are plane curves with a genus and a point count.
"Geometrically irreducible fibers" is a real hypothesis: a fiber with `r`
components carries about `r·q` points, which swamps the Weil-scale corrections
(`1 − C ≍ log r/log q` instead of `≍ 1/(√q log q)`).

**`σ(f)`** — the fiber signature, the multiset `{N_c}` of fiber cardinalities,
`N_c = #f⁻¹(c)`. Always `Σ_c N_c = q²`.

**`L`** — the flat signature `(q, q, …, q)` with `q` entries: the signature of a
linear map `f(x,y) = x`, whose fibers are `q` lines of `q` points each. It is the
reference resource — same fiber count and same total as any `f : A² → A¹`, but
perfectly evenly distributed.

**`C(g→f)`** — the implementer is written **first**: the number of copies of the
target `f` obtainable per copy of the resource `g`, `inf_β log Z_g/log Z_f`. So
`C(L→f)` and `C(f→L)` read different things, and that asymmetry is the content
of T2.1 and T2.2.

**`a_c`** — the trace of Frobenius of the fiber over `c`:

```
a_c = q − N_c
```

the deviation of the fiber's point count from the generic value `q`. The name is
Weil's: a smooth projective curve of genus `g` over `F_q` has `q + 1 − a` points
with `|a| ≤ 2g√q`. Since there are `q` fibers carrying `q²` points in total,

```
Σ_c a_c = q·q − q² = 0     identically, for every f
```

which is why the first moment drops out of the expansion and the linear-in-β
term vanishes.

**`m₂`** — the normalised second moment:

```
m₂ = (1/q²) Σ_c a_c²
```

The normalisation makes it `O(1)`: Weil gives `|a_c| ≲ 2g√q`, so `a_c² ≲ 4g²q`
and the sum over `q` fibers is `≲ 4g²q²`. For a family with large monodromy the
`a_c/√q` equidistribute by Sato–Tate, `E[a_c²] ≈ q`, hence **`m₂ → 1`**; the CM
family `y² = x³ + c` gives `m₂ → 2` at `q ≡ 1 mod 3` and `m₂ = 0` at
`q ≡ 2 mod 3`, where every fiber is supersingular.

Equivalently `Σ_c a_c² = Σ_c N_c² − q³`, that is `m₂ = Z_f(2)/q² − q`, where
`Z_f(2) = Σ_c N_c²` counts the points of the **fiber square** `X ×_Y X`. So `m₂`
is not an auxiliary statistic but the partition function itself at the integer
point `β = 2`.

*Worked example*, `f = y² − x³ − x` at `q = 11`:

```
N_c = [17, 15, 14, 13, 12, 11, 10, 9, 8, 7, 5]     sum = 121 = q²
a_c = [−6, −4, −3, −2, −1,  0,  1, 2, 3, 4, 6]     sum = 0
Σ a_c² = 132        m₂ = 132/121 = 1.0909
```

| q | m₂ | `1 − C(f→L)` | `(3−2√2)·m₂/(2q log q)` |
|---:|---:|---:|---:|
| 11 | 1.090909 | 3.869·10⁻³ | 3.548·10⁻³ |
| 1009 | 0.997027 | 1.227·10⁻⁵ | 1.226·10⁻⁵ |

At `q = 11` the discrepancy is ~9%, since the dropped term is of order
`q^{−3/2} ≈ 0.027`; by `q = 1009` the two agree to three significant figures.

---

## Track 1 — the M ↔ E connection

### T1.4 The Weil pairing evaluated on signatures *(theorem + numerics)*

A signature `a = (a_1,…,a_r)` **is** a Weil test measure: the atomic measure
`Σ_i δ_{a_i}` has Mellin transform `Z_a(s) = Σ_i a_i^s`, exactly the partition
function. So both matrices live on the same objects,

```
M_ab = C(a→b) = inf_β log Z_a(β)/log Z_b(β)        (real spectrum)
E_ab = Σ_ρ Z_a(ρ) conj(Z_b(ρ))                     (complex spectrum)
```

and `E`, truncated to the first `N` zeros, is a Gram matrix, hence PSD by
construction (min eigenvalue `−1.7e−11`).

Writing `ρ = ½ + iγ` gives `Z_a(ρ) = Σ_i √a_i · e^{iγ log a_i}`, so
`E_ab = Σ_{i,j} √(a_i b_j) · S(a_i/b_j)` with `S(x) = Σ_{0<γ≤T} x^{iγ}`.
Landau's theorem evaluates `S`, and since `√(a_i b_j)/√y = min(a_i,b_j)`:

```
E_ab = N·O(a,b) − (T/2π)·A(a,b) + O(rs log T)
O(a,b) = Σ_{a_i = b_j} a_i                             multiset overlap
A(a,b) = Σ_{a_i ≠ b_j} min(a_i,b_j)·Λ(max/min)         prime-power ratios
```

Confirmed on 1200 zeros over 136 pairs: **max error 49.9 against a predicted
Landau remainder of ~67**, mean 15.5, diagonal scale ~24000.

**Consequence.** `E` is driven entirely by *multiplicative coincidences* —
equality of entries, and ratios that are exactly prime powers — and vanishes for
a generic pair. `M` reads the same entries through their logarithms and is
finite for every pair. Weil-orthogonality is generic: `(4,2)` and `(5,3)` are
among the closest pairs in the exchange metric (`d = 0.163`) and exactly
orthogonal in the Weil geometry (`E = 3.14` against a scale of 9600). Measured
correlation between exchange distance and Weil angle: **+0.19**.

The appearance of `Λ` is the explicit formula in its smallest instance, one
delta per fiber.

*Caveat:* atomic measures are not admissible Weil test functions, so this is the
finite-rank truncation; `|Z_a(½+iγ)|` does not decay and the full sum does not
converge. Smoothing would restore admissibility at the cost of blurring `Λ`.

### T1.5 No functional relation is possible *(two exact theorems)*

Designed families *do* reach high correlation — **+0.86** on geometric-ladder
subsets, **−0.98** on the receding staircase `a_t = (p^t, 1)` — so the `+0.19`
above is an artefact of generic signatures. **But multiplicative design is not
what buys it.** The control ladder `{3,7,17,37,67,131,257,521}` has the same
log-spacing but no prime-power ratio, so `A ≡ 0` identically and its Weil matrix
is exactly `N·O` with no zeta content; it matches the geometric ladders cell for
cell at every subset size and every `N`. The climb is bought entirely by shared
entries. Where `Λ` does dominate (`O = 0`) the correlation flips sign.

**Theorem A (Weil angles are scale-invariant).** `Z_{λa}(ρ) = λ^ρ Z_a(ρ)`, so
`E_{λa,λb} = λ^{2 Re ρ} E_{a,b}`. Verified to `6e−14` while mean `d` falls
0.393 → 0.151 as λ: 1 → 210.

> **This step uses RH.** The factor is constant across zeros *only because every
> `Re ρ = ½`*. Counterfactual: moving 60 of 1200 zeros to `Re = 0.7` breaks the
> invariance from `6e−14` to **`1.3e−1`**.

**Theorem B (d is Cartesian-power invariant).** `log Z_{a^⊗k} = k log Z_a` gives
`C(a^⊗k→b) = k C(a→b)` and `C(b→a^⊗k) = C(b→a)/k`, so the product — hence `d` —
is unchanged. Verified to `2.5e−16`; immediate from `u_{a^⊗k} = u_a + log k`
being a constant shift while `d` is an oscillation. The angles meanwhile move
0.224 rad on average, 0.821 at worst.

**Conclusion.** Transverse invariance groups ⇒ neither geometry is a function of
the other, and any measured correlation is a property of the family chosen. The
staircase makes it concrete: `d(a_s,a_t) = log(t/s)` exactly while
`angle − π/2 = K·p^{−(t−s)/2}` — **d sees the ratio `t/s`, the angle sees the
difference `t−s`.**

### T1.1 The minimal negative-type violation has exactly 5 points *(proved)*

```
a₁=(12,10,8,8,2,1)  a₂=(11,9,7,7,4,1)  a₃=(12,12,6,5,4,4)
a₄=(12,10,7,4,3,3)  a₅=(11,11,7,7,4,3)
x = (302626, −510642, −576418, 330027, 454407)/10⁶      Σx = 0 exactly
xᵀDx = +9.8126948851·10⁻⁴
```

Verified four ways: package solver, dense `2e6`-point β-grid on `[0,600]`
(`1.5e−9`), the same on `[0,2000]`, and 40-digit mpmath (`3.8e−16`). Margin
`~10⁻³`, six orders above solver tolerance. Triangle slack `+4.87e−5`, so it is
a genuine metric; the twenty infima are attained at `β = 0` (4 pairs), interior
β up to 12.73, and `β = ∞` (3 pairs) — the violation uses the whole temperature
range.

**Minimality is proved, not searched.** Every cut semimetric satisfies
`xᵀδ_S x = −2(Σ_{i∈S} x_i)² ≤ 0`, so `CUT_n ⊆ NEG_n`; and `MET₄ = CUT₄`, so no
four-point metric can violate negative type at all.

**Position in the hierarchy, exact.** Triangle holds (min slack `−4.4e−16` over
`8e9` triples). The **pentagonal** inequality `b = (1,1,1,−1,−1)` fails by
`+5.39e−2`. Since `CUT₅ = HYP₅`, ℓ₁-embeddability also breaks at exactly five
points. **The exchange metric lies in `MET \ HYP` and leaves the hierarchy at
the first opportunity.**

*Methodological trap, worth propagating:* `−½JDJ` always has the constant vector
in its kernel, so its smallest eigenvalue is capped at 0 and gives no search
gradient — an earlier run reported "no 5-point violation" purely because of
that. One must work in an orthonormal basis of `{Σx = 0}`.

### T1.2 The PSD ray *(mixed: one proof, one conjecture)*

`S(d) = {t : exp(−t·d) ⪰ 0}` was a **single closed ray `[t*, ∞)`** in every case
tested — 26 signature families, 55 bipartite metrics, 6805 random
non-negative-type metrics including an adversarial two-scale ensemble built to
force a gap. One sign change of `λ_min`, never more.

**Proved:** `S` is closed, additive (Schur product theorem) and contains a ray,
hence has finitely many gaps. The interval claim itself remains a conjecture.
**Proved exactly:** `t*(K_{m,n}) = ½ log((m−1)(n−1))`, so `t*` is **unbounded at
fixed diameter** and no bound `t* ≤ f(diam)` exists.

**Positivity is worst in the middle.** Via `D_t = t·κ₂ + t³·κ₄/12 + …` (odd
cumulants cancel), `D_t/t → Var_μ(u_i − u_j)` at the Hilbert end and `D_t → d`
at the tropical end; the threshold peaks `4×`–`199×` above `t*(d)` near
`t ≈ 15–60`. The Laplace limit `P_t → 1/C` was confirmed to `1e−5` at `t = 1e6`.

**Caveat with consequences.** *Every bounded spectral window destroyed the
13-signature violation*, including `[0.001, 200]` — truncation error `3.5e−3`
against a violation of `5.45e−4`; it reappears only at `β_max ≳ 500`. Any
windowed computation on these objects needs `β ~ 10³`.


### T1.3 The failure of negative type is quantitatively mild *(computed, bracketed)*

Since `d` is not of negative type it does not embed isometrically in Hilbert
space. The `l2`-distortion measures how badly:

```
c2(d) = min over f : X -> l2 of (max stretch) x (max shrink)
```

and Bourgain's theorem gives `c2 = O(log n)` for every n-point metric, attained
on expanders. Computing `c2` exactly is an SDP; with no solver available it was
bracketed from both sides.

* **Upper bound** — any embedding gives one. Parametrise `G = X Xᵀ`, so
  `Q_ij = ||x_i − x_j||²`, and minimise the spread
  `max_ij(Q_ij/d_ij²) / min_ij(Q_ij/d_ij²)`.
* **Lower bound** — a Poincaré certificate. For symmetric `Δ` with zero row sums
  and `Δ ⪰ 0`, and any `f`, `Σ_ij Δ_ij ||f_i − f_j||² = −2 tr(ΔF) ≤ 0` because
  the Gram matrix `F` is PSD. Splitting `Δ` by sign off the diagonal and using
  `d_ij² ≤ ||f_i−f_j||² ≤ D² d_ij²` gives
  `c2² ≥ Σ_{Δ>0} Δ d² / Σ_{Δ<0} |Δ| d²`. Writing `Δ = J Y Yᵀ J` makes both
  constraints automatic.

Solver check: `c2(C₄) = √2 = 1.414214` is reproduced from both sides.

| family | n | c2 lower | c2 upper | c2 / log n |
|---|---:|---:|---:|---:|
| **minimal certificate** | 5 | **1.337503** | **1.337774** | 0.83 |
| random | 5 | 1.008199 | 1.008385 | 0.63 |
| random | 8 | 1.021893 | 1.025293 | 0.49 |
| random | 12 | 1.114059 | 1.129574 | 0.45 |
| random | 16 | 1.102484 | 1.134143 | 0.41 |
| random | 20 | 1.052543 | 1.126801 | 0.38 |
| random | 25 | 1.070405 | 1.143308 | 0.36 |

**The metric is almost Hilbertian.** On the minimal five-signature witness
`c2 = 1.3375`, bracketed to `3·10⁻⁴`; on random families `c2 ≈ 1.1` and it does
**not grow** — `c2/log n` falls monotonically from 0.63 to 0.36, so Bourgain's
bound is nowhere near attained (at `n = 25` it would allow ≈ 3.2 against an
observed 1.14). Greedily extending the certificate improved nothing over 60
candidate additions, so the largest distortion seen anywhere belongs to the
*smallest possible* witness.

This agrees with T1.2 from the other direction: the certificate's PSD threshold
is only `t* = 0.124`, while artificially hill-climbed metrics reach `t* ≈ 12–17`.
By both measures the exchange metric is a weak violator.

*Open:* whether `c2` is bounded uniformly in `n` or grows very slowly. A greedy
search is weak evidence and `n ≤ 25` is small, though the downward trend in
`c2/log n` is consistent.

Reproduce: `analysis/l2_distortion.py`.

---

## Track 2 — what M knows about {a_c}

### T2.1 The rate against L is an endpoint, and the convergence exponent is the invariant *(proved + numerics)*

**Proposition.** For every `f : A² → A¹` over `F_q`,

```
C(L→f) = log q / log(max_c N_c),    attained at β = ∞.
```

*Proof.* `Z_f(β) = Σ N_c^β ≤ q(max N)^β`, so with `A = log q`, `B = log max N`,
`R(β) ≥ (1+β)A/(A+βB) > A/B` for finite β when `B > A`, and `R(∞) = A/B`. ∎

Verified to ten decimals for genus 1, 2, 3. **So this rate sees the largest
fiber and nothing else.**

Hence `(1−C)√q log q = μ − μ²(½ + 1/log q)/√q + O(1/q)` with
`μ = max_c(−a_c)/√q`, and the `2g` law is entirely a question of extreme-value
statistics. The `USp(2g)` lower tail `P(2g − T < ε) ~ K_g ε^{d/2}` with
`d = dim USp(2g) = 2g² + g` gives `2g − μ ≈ Γ(1+2/d)(qK_g)^{−2/d}`:

| g | fitted exponent | predicted |
|---|---:|---:|
| 2 | −0.1969 | −0.2000 |
| 3 | −0.0931 | −0.0952 |
| 4 | −0.0583 | −0.0556 |

> **The limit is `2g`; the *rate of approach* encodes `dim USp(2g)` — a finer
> invariant than the genus, read off the same curve.**

At `q = 10⁶`, `μ` = 3.684, 4.663, 4.760 against targets 4, 6, 8; reaching within
10% of `2g` needs `q ≈ 6e4`, `6e9`, `4e16`. For `g ≥ 3` at reachable `q` the
observed maximum is a **Gaussian** extreme, not a Weil edge, with crossover at
`log q ≈ 2g²`.

### T2.2 A universal bottleneck temperature *(derived analytically, confirmed)*

Since `m₁ = Σ_c a_c/q = 0` identically,
`1 − C(f→L) = (m₂/(2q log q))·min_β β(β−1)/(β+1) + O(q^{−3/2})`, and
`d/dβ[β(β−1)/(β+1)] = (β²+2β−1)/(β+1)²`, so the minimiser is the positive root
of `β² + 2β − 1 = 0`:

```
β* = √2 − 1 = 0.414213562…        1 − C(f→L) = (3−2√2)·m₂/(2q log q)
```

**Independent of family and genus.** At `q = 1009` the measured `β*` is
0.414089, 0.414127, 0.413946 for g = 1, 2, 3, and `(1−C)·2q log q/m₂` is
0.171802, 0.171671, 0.172101 against `3−2√2 = 0.1715729`.

**The interior carries the whole moment ladder**, each order damped by ≈`0.6/√q`.
Subtracting moments 2..K drops the residual rms `2.75e−5 → 8.5e−8 → 5.7e−8 →
3.7e−9`, then plateaus at exactly the first dropped term (predicted `3.30e−9`,
observed `3.36e−9`). `m₂` is recoverable from one rate to relative `O(q^{−1/2})`.

**Sharp separation.** Two maps at `q = 211` with identical image size, identical
largest fiber and identical `m₂`, smallest fibers 173 vs 167: all endpoint
probes agree to `0.000e+00`, interior probes separate by `2.42e−7`, predicted in
advance by `c₃Δm₃ + c₄Δm₄ = 2.32e−7` (4%). **The interior is not redundant.**

**What is invisible, and the price of a convention.** The *smallest* fiber
(`max_c a_c`) is not determined — R² = 0.80 against `1.0000000` for the largest
— because isolating it needs `β < 0`, precisely the region the first paper
excludes since `Z_a` is order-preserving only for `β ≥ 0`. **The exchange rate
sees the largest fiber exactly and the smallest not at all.**

**Structural collisions.** 400 random elliptic fibrations give only 5/3/3
distinct signatures at `q = 101/211/503` — exactly `gcd(4,q−1)+1`, the quartic
twist classes. Genus ≥ 2: 398–400 of 400 distinct. And `m₂ = K_P/q − 1` with
`K_P = #{(x,x'): P(x)=P(x')}`, the fiber-square count.

### T2.3 Flatness is a permutation-polynomial condition *(theorem)*

**Theorem.** For `f = y² − P(x)` over `F_q`, the signature is flat — hence
`C(f→L) = C(L→f) = 1` exactly — **iff `P` is a permutation polynomial of `F_q`**.
Two lines from `Σ_c a_c² = q(ν(P) − q)`, `ν(P) = #{(x,x′) : P(x) = P(x′)}`:
`Σ n_u = Σ n_u² = q` forces `n_u ≡ 1`. Verified on 2093 `(P,q)` pairs, zero
violations; independently re-verified across five polynomials and eight primes.

The known `q ≡ 2 mod 3` result is exactly the case `P = x³`, since `x³` permutes
`F_q` iff `gcd(3,q−1) = 1`. Extensions: `y^r − x^d` flat iff
`min(gcd(r,q−1), gcd(d,q−1)) = 1` (4368 triples, 0 violations); **Dickson
polynomials `y² − D_n(x)` flat iff `gcd(n, q²−1) = 1`** — a congruence involving
`q+1`, which no monomial family can produce.

**Beyond congruences.** The quadratic-twist family `f = P(x)y²` equals the split
conic iff `a_E = 0`; for the non-CM curve `y² = x³+x+1` that set is
`{17, 179, 227}` below 500 — the **supersingular primes**, provably not a union
of residue classes.

**Symmetry type is invisible, by construction.** `Σ_c a_c = 0` holds identically
for every fibration of `A²`, and that first moment is exactly the Katz–Sarnak
statistic separating orthogonal from symplectic. What is visible beyond genus is
`m₂ = ν(P)/q − 1` = monodromy rank of the branch map minus one: cyclic `e−1`,
dihedral `≈⌊n/2⌋`, 2-transitive `1`, at every genus. Monodromy at fixed genus
separates clearly (`ΔC(f→L) = 2.2e−5` against a `1e−13` floor); symmetry type at
a single `q` does not, reappearing only in the fluctuation over `q`.

**No exotic exact values.** A scan of 69426 ordered rates finds only `1`, `1/2`
and `0` persisting. The arithmetic content is carried entirely by *which
signatures coincide*.

---

## Corrections to earlier claims

The research overturned five things asserted before it ran. All have been
propagated into the paper notes.

1. **The 13-signature negative-type certificate was a search artefact.** The
   minimum is exactly 5, and minimality is provable (T1.1). Theorem 2 in
   `exchange_positivity_and_weil.md` now carries the minimal certificate.
2. **The M ↔ E "dictionary" is a proven non-relation**, not a hedge. Two
   transverse invariances rule out any functional relation (T1.5).
3. **`1 − C(L→f) ≍ 2g/(√q log q)` needs geometrically irreducible fibers.** A
   fiber with `r` components gives `1 − C ≍ log r/log q`, swamping the Weil
   scale. The correct genus-0 control is `x²+y²` at `q ≡ 3 mod 4`, not the split
   conic (T2.1).
4. **`max_c(−a_c) ≠ max_c |a_c|`** at finite `q` — they agree only 50–68% of the
   time for g = 2,3,4. The rate tracks the former only (T2.1).
5. **The CM-by-`Z[i]` prediction was false as stated.** The quartic-twist family
   does not go flat at `q ≡ 3 mod 4`; it becomes exactly the split conic, and at
   `q ≡ 1 mod 4` its signature records the factorisation `a² + b² = 4q` in
   `Z[i]` rather than merely the congruence (T2.3).

---

## Open

* Whether `S(d)` is always an interval — proved only to have finitely many gaps.
* Smoothed test measures, to make the Weil side admissible rather than
  finite-rank truncated, and to see what survives of the `Λ` term.
* Whether the convergence exponent of T2.1 gives a *practical* symmetry-type
  detector, given that T2.3 shows the first moment is unavailable.

## Files

Per-thread notes `T1_1`, `T1_2`, `T1_4`, `T1_5`, `T2_1`, `T2_2`, `T2_3`, with
scripts and CSVs alongside. Zeta zeros are cached in `zeta_zeros_1200.npy` and
`zeta_zeros_2400.npy` (they take minutes to regenerate). The 48 MB search cache
under `cache/` is gitignored and regenerated on demand.
