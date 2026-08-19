# Findings — the genus-four witnesses, and a wider cycle that needs none

Answer to session brief J. Three results, one of which was not asked for and is
the largest.

1. **`SU(2) × 3·SU(2)` at genus four is arithmetic.** The route
   `REPEATED_FACTOR.md` wrote out works, and the locus is not merely a curve but
   a **rational** one, parametrised by a single `n`:

   ```
   B₃ = {i, −i, n, −n, 1/n, −1/n},   p = n + i − 1/n
   f₁ = (x+i)(x+n)(x−1/n)(x−p)       f₂ = (x−i)(x−n)(x+1/n)(x−p)
   f₁f₂ = (x−p)² (x²+1)(x²−n²)(x²−n⁻²)
   ```

   `C_n` = the `(ℤ/2)²`-cover `y₁² = f₁`, `y₂² = f₂`; genus 4, non-hyperelliptic,
   `Jac(C_n) ∼ J₁ × E³` with `E` varying. `J₂ ≅ E` **over the base field** — the
   pullback constant is `(n²+1)⁴`, a perfect square, so there is no quadratic
   twist to fight. Verified: **`a(C_n) = a(J₁) + 3a(E)` on all 50548 fibres of
   two pencils over 32 primes `q ≡ 1 (mod 4)`, zero mismatches**, with
   `(m₂,m₄,m₆) = (10.059, 218.30, 6299)` at `q = 16001` against the
   `SU(2)×3·SU(2)` values `(10, 218, 6350)`. With a second new pencil for
   `SU(2)×USp(4)×USp(4)` at genus five, **cycles 8 and 9 of the nine become
   theorems** for `q ≡ 1 (mod 4)`.

2. **`4·SU(2)` is zero-dimensional on this route — a negative result.** Every
   one of the 10 splittings `S ⊔ T` and 6 `S₃`-branches gives `p = p(n)`
   rational, and on each of the 60 the extra condition `j(J₁) = j(E)` is a
   non-zero polynomial in `n` of degree ≤ 28. Counted over `F_q`, the number of
   solutions stays in `[0, 18]` while `q` runs over a **100-fold** range
   `401 … 40009`: a finite set, never a pencil. Two structural facts explain
   it: a `(ℤ/2)²`-cover of `P¹` of genus 4 has quotient genera `(1,1,2)` or
   `(0,2,2)` **and nothing else**, and **no `(ℤ/2)³`-cover of `P¹` has genus 4
   at all** — so the mechanism that produces `3·SU(2)` at genus 3, four
   elliptic quotients cyclically permuted, is unavailable. The six widest
   cycles stay unwitnessed. This is evidence, **not a proof**, of
   non-existence.

3. **The reason the six widest cycles were the target is gone: the cutoff that
   made them so was an artefact.** `FINDINGS.md` searched the symplectic cone
   only to `α_max ≤ 12` and concluded that *every* limiting 3-cycle needs a
   repeated isogeny factor, the multiplicity-free sub-cone being transitive
   (7308 oriented triangles, 0 cycles). Pushing the identical search to
   `α_max ≤ 14` breaks that:

   ```
   SU(2)⁵   ≺   USp(14)   ≺   USp(6) × USp(6)   ≺   SU(2)⁵
    genus 5      genus 7        genus 6
   ```

   is a strict 3-cycle of **multiplicity-free** measures with margin
   `4.0116·10⁻²`, and all three vertices have explicit pencils, over **every**
   odd `q` past a small bound, with no congruence condition. `SU(2)⁵` is the
   `(ℤ/2)³`-cover of `P¹` branched at `{0, ∞, 1, 2, t, t+1}`, whose five
   elliptic quotients include the Legendre pencil itself. A second cycle at the same `α_max`,
   needing only the multiplicity two that Theorem A already supplies, is
   **five times wider than anything witnessed before**:

   ```
   SU(2) × 2·SU(2) × USp(4)   ≺   USp(14)   ≺   USp(6) × USp(6)
        genus 5                     genus 7        genus 6
   ```

   margin `6.1583·10⁻²`. So `FINDINGS.md`'s reduction of the whole `q → ∞`
   question to "one arithmetic input, a Jacobian isogenous to `A^k`" is
   **withdrawn**: the limiting comparison is non-transitive on families with no
   repeated isogeny factor whatever.

Everything below is either proved, or computed and independently re-verified;
each statement says which. One claim of `FINDINGS.md` is corrected in
[Corrections](#corrections).

**The standing caveat, unchanged.** These are theorems about the `q → ∞`
**limit**. `FINDINGS.md` measures `sup_τ|Ψ_f − Ψ_μ| = 0.28` at genus 2 and
`1.52` at genus 3 at `q = 4·10⁵`, and extrapolates `q ~ 10²⁰` before a genus-6
pencil is within `10⁻²` of its limit. The widest margin here is `6.2·10⁻²`,
five times the old one, so the implied `q₀` is much smaller — but it is still
astronomically beyond any census, and none of this predicts anything about one.

---

## Notation

As `FINDINGS.md` and `REPEATED_FACTOR.md`. `a_c = q − N_c`, `α_c = −a_c/√q`,
`m_j = E[α^j]` in the limit; `Ψ_μ(τ) = K_μ(τ)/τ`, `μ ≺ ν` iff
`mid_τ(Ψ_μ − Ψ_ν) < 0`, `mid = ½(sup+inf)` over `[0,∞]`. `k·G` is `k` isogenous
copies of a block: `α = k·tr(g)`.

Moment targets used as detectors:

| measure | `α_max` | `m₂` | `m₄` | `m₆` |
|---|---:|---:|---:|---:|
| `USp(8)` (generic genus 4) | 8 | 1 | 3 | 15 |
| **`SU(2) × 3·SU(2)`** | 8 | **10** | **218** | **6350** |
| **`4·SU(2)`** | 8 | **16** | **512** | **20480** |
| half-twisted mixture (the trap below) | 8 | 6 | 114 | 3210 |
| `SU(2)⁵` | 10 | 5 | 70 | 1525 |
| `SU(2) × 2·SU(2) × USp(4)` | 10 | 6 | 91 | 2034 |
| `SU(2) × USp(4) × USp(4)` | 10 | 3 | 26 | 363 |
| `SU(2)² × USp(6)` | 10 | 3 | 25 | 325 |
| `SU(2)³ × 2·SU(2)` (the `SU(2)⁵` trap) | 10 | 7 | — | — |
| `USp(6) × USp(6)` | 12 | 2 | 12 | 120 |
| `USp(14)` | 14 | 1 | 3 | 15 |

---

## The shapes available at genus four *(proved)*

**Lemma 1.** Let `G = (ℤ/2)^r` act on a smooth projective curve `C` over a
field of characteristic `≠ 2` with `C/G = P¹`. Then

```
2g(C) − 2 = 2^r(−2) + B·2^{r−1},        g(C) = Σ_{χ ≠ 1} g_χ,
```

`B` the number of branch points, `g_χ` the genus of the double cover
`C/ker χ`, branched at `B_χ = {P : χ(v_P) = −1}` where `v_P` generates the
inertia at `P`. The inertia vectors satisfy `Σ_P v_P = 0` (the cover exists)
and `⟨v_P⟩ = G` (it is connected).

*Proof.* Riemann–Hurwitz with all inertia of order 2 (cyclic, hence of order 2
in an elementary abelian 2-group); and `H¹(C) = ⊕_{χ≠1} H¹(C)_χ` because
`H¹(C)^G = H¹(P¹) = 0`, with `H¹(C/ker χ) = H¹(C)_χ` for `χ ≠ 1`. ∎

**Corollary (enumerated in `genus4_symbolic.py`, part 0).**

* `r = 2`, genus 4: the quotient-genus multiset is `(1,1,2)` or `(0,2,2)`, and
  nothing else.
* `r = 3`: `2g − 2 = −16 + 4B`, so `g` is **odd**. **There is no `(ℤ/2)³`-cover
  of `P¹` of genus 4.** At genus 3 the shapes are `(0,0,0,0,1,1,1)`, at genus 5
  they are `(0,0,0,1,1,1,2)` and `(0,0,1,1,1,1,1)`.

This is why `4·SU(2)` at genus 4 cannot be produced the way `3·SU(2)` was at
genus 3. There, three elliptic quotients of a `(ℤ/2)²`-cover are cyclically
permuted by `u ↦ ζu` and hence isomorphic. Four elliptic quotients would need
`(ℤ/2)³`, and Lemma 1 says the genus is then odd. The `(0,0,1,1,1,1,1)` shape
at genus **5** is real — it is what produces `SU(2)⁵` below — but it is one
genus too high for `4·SU(2)`.

## The `SU(2) × 3·SU(2)` pencil *(proved)*

**Theorem C.** Let `K` be a field of characteristic `≠ 2` containing `i` with
`i² = −1`, let `e ∈ {±1}` and let `n ∈ K` with

```
n ≠ 0,   n⁴ ≠ 1,   and   p := n + i + e/n  ∉ {±i, ±n, ±1/n}
```

(explicitly: `e = −1` excludes `n = ±1` and the roots of `n²+in−2` and of
`2n²+in−1`; `e = +1` excludes `n = −2i` and `n = i/2` and the roots of
`n²+2in+1`). Put

```
S = {i, n, e/n},   T = −S,   B₃ = S ⊔ T = {±i, ±n, ±1/n}
f₁ = ∏_{r ∈ T∪{p}} (x−r),      f₂ = ∏_{r ∈ S∪{p}} (x−r)
f₃ = f₁f₂/(x−p)² = (x²+1)(x²−n²)(x²−n⁻²)
```

and let `C_n` be the smooth `(ℤ/2)²`-cover of `P¹_x` with `y₁² = f₁`,
`y₂² = f₂`. Then

1. `C_n` has genus 4 and is **non-hyperelliptic**. (No element of `G` is the
   hyperelliptic involution, since no quotient has genus 0; and if the
   hyperelliptic involution `ι_h` lay outside `G` it would commute with `G`,
   being central in `Aut(C_n)`, so `⟨G, ι_h⟩ ≅ (ℤ/2)³` would act with rational
   quotient — impossible at genus 4 by Lemma 1.)
2. `Jac(C_n) ∼ J₁ × J₂ × J₃` (Lemma 1), with `J₁, J₂` elliptic and `J₃` of
   genus 2;
3. `J₃ ∼ E²` where `E : v² = (u+1)(u−n²)(u−n⁻²)`, `u = x²` — this is
   Theorem A of `REPEATED_FACTOR.md`, `f₃` being even and palindromic of degree
   `2m` with `m = 3` odd;
4. **`J₂ ≅ E` over `K`**, not merely over `K̄`;

hence `Jac(C_n) ∼ J₁ × E³` and `a(C_n) = a(J₁) + 3a(E)` over every finite
extension.

*Proof of (4).* The λ-invariant of the four-point set `S ∪ {p}`, computed with
the Möbius map sending `(n, e/n, i) ↦ (0, ∞, 1)`, is `p`-linear-fractional and
equals `1/n²` exactly when `p = n + i + e/n`; and `1/n²` is the λ-invariant of
`E`'s branch set `{−1, n², n⁻², ∞}`. So the two four-point sets are
`PGL₂(K)`-equivalent. Explicitly for `e = −1`, with

```
M(x) = ( (−i n² − n − i)x + (i n³ + i n − 1) ) / ( n x + 1 )
```

one has `M(n) = −1`, `M(i) = n²`, `M(−1/n) = ∞`, `M(p) = n⁻²`, and

```
(u+1)(u−n²)(u−n⁻²) |_{u = M(x)} · (nx+1)⁴  =  (n²+1)⁴ · f₂(x).
```

The constant `(n²+1)⁴` is a **square** in `K(n)`, so the pullback of `E` along
`M` is `w² = f₂(x)` on the nose and `J₂ ≅ E` over `K`. For `e = +1` the same
computation gives the constant `(n−1)²(n+1)²(n+i)⁴ = ((n²−1)(n+i)²)²`, again a
square. ∎ *(verified symbolically in `genus4_symbolic.py`)*

**Proposition C1 (non-isotriviality and independence).** Over `K = F̄_q(n)`:

```
j(E)  = 256 (n⁴−n²+1)³ / ( n⁴ (n−1)²(n+1)² )                    poles {0, ±1, ∞}
j(J₁) = −256 (n²+in−1)³ (n⁴+3in³−8n²−3in+1)³
        / ( n² (n+i)⁶ (n²+in−2)² (2n²+in−1)² )      [e = −1]
        poles {0, −i, roots(n²+in−2), roots(2n²+in−1), ∞}
```

Both are non-constant, so neither pencil is isotrivial; both have poles, so both
have places of potentially multiplicative reduction, so — by Deligne's
semisimplicity plus the unipotent argument of `FINDINGS.md` — each has
`G_geom = SL₂`. **The two pole sets differ** (`−i` is a pole of `j(J₁)` and not
of `j(E)`), and isogenous elliptic surfaces over the same base have the *same*
places of potentially multiplicative reduction; so `J₁` and `E` are **not
isogenous** over `K`. By Goursat, `G_geom ⊆ SL₂ × SL₂` surjecting onto both
factors is either the full product or the graph of an isomorphism, and the graph
would force an isogeny. Hence `G_geom = SL₂ × SL₂`, Deligne equidistribution
gives `α = t₁ + 3t` with `t₁, t` independent, and the measure is
**`SU(2) × 3·SU(2)`**, `α_max = 8`, `m₂ = 10`. For `e = +1` the pole sets are
`{0,±1,∞}` against `{0, ±1, ±i, −2i, i/2, i(−1±√2), ∞}` — again different. ∎

## Verification of the genus-four pencil *(computed)*

`genus4_witness.py`, part 1. Traces by exact character sums with the correct
count of points at infinity; the `(ℤ/2)²` point count on `C_n` itself is done
independently of Lemma 1, from `1 + χ(f₁) + χ(f₂) + χ(f₃)` per affine `x` plus
four unramified points over `∞`, on a sample of 16 fibres per prime.

| `q` | fibres | `a_C ≠ a₁+3a_E` | `a₂ ≠ a_E` | `a₃ ≠ 2a_E` | direct count fail | distinct `j(E)` | `m₂` | `m₄` | `m₆` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 401 | 392 | **0** | **0** | **0** | **0** | 67 | 10.928 | 223.14 | 5601 |
| 1009 | 1000 | **0** | **0** | **0** | **0** | 169 | 10.021 | 219.66 | 6317 |
| 4001 | 3992 | **0** | **0** | **0** | **0** | 667 | 10.040 | 220.76 | 6502 |
| 16001 | 15996 | **0** | **0** | **0** | **0** | 2667 | **10.059** | **218.30** | **6299** |
| target | | | | | | `≈ q/6` | **10** | **218** | **6350** |

(the `e = −1` pencil; the `e = +1` pencil is in `genus4_witness.csv` and
behaves identically, `m₂ = 9.925` at `q = 16001`.)

Over all 32 primes `13 … 16001` with `q ≡ 1 (mod 4)`: **25292 fibres for
`e = −1` and 25256 for `e = +1`, 50548 in all, 0 mismatches** of
`a(C_n) = a(J₁) + 3a(E)`. The count of distinct `j(E)` values is `≈ q/6`, as it
should be for a degree-12 `j`-map whose six-element `S₃`-fibre meets the squares
about half the time — so the family is **not isotrivial**.

Independence of `J₁` from `E` is not only Proposition C1: `#{n : a(J₁) = ±a(E)}`
is 362 of 15996 at `q = 16001`, against 320–380 for eight random reshufflings of
`a(J₁)` against `a(E)` — indistinguishable from independence. (The naive
estimate `2∫ρ²/√q` is much smaller because `a(E)` takes only 126 distinct values
on the pencil; the shuffle control is the right null.)

## Two traps on this construction *(computed)*

Both are the failure modes brief J named, and one of them is new.

**Trap 3 — the branch that twists half the fibres.** Requiring
`j(J₂) = j(E)` gives six values of `p` per splitting, one for each element of
the `S₃`-orbit of `λ_E = 1/n²`. Only the branch `λ(J₂) = λ_E` itself, and only
on the four splittings with `T = −S`, makes `J₂ ≅ E` over `F_q`; on the other 56
of the 60 branches `J₂` is the quadratic **twist** of `E` on about half the
fibres, so `α` is `t₁ + 3t` on half and `t₁ + t` on the other half, and the
measure is the mixture `½(SU(2)×3·SU(2)) + ½(SU(2)×SU(2))` with

```
m₂ = ½·10 + ½·2 = 6,       m₄ = ½·218 + ½·10 = 114.
```

Measured at `q = 4001`: `a₂ = a_E` on 1993 fibres and `a₂ = −a_E` on 2000 of
3993, `m₂ = 6.06`, `m₄ = 116.1`. Classifying all 60 splitting/branch pairs at
`q = 1009` gives exactly **4 with `J₂ ≅ E` on every fibre and 56 half-twisted**,
none with `J₂` the twist of `E` on every fibre. A geometric `j`-match is not
enough; the pullback constant has to be a square, which is what makes
`(n²+1)⁴` load-bearing.

**Trap 4 — an `S₃`-coincidence between two quotients.** In the `SU(2)⁵`
construction below, the naive parametrisation `a = t, b = t+1, c = t+2` of the
six branch points has two of the five λ-invariants equal to `t` and `(t+1)/t`,
which lie in the **same `S₃`-orbit** (`μ ↦ μ/(μ−1)` carries `t+1` to
`(t+1)/t`). Two of the five elliptic quotients are then isomorphic and the
measure is `SU(2)³ × 2·SU(2)` with `m₂ = 7`, not `SU(2)⁵` with `m₂ = 5`.
Measured `7.040` at `q = 16001`, and the offending pair is identified exactly.
The cure is `a = 2, b = t, c = t+1`.

## `4·SU(2)`: the locus is zero-dimensional *(computed + a moduli count)*

**What the construction gives.** With `B₃` the `Jac ∼ E²` sextic as above, the
free parameters are `(n, p)`. Requiring `j(J₂) = j(E)` — the condition for the
third and fourth copies of `E` — puts `λ(J₂)` in the six-element `S₃`-orbit of
`λ_E`, and since `λ(J₂)` is a Möbius function of `p`, **each of the 10
splittings × 6 branches gives `p` as an explicit rational function of `n`**. So
the locus `j(J₂) = j(E)` is a union of 60 rational curves — four of them
carrying `SU(2)×3·SU(2)` and the other 56 the twisted mixture of Trap 3 — and on
each of them the remaining condition for `4·SU(2)` is `j(J₁) = j(E)`, one
equation in `n`.

**It is never an identity.** Symbolically, for the two splittings that matter
(`genus4_symbolic.py`):

```
e = −1:   num( j(J₁) − j(E) ) = −512 (n−i)² (n²+in−1)³ (n²+(i−1)n−1) (n²+(i+1)n−1)
                                 (n³−2n−i)(n³+2in²−i)(n⁴+3in³−4n²−2in+2)
                                 (2n⁴+2in³−4n²−3in+1)          degree 26
e = +1:   degree 28, likewise not identically zero.
```

Numerically over all 60 branches, `genus4_witness.py` part 3 counts
`#{n ∈ F_q : j(J₁) = j(E)}`:

| `q` | 401 | 1009 | 4001 | 16001 | 40009 |
|---|---:|---:|---:|---:|---:|
| minimum over the 60 branches | 2 | 4 | 4 | 0 | 2 |
| maximum over the 60 branches | 10 | 13 | 18 | 10 | 10 |

**bounded by 18 over a hundredfold range of `q`.** A curve in the `(n,p)`-plane
would give `≈ q` solutions. The locus is a finite set of fibres, and a finite
set of fibres has an atomic vertical Sato–Tate measure, not `4·SU(2)`.

**Why, and how far it generalises.** Over `ℂ`, `Jac ∼ E⁴` with `End(E) = ℤ`
means `End⁰(Jac) ⊇ M₄(ℚ)`, and the locus of such principally polarised abelian
fourfolds is a countable union of **1-dimensional** special (PEL) subvarieties
of `A₄` — modular curves, parametrised by `j(E)` alone. But `dim A₄ = 10` and
the Jacobian locus `M₄ ⊂ A₄` is the **Schottky divisor**, of dimension 9. A
curve in a tenfold meets a divisor in a finite set unless it is contained in it,
so a one-parameter family of genus-4 curves with `Jac ∼ E⁴` and `E` varying
requires one of countably many modular curves to lie inside the Schottky locus —
a Coleman–Oort-type coincidence. The same count run on
`End⁰ ⊇ ℚ × M₃(ℚ)` — the condition for `SU(2)×3·SU(2)` — gives a special locus
of dimension `1 + 1 = 2`, and a **surface** meets a divisor in a curve. That is
exactly the dichotomy observed: `SU(2)×3·SU(2)` came out one-dimensional and
`4·SU(2)` zero-dimensional. It also explains why `2·SU(2)` (genus 2) and
`3·SU(2)` (genus 3) were unobstructed: there `M_g` is dense in `A_g` and there
is no divisor to be inside.

**This is a heuristic, not a proof.** The count is a generic-intersection
expectation, and the whole point of Coleman–Oort is that exceptional special
subvarieties inside the Torelli locus do exist in low genus. What is *proved*
here is only Lemma 1 (no `(ℤ/2)³` route at genus 4) and what is *computed* is
the zero-dimensionality of the 60 branches. Combined with the recorded fact
that the even+palindromic route also gives a zero-dimensional locus, three
independent constructions now fail, which is real evidence and nothing more.

The remaining `(ℤ/2)²` shape at genus 4, `(0,2,2)`, needs both genus-2
quotients to be `∼ E²` with the *same* `E` while sharing five of their six
branch points: 4 moduli against `2 + 2 + 1 = 5` conditions, so the same count
says empty. Non-hyperelliptic genus-4 curves are **not** an untouched second
front — every curve produced by the `(1,1,2)` shape is already non-hyperelliptic
— but genus-4 curves with no `(ℤ/2)²` action at all remain unexamined.

## Genus five, `SU(2) × USp(4) × USp(4)` *(split proved, independence computed)*

The third vertex of cycle 8. A `(ℤ/2)²`-cover with `|D₁| = 4`,
`|D₂| = |D₃| = 2`, quotient genera `(1,2,2)`:

```
u₁ = x⁴ + x + c    u₂ = (x−1)(x−2)    u₃ = (x−3)(x−c)
f₁ = u₂u₃  (genus 1)   f₂ = u₃u₁  (genus 2)   f₃ = u₁u₂  (genus 2)
```

| `q` | fibres | `m₂(C)` | `m₄(C)` | `m₆(C)` | `m₂(E)` | `m₂(J₂)` | `m₂(J₃)` | `m₂(J₂+J₃)` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 401 | 392 | 3.038 | 27.35 | 379 | 0.999 | 0.964 | 0.904 | 1.941 |
| 1009 | 1002 | 3.016 | 26.47 | 368 | 1.001 | 1.075 | 1.016 | 2.043 |
| 4001 | 3992 | 3.048 | 26.59 | 385 | 1.001 | 1.008 | 1.009 | 2.063 |
| 16001 | 15992 | **2.965** | **25.91** | **374** | 1.000 | 1.006 | 0.990 | **1.996** |
| target | | 3 | 26 | 363 | 1 | 1 | 1 | 2 |

Goursat: `G_geom ⊆ SL₂ × Sp₄ × Sp₄` surjects onto each factor; `SL₂ ≇ Sp₄`, so
the elliptic block cannot be linked, and an `Sp₄`–`Sp₄` linkage would give
`m₂(J₂+J₃) = 4`, measured `1.996`. The genus profile `(1,2,2)` is forced by the
branch data (Lemma 1), and `m₄(C) = 25.91` separates `SU(2)×USp(4)×USp(4)`
(26) from `SU(2)²×USp(6)` (25) only weakly — the structural argument, not the
fourth moment, is what identifies it.

## Cycles 8 and 9, at 40 digits *(computed)*

`genus4_cycles.py`, with the **independent** Bessel Toeplitz-minus-Hankel
determinant `E[e^{τ tr}]_{USp(2N)} = det(I_{i−j}(2τ) − I_{i+j}(2τ))` carried in
`mpmath` and both extrema polished by golden section on `log τ` to `10⁻³⁰`, so
no float64 touches the headline numbers. The cross-check against `st_lib.py`'s
Andreief/Chebyshev construction agrees to all 12 printed digits.

```
cycle 8:   SU2 × 3·SU2  ≺  USp12  ≺  SU2 × USp4 × USp4  ≺  SU2 × 3·SU2
           −0.0104701004980744175618332104908580629417
           −0.3193781676625367999674667509788347439026
           −0.0694175511594218671257637055500777423878

cycle 9:   SU2 × 3·SU2  ≺  USp12  ≺  SU2² × USp6  ≺  SU2 × 3·SU2
           −0.0104701004980744175618332104908580629417
           −0.1500338202737787826586997664872715993093
           −0.1187123632823770305842286959426434046775
```

smallest margin, both cycles, `1.04701004980744175618332104909·10⁻²`
(`FINDINGS.md` reported `1.050·10⁻²` on its 1201-point search grid). The
control — cycle 7, the one `REPEATED_FACTOR.md` witnessed — reproduces its three
midranges `−0.07502532855469174755865859708565165487046`,
`−0.15003382027377878265869976648727159930930`,
`−0.01206057421855458525935346125780578377559` to all 40 digits.

> **Theorem 1.** For every odd prime power `q ≡ 1 (mod 4)` with `q ≥ q₀`, the
> genus-4 pencil `C_n` of Theorem C, the genus-6 pencil `y² = x¹³+x+c`, and the
> genus-5 `(ℤ/2)²`-cover with `f₁ = u₂u₃`, `f₂ = u₃u₁`, `f₃ = u₁u₂` for
> `u₁ = x⁴+x+c`, `u₂ = (x−1)(x−2)`, `u₃ = (x−3)(x−c)` have vertical Sato–Tate
> measures `SU(2)×3·SU(2)`, `USp(12)` and `SU(2)×USp(4)×USp(4)`, and the
> limiting midrange comparison on them is a strict 3-cycle with margin
> `1.04701·10⁻²`. The same holds with the genus-5 vertex replaced by
> `REPEATED_FACTOR.md`'s `SU(2)²×USp(6)` pencil
> `y² = f_c(x²)`, `f_c(u) = ((u−1)²−1)((u−1)²−2)((u−1)²−c)`.

Three of `FINDINGS.md`'s nine cycles are now witnessed. **The congruence
`q ≡ 1 (mod 4)` is real**: the splitting `B₃ = S ⊔ T` requires the individual
points `±i` to be rational, and for `q ≡ 3 (mod 4)` Frobenius interchanges `S`
and `T`, giving a swap-coset measure instead. (The `3·SU(2)` construction of
`REPEATED_FACTOR.md` carries the same kind of restriction, `q ≡ 1 (mod 3)`.)

---

## The wider cycles, and the correction they force *(computed)*

`genus4_cycles.py --wide` repeats `symplectic_search.py`'s search with the
isogeny multiplicity capped at `k` and `α_max` pushed past 12. The Bessel
determinants for ranks 7 and 8 are computed from a single table of
`I_0(2τ) … I_{16}(2τ)` per `τ` at `40 + r(r−1)log₁₀(2τ)` digits, and agree with
`witness_search.py`'s `K_usp` to the last printed digit on 16 probes.

| `α_max ≤` | cap `k` | measures | oriented triangles | distinct 3-cycles | widest margin |
|---:|---:|---:|---:|---:|---:|
| 12 | **1** | 29 | 7308 | **0** | — |
| 12 | 2 | 52 | 44200 | 1 | `1.2041·10⁻²` |
| 12 | 3 | 63 | 79422 | 3 | `1.2041·10⁻²` |
| 14 | **1** | 44 | 26488 | **1** | **`4.0050·10⁻²`** |
| 14 | 2 | 83 | 183762 | 7 | `6.1513·10⁻²` |
| 14 | 3 | 105 | 374920 | 18 | `9.5196·10⁻²` |
| 16 | **1** | 66 | 91520 | **3** | `4.0050·10⁻²` |
| 16 | 2 | 137 | 838440 | 31 | `7.5069·10⁻²` |
| 16 | 3 | 177 | 1817200 | 126 | `1.6403·10⁻¹` |

The `α_max ≤ 12` rows reproduce `FINDINGS.md` and `REPEATED_FACTOR.md`
**exactly** — 29 multiplicity-free measures, 7308 oriented triangles, 0 cycles;
then 52/44200/1 and 63/79422/3 at margin `1.2041·10⁻²`. This is an independent
re-implementation (a single Bessel table per `τ` rather than `st_lib`'s
Chebyshev/Andreief route), so the agreement is a check on both. Everything
above `α_max = 12` is new, and the first row of it — one multiplicity-free
3-cycle at `α_max ≤ 14` — is the result that changes the picture.

### The multiplicity-free cycle *(computed; all three vertices witnessed)*

```
SU(2)⁵   ≺   USp(14)   ≺   USp(6) × USp(6)   ≺   SU(2)⁵
 genus 5      genus 7        genus 6
```

40-digit midranges

```
−0.1051934832886848774201778453478543527828      SU2⁵ → USp14
−0.1106336225035651869285581462721661558688      USp14 → USp6×USp6
−0.0401164121961906983999765284263606551578      USp6×USp6 → SU2⁵
```

margin `4.01164121961906983999765284264·10⁻²`. The `α_max ≤ 12` near-miss that
`FINDINGS.md` records, `SU2⁵ → USp4×USp8 → USp6×USp6` at `−7.936·10⁻²`, closes
as soon as the middle vertex is allowed to be `USp(14)` instead of a genus-6
measure: the endpoint gap `Ψ(∞) = α_max` is the third level a midrange 3-cycle
needs, and `14 > 12` supplies it.

**`SU(2)⁵` at genus five** *(split proved, independence proved + computed)*.
By Lemma 1 the shape `(0,0,1,1,1,1,1)` exists at genus 5; take the
`(ℤ/2)³`-cover of `P¹` branched at the six points

```
0 ↦ e₁,   ∞ ↦ e₁,   1 ↦ e₂,   2 ↦ e₁+e₂,   t ↦ e₃,   t+1 ↦ e₁+e₃
```

(sum zero, spanning). Five of the seven non-trivial characters have a
four-point branch locus and two have a two-point one, so

```
Jac(C_t) ∼ E_a × E_b × E_c × E_d × E_e
E_a : y² = x(x−2)(x−t−1)        E_b : y² = x(x−1)(x−t−1)
E_c : y² = x(x−2)(x−t)          E_d : y² = (x−1)(x−2)(x−t)(x−t−1)
E_e : y² = x(x−1)(x−t)          — the Legendre pencil itself
```

λ-invariants `(t+1)/2, t+1, t/2, (1−t)²/(t(t−2)), t`, whose `j`-maps have the
five **distinct** pole divisors

```
{−1, 1, ∞}   {−1, 0, ∞}   {0, 2, ∞}   {0, 1, 2, ∞}   {0, 1, ∞}
```

so no two of the five are isogenous (different places of potentially
multiplicative reduction), each is non-isotrivial with `G_geom = SL₂`, and
Goursat gives `G_geom = SL₂⁵`: the measure is `SU(2)⁵`, `m₂ = 5`. Verified:

| `q` | fibres | `m₂` | `m₄` | `m₆` | linked pairs | direct count fail |
|---:|---:|---:|---:|---:|---:|---:|
| 401 | 397 | 5.301 | 84.11 | 2159 | none | **0** |
| 1009 | 1005 | 4.989 | 67.30 | 1359 | none | **0** |
| 4001 | 3997 | 4.958 | 71.47 | 1624 | none | **0** |
| 16001 | 15997 | **5.031** | **70.36** | 1520 | none | **0** |
| target | | 5 | 70 | 1525 | | |

"direct count fail" compares `a(C_t) = Σ_χ a(C/ker χ)` against a point count of
the smooth `(ℤ/2)³`-cover done from the local structure of the normalisation —
the product rule away from the branch locus, and at each branch point the rule
that the four geometric points above it are rational exactly when `f_ψ` is a
square there for every `ψ` in the order-4 annihilator of the inertia. Zero
mismatches on 12 sampled fibres per prime. (Getting this right matters: the
naive product-over-a-basis formula undercounts at branch points and was wrong on
5 of 8 samples before it was fixed.)

**`USp(14)` at genus seven**: `y² = x¹⁵ + x + c`, the generic hyperelliptic
pencil, big monodromy (Katz–Sarnak); measured `(m₂, m₄) = (1.0071, 3.032)` at
`q = 16001` against `(1, 3)`.
**`USp(6) × USp(6)` at genus six**: `y² = (x²)⁷ + x² + c`; the involution
`x ↦ −x` splits it as `v² = f(u)+c` (genus 3) and `w² = u(f(u)+c)` (genus 3),
`a_C = a_v + a_w` on every fibre with **0 mismatches at every `q` tested**, each
half with `m₂ → 1`, and `(m₂, m₄) = (1.9879, 12.173)` at `q = 16001` against
`(2, 12)`. Goursat: `Sp₆ × Sp₆` surjecting onto each factor is the product or a
graph, and a graph gives `m₂ = 4`.

**No congruence condition anywhere**: `{0, ∞, 1, 2, t, t+1}`, `x¹⁵+x+c` and
`(x²)⁷+x²+c` are defined over the prime field for every odd `q`.

### The wider cycle at multiplicity two *(computed; all three vertices witnessed)*

```
SU(2) × 2·SU(2) × USp(4)   ≺   USp(14)   ≺   USp(6) × USp(6)
      genus 5                   genus 7        genus 6
```

40-digit midranges `−0.1196695279279056449131927908347005053983`,
`−0.1106336225035651869285581462721661558688`,
`−0.0615830920166186043462861218790285729408`; margin
`6.15830920166186043462861218790·10⁻²` — **5.1 times** the `1.206·10⁻²` of the
cycle `REPEATED_FACTOR.md` witnessed, and the widest witnessed anywhere in this
project.

The genus-five vertex is a `(ℤ/2)²`-cover with `|D₁| = 4`, `|D₂| = |D₃| = 2` in
which the genus-2 quotient `J₂` is made to be the `Jac ∼ E²` sextic:

```
D₁ = {i, −i, 1/m, −1/m}    D₂ = {1, 2}    D₃ = {m, −m}
f₁ = (x−1)(x−2)(x²−m²)                  genus 1
f₂ = (x²+1)(x²−m²)(x²−m⁻²)              genus 2, Jac ∼ E²  (Theorem A)
f₃ = (x²+1)(x²−m⁻²)(x−1)(x−2)           genus 2
```

All three are polynomials over the prime field for `m ∈ F_q` — the individual
points `±i` are never separated, so unlike Theorem C this needs **no congruence
on `q`**.

| `q` | fibres | `m₂(C)` | `m₄(C)` | `m₆(C)` | `a(J₂) ≠ 2a(E)` | `m₂(E₁)` | `m₂(J₂)` | `m₂(J₃)` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 401 | 392 | 6.494 | 107.6 | 2814 | **0** | 1.130 | 4.304 | 0.984 |
| 1009 | 1000 | 6.091 | 91.6 | 1914 | **0** | 1.075 | 3.971 | 1.077 |
| 4001 | 3992 | 5.848 | 87.0 | 1963 | **0** | 0.962 | 4.004 | 0.974 |
| 16001 | 15992 | **6.130** | **95.4** | 2200 | **0** | 0.980 | 4.040 | 0.990 |
| target | | 6 | 91 | 2034 | | 1 | 4 | 1 |

> **Theorem 2.** For every odd prime power `q ≥ q₀`, the three pencils
> `C_t` (genus 5, the `(ℤ/2)³`-cover above), `y² = x¹⁵+x+c` (genus 7) and
> `y² = (x²)⁷+x²+c` (genus 6) have vertical Sato–Tate measures `SU(2)⁵`,
> `USp(14)` and `USp(6)×USp(6)`, and the limiting midrange comparison on them
> is a strict 3-cycle with margin `4.0116·10⁻²`. **None of the three has a
> repeated isogeny factor.** Replacing the genus-5 vertex by the
> `SU(2)×2·SU(2)×USp(4)` pencil gives a strict 3-cycle of margin
> `6.1583·10⁻²`.

---

## Addendum — brief K's genus-13 cycle, and which of its vertices has a curve

`TRANSITIVITY.md` (brief K, finished while this session was running) settles
same-genus transitivity: it holds to genus 6, fails from genus 13, and its
genus-13 counterexample

```
USp14×USp4²×SU2²   ≺   USp10×USp8×USp4²   ≺   USp12×USp8×SU2³
   (7,2,2,1,1)             (5,4,2,2)             (6,4,1,1,1)
```

margin `7.998·10⁻³`, is **multiplicity-free** — the same conclusion this session
reaches independently and one genus span lower. The coordinator asked whether to
switch to realising those vertices. **This session did not switch**, for a
reason that the computation below makes concrete, and instead spent a bounded
effort answering the realisability question directly.

**Two elementary constraints** *(proved)*.

**Lemma 2.** For a `(ℤ/2)^r`-cover of `P¹` of genus `g` (Lemma 1) the number of
branch points is `B = (2g − 2 + 2^{r+1})/2^{r−1}`, there are at most `2^r − 1`
non-zero blocks, and every block has genus at most `⌊B/2⌋ − 1`. `B` *shrinks*
as `r` grows, so many blocks and one large block are incompatible:

| genus | `r = 2` | `r = 3` | `r = 4` |
|---:|---|---|---|
| 13 | `B = 16`, ≤ 3 blocks, genus ≤ 7 | `B = 10`, ≤ 7 blocks, genus ≤ 4 | `B = 7`, ≤ 15 blocks, genus ≤ 2 |
| 14 | `B = 17`, ≤ 3 blocks, genus ≤ 7 | — | — |
| 15 | `B = 18`, ≤ 3 blocks, genus ≤ 8 | `B = 11`, ≤ 7 blocks, genus ≤ 4 | — |

**Lemma 3.** If `C → X` is a double cover with `g(X) = h`, then
`Jac(C) ∼ Jac(X) × Prym`, `dim Prym = g − h`, and Riemann–Hurwitz gives
`2g − 2 = 2(2h − 2) + #branch ≥ 4h − 4`, i.e. `g ≥ 2h − 1`. **The new block is
at least the genus of the base minus one.** Recursively: a partition is
reachable by a tower of double covers only if its largest part is at least the
sum of the others minus one, and that condition then has to hold again on the
remainder.

**Lemma 4.** In a `(ℤ/2)^r`-cover of `P¹`, **at most one character can have
`B_χ` equal to the whole branch locus.** If `χ₁ ≠ χ₂` both did, then
`χ₁χ₂(v_P) = 1` for every `P`, and the `v_P` span `G`, so `χ₁χ₂ = 1`. ∎

**The verdict on the six vertices** *(computed, `realisable_partitions.py`)*.
Closing up the three rules — generic hyperelliptic pencil, `(ℤ/2)^r`-cover of
`P¹` for `r ≤ 4`, and the Prym of a double cover of an already-realised base —
gives, at genus 13, 57 realisable partitions out of 101. `TRANSITIVITY.md` lists
two genus-13 cycles; of their six vertices:

| cycle | vertex | blocks | largest | realisable? |
|---|---|---:|---:|---|
| A, margin `8.0·10⁻³` | `(7,2,2,1,1)` | 5 | 7 | **no** |
| | `(5,4,2,2)` | 4 | 5 | **no** |
| | `(6,4,1,1,1)` | 5 | 6 | yes — Prym tower `13 → 7 → 3` |
| B, margin `4.4·10⁻⁴` | `(8,1,1,1,1,1)` | 6 | 8 | yes — Prym `13 → 5` over the `SU(2)⁵` curve above |
| | `(4,4,4,1)` | 4 | 4 | **no** |
| | `(6,3,2,2)` | 4 | 6 | yes — Prym `13 → 7` over a `(ℤ/2)²`-cover with genera `(3,2,2)` |

The three failures each have a one-line reason.

* `(5,4,2,2)`: `5 < 4+2+2−1`, so no Prym step; four blocks needs `r ≥ 3`, and at
  genus 13 an `r = 3` cover has `B = 10`, so no block can exceed genus 4.
* `(7,2,2,1,1)`: the Prym step is allowed (`7 ≥ 6−1`) but leaves `(2,2,1,1)` at
  genus 6, and genus 6 admits **no** `(ℤ/2)³`- or `(ℤ/2)⁴`-cover of `P¹` at all
  (`B` is not an integer), its `(ℤ/2)²`-covers give three blocks, and the Prym
  rule there would need `2 ≥ 3`.
* `(4,4,4,1)`: no Prym step (`4 < 4+4+1−1`), and an `r = 3` cover at genus 13
  has `B = 10`, so three blocks of genus 4 would need `|B_χ| = 10 = B` for three
  different characters — impossible by **Lemma 4**.

**So each of the two genus-13 cycles is missing at least one vertex, and each
missing vertex is missing for a structural reason, not for want of trying.**
Realising them needs a genuinely different mechanism — an Ekedahl–Serre style
curve with several independent low-genus quotients that do *not* come from an
abelian 2-group action, or a non-abelian group action with the right Kani–Rosen
idempotents. That is a research project, not a session. Encouragingly, two of
the three vertices of cycle B *are* reachable, and one of them is built directly
on the genus-5 `SU(2)⁵` curve constructed above — so cycle B, at margin
`4.4·10⁻⁴`, is the one to attack, and the single obstruction is `(4,4,4,1)`.

**Why this session stayed with its own cycle.** The
`SU(2)⁵ ≺ USp(14) ≺ USp(6)×USp(6)` cycle found above already delivers the
headline that brief K's genus-13 cycle delivers — *the `q → ∞` non-transitivity
needs no repeated isogeny factor* — and delivers it strictly better on every
axis that matters:

| | this session's cycle | brief K's genus-13 cycle |
|---|---|---|
| repeated factor | none | none |
| margin | `4.0116·10⁻²` | `7.998·10⁻³` |
| genus span | 5 / 7 / 6 | 13 / 13 / 13 |
| vertices with an explicit pencil | **3 of 3** | 1 of 3 |
| congruence on `q` | none | — |

What brief K's cycle has and this one does not is that it is **same-genus**,
which is a genuinely different and important statement. But a *witnessed*
same-genus cycle is out of reach of every construction enumerated here, and the
cheapest honest next step is not to attack genus 13 head-on but to ask the
restricted question: **is there a same-genus 3-cycle among the realisable
partitions at all?**

> **Corollary (proved, given `TRANSITIVITY.md`'s exhaustive search).** There is
> **no** same-genus 3-cycle all of whose vertices are realisable, at any genus
> `≤ 13`. `TRANSITIVITY.md` §4.3 is exhaustive over every partition of every
> genus to 15 and finds 0 cycles at genus `≤ 12` and exactly 2 at genus 13; the
> table above shows each of those 2 has an unrealisable vertex. ∎

So the arithmetic input `TRANSITIVITY.md` asks for is not merely unfound: at the
first genus where a same-genus cycle exists at all, **no** cycle is within reach
of the constructions catalogued here. The `α_max ≤ 14` cross-genus cycle above
remains the only witnessed non-transitivity that needs no repeated factor.

`realisable_partitions.py` also runs the restricted search directly, genus by
genus, which both re-derives `TRANSITIVITY.md`'s zero counts and measures how
much of the cone a cover construction actually reaches:

| genus | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| partitions | 2 | 3 | 5 | 7 | 11 | 15 | 22 | 30 | 42 | 56 |
| realisable | 2 | 3 | 4 | 7 | 8 | 14 | 12 | 27 | 20 | 33 |
| 3-cycles, all | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3-cycles, realisable | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

(the realisable fraction is much larger at odd genus, where a `(ℤ/2)³`-cover of
`P¹` exists and a `(ℤ/2)²`-cover's three blocks can be broken up further). Genus
14 and 15, where `TRANSITIVITY.md` reports 3 and 6 cycles without listing them
individually, are not settled here.

---

## Corrections

* **To `FINDINGS.md`, and to `REPEATED_FACTOR.md` which inherits it.**
  `FINDINGS.md` states: "**Every one of the nine symplectic cycles uses a
  measure with a repeated isogeny factor**", and "On these — the
  multiplicity-free sub-cone and its swap cosets — the comparison is
  transitive. 29 multiplicity-free measures to genus 6: 7308 oriented
  triangles, **0 cycles**", and concludes "**The `q → ∞` question is now a
  single arithmetic question** … is there a one-parameter family whose
  Jacobian has an isogeny factor of multiplicity `k ≥ 2`". The computation is
  right and reproduces; **the conclusion is an artefact of the `α_max ≤ 12`
  cutoff**. At `α_max ≤ 14` the multiplicity-free sub-cone is *not* transitive:
  `SU(2)⁵ ≺ USp(14) ≺ USp(6)×USp(6)` is a strict 3-cycle with margin
  `4.0116·10⁻²`, and all three vertices have explicit pencils over every odd
  `q`. The repeated isogeny factor is not needed at all — it was only needed
  under `α_max ≤ 12`. `FINDINGS.md`'s own remark that "cycle margins grow with
  the genus span, so genus 7–10 should give wider cycles and possibly a
  multiplicity-free one" was exactly right and is now confirmed.
* **A refinement of brief J's framing, not a correction of the repo.** Brief J
  says the six widest cycles "need `4·SU(2)` or `SU(2)×3·SU(2)`". The grading
  in `REPEATED_FACTOR.md` is the accurate one: the **six widest** need
  `4·SU(2)` alone, and `SU(2)×3·SU(2)` unlocks cycles 8 and 9, which are the
  two **narrowest** (`1.047·10⁻²` against cycle 7's `1.206·10⁻²`). Witnessing
  `SU(2)×3·SU(2)` therefore takes the count from 1 of 9 to 3 of 9 and leaves
  the six widest exactly where they were.
* **No numerical claim of `FINDINGS.md` or `REPEATED_FACTOR.md` is wrong.**
  Cycle 7's three midranges reproduce to all 40 digits; the capped-cone table at
  `α_max ≤ 12` reproduces exactly (`52/44200/1`, `63/79422/3`, margin
  `1.2041·10⁻²`); `FINDINGS.md`'s grid margin `1.050·10⁻²` for cycles 8 and 9
  polishes to `1.04701·10⁻²`.
* **To this session's own first attempt.** The direct point count used to
  check the `(ℤ/2)³` decomposition was written as a product over a basis of
  characters, `∏_i (1+χ(f_i(x)))`. That is correct away from the branch locus
  and **wrong at a branch point**, where two of the three basis polynomials can
  vanish while the fibre still has four rational points; it reported 5 of 8
  sampled fibres as mismatches. The correct local rule — four points, rational
  iff `f_ψ(x)` is a square for every `ψ` in the order-4 annihilator of the
  inertia — gives 0 mismatches everywhere. The `(ℤ/2)²` version used in
  `genus4_witness.py` does not have this problem (the annihilator has order 2
  and the naive sum is right), and reported 0 throughout.

## Open

* **Is `4·SU(2)` realisable at all?** Still open, now with three failed
  constructions (even+palindromic, the 60 branches of the `(1,1,2)` route, and
  by a moduli count the `(0,2,2)` route) and a clean structural reason for the
  difficulty: `M₄` is a divisor in `A₄` and the `E⁴`-locus is a curve. A
  genuine non-existence proof would need a Coleman–Oort-type statement for the
  `M₄(ℚ)`-PEL curves in `A₄`, and would permanently cap the six widest of
  `FINDINGS.md`'s nine cycles. A genuine construction would need a genus-4
  curve with `Jac ∼ E⁴` and no `(ℤ/2)²`-action — none is known to this session.
* **Wider still.** The `k ≤ 3` search at `α_max ≤ 16` finds 126 cycles with a
  widest margin of `1.6403·10⁻¹`, the widest being
  `2·SU(2)×3·SU(2) ≺ USp(16) ≺ 2·USp(6)` at genus `5/8/6`. Two of its three
  vertices are already available — `USp(16)` is `y² = x¹⁷+x+c` and `2·USp(6)` is
  Theorem A at `m = 7` — and only the genus-5 vertex `Jac ∼ E₁² × E₂³` is
  missing. The same moduli count as above says it should not exist:
  `End⁰ ⊇ M₂(ℚ) × M₃(ℚ)` cuts out a **surface** in `A₅` (dim 15) while `M₅` has
  codimension 3, so `2 + 12 − 15 < 0`. **A margin of `0.16` is worth another
  look nonetheless**, because it would shrink `q₀` by another factor. Whether
  the multiplicity-free sub-cone keeps producing wider cycles as `α_max` grows
  is the cheapest remaining question — at `α_max ≤ 16` it is still stuck at
  `4.0050·10⁻²` (3 cycles, but the widest is the `α_max = 14` one padded with an
  extra `SU(2)` on all three vertices).
* **How large is `q₀`, now that the margin is 5× bigger?** The Monte-Carlo
  experiment `FINDINGS.md` proposes — draw `q` samples from each limiting
  measure and run the exact comparison — is now much more likely to land
  somewhere reachable, and is the obvious next computation. The margin went
  from `1.2·10⁻²` to `6.2·10⁻²`; the edge law `q^{−1/t}` with `t = N²+N/2 = 52.5`
  for `USp(14)` is what will fight back.
* **Removing the congruence from Theorem 1.** `C_n` needs `q ≡ 1 (mod 4)`. A
  genus-2 family with `Jac ∼ E²` whose six branch points are individually
  rational for every `q` would remove it; the `V`-symmetric sextics cannot do
  this (the special orbit that must lie in the branch set is
  `Fix(x ↦ −1/x) = {±i}`), so it would need elliptic subcovers of degree > 2.
* **A witnessed same-genus cycle.** By the Corollary above none exists at genus
  `≤ 13` within the constructions catalogued here, and the single obstruction at
  genus 13 is the vertex `(4,4,4,1)` of `TRANSITIVITY.md`'s second cycle — the
  other two vertices of that cycle *are* reachable, one of them directly on top
  of this session's genus-5 `SU(2)⁵` curve. A genus-13 curve family with
  `Jac ∼ A₄ × B₄ × C₄ × E` and independent monodromy would close it and would be
  the strongest result in this line: it removes the genus constraint entirely.
  Lemma 4 says no elementary abelian 2-cover of `P¹` can supply it, so the
  construction has to come from somewhere else — a non-abelian group action, or
  an Ekedahl–Serre style curve with three independent genus-4 quotients.
  Genus 14 and 15 (3 and 6 cycles in `TRANSITIVITY.md`, not listed there
  individually) are not settled here and may well contain a realisable one:
  appending a common block to `(7,2,2,1,1)` gives `(7,2,2,1,1,1)`, whose Prym
  step leaves `(2,2,1,1,1)` at genus 7 — and that *is* within the `r = 3` range,
  `B = 7`, max block genus 2. **That is the cheapest next computation in this
  direction.**
* **Same-genus transitivity below genus 13**, settled by `TRANSITIVITY.md`.

---

## Files

| file | what |
|---|---|
| `genus4_symbolic.py` | the shapes of elementary abelian covers of `P¹` (Lemma 1, enumerated); `p = n + i + e/n` and `λ(J₂) = λ(E)`; the Möbius map `M` and the perfect-square pullback constant proving `J₂ ≅ E` over the base field; the `j`-invariants and their pole divisors; the degree of `j(J₁) − j(E)`; the five `SU(2)⁵` λ-invariants and their distinct pole divisors |
| `genus4_witness.py` | the `SU(2)×3·SU(2)` pencil verified fibre by fibre over 32 primes with an independent `(ℤ/2)²` point count; the branch controls showing the half-twisted mixture; the 60-branch `4·SU(2)` locus count over five primes; the genus-5 `SU(2)×USp(4)×USp(4)` pencil. Writes `genus4_witness.csv` |
| `genus4_cycles.py` | cycles 8 and 9 at 40 digits through the independent Bessel determinant, with cycle 7 as control; `--wide` runs the capped cone search at `α_max ≤ 12, 14, 16` and `k ≤ 1, 2, 3`. Writes `genus4_cycles.csv` |
| `wide_cycle.py` | the four vertices of the `α_max = 14` cycles: `SU(2)⁵` (with the direct `(ℤ/2)³` point count and the `S₃`-coincidence trap), `SU(2)×2·SU(2)×USp(4)`, `USp(14)`, `USp(6)×USp(6)`. Writes `wide_cycle.csv` |
| `realisable_partitions.py` | which multiplicity-free partitions an elementary abelian 2-cover of `P¹` or a Prym tower can supply (Lemmas 2 and 3, enumerated to closure), the verdict on the three vertices of `TRANSITIVITY.md`'s genus-13 cycle, and the same-genus cycle search restricted to the realisable set. Writes `realisable_partitions.csv` |
| `*_output.txt` | console output of the five scripts, kept verbatim |
| `REPEATED_FACTOR.md` | the session this answers |
| `TRANSITIVITY.md` | brief K, whose genus-13 cycle the addendum addresses |
| `FINDINGS.md` | brief F, corrected above |
| `curve_lib.py`, `st_lib.py`, `witness_search.py` | libraries used unchanged |
