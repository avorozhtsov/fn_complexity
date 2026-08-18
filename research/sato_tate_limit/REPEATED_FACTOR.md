# Findings — the repeated isogeny factor, and the `q → ∞` theorem

Answer to the single question `FINDINGS.md` left open. **Yes: the families
exist, they are elementary, and the `q → ∞` non-transitivity is now a
theorem.**

`FINDINGS.md` reduced the whole limiting question to one arithmetic input — a
one-parameter family of curves over `F_q` whose Jacobian is isogenous to `A^k`
for `k ≥ 2` with `A` varying. The construction is one line:

> **`y² = f(x)` with `f` even *and* palindromic of degree `2m`, `m` odd, has
> `Jac ∼ J × J`**, where `J` is the Jacobian of the even half `v² = g(u)`,
> `f(x) = g(x²)`.

The mechanism is that `ι : (x,y) ↦ (1/x, y/x^m)` conjugates the even involution
`σ : (x,y) ↦ (−x,y)` to `σ·(hyperelliptic)` **exactly when `m` is odd**, so the
two elliptic-type quotients — which are otherwise unrelated abelian varieties —
are isomorphic over the base field. `m = 3` is

```
C_c :  y² = (x²+1)(x⁴ + c x² + 1) = x⁶ + (c+1)x⁴ + (c+1)x² + 1
E_c :  v² = (u+1)(u² + c u + 1),        j(E_c) = 256 (c+1)³ / (c+2)
```

with `Jac(C_c) ∼ E_c²`, hence `a_c(C) = 2 a_c(E)` and the vertical Sato–Tate
measure `2·SU(2)`. `j` is non-constant, so the family is **not isotrivial** —
the trap the brief warned about. Verified by exact point counting: **32721
fibres over 29 primes from `q = 11` to `q = 16001`, `a_c(C) = 2 a_c(E)` on every
single one, zero mismatches**, with `(m₂, m₄, m₆) = (3.9997, 31.995, 319.30)` at
`q = 16001` against the `2·SU(2)` values `(4, 32, 320)`.

The same construction at `m = 5, 7` gives `2·USp(4)` (genus 4) and `2·USp(6)`
(genus 6); choosing the half `J` to be a **rationally** split abelian surface
gives `2·SU(2) × 2·SU(2)` at genus 4, which is the vertex `FINDINGS.md` needed.
With two more pencils — one new, one already in the library — all three vertices
of the smallest limiting cycle are witnessed:

| genus | measure | pencil |
|---:|---|---|
| 4 | `2·SU(2) × 2·SU(2)` | `y² = g_w(x²)`, `g_w(u) = (u+1)(u²−Ru+1)(u²−Su+1)`, `R = r+1/r`, `S = s+1/s`, `r = −1/(w⁴−3w²+1)`, `s = 2−w²` |
| 5 | `SU(2)² × USp(6)` | `y² = ((x²−1)²−1)((x²−1)²−2)((x²−1)²−c)` |
| 6 | `USp(12)` | `y² = x¹³ + x + c` |

and the limiting comparison on them is a strict 3-cycle,

```
2·SU(2)×2·SU(2)  ≺  USp(12)  ≺  SU(2)²×USp(6)  ≺  2·SU(2)×2·SU(2)
   genus 4           genus 6        genus 5
```

with midranges, at 40 working digits through an **independent** implementation
of the Sato–Tate cumulant generating function,

```
−0.07502532855469174755865859708565165487046
−0.15003382027377878265869976648727159930930
−0.01206057421855458525935346125780578377559
```

smallest margin `1.20605742185545852593534612578·10⁻²`.

Two further results sharpen the picture. **Multiplicity two is exactly the
threshold**: capping the isogeny multiplicity at `k` in the symplectic library
gives `0, 1, 3, 9, 9` distinct 3-cycles for `k = 1, 2, 3, 4, 6`. And
**multiplicity three is also arithmetic** — a `(ℤ/2)²`-cover of `P¹` branched at
`{a, 1, ζa, ζ, ζ²a, ζ²}` has three elliptic quotients cyclically permuted by
`u ↦ ζu`, hence isomorphic over `F_q` whenever `3 | q−1`, giving `Jac ∼ E³` and
the measure `3·SU(2)` at genus 3 (`a₁ = a₂ = a₃` on every fibre tested,
`m₂ = 8.96` against `9`). What is still missing is `4·SU(2)` and
`SU(2) × 3·SU(2)`, both at genus 4, which is why six of the nine cycles remain
unwitnessed.

Everything below is either proved, or computed and independently re-verified;
each statement says which. Two constructions that look right and are not are
recorded in [Traps](#traps), because both cost this session real time and both
are exactly the failure modes the brief named.

---

## Notation

As `FINDINGS.md`. `a_c = q − N_c`, `α_c = −a_c/√q`, `m_j = E[α^j]` in the limit;
`Ψ_μ(τ) = K_μ(τ)/τ` with `K_μ(τ) = log E_μ[e^{τα}]`; `μ ≺ ν` iff
`mid_τ(Ψ_μ − Ψ_ν) < 0` with `mid = ½(sup + inf)` over `[0,∞]`. `k·G` is `k`
isogenous copies of a block: `α = k·tr(g)`, so `α_max` scales by `k`, `m₂` by
`k²`, and the edge exponent `t` is unchanged.

For a curve `C` over `F_q`, `a(C) = q + 1 − #C(F_q)` is the Frobenius trace of
the smooth projective model. Moment targets used as detectors throughout:

| measure | `m₂` | `m₄` | `m₆` |
|---|---:|---:|---:|
| `SU(2)` | 1 | 2 | 5 |
| `USp(4)` | 1 | 3 | 14 |
| `USp(2g)`, `g ≥ 3` | 1 | 3 | 15 |
| `SU(2)²` | 2 | 10 | 70 |
| `SU(2)² × USp(6)` | 3 | 25 | 325 |
| `2·SU(2)` | 4 | 32 | 320 |
| `2·USp(4)` | 4 | 48 | 896 |
| `2·USp(6)` | 4 | 48 | 960 |
| `2·SU(2) × 2·SU(2)` | 8 | 160 | 4480 |
| `3·SU(2)` | 9 | 162 | 3645 |
| `½(SU2 ∗ SU2) + ½δ₀` (swap) | 1 | 5 | 35 |

---

## The construction *(proved)*

**Theorem A.** Let `K` be a field of characteristic `≠ 2`, let `m` be **odd**,
and let `f ∈ K[x]` be squarefree of degree `2m` with

1. `f(−x) = f(x)`  (**even**), so `f(x) = g(x²)` with `deg g = m`; and
2. `x^{2m} f(1/x) = f(x)`  (**palindromic**).

Let `C : y² = f(x)`, of genus `m − 1`, and `J : v² = g(u)`, of genus
`(m−1)/2`. Then

```
Jac(C)  ∼  Jac(J) × Jac(J)     over K,
```

and in particular `a(C) = 2 a(J)` over every finite extension.

*Proof.* `σ(x,y) = (−x,y)` and `σ' = σ∘ι_h : (x,y) ↦ (−x,−y)` are commuting
involutions of `C`, and `⟨σ,σ'⟩ ≅ (ℤ/2)²` contains the hyperelliptic
involution `ι_h = σσ'`. The quotients are

```
C/σ  :  v² = g(u),        u = x², v = y            genus ⌊(m−1)/2⌋
C/σ' :  w² = u g(u),      u = x², w = xy           genus ⌊m/2⌋
C/ι_h = P¹,               C/⟨σ,σ'⟩ = P¹
```

so the Kani–Rosen idempotent relation for a `(ℤ/2)²`-action with rational
quotient gives `Jac(C) ∼ Jac(C/σ) × Jac(C/σ')`. Because `m` is odd both
quotients have genus `(m−1)/2`.

Now set `ι(x,y) = (1/x, y/x^m)`. Palindromy gives
`(y/x^m)² = f(x)/x^{2m} = f(1/x)`, so `ι ∈ Aut(C)`, and `ι² = id`. Then

```
ι σ ι (x,y) = ι(σ(1/x, y/x^m)) = ι(−1/x, y/x^m)
            = (−x, (y/x^m)·(−x)^m) = (−x, (−1)^m y)
```

which is `σ'(x,y)` **because `m` is odd**. So `ι` conjugates `σ'` to `σ` and
descends to an isomorphism `C/σ' ≅ C/σ` defined over `K`. ∎

The parity is sharp: for `m` even `ισι = σ`, and the two quotients have
different genera `m/2 − 1` and `m/2` and are not isomorphic.

**Corollary A1 (`m = 3`, the genus-two family).** `f_c(x) = (x²+1)(x⁴+cx²+1)`
is even and palindromic of degree 6, squarefree iff `c ≠ ±2`, and

```
C_c : y² = f_c(x)        genus 2,      Jac(C_c) ∼ E_c²
E_c : v² = (u+1)(u²+cu+1) = u³ + (c+1)u² + (c+1)u + 1
```

Writing `a = c+1` and shifting `u = t−1` gives `E_c : v² = t³ + (a−3)t² +
(3−a)t`, whence

```
j(E_c) = 256 a³/(a+1) = 256 (c+1)³ / (c+2)
```

(verified symbolically in `sympy`). `j` is a degree-three rational function of
`c`, so it is **non-constant** and the family is not isotrivial; on `F_q` it
takes about `2q/3` values, the fraction expected of a degree-three map with
Galois group `S₃` (measured 10667 of 15999 at `q = 16001`, `0.66673` against
`2/3`).

**Corollary A2 (higher `m`).** `m = 5` gives genus-four curves with
`Jac ∼ A₂ × A₂`, `A₂` an abelian surface — measure `2·USp(4)`; `m = 7` gives
genus six and `2·USp(6)`. The even+palindromic sextic pencils have
`⌈(m+1)/2⌉ − 1` parameters after normalising `c₀ = c_m = 1`, i.e. 1, 2, 3
parameters for `m = 3, 5, 7`.

## The monodromy, and why the limit is `2·SU(2)` *(proved, modulo standard theorems)*

For the pencil `{E_c}` of Corollary A1: `j` is non-constant, so at a place of
`P¹_c` where `ord(j) < 0` the reduction is potentially multiplicative; local
monodromy there is (quadratic character) ⊗ (non-trivial unipotent), so its
square is a non-trivial unipotent in `G_geom ⊆ SL₂`. By Deligne (Weil II,
3.4.1(iii)) `G_geom` is semisimple, hence not contained in a Borel; a closed
subgroup of `SL₂` containing a non-trivial unipotent and acting irreducibly is
`SL₂`. So `G_geom = G_arith = SL₂`.

Theorem A gives `H¹(C_c) ≅ H¹(E_c) ⊕ H¹(E_c)` as lisse sheaves on the
`c`-line, so the monodromy of the `C`-family is the **diagonal** `SL₂ ⊂ Sp₄`
and Deligne equidistribution gives, as `q → ∞`, the trace measure of the
diagonal `SU(2)`: `α = 2·tr`, i.e. **`2·SU(2)`**, with `α_max = 4`, `m₂ = 4`,
`t = 3/2`. This is exactly the seed measure `FINDINGS.md` could not realise.

## Verification of the genus-two family *(computed)*

`repeated_factor.py`, part 1. `a` computed by exact character sums over `F_q`
with the correct count of points at infinity (`deg 6`, leading coefficient a
square: two points; `deg 3`: one point).

| `q` | fibres | `a_C ≠ 2a_E` | distinct `j` | `m₂` | `m₄` | `m₆` |
|---:|---:|---:|---:|---:|---:|---:|
| 101 | 99 | **0** | 67 | 3.9588 | 31.046 | 309.36 |
| 1009 | 1007 | **0** | 673 | 3.9960 | 31.949 | 317.78 |
| 4001 | 3999 | **0** | 2667 | 3.9990 | 31.991 | 318.42 |
| 8009 | 8007 | **0** | 5339 | 3.9995 | 31.991 | 320.42 |
| 16001 | 15999 | **0** | 10667 | 3.9997 | 31.995 | 319.30 |
| target | | | `2q/3` | **4** | **32** | **320** |

Over all 29 primes `11 … 16001`: **32721 fibres, 0 mismatches**. The identity
`a_c(C) = 2 a_c(E)` is what the brief asked for — the moment alone could be a
coincidence, an exact identity on every fibre of 29 primes is not.

The same table for `m = 5` and `m = 7`, again with `a_c(C) = 2 a_c(J)` on every
fibre and moments converging to `(4, 48, 896)` and `(4, 48, 960)`, is in
`repeated_factor.csv`.

## The three vertices of the cycle

### Genus 4: `2·SU(2) × 2·SU(2)` *(construction derived, split proved, independence computed)*

By Theorem A it suffices to make the half `J` an abelian surface that splits
**rationally** into two independent elliptic curves. `J : v² = g(u)` with `g`
the palindromic quintic with roots `−1, r, 1/r, s, 1/s`, so the branch set is
`B = {∞, −1, r, 1/r, s, 1/s}`. `J` is bielliptic iff some Möbius involution `θ`
preserves `B`. Enumerating the pairings of `B`:

* `(∞,−1)(r,s)(1/r,1/s)` and `(∞,−1)(r,1/s)(1/r,s)` force `s ∈ {r, 1/r, −1}` —
  degenerate;
* `(∞,r)(−1,1/r)(s,1/s)` forces `r²s² + (r²−2r−1)s + r² = 0`, whose
  discriminant is `−(r+1)²(3r+1)(r−1)`;
* `(∞,r)(−1,s)(1/r,1/s)` forces, with `θ(u) = (ru+β)/(u−r)`, `β = r − s(1+r)`,

  ```
  r s² − r s + (1 − r) = 0,      discriminant  r(5r − 4).
  ```

The last is the usable one. `r(5r−4) = □` is a conic, parametrised by
`r = 4/(5−n²)`, `s = (1+n)/2`. But a bielliptic `θ` splits `J` **over the base
field** only when the lift of `θ` to `J` is rational, which happens exactly when
the fixed points of `θ` are — otherwise the two elliptic quotients are conjugate
over a quadratic extension and the geometric monodromy contains the swap (see
[Traps](#traps)). The fixed points satisfy `u² − 2ru − β = 0`, so the condition
is `(r+1)(r−s) = □`; substituting,

```
(r+1)(r−s) = [ (n−1)(n+3) / (5−n²) ]² · (3−n)/2
```

so it reduces to `(3−n)/2 = □`, i.e. `n = 3 − 2w²`. Both conics are rational
and the family is a genuine pencil:

```
r = −1/(w⁴ − 3w² + 1),   s = 2 − w²,   R = r + 1/r,   S = s + 1/s
g_w(u) = (u+1)(u² − R u + 1)(u² − S u + 1)
C_w  :  y² = g_w(x²)      degree 10, genus 4
```

Theorem A **proves** `Jac(C_w) ∼ J_w²` (and the point counts confirm
`a(C) = 2a(J)` on every fibre, 0 mismatches). What is *computed* is that
`J_w ∼ E₁ × E₂` with `E₁, E₂` non-isogenous, hence `Jac(C_w) ∼ E₁² × E₂²`:

| `q` | fibres | `m₂(C)` | `m₄(C)` | `m₆(C)` | `m₂(J)` | `m₄(J)` | `m₆(J)` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 401 | 392 | 7.508 | 142.80 | 3872.6 | 1.877 | 8.925 | 60.5 |
| 1009 | 998 | 7.596 | 151.97 | 4300.3 | 1.899 | 9.498 | 67.2 |
| 4001 | 3992 | 7.523 | 146.08 | 4017.9 | 1.881 | 9.130 | 62.8 |
| 16001 | 15992 | **7.9965** | **162.06** | **4589.1** | **1.9991** | **10.129** | **71.70** |
| target | | 8 | 160 | 4480 | 2 | 10 | 70 |

The identification is not just a moment match. `G_geom` of the `J`-family is a
closed subgroup of `SL₂ × SL₂` surjecting onto each factor (each elliptic
sub-family is non-isotrivial with full `SL₂`), so by Goursat it is either the
full product or the graph of an isomorphism `SL₂ → SL₂` — and the graph would
give `α = 2t` and `m₂ = 4`, not `2`. `m₂(J) = 1.9991` at `q = 16001` excludes
it. So `G_geom = SL₂ × SL₂` and the measure is `2·SU(2) × 2·SU(2)`.

### Genus 5: `SU(2)² × USp(6)` *(split proved, independence computed)*

```
f_c(u) = ((u−1)² − 1)((u−1)² − 2)((u−1)² − c),      C_c : y² = f_c(x²)
```

degree 12, genus 5. The two halves are `v² = f_c(u)` (genus 2) and
`w² = u f_c(u)` (genus 3), independent by the even-split construction of
`FINDINGS.md`. The genus-two half is **even in `u−1`**, so it splits rationally
into `E₁ : v² = (z−1)(z−2)(z−c)` and `E₂ : v² = z(z−1)(z−2)(z−c)`, `z = (u−1)²`.
The genus-three half does **not** inherit the symmetry, because `u ↦ 2−u` does
not fix `{0, ∞}` and so does not preserve `{0, ∞} ∪ roots(f_c)` — which is
precisely what goes wrong in the palindromic variant (see [Traps](#traps)).

| `q` | `m₂` | `m₄` | `m₆` | genus-2 half `m₂` | genus-3 half `m₂` |
|---:|---:|---:|---:|---:|---:|
| 401 | 2.9657 | 24.340 | 315.1 | 1.9948 | 0.9690 |
| 1009 | 2.9866 | 23.316 | 272.4 | 1.9963 | 0.9919 |
| 4001 | 2.9924 | 24.498 | 303.8 | 1.9983 | 0.9985 |
| 16001 | **2.9997** | **25.121** | **335.3** | **1.9999** | **0.9996** |
| target | 3 | 25 | 325 | 2 | 1 |

Goursat again: `G_geom ⊆ SL₂ × SL₂ × Sp₆` surjects onto each factor; `Sp₆` is
not isomorphic to `SL₂`, so it cannot be linked to either; and an `SL₂`–`SL₂`
linkage would give `m₂ = 5`. Measured `3.0000`.

### Genus 6: `USp(12)` *(existing)*

`y² = x¹³ + x + c`, the generic hyperelliptic pencil, big monodromy
(Katz–Sarnak). Measured `(0.950, 2.62, 11.7)` at `q = 16001` against `(1,3,15)`
— the slow convergence at genus 6 that `FINDINGS.md` measures and explains
(`q^{−1/t}` with `t = 39`). The single fibre `c = 0` is genuinely special
(`y² = x(x¹²+1)`) and carries `α = 8.69` at `q = 1009`, which is what inflates
`m₆` in that row; it is an honest member of the pencil, not a counting error.

## The theorem *(proved + computed)*

> **Theorem.** Let `q` be an odd prime power, `q > 5`. The three pencils above
> — `C_w` (genus 4), `C_c` (genus 5), `y² = x¹³+x+c` (genus 6) — have vertical
> Sato–Tate measures `2·SU(2)×2·SU(2)`, `SU(2)²×USp(6)` and `USp(12)`, and the
> limiting midrange comparison on them is a **strict 3-cycle**
>
> ```
> 2·SU(2)×2·SU(2)  ≺  USp(12)  ≺  SU(2)²×USp(6)  ≺  2·SU(2)×2·SU(2)
> ```
>
> with midranges `−0.0750253285546917`, `−0.1500338202737788`,
> `−0.0120605742185546` and smallest margin `1.2060574218554585·10⁻²`.
> Consequently the midrange comparison `≺` on one-parameter families of curves
> over `F_q` is **not transitive** in the `q → ∞` limit; and there is a `q₀`
> such that for every `q ≥ q₀` the finite-`q` comparison on these three pencils
> is not transitive either.

**Hypotheses, stated honestly.**

* The measure identifications are proved for the `Jac ∼ J²` structure
  (Theorem A) and for the `SL₂` monodromy of the elliptic pencils, and
  *computed* for the independence of the factors (via Goursat plus the `m₂`
  measurements above). A fully proved version needs the standard big-monodromy
  input for the individual blocks, which is available in the literature for
  generic hyperelliptic pencils but which this session did not re-derive.
* The finite-`q` half of the conclusion needs `sup_τ |Ψ_f − Ψ_μ| → 0` for each
  of the three families; `mid` is `1`-Lipschitz for the sup norm, so the signs
  then follow. That convergence is what `FINDINGS.md` establishes — bulk error
  `O(q^{−1/2})`, edge error `q^{−1/t}` — measured rather than proved, and it is
  the weakest link in the finite-`q` statement.
* This is a statement about the **limit**. `FINDINGS.md` measures
  `sup_τ |Ψ_f − Ψ_μ| = 0.28` at genus 2 and `1.52` at genus 3 at `q = 4·10⁵`,
  and extrapolates `q ~ 10²⁰` before a genus-6 pencil is within `10⁻²` of its
  limit. The cycle margin here is `1.2·10⁻²`. **`q₀` is astronomically beyond
  any census** and the theorem makes no prediction about one. That is not a
  weakness of the theorem; it is the same separation `FINDINGS.md` already
  established between the limiting comparison and the finite-`q` census, and it
  is why the certified `F_11` and `F_101` cycles are a genuinely different
  phenomenon.
* `≺` here is the `Ψ`-midrange comparison of brief D Part 0, which
  `FINDINGS.md` shows gets all three signs of the certified `F_101` cycle right
  but its magnitudes only to a factor of about two at `q = 101`. The theorem is
  about `≺`, not directly about the exchange rate `C(g → f)`.

The three midranges were recomputed at 40 working digits by
`witness_search.py` with a **second, independent** implementation of the
cumulant generating function — the Toeplitz-minus-Hankel Bessel determinant

```
E[e^{τ tr}]_{USp(2N)} = det( I_{i−j}(2τ) − I_{i+j}(2τ) )_{1 ≤ i,j ≤ N}
```

against `st_lib.py`'s Andreief/Chebyshev construction (agreement to all 12
printed digits at `τ = 0.7, 3, 40` for `SU2, USp4, USp6, USp12`), with both
extrema located by golden section on `log τ` in `mpf` arithmetic to `10⁻³⁰`, so
no float64 touches the headline numbers.

## Multiplicity two is exactly the threshold *(computed)*

`witness_search.py`, the symplectic library of `FINDINGS.md` (`α_max ≤ 12`,
genus 1–6) with the isogeny multiplicity capped at `k`:

| cap `k` | measures | oriented triangles | distinct 3-cycles | widest margin |
|---:|---:|---:|---:|---:|
| 1 | 29 | 7308 | **0** | — |
| 2 | 52 | 44200 | **1** | `1.2041·10⁻²` |
| 3 | 63 | 79422 | **3** | `1.2041·10⁻²` |
| 4 | 68 | 100232 | **9** | `7.4507·10⁻²` |
| 6 | 71 | 114310 | 9 | `7.4507·10⁻²` |

The single `k ≤ 2` cycle is the one witnessed above. The two extra cycles at
`k = 3` are `SU2 × 3·SU2 ≺ USp12 ≺ SU2 × USp4 × USp4` and
`SU2 × 3·SU2 ≺ USp12 ≺ SU2² × USp6`, both needing `E₁ × E₂³` at genus 4; the six
extra at `k = 4` all need `4·SU(2)`, i.e. `Jac ∼ E⁴` at genus 4. So the
arithmetic content of `FINDINGS.md`'s nine cycles is graded exactly by
multiplicity, and multiplicity two — which is now constructed — already
suffices.

## Multiplicity three is arithmetic too *(proved + computed)*

**Proposition.** Let `3 | q − 1`, let `ζ ∈ F_q` be a primitive cube root of
unity and let `a ∈ F_q` with the six points `{a, 1, ζa, ζ, ζ²a, ζ²}` distinct.
The `(ℤ/2)²`-cover of `P¹` with

```
B₁ = {a, 1, ζa, ζ},   B₂ = {ζa, ζ, ζ²a, ζ²},   B₃ = B₁ Δ B₂ = {a, 1, ζ²a, ζ²}
```

is a curve of genus 3 with `Jac ∼ E³`, `E : y² = ∏_{r∈B₁}(u − r)`.

*Proof.* Kani–Rosen gives `Jac ∼ J₁ × J₂ × J₃` with `J_i` the Jacobian of
`y² = f_i`, `f_i = ∏_{r ∈ B_i}(u − r)`, each of genus 1, so the total genus is
3. Substituting `u ↦ ζu` in `f₁` gives `f₁(ζu) = ζ⁴ ∏(u − ζ^{-1}r) = ζ f₃(u)`,
and likewise `f₃(ζu) = ζ f₂(u)`. So each `J_{i+1}` is the twist of `J_i` by
`ζ`. Now `ζ = g^{(q−1)/3}` for a generator `g`, so
`ζ^{(q−1)/2} = g^{(q−1)²/6} = 1` whenever `6 | q − 1`, which holds for odd `q`
with `3 | q−1`. Hence `ζ` is a square, the twists are trivial, and
`J₁ ≅ J₂ ≅ J₃`. ∎

Computed (`repeated_factor.py`, part 3b): `a₁ = a₂ = a₃` on **every** fibre at
`q = 211, 601, 1009, 4021, 16033`, and `m₂ = 8.96, m₄ = 161.0, m₆ = 3624` at
`q = 16033` against the `3·SU(2)` targets `9, 162, 3645`.

This does not yet unlock cycles 8 and 9, which need `SU(2) × 3·SU(2)` at genus
4 — an independent elliptic factor **alongside** the `E³` — and the genus-3
`E³` curve admits no degree-two cover of genus 4 (Riemann–Hurwitz forces genus
`≥ 5`).

## The pencil as a map `A² → A¹` *(proved + computed)*

`FINDINGS.md` notes that the framework wants a family presented as a map
`A² → A¹` so that `Σ_c N_c = q²` and `Σ_c a_c = 0` hold exactly. The
`Jac ∼ E²` pencil is not additive — its base locus is `x²(x²+1) = 0`, so

```
Φ(x,y) = ( y² − (x²+1)(x⁴+1) ) / ( x²(x²+1) )
```

is a rational, not polynomial, map. This costs nothing, because the framework of
this repository is about **maps of finite sets** (`README.md`: "for an onto map
between finite sets, its signature is the decreasing multiset of its non-empty
fibre sizes"), and the exceptional locus is exactly `q` points (if `−1` is a
non-square in `F_q`) or `3q` points (if it is a square) — in either case a
multiple of `q`. Distributing them one or three to a fibre gives an honest map
`Φ : F_q² → F_q` with

```
Σ_c N_c = q²  exactly,        Σ_c a_c = 0  exactly,
a_c(Φ) − a_c(C_c) = constant  (independent of c),
```

so `α_c` shifts by `O(q^{−1/2})` uniformly and the limit measure is unchanged.
Verified in `repeated_factor.py` part 6: the shift set has exactly one element
at every `q` tested, and `max|Δα| ≤ 2/√q`.

## Traps

Both are the failure modes the brief named, and both were hit in this session
before being diagnosed. They are recorded because a reader repeating the work
will hit them too.

**Trap 1 — the split defined only over a quadratic extension.** The first
bielliptic parametrisation tried here used the pairing
`(∞, ρ)(1, −ρ)(τ, −τ)` and the conic `τ² = ρ(2−ρ)`, giving the rational pencil
`ρ = 2/(1+k²)`, `τ = 2k/(1+k²)`. `Jac(C) ∼ J²` holds (`a_C = 2a_J` on every
fibre), and `J` is bielliptic — but the fixed points of `θ` are irrational, the
lift of `θ` to `J` needs a square root, the two elliptic quotients are conjugate
over a quadratic extension, and the geometric monodromy contains the swap. The
measured moments at `q = 16001` are

```
J : (m₂, m₄, m₆) = (1.0169, 5.188, 37.3)   -> 1/2 (SU2 * SU2) + 1/2 delta_0
C : (m₂, m₄, m₆) = (4.0675, 83.01, 2384)   -> 2 x that
```

against `(1, 5, 35)` and `(4, 80, 2240)` — not `(2,10,70)` and `(8,160,4480)`.
This is exactly `FINDINGS.md`'s "swap coset" construction arising unbidden, and
it is why the rationality condition `(r+1)(r−s) = □` is part of the genus-four
pencil and not decoration.

**Trap 2 — a symmetry that propagates to the other half.** For the genus-5
vertex the natural first try is `y² = F(x²)` with `F` a *palindromic* sextic:
the genus-two half then does split rationally into `E₁ × E₂` via `u ↦ 1/u`. But
the branch set of the genus-three half is `{0, ∞} ∪ roots(F)`, which `u ↦ 1/u`
**also** preserves (it swaps `0` and `∞`), so that half splits too and the
measure is `SU(2)³ × USp(4)` with `m₂ = 4`, not `SU(2)² × USp(6)` with `m₂ = 3`.
Measured `(4.031, 45.24, 783.6)` at `q = 16001` against the `SU(2)³ × USp(4)`
values `(4, 45, 794)`. The same thing kills `y² = g(x⁴)`: measured `(3.97, 45.5)`, again
four independent factors of genus `1,1,1,2`, not three of genus `1,1,3`. The
cure is an involution that does **not** normalise `{0, ∞}` — here `u ↦ 2 − u`.

## Same-genus transitivity: the exceptions, characterised *(computed)*

`FINDINGS.md` conjectures that inside a genus the comparison is lexicographic in
`(m₂ ↑, t ↓)` and reports 745 of 765 pairs. `lex_exceptions.py` isolates the 20
failures. **Every one of them has `Δm₂ < 0` and `Δt < 0` simultaneously** — that
is, every exception is a *comonotone* pair, exactly the configuration in which
`FINDINGS.md`'s own criterion says the two `Ψ`s must cross, so the exceptions are
not scattered: they are the crossing pairs on which the edge term outweighs a
small `m₂` gap. `|Δm₂| ∈ {1, 2, 3}` on all 20 and `|Δt| ∈ [1, 9]`.

Replacing the lexicographic rule by a single trade-off,

```
sign mid(Ψ_μ − Ψ_ν)  =  sign( Δm₂ − κ Δt ),      κ ≈ 0.668,
```

raises the agreement from 715 to **731 of the 735 pairs with `Δm₂ ≠ 0`**
(97.28% → 99.46%), and all 30 pairs tied in `m₂` are decided by `t` with
larger `t` preceding, 30 of 30. So the two-coordinate description is much
sharper than `FINDINGS.md` records — but four pairs still escape it, and **this
is not a proof of same-genus transitivity**, which remains open.

---

## Corrections

* **No claim of `FINDINGS.md` is corrected here.** Its cycle margin for the
  smallest cycle, `1.204·10⁻²` on the search grid and `1.206·10⁻²` polished,
  reproduces at 40 digits as `1.20605742185545852593534612578·10⁻²`; its
  multiplicity-free sub-cone is confirmed transitive (0 cycles among 7308
  oriented triangles); and its reading of the required arithmetic input —
  "a genus-4 family with `Jac ∼ E₁²×E₂²`, a generic genus-6 pencil, and a
  genus-5 family splitting as `E₃ × E₄ × A₃`" — is exactly what is built above.
  Its guess that the `k = 2` case "looks like a matter of writing down a pencil
  rather than of existence" was right.
* **Refinement, not correction.** `FINDINGS.md` frames the requirement as
  "multiplicity `k ≥ 2`". The grading is sharper: `k = 2` gives one cycle,
  `k = 3` gives three, `k = 4` gives all nine. So the arithmetic difficulty is
  not uniform across the nine, and the one that is now a theorem is the one
  `FINDINGS.md` singled out.
* **To this session's own first attempt.** The bielliptic parametrisation of
  Trap 1 was believed for one round of computation to realise
  `2·SU(2) × 2·SU(2)`; the moment `m₂ = 4.12` rather than `8` exposed it as the
  swap. Recorded above in full so that the rationality condition is not lost.
* **A moment target in an earlier draft of this session's script was wrong**:
  `E[(t₁+t₂)⁴]` for two independent `SU(2)` traces is `10`, not `14`
  (`2 + 6 + 2`), and `E[(t₁+t₂+t₃)⁶]` for `SU2, SU2, USp6` is `325`, not `245`.
  Both are corrected in `repeated_factor.py`; the conclusions did not depend on
  them.

## Open

* **`4·SU(2)` and `SU(2) × 3·SU(2)`, both at genus 4.** These are the only
  measures still missing from the nine cycles. Within Theorem A, `4·SU(2)`
  needs the genus-two half `J` to itself satisfy `Jac(J) ∼ E²`; the palindromic
  quintic family has two moduli, the bielliptic condition consumes one, and
  "the two elliptic quotients isomorphic" consumes the second — so the
  even+palindromic route gives a zero-dimensional locus, not a pencil. That is
  an obstruction *for this construction only*, not a proof of non-existence.
  The concrete route not taken: a `(ℤ/2)²`-cover with `|B₁| = |B₂| = 4`,
  `|B₁ ∩ B₂| = 1`, hence `|B₃| = 6` and quotient genera `(1, 1, 2)`. Take `B₃`
  to be the branch set of the `Jac ∼ E²` sextic `(x²+1)(x⁴+cx²+1)`, split it
  into two triples `S ⊔ T` and put `B₁ = S ∪ {p}`, `B₂ = T ∪ {p}` for a seventh
  point `p`. Then `Jac ∼ J₁ × J₂ × E²` with `J₁, J₂` elliptic and two free
  parameters `(c, p)`; imposing `J₂ ≅ E` (one condition) gives
  `SU(2) × 3·SU(2)`, and imposing `J₁ ≅ J₂ ≅ E` (two) gives `4·SU(2)`. The
  first is a curve in the `(c,p)`-plane over `ℚ`; whether it is rational, and
  whether `E` varies on it, is a finite computation this session did not run.
  Non-hyperelliptic genus-4 curves were not examined at all. **This is the
  natural next session.**
* **Whether same-genus transitivity is a theorem** — still open; the
  characterisation above narrows the target to four pairs but is not a proof.
* **How large is `q₀`?** Unchanged from `FINDINGS.md`: the honest way to find
  out is the Monte-Carlo experiment it proposes — draw `q` samples from each of
  the three limiting measures and run the exact comparison.
* **Higher genus with multiplicity two.** The `k ≤ 2` sub-cone was searched only
  to `α_max = 12`. `FINDINGS.md` expects wider cycles as the genus span grows,
  so a `k ≤ 2` search to genus 8 could produce a cycle with a much larger margin
  — and therefore a much smaller `q₀`. The cost is `USp(16)`, a few hundred
  digits per `τ`.
* **Multiplicity `k ≥ 4` at all.** `3·SU(2)` is now arithmetic; nothing here
  says whether `4·SU(2)` is, in any construction. A negative answer would be as
  interesting as a positive one, since it would cap the six widest cycles
  permanently.

---

## Files

| file | what |
|---|---|
| `curve_lib.py` | exact point counting for one-parameter families of hyperelliptic curves over prime fields: quadratic character, vectorised Horner over a `(parameters × F_q)` grid, points at infinity, squarefree test, moments |
| `repeated_factor.py` | the constructions and their verification: the `Jac ∼ E²` pencil and the fibre-by-fibre identity over 29 primes; `2·USp(4)`, `2·USp(6)`; the three cycle vertices; `Jac ∼ E³` at genus 3; the two traps; controls; the map-`A²→A¹` presentation. Writes `repeated_factor.csv` |
| `witness_search.py` | the cone search with the multiplicity cap, and the 40-digit recomputation of the cycle through an independent Bessel-determinant implementation. Writes `witness_search.csv` |
| `lex_exceptions.py` | the 20 exceptions to the same-genus lexicographic rule and the trade-off that reduces them to four. Writes `lex_exceptions.csv` |
| `repeated_factor_output.txt`, `witness_search_output.txt`, `lex_exceptions_output.txt` | the console output of the three scripts, kept verbatim |
| `FINDINGS.md` | brief F, the session this answers |
| `verify_cycle.py`, `st_lib.py`, `symplectic_search.py` | brief F's library and cycle verification, used unchanged |
