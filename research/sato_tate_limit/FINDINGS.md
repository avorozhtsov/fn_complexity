# Findings — do the cycles survive `q → ∞`?

Answer to session brief F. **Yes — but not where the brief and its seed looked,
and not for the reason either expected.**

The limiting comparison is **transitive inside a genus** and **non-transitive
across genera**. That is the exact reverse of the seed's diagnosis, which
concluded "cycles cannot come from mixing genera — fix `α_max` and vary the
measure". Every `α_max` class tested is a total order: 0 cycles among 330, 8990
and 190568 oriented triangles at `α_max = 4, 6, 8` in a 149-measure library, and
0 among 11968 at `α_max = 12` in the symplectic one. Every cycle found — 9 among
the symplectic measures to genus six, and, in the separate `α_max ≤ 8` library,
2 once torus factors are admitted and 42 more once orthogonal ones are (44 in
all) — crosses genus classes, and does so because the endpoint gap
`Ψ_μ(∞) − Ψ_ν(∞) = α_max(μ) − α_max(ν)` is the third level a midrange 3-cycle
needs and a fixed genus cannot supply.

The widest is

```
4·SU(2)   ≺   USp(6) × USp(6)   ≺   SU(2) × 2·USp(4)   ≺   4·SU(2)
 genus 4         genus 6              genus 5
```

with midranges `−7.447·10⁻², −1.605·10⁻¹, −1.044·10⁻¹` in `Ψ` units — three
orders above any grid artefact, agreeing to five digits on two independent
`τ`-grids and on golden-section-polished extrema. Two of its six contacts are at
`τ = ∞` and four at interior `τ = 1.63, 2.84, 0.89`, the same mixed pattern the
certified `F_11` and `F_101` cycles show.

**All nine symplectic cycles need a Jacobian with a repeated isogeny factor.**
The multiplicity-free sub-cone — the measures with an explicit pencil witness —
is transitive and misses by `0.0794`; adding the one coset construction that is
demonstrably realisable (a Jacobian whose two halves are exchanged by the
geometric monodromy) leaves it transitive. **So the whole `q → ∞` question
reduces to one arithmetic input**: is there a one-parameter family of curves over
`F_q` whose Jacobian is isogenous to `A^k` for `k ≥ 2` with `A` varying? For
`k = 2` and `A` elliptic that is the classical locus of genus-two curves with
`Jac ∼ E²`, and it supplies a cycle on its own,

```
2·SU(2) × 2·SU(2)   ≺   USp(12)   ≺   SU(2) × SU(2) × USp(6)   ≺   …
   genus 4              genus 6             genus 5
```

margin `1.206·10⁻²` (polished; `1.204·10⁻²` on the search grid).

**And the limit is silent about the certified finite-`q` cycles.** All three
pencils of the `F_101` witness are generic genus-two hyperelliptic pencils, so
all three have the *same* limiting measure `USp(4)`; `Ψ_u − Ψ_v → 0` and the
limiting comparison makes no prediction about them at all. A transitive limit
could never have killed them. Measured: the same-genus cycling fraction shows no
decline over a 260-fold range of `q` — it runs `3.2·10⁻⁵, 2.1·10⁻⁴, 6.3·10⁻⁵,
3.0·10⁻⁴, 1.8·10⁻⁴, 1.6·10⁻⁴, 2.2·10⁻⁴, 2.4·10⁻⁴` at `q = 31…8009`, noisy and if
anything rising.
**There is no `q` at which the census collapses**, and brief F's step 5, as
posed, has no answer because its premise is false. What does die is the *mixing*
of genera, and that is predicted and tested below: genus 2 and genus 3 stop
sharing a largest fibre at `q ≈ 5.9·10⁵`, and the measured sharing probability
falls `0.127 → 3.0·10⁻⁴` from `q = 31` to `q = 50021` on exactly that law.

Everything below is either proved, or computed and independently re-verified;
each statement says which. Five earlier claims are corrected in
[Corrections](#corrections).

---

## Notation

Conventions of `research/m_and_e_and_a_c/PLAN.md` and brief D Part 0.
`a_c = q − N_c`, `α_c = −a_c/√q`, and for a map `f : A² → A¹` over `F_q`

```
Λ_f(β) = log( (1/q) Σ_c (N_c/q)^β ),     log Z_f(β) = (1+β) log q + Λ_f(β)
Ψ_f(τ) = Λ_f(τ√q)/τ,                     Ψ_f(0) = 0,  Ψ_f(∞) = α_max(f)
```

exactly. For a limiting measure `μ` of the `α_c`,

```
K_μ(τ) = log E_μ[e^{τα}],    Ψ_μ(τ) = K_μ(τ)/τ,
Ψ_μ(0) = E_μ[α] = 0,         Ψ_μ(∞) = ess sup supp μ = α_max
```

and `Ψ_μ` is non-decreasing (`K_μ` is convex with `K_μ(0) = 0`, so `K/τ` is a
chord slope from the origin). The comparison is the midrange, brief D Part 0:

```
μ ≺ ν   ⟺   mid_τ (Ψ_μ − Ψ_ν) < 0,     mid = ½(sup + inf) over [0, ∞]
```

Three coordinates of a measure are used throughout, all **additive over
independent products**:

* `α_max` — the endpoint, `= 2g` for any family with infinite monodromy;
* `m₂ = E[α²]` — the variance, which sets the small-`τ` slope `Ψ ≈ m₂τ/2`;
* `t` — the **edge exponent**, `P(α_max − α < ε) ∼ Kε^t`, which sets the
  large-`τ` approach `Ψ = α_max − (t log τ − log A)/τ + O(τ^{−2})`.

`k·G` means `k` isogenous copies of a block: `α = k·tr(g)`, so `α_max` is
multiplied by `k`, `m₂` by `k²`, and `t` is unchanged.

---

## The library, and what validates it *(computed)*

`st_lib.py` builds `K_μ` from the **Weyl integration formula** in Andreief
(determinantal) form. For a classical group of rank `N` with eigenvalues
`e^{±iθ_j}` plus fixed eigenvalues contributing `ε` to the trace, in the
variables `x_j = cos θ_j`,

```
E[e^{τ tr}] = e^{τε} · det( L_τ[x^{i+j}] ) / det( L_0[x^{i+j}] ),   0 ≤ i,j < N
L_τ[p] = ∫ p(x) e^{2τx} w(x) dx
```

with the three Jacobi weights handled in the Chebyshev basis, where every moment
is a Bessel function of argument `2τ`:

| group | `w(x)` | `L_τ[T_k]` |
|---|---|---|
| `Sp(2N)` | `(2/π)√(1−x²)` | `I_k − (I_{k−2}+I_{k+2})/2` |
| `SO(2N)` | `(1/π)/√(1−x²)` | `I_k` |
| `SO(2N+1)` | `(1/π)(1−x)/√(1−x²)` | `I_k − (I_{k−1}+I_{k+1})/2` |

computed at `mpmath` precision set from the cancellation budget
`rank(rank−1)·log₁₀(2τ)`, which is 159 digits for `USp(12)` at `τ = 10⁵`.

Three independent checks, `validate_library.py`:

* **against brute-force Weyl quadrature** in `g` dimensions (tensor Gauss–
  Legendre in the angles, up to 6000 nodes per axis) for `USp(2), USp(4),
  USp(6), SO(3..6), U(1)` at `τ = 0.25, 1, 3, 6`: worst relative deviation
  **`9.3·10⁻¹³`**;
* **against exactly known moments** — `E[tr^{2k}]` over `USp(2g)` is the number
  of perfect matchings of `2k` points, cut down by the Brauer relations when
  `2g < 2k`:

  | | `E[tr²]` | `E[tr⁴]` | `E[tr⁶]` | `E[tr⁸]` |
  |---|---:|---:|---:|---:|
  | `USp(2)` | 1 | 2 | 5 | 14 |
  | `USp(4)` | 1 | 3 | 14 | 84 |
  | `USp(6)` | 1 | 3 | 15 | 104 |
  | `USp(8)` | 1 | 3 | 15 | 105 |
  | `U(1)` | 2 | 6 | 20 | 70 |

  reproduced to `≤ 1.1·10⁻⁹` (the finite-difference scheme is second order on the
  symmetric weights; on `SO(3)`, whose weight is asymmetric, it is first order
  and the error is `10⁻⁶`), with `E[tr²] = 1` for `SO(3..7)`. The mean-zero
  constraint `E[α] = 0`, which every limit measure of a fibration of `A²`
  satisfies exactly because `Σ_c a_c = 0`, holds to `≤ 1.7·10⁻²¹` by central
  difference for all ten groups;
* **the edge exponent**, fitted from `τ(α_max − Ψ(τ))` at `τ = 2·10³…1.6·10⁴`
  against the closed form derived here — for a weight behaving like `(1−x)^a` at
  the top edge the density `∏(u_i−u_j)²∏u_i^a` is homogeneous of degree
  `N(N−1)+aN` in `u_i = 1−x_i` and there are `N` differentials, so

  ```
  t = N² + aN
  ```

  which is `g(2g+1)/2 = dim USp(2g)/2` for `Sp(2g)`, `N²−N/2` for `SO(2N)` and
  `N²+N/2` for `SO(2N+1)`: agreement to `≤ 4.7·10⁻⁴` on nine groups.

Endpoints check out: `Ψ(10⁻⁶) ≤ 10⁻⁶` and `Ψ(10⁵)` is within `6·10⁻⁴` of
`α_max` for every group in the library.

---

## `α_max` is the genus, and nothing else *(proved)*

**Proposition.** Let `G ⊆ USp(2g)` be a compact subgroup of positive dimension.
Then its trace measure has `ess sup = 2g`.

*Proof.* `tr` is continuous and `tr(1) = 2g`; a neighbourhood of the identity in
`G` has positive Haar measure, so `tr > 2g − δ` on a set of positive measure for
every `δ > 0`. ∎

So for every family of curves with infinite monodromy `α_max = 2g` **whatever
the monodromy group is**, and fixing `α_max` is the same as fixing the genus.
The seed brief's instruction — "fix `α_max` and vary the measure", on the grounds
that "the endpoint gap between genus classes is `O(1)` and so is the interior, so
the endpoint wins" — throws away the one degree of freedom a cycle needs. The
endpoint gap and the interior excursion being *the same order* is exactly what
lets an interior excursion overturn an endpoint gap; that is brief B's mechanism
at `q < ∞`, and it survives the limit.

---

## What decides the comparison: `(α_max, m₂, t)` *(derived + computed)*

From the two asymptotics of `Ψ`, at fixed `α_max`

```
Ψ_μ − Ψ_ν  ≈  (m₂(μ) − m₂(ν))·τ/2                        (τ → 0)
Ψ_μ − Ψ_ν  ≈  (t(ν) − t(μ))·log τ/τ                      (τ → ∞)
```

so the small-`τ` order is by `m₂` and the large-`τ` order is by `−t`, and
**`Ψ_μ − Ψ_ν` changes sign iff `m₂` and `t` are ordered the same way** — the
trade the cycles session identified independently from the finite-`q` side ("a
signature whose `Ψ` rises faster at small `τ` (larger `m₂`) must flatten before
the common endpoint `α_max`"). Tested on every same-`α_max` pair of the
149-measure library, the criterion `(m₂(μ)−m₂(ν))(t(μ)−t(ν)) > 0` agrees with the
computed sign-change count on **3954 of 4009 pairs (98.63%)**; all 55
disagreements are exact `m₂` ties, where the product criterion is vacuous and the
crossing is decided at higher order.

A pair whose `Ψ`s do not cross is ordered by pointwise domination, and a
tournament all of whose edges are pointwise dominations is transitive, so **a
3-cycle needs all three pairs to cross**, hence needs `(m₂, t)` comonotone on the
triple.

**Inside a genus the comparison is a lexicographic order** — `m₂` ascending,
ties broken by `t` descending — on 745 of 765 same-genus pairs (97.39%) of the
symplectic library to genus six. The 20 exceptions all have `m₂` differing by 1
and `t` differing by a large amount: the edge exponent overturning a small `m₂`
gap. In particular:

* **335 pairs of the full library are matched in `α_max` *and* `m₂`** — the
  configuration of the certified `F_101` witness, where two of the three pencils
  share max fibre, `m₂` and `ν(P)`. 16 of them also tie in `t`; of the remaining
  319, **302 (94.7%) follow the edge exponent**, larger `t` preceding, and 42 of
  the 335 cross. So `t` is the decider on a matched pair but not an infallible
  one — which is the limiting form of brief E's finding that `log μ`, the
  multiplicity of the largest fibre, predicts the certified edge where the moment
  ladder does not. The identification of `log μ` with `t` is made exact below.
* The limit reproduces brief E's **sign**: larger `m₂` *follows*. Brief B's
  addendum 1 has it the other way, and is wrong.

`level_lemma.py` says how much freedom a cycle needs. Sampling the three
differences at `n` interior levels and solving the resulting LP exactly over all
argmax/argmin patterns, with `|D| ≤ 1`:

| interior levels `n` | both ends pinned (one genus) | one end free (two genera) |
|---:|---:|---:|
| 1 | `0` | `0` |
| 2 | **`0`** | **`+1/4`** |
| 3 | `+1/4` | `+1/4` |

**The endpoint gap is worth a whole level.** This is an explanation and not a
proof — three continuous functions attain six extrema at up to six distinct `τ`,
so `n = 2` is a hypothesis about how aligned those extrema are — but it matches
the library exactly.

---

## Inside a genus: a total order *(computed)*

| library | `α_max` | measures | oriented triangles | 3-cycles |
|---|---:|---:|---:|---:|
| full (products of `SU(2)`, `U(1)`, `CM`, `USp(2g)`, `SO(n)`, mult ≤ 4) | 4 | 11 | 330 | **0** |
| | 6 | 31 | 8990 | **0** |
| | 8 | 84 | 190568 | **0** |
| symplectic only, to genus 6 | 4 | 3 | 2 | **0** |
| | 6 | 5 | 20 | **0** |
| | 8 | 11 | 330 | **0** |
| | 10 | 17 | 1360 | **0** |
| | 12 | 34 | 11968 | **0** |

The closest any same-genus triangle comes to closing is `−8.05·10⁻²` at genus 6.

**The seed's ten genus-two measures reproduce exactly** at proper resolution
(3001-point `τ` grid to `10⁵`, against the seed's 4000-bin measures and 400-point
grid to `τ = 400`):

```
0 cycles of 240 oriented triangles           (seed: 0 of 120 unordered triples)
curl residual ‖mid − grad ψ‖/‖mid‖ = 6.779222·10⁻²        (seed: 6.8·10⁻²)
crossing pairs 3 of 45                                    (seed: 3 of 45)
USp4 ≺ SU2² ≺ SU2×CM ≺ CM² ≺ SU2×U1 ≺ U1×CM ≺ 2·SU2 ≺ U1² ≺ 2·CM ≺ 2·U1
```

identical to the seed's order, and the three crossing pairs are the same three.
**The seed's numbers were right; the conclusion drawn from them was not.**

---

## Across genera: the cycles *(computed, independently verified)*

`symplectic_search.py`, over all 71 products of `Sp(2g_i)` blocks with
multiplicities and `α_max ≤ 12` (genus 1–6): **114310 oriented triangles, 9
3-cycles**, every one with genus profile `4 / 6 / 5`.

| cycle | margin |
|---|---:|
| `4·SU2 ≺ USp6×USp6 ≺ SU2×2·USp4` | `7.451·10⁻²` |
| `4·SU2 ≺ USp4×USp8 ≺ SU2×2·USp4` | `6.912·10⁻²` |
| `4·SU2 ≺ USp4×USp8 ≺ SU2⁵` | `6.912·10⁻²` |
| `4·SU2 ≺ USp6×USp6 ≺ SU2⁵` | `4.009·10⁻²` |
| `4·SU2 ≺ USp6×USp6 ≺ SU2×2·SU2×USp4` | `2.841·10⁻²` |
| `4·SU2 ≺ USp4×USp8 ≺ SU2×2·SU2×USp4` | `2.841·10⁻²` |
| `2·SU2×2·SU2 ≺ USp12 ≺ SU2²×USp6` | `1.204·10⁻²` |
| `SU2×3·SU2 ≺ USp12 ≺ SU2×USp4×USp4` | `1.050·10⁻²` |
| `SU2×3·SU2 ≺ USp12 ≺ SU2²×USp6` | `1.050·10⁻²` |

(margins on the 1201-point search grid of `symplectic_search.py`; the two
verified in `verify_cycle.py` move by `2·10⁻⁵` and `2·10⁻⁵` when the extrema are
polished.)

The widest, verified three ways in `verify_cycle.py` (search grid `10⁻⁴…10⁵`
3001 points; independent grid `10⁻⁵…10⁶` 1501 points; extrema polished by golden
section on `log τ` to `10⁻¹³`, so the grid discretisation drops out):

| edge | grid A | grid B | polished | `sup` at `τ` | `inf` at `τ` |
|---|---:|---:|---:|---|---|
| `4·SU2 → USp6×USp6` | `−0.07447040` | `−0.07448845` | `−0.07446973` | `+3.8511` at `1.63` | `−4.0000` at `∞` |
| `USp6×USp6 → SU2×2·USp4` | `−0.16054317` | `−0.16051889` | `−0.16054706` | `+2.0000` at `∞` | `−2.3211` at `2.84` |
| `SU2×2·USp4 → 4·SU2` | `−0.10437566` | `−0.10438033` | `−0.10438075` | `+2.0000` at `∞` | `−2.2088` at `0.89` |

Four interior contacts and two endpoint contacts — the same signature as the
certified `F_11` cycle, one regime up.

Two more are verified the same way in `verify_cycle.py`. The one that puts the
least weight on the unproved arithmetic is

```
4·SU2   ≺   USp(4) × USp(8)   ≺   SU(2)⁵   ≺   4·SU2
```

midranges `−0.06910, −0.07936, −0.14261`, smallest margin `6.910·10⁻²`: **two of
its three vertices are multiplicity-free**, so it needs a repeated isogeny factor
at only one vertex. And the one that needs only multiplicity two,

```
2·SU2 × 2·SU2   ≺   USp(12)   ≺   SU(2)² × USp(6)
```

midranges `−0.07503, −0.15003, −0.01206`, smallest margin `1.206·10⁻²`.

**The mechanism.** `4·SU2` has `α_max = 8` but `m₂ = 16`, so its `Ψ` rises
fastest and then stalls; `USp6×USp6` has `α_max = 12` but `m₂ = 2` and `t = 21`,
so it rises slowest and arrives highest; `SU2×2·USp4` is intermediate in both.
The `4/6/5` genus profile is forced: the cycle needs the endpoint order and the
interior order to disagree, and only a genus gap can put a large-`m₂` measure
below a small-`m₂` one at `τ = ∞`.

---

## Realisability — where the arithmetic re-enters *(proved + open)*

Brief F: "this is where the arithmetic re-enters, and it is the only place it
does". Three tiers, and only the third survives.

**Orthogonal measures are not `H¹` of a curve family** *(proved)*. `H¹` of a
curve is symplectically self-dual, so `G_geom ⊆ Sp(2g)`. The `SO(n)` trace
measures are in the library for the contrast the brief asked for; they take the
`α_max ≤ 8` library from 2 cycles to 44, and must be discarded.

**Torus measures — `U(1)` (arcsine) and `CM = ½δ₀ + ½·arcsine` — are not
vertical Sato–Tate measures** *(proved, modulo the cited theorem)*. By Deligne
(Weil II, 3.4.1(iii)) the geometric monodromy group of a pure lisse sheaf on a
normal variety over a finite field is **semisimple**; a positive-dimensional
torus is not, so it cannot be the identity component of `G_geom` for
`R¹π_*Q_ℓ`. Independently, for an elliptic pencil: if `j` is non-constant then
some place has potentially multiplicative reduction, so `G_geom` contains a
non-trivial unipotent; a unipotent lies neither in a torus nor in the non-trivial
coset of its normaliser (whose elements have trace `0` and are semisimple), so
`G_geom ⊄ N(T)`; and a closed subgroup of `SL₂` containing a non-trivial
unipotent acts irreducibly and is `SL₂`. If `j` is constant the family is
isotrivial and `G_geom` is finite, giving an atomic measure, not the arcsine.
**The arcsine is a horizontal (vary the prime) Sato–Tate measure, not a vertical
one.** Discarding it and the `CM` mixture removes the remaining 2 cycles of that
library and 7 of the seed's 10 genus-two measures.

**Symplectic measures.** Explicit witnesses:

| measure | witness |
|---|---|
| `USp(2g)` | `y² = h(x) + c`, `deg h = 2g+1` generic — big monodromy |
| `USp(2⌊(n−1)/2⌋) × USp(2⌊n/2⌋)` | `y² = f(x²) + c`, `deg f = n`: the involution `x ↦ −x` gives quotients `v² = f(u)+c` and `w² = u(f(u)+c)` with `u = x²`, of genus `⌊(n−1)/2⌋` and `⌊n/2⌋` with independent monodromy. `n = 3` gives `SU(2)²`, `n = 4` gives `SU(2)×USp(4)`, `n = 5` gives `USp(4)²`, `n = 7` gives `USp(6)²`. Verified numerically: `y² = h(x²)+c` has `m₂ → 2` (table below), the `SU(2)²` value. |
| `½(μ ∗ μ) + ½δ₀` | a Jacobian whose two halves are exchanged by the *geometric* monodromy: the non-identity coset consists of `[[0,X],[Y,0]]`, whose trace vanishes identically |

**On these — the multiplicity-free sub-cone and its swap cosets — the comparison
is transitive.** 29 multiplicity-free measures to genus 6: 7308 oriented
triangles, **0 cycles**, closest approach `−7.936·10⁻²`, the near-miss being
`SU2⁵ → USp4×USp8 → USp6×USp6`. Adding the swap cosets that fit under
`α_max = 12` (35 measures, 13090 triangles): still **0 cycles**, closest approach
unchanged at `−7.936·10⁻²`.

**Every one of the nine symplectic cycles uses a measure with a repeated isogeny
factor**, `k·G` with `k ≥ 2` — a family whose Jacobian is isogenous to `A^k` with
`A` varying. Hence:

> **The `q → ∞` question is now a single arithmetic question.** Cycles persist
> for fixed genus profile as `q → ∞` **iff** there are one-parameter families of
> curves over `F_q`, of the genera listed, whose Jacobians have an isogeny factor
> of multiplicity `k ≥ 2` with non-isotrivial monodromy.

For the mildest case — `k = 2` with `A` elliptic, i.e. `Jac ∼ E²` for a genus-two
family, the classical locus in `M₂` — a cycle already exists,
`2·SU2 × 2·SU2 ≺ USp(12) ≺ SU(2)²×USp(6)`, verified the same three ways at
midranges `−0.07503, −0.15003, −0.01206` (smallest margin `1.206·10⁻²`), needing a
genus-4 family with `Jac ∼ E₁²×E₂²`, a generic genus-6 pencil, and a genus-5
family splitting as `E₃ × E₄ × A₃`. **I have not exhibited any multiplicity-`k`
family explicitly, and that is the single open point of this session.** Two
further caveats, both real: the family must have monodromy exactly the diagonal
`SL₂ ⊂ Sp(2k)` and not something larger, and the brief's framework wants it
presented as a map `A² → A¹` so that `Σ_c N_c = q²` and `Σ_c a_c = 0` hold on the
nose.

Convex combinations of *arbitrary* pairs of library measures do cycle, and
abundantly — **3091 3-cycles** among the 435 measures obtained by adding all
uniform averages of pairs of multiplicity-free measures, widest margin
`1.686·10⁻¹`, the widest being
`SU2⁵ ≺ ½(SU2×USp4×USp4) + ½USp12 ≺ ½SU2 + ½(USp6×USp6)`; free (non-uniform)
weights do the same, `−5.306·10⁻²` from a six-restart Powell spot check. But a
prescribed average of two prescribed measures corresponds to a monodromy group
whose cosets have exactly those trace measures with exactly those weights, and I
can exhibit no such group beyond the swap construction above. **These are a
statement about the convex hull, not about curves** — and the contrast is sharp:
mixtures of multiplicity-free measures cycle, while *products* of them do not.

---

## The finite-`q` regime, and how far it is from the limit *(computed)*

Brief F requires `Ψ_f → Ψ_μ` be demonstrated before anything is concluded from
`Ψ_μ`. `finite_q.py`, 12 random pencils per point, traces computed exactly by FFT
from `N_c = q + (m_h ⋆ χ)[c]` and checked against brute force at `q = 101, 211`:

| family | limit | `m₂` at `q = 4·10⁵` | limit `m₂` | `|ΔΨ|` at `τ = 2` | edge deficit |
|---|---|---:|---:|---:|---:|
| `y² = h(x)+c`, `deg h = 3` | `SU(2)` | 1.0000 | 1 | `5.8·10⁻⁴` | `−1.5·10⁻⁴` |
| `deg h = 5` | `USp(4)` | 0.9989 | 1 | `1.8·10⁻³` | `−0.2797` |
| `deg h = 7` | `USp(6)` | 0.9994 | 1 | `2.9·10⁻³` | `−1.5171` |
| `y² = h(x²)+c`, `deg h = 3` | `SU(2)×SU(2)` | 2.0012 | 2 | `1.7·10⁻³` | `−0.0433` |

`m₂` identifies the symmetry type unambiguously, including the split family at
`m₂ → 2`. **The deviation has two scales with different exponents**, fitted over
`q = 101…409609`:

| family | bulk exponent (`τ = 2`) | edge exponent | `1/t` |
|---|---:|---:|---:|
| `SU(2)` | 0.641 | 0.797 | 0.667 |
| `USp(4)` | 0.480 | 0.210 | 0.200 |
| `USp(6)` | 0.444 | 0.102 | 0.095 |
| `SU(2)×SU(2)` | 0.628 | 0.305 | 0.333 |

The bulk is the `O(q^{−1/2})` error of a `q`-sample empirical mean. The edge is
an extreme value, and **the general law is `q^{−1/t}` with `t` the edge
exponent** — T2.1's `2/dim USp(2g)` is the special case `t = dim/2`, and the
split family, which T2.1's form does not cover, comes out at `1/3` against a
fitted `0.305`. This is an independent *arithmetic* confirmation of the `t`
values that decide the limiting comparison.

**The edge dominates by `q^{3/10}` at genus two, and it is enormous at every
reachable `q`.** At `q = 409609` the genus-two `Ψ_f` is still `0.28` below
`Ψ_{USp(4)}` and the genus-three one `1.52`; `sup_τ |Ψ_f − Ψ_μ|` equals the edge
deficit in almost every row. **The finite-`q` discrepancy is concentrated at
`τ = ∞`, which is exactly where `φ` reads.** The limiting cycles above, with
margins `10⁻²`, are statements about `q` far beyond any census.

### The edge datum at finite `q` is `log(multiplicity of the largest fibre)`

Exactly, for large `τ`,

```
Ψ_f(τ) = α_max(f) + log(mult_f/q)/τ + O(τ^{−2})
Ψ_μ(τ) = α_max     − (t log τ − log A)/τ + O(τ^{−2})
```

so `log mult` and `−t log τ + log A` occupy the same slot — the coefficient of
`1/τ` at the top end. Measured on the `F_101` witness at `τ = 200`:

| | mult | `log(mult/q)` | `τ(Ψ_f − α_max)` |
|---|---:|---:|---:|
| `f₁` | 2 | `−3.921973` | `−3.921973` |
| `f₂` | 1 | `−4.615121` | `−4.615121` |
| `f₃` | 3 | `−3.516508` | `−3.516508` |

to six decimals. **This identifies brief E's `log μ` — the one variable that
predicts the certified edge where the moment ladder does not — with the limiting
edge exponent `t`**, and the two signs agree: larger `t` precedes in the limit,
larger `mult` (effectively smaller `t`) follows at finite `q`, and on the
`(M, m₂, ν)`-tied pair of the `F_101` witness `f₂ ≺ f₁` with `mult = 1` against
`mult = 2`.

### The `F_101` cycle, recomputed and read in the limit

Signatures recomputed from the three polynomials by FFT point count (each sums to
`q² = 10201`); comparison on a `β` grid to `360q`:

| edge | `C(a→b)` | `C(b→a)` | margin |
|---|---:|---:|---:|
| `f₁ → f₂` | `0.999320082239` | `0.998996411149` | `3.237·10⁻⁴` |
| `f₂ → f₃` | `0.998563509180` | `0.998303616396` | `2.599·10⁻⁴` |
| `f₃ → f₁` | `0.998303616396` | `0.998071094651` | `2.325·10⁻⁴` |

All three edges agree in sign — `f₁ ≻ f₂ ≻ f₃ ≻ f₁`, smallest margin
`2.325·10⁻⁴` — confirming addendum 2 independently. And the **exact `Ψ`-midrange
law gets all three signs right**:

```
mid(Ψ_{f1} − Ψ_{f2}) = +0.012025   →  margin −5.185·10⁻⁴  (exact −3.237·10⁻⁴)
mid(Ψ_{f2} − Ψ_{f3}) = +0.009220   →  margin −3.976·10⁻⁴  (exact −2.599·10⁻⁴)
mid(Ψ_{f3} − Ψ_{f1}) = +0.002425   →  margin −1.046·10⁻⁴  (exact −2.325·10⁻⁴)
```

This is *not* brief B's leading-order law, which gets the sign of a certified
`F_11` edge wrong: `Ψ_f` here is the exact `Λ_f(τ√q)/τ`, with no expansion of
`log(1+α/√q)`; only the denominator `(1+β)log q + Λ_v` is truncated. Magnitudes
are right to a factor of about two at `q = 101`. **The two regimes still have to
be kept apart**: the limit is `Ψ_μ`, and `Ψ_{f_i}` at `q = 101` is `1.8` away
from it.

---

## Step 5: what dies, what does not, and the tested prediction *(computed)*

**What does not die.** Two families with the same limiting measure have
`Ψ_u − Ψ_v → 0`, so the limiting comparison is silent and their sign is pure
fluctuation. All three `F_101` pencils are of that kind. Measured on pools of
3000 genus-two pencils, restricted to the largest class of equal largest fibre —
the only pairs that can carry a cycle:

| `q` | members | oriented triangles | 3-cycles | fraction | median margin `×√q log q` |
|---:|---:|---:|---:|---:|---:|
| 31 | 80 | 158080 | 5 | `3.2·10⁻⁵` | 0.0728 |
| 101 | 80 | 158080 | 33 | `2.1·10⁻⁴` | 0.0678 |
| 211 | 80 | 158080 | 10 | `6.3·10⁻⁵` | 0.0646 |
| 401 | 80 | 158080 | 47 | `3.0·10⁻⁴` | 0.0520 |
| 1009 | 80 | 158080 | 28 | `1.8·10⁻⁴` | 0.0621 |
| 2003 | 80 | 158080 | 26 | `1.6·10⁻⁴` | 0.0524 |
| 4001 | 80 | 158080 | 34 | `2.2·10⁻⁴` | 0.0391 |
| 8009 | 80 | 158080 | 38 | `2.4·10⁻⁴` | 0.0353 |

Over a 260-fold range of `q` the cycling fraction does not decline; it is noisy
(the counts are 5 to 47 out of 158080) and if anything rises. The scaled margin
does drift down, like `q^{−0.13}`, consistent with an edge-driven rather than
bulk-driven fluctuation — but the *count* does not fall. **Brief F's step 5
presumes a collapse that does not happen**, because it implicitly assumed the
three families would have different limits.

**What does die: mixing genera.** A cycle needs all three largest fibres equal
(otherwise `φ` orders the triple, and `φ` is a total order), so a genus-2/genus-3
cycle needs the two `max_c α_c` distributions to overlap. Fitting

```
2g − E[max_c α_c] = Γ(1 + 1/t)·(c_g q)^{−1/t},   t = g(2g+1)/2
```

on 300 pencils per point over `q = 31…8009` gives fitted exponents `−0.2121`
(predicted `−0.2000`) and `−0.1011` (predicted `−0.0952`), with
`c₂ = 5.82·10⁻⁴` and `c₃ = 9.05·10⁻⁹`; the spreads scale the same way,
`0.826·q^{−0.192}` and `0.648·q^{−0.099}`.

> **Prediction.** Genus 2 and genus 3 stop sharing a largest fibre at
> **`q ≈ 5.9·10⁵`** (criterion: mean gap `> 3(sd₂ + sd₃)`). Beyond it no
> mixed-genus class exists and no cycle can mix those two genera, while
> same-genus cycles continue at the constant rate above.

**Test** — the probability that a random genus-2 and a random genus-3 pencil
share a largest fibre, 600 of each per point:

| `q` | 31 | 101 | 401 | 1009 | 4001 | 8009 | 20011 | 50021 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `P(max₂ = max₃)` | 0.1267 | 0.0744 | 0.0453 | 0.0318 | 0.0132 | 0.0078 | 0.0024 | 0.00030 |

a clean power law `≈ 0.127·(q/31)^{−0.82}`, extrapolating to `4·10⁻⁵` at the
predicted `q = 5.9·10⁵` — one shared class in 25000 random pairs, i.e. gone for
any pool a census can build. The two independent routes to the same `q` agree.
Mixed-genus cycles are still found at every `q` tested (5, 10, 12, 27, 16, 17,
26, 10 at `q = 31…8009`), as they must be while the classes still overlap.

**Brief E's curl prediction, confirmed.** Brief E finds `‖curl A‖/‖A‖` is a
function of the trace spread; the trace spread converges to its Sato–Tate value,
so the curl fraction should converge. Inside the largest class, 80 members:

| `q` | 31 | 101 | 211 | 401 | 1009 | 2003 | 4001 | 8009 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| trace spread | 0.9814 | 0.9737 | 0.9987 | 0.9994 | 1.0026 | 1.0023 | 1.0007 | 0.9999 |
| `‖curl‖/‖A‖` | 0.0857 | 0.0905 | 0.0860 | 0.0875 | 0.0948 | 0.0903 | 0.0993 | 0.0998 |

The spread converges to `√m₂ = 1`, the `USp(4)` value, and the curl fraction
converges to `≈ 0.09` with no trend. **The flux does not vanish in the limit** —
the same conclusion the limiting library reaches from the other side, where the
curl residual is `6.78·10⁻²` on the seed's ten genus-two measures and `0.20–0.22`
on the larger libraries.

---

## Corrections

* **To the seed of this brief.** "The `mid` matrix … is a total order, dominated
  by `α_max = 2g`: … the endpoint wins and `φ` decides. **Cycles therefore cannot
  come from mixing genera. Fix `α_max` and vary the measure.**" The numerical
  content of the seed run is correct and reproduces exactly, but this conclusion
  is backwards. Fixing `α_max` is provably the case with *no* cycles (0 among
  330, 8990, 190568, 11968 oriented triangles across four `α_max` classes in two
  libraries); mixing genera is where every cycle found lives.
* **To the seed's library.** `U(1)` (arcsine) and `CM = ½δ₀ + ½·arcsine` are
  *horizontal* Sato–Tate measures. As vertical limits — fixed family, `q → ∞` —
  they are excluded by Deligne's semisimplicity theorem. Seven of the ten
  genus-two measures in the seed table (`SU2×U1`, `SU2×CM`, `U1×CM`, `CM²`,
  `U1²`, `2·U1`, `2·CM`) are therefore not realisable; the realisable genus-two
  library is `{USp(4), SU(2)², 2·SU(2)}`.
* **To brief F's step 5, as posed.** "If no cycle is found … cycles die once
  `G ≫ q^{−1/2}`, i.e. beyond `q ≈ …`. Predict the `q` at which the census counts
  should collapse." The census counts do not collapse, because the census
  compares families with the *same* limiting measure, about which the limit says
  nothing. And the exponent is not `q^{−1/2}`: the dominant deviation is the
  edge, `q^{−1/t}`, i.e. `q^{−1/5}` at genus two and `q^{−2/21}` at genus three.
* **To brief B's addendum 1, independently.** The addendum's rule "the larger
  `m₂` precedes" is wrong; the limiting order is `m₂` *ascending*, agreeing with
  brief E's within-class regression. A third independent confirmation of the same
  sign error.
* **Generalising T2.1.** The extreme-value exponent is `1/t`, not
  `2/dim USp(2g)`; the two agree for `USp(2g)` because `t = dim/2`, and only the
  general form covers split Jacobians (`SU(2)²`: `t = 3`, exponent `1/3`, fitted
  `0.305`).

---

## Open

* **The one arithmetic input.** An explicit one-parameter family of curves over
  `F_q` whose Jacobian is isogenous to `A^k`, `k ≥ 2`, with `A` varying and
  monodromy exactly the diagonal `SL₂` (resp. `Sp(4)`) — ideally presented as a
  map `A² → A¹`. With it, the nine cycles above become a theorem: *the exchange
  comparison on curve families over `F_q` is non-transitive for all sufficiently
  large `q`, at genus profile `(4, 6, 5)`.* Without it, the realisable sub-cone is
  transitive and the answer to brief F is negative. For `k = 2` with `A` elliptic
  the relevant locus in `M₂` is classical, so this looks like a matter of writing
  down a pencil rather than of existence.

  > **Settled — see [`REPEATED_FACTOR.md`](REPEATED_FACTOR.md).** `y² = f(x)`
  > with `f` even *and* palindromic of degree `2m`, `m` odd, has `Jac ∼ J²`;
  > `m = 3` gives `y² = (x²+1)(x⁴+cx²+1)` with `Jac ∼ E_c²`,
  > `j(E_c) = 256(c+1)³/(c+2)` non-constant, and `a_c(C) = 2a_c(E)` verified on
  > 32721 fibres over 29 primes with zero exceptions. The genus-4/5/6 vertices
  > of the `1.206·10⁻²` cycle are all witnessed, so that cycle is now a theorem.
  > `4·SU(2)` and `SU(2)×3·SU(2)`, hence the other six cycles, are still open.
* **How large is "sufficiently large"?** The cycle margins are `10⁻²` in `Ψ`
  units while `sup|Ψ_f − Ψ_μ|` is `0.28` at genus two and `1.52` at genus three at
  `q = 4·10⁵`. Extrapolating the fitted edge law, genus-6 families need `q` of
  order `10²⁰` before their `Ψ` is within `10⁻²` of the limit. A Monte Carlo test
  — draw `q` samples from each of the three limiting measures and run the exact
  comparison — would locate the crossover honestly, and is the obvious next
  computation.
* **Non-connected monodromy beyond the swap.** Which convex combinations of
  library measures are actually coset measures of a compact group? The hull
  cycles (3091, margin up to `0.169`) are much wider than the product cycles, so
  after the multiplicity question this is the highest-value direction.
* **Whether same-genus transitivity is a theorem.** It is exhaustive over every
  library built here and explained by the lexicographic `(m₂, −t)` rule and the
  level lemma, but the rule has 20 exceptions in 765 and the level lemma is not a
  proof. A genuine theorem would need a total-positivity argument controlling the
  sign changes of `∫e^{τα}d(μ−ν)` together with a bound on how far the extrema of
  the three differences can separate.
* **Higher genus.** The library stops at genus 6 because `USp(12)` costs 200
  digits and half a second per `τ`. Cycle margins grow with the genus span, so
  genus 7–10 should give wider cycles and possibly a multiplicity-free one.

---

## Files

| file | what |
|---|---|
| `st_lib.py` | the measure library: Weyl integration as a Bessel determinant, products, multiplicities, atoms, mixtures |
| `validate_library.py` | quadrature, moment, mean-zero, edge-exponent and endpoint checks; writes `validation.csv` |
| `cone_search.py` | the 149-measure search by realisability tier and `α_max` class, and the seed table recomputed; writes `cone_summary.csv`, `cycles.csv`, `library.csv` |
| `symplectic_search.py` | the arithmetically decisive search: 71 symplectic measures to genus 6, the multiplicity-free sub-cone, coset mixtures; writes `symplectic_search.csv`, `symplectic_library.csv` |
| `verify_cycle.py` | the cycles on two independent grids plus polished extrema; writes `verify_cycle.csv` |
| `level_lemma.py` | the LP saying how many levels a midrange 3-cycle needs; writes `level_lemma.csv` |
| `finite_q.py` | `Ψ_f → Ψ_μ` and the two deviation exponents; writes `convergence.csv` |
| `census_test.py` | the `F_101` cycle recomputed, the edge-datum identity, the extreme-value separation prediction and its test, the curl fraction; writes `census_test.csv` |
| `psi_library.py` | the seed script quoted in the brief, kept for reference |
