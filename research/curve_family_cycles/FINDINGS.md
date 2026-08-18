# Findings — do arithmetic families of curves ever cycle?

Answer to session brief B. **Yes, and abundantly.** Three pencils of genus-two
curves over `F_11`, every fiber of every one of them a smooth curve, beat each
other in a circle; the cycle is proved by interval arithmetic in `certify.py`.
Over `F_11` the search is exhaustive: among the 296 fiber signatures realised by
genus-two pencils `y² = P(x) + c` with `P` monic of degree 5 or 6 there are
exactly **132 strict three-cycles**, and over `F_13` there are **1475**.

The brief expected the opposite, and the reason it expected the opposite is
worth recording, because the argument inverts cleanly:

> For curve families `C(L→f)` is always the `β = ∞` endpoint and the signatures
> are close to flat, which is exactly the regime the theorem forbids cycles in.

Flatness does not put these families *into* the endpoint regime. It **destroys
the endpoint data**, which is the opposite thing. Every `f : A² → A¹` that is
onto has exactly `q` fibers, so the `β = 0` endpoint is a tie for the whole pool
and the index `φ = log(#fibers)·log(max fiber)` degenerates to a function of the
single integer `max_c N_c`, which Weil confines to `[q, q + 2g√q]`. A pool of
thousands of families therefore falls into `O(g√q)` `φ`-classes — measured: 49,
75 and 117 classes for pools of about 6000 mixed-genus signatures at
`q = 31, 101, 211`. Inside a class **both** endpoints are exact ties, so the
endpoint-regime theorem has no content there and every strict comparison is
decided by an interior tangency. The flat regime is the most favourable one for
cycles, not the least.

---

## Notation

Conventions of `research/m_and_e_and_a_c/PLAN.md` throughout: `Z_a(β) = Σ a_i^β`,
`C(g→f) = inf_β log Z_g/log Z_f` with the implementer written first, `a ≺ b` iff
`C(a→b) < C(b→a)`, and `φ(a) = log(#fibers)·log(max fiber)`.

Fiber signatures come from maps `f : A² → A¹` over `F_q`, so `Σ_c N_c = q²`; only
onto maps are used, so every signature has exactly `q` entries.

---

## The comparison is a midrange, exactly

Put

```
u_a(β) = log log Z_a(β)          (the isometry of PLAN.md)
Δ(β)   = u_a(β) − u_b(β)
```

Then `C(a→b) = inf_β exp Δ = exp(min Δ)` and `C(b→a) = inf_β exp(−Δ) =
exp(−max Δ)`, so

```
d(a,b) = −log(C(a→b)·C(b→a)) = max Δ − min Δ  =  osc Δ        (known)
a ≺ b  ⟺  min Δ + max Δ < 0                   =  mid Δ < 0    (used here)
```

Two exact consequences.

**The endpoint-regime theorem, restated.** The two endpoint values are
`Δ(0) = log(x_a/x_b)` and `Δ(∞) = log(y_a/y_b)` with `x = log(#fibers)` and
`y = log(max fiber)`. If both extremes of `Δ` are attained at the endpoints then
`mid Δ = ½ log(x_a y_a / x_b y_b)`, whose sign is the sign of `φ(a) − φ(b)`.
That is exactly the theorem of `finite_field_exchange_matrix.md`, and it makes
the mechanism visible: **`φ` reads the two ends of `Δ`; `≺` reads the midrange of
the whole curve.** A `φ`-violation is an interior excursion deeper than the
endpoint values.

**Why three of them can close a loop.** `Δ_ab + Δ_bc + Δ_ca = 0` pointwise, and
`mid` is odd and positively homogeneous but *not* additive, so three functions
summing to zero can each have a negative midrange. The model is
`(−2,1,1), (1,−2,1), (1,1,−2)`: each pair dips deeply at one temperature and
rises mildly at the other two. That is precisely what the picture of the
certified cycle shows (`curve_family_cycle_f11.svg`).

The margin follows: `C(b→a) − C(a→b) = e^{−max Δ} − e^{min Δ} ≈ −2·mid Δ`, which
reproduces the three certified margins to two digits: `−2 mid` is
`4.73, 1.75, 2.14 ·10⁻³` against computed `4.70, 1.71, 2.10 ·10⁻³`.

---

## The certified cycle

Three pencils over `F_11`, each `P(x) + c` squarefree for all eleven `c`, so
every fiber is a smooth genus-two curve:

```
A:  y² = x⁵ +  3x⁴ + 4x³ + x² +  x       + c
B:  y² = x⁶ +  9x⁵ + 7x⁴ + 2x³ +  x      + c
C:  y² = x⁶ + 10x⁵ + 8x⁴ + 8x³ + x² + 2x + c
```

| | signature | `φ` | `m₂` | max fiber |
|---|---|---:|---:|---:|
| `A` | `{18,16,15,15,14,12,9,6,6,5,5}` | 6.930809 | 2.000000 | 18 |
| `B` | `{18,18,14,13,12,9,9,9,8,7,4}` | 6.930809 | 1.636364 | 18 |
| `C` | `{19,14,12,11,11,10,10,10,9,9,6}` | 7.060456 | 0.909091 | 19 |

`A ≺ B ≺ C ≺ A`:

| edge | `C(a→b)` | contact | `C(b→a)` | contact | margin | vs `φ` |
|---|---:|---|---:|---|---:|---|
| `A ≺ B` | `0.990213498322` | `β* = 16.508` | `0.994908823027` | `β* = 2.956` | `4.70·10⁻³` | **blind** (`φ` equal) |
| `B ≺ C` | `0.981637513410` | `β = ∞` | `0.983352120195` | `β* = 4.455` | `1.71·10⁻³` | consistent |
| `C ≺ A` | `0.979537663762` | `β* = 3.830` | `0.981637513410` | `β = ∞` | `2.10·10⁻³` | **violating** |

Four of the six rates are attained at an interior temperature. The two endpoint
contacts are the two rates *into* `C`, both equal to `log 18 / log 19` because
`A` and `B` share the largest fiber 18; and the violating edge is exactly the one
where an interior minimum, `C(C→A) = 0.9795` at `β* = 3.83`, dips below that
endpoint value `0.9816` which is all `φ` reads. That is the pattern the companion
paper predicts, with one addition specific to this regime: a second edge is not
*violating* but **blind**, because `φ` is exactly equal on `A` and `B`.

`certify.py` proves all three edges with interval arithmetic (branch and bound
on `[0,60]` plus a closed-form tail), in under a second:

```
edge A < B   C(A->B) <= 0.990213500577223071   C(B->A) > 0.99027851653160887
edge B < C   C(B->C) <= 0.981637513409912134   C(C->B) > 0.981638969617049568
edge C < A   C(C->A) <= 0.979537680325369929   C(A->C) > 0.979543642232955923
cycle certified
```

**Conditioning, worth propagating.** A naive interval enclosure of
`log Z = log Σ N_c^β` has width `≈ w log q` on a box of width `w`, and that width
enters the quotient twice, which on these nearly flat signatures forces boxes
around `10⁻⁵`. Factoring the common `q^β` out first,
`log Z = β log q + log S`, `S(β) = Σ_c (N_c/q)^β`, lets the two `β log q` terms
cancel *before* the interval division, and `log S` moves an order of magnitude
more slowly. The branch and bound then closes in a few hundred boxes.

The widest cycle found at `F_11` is slightly wider still — minimum margin
`1.843·10⁻³` — but one of its three pencils has a singular member, so the
all-smooth cycle above is the one reported.

---

## The census

### Exhaustive: genus-two pencils `y² = P(x) + c`, `P` monic of degree 5 or 6

`P(0)` may be set to zero because shifting it only permutes the fibers, so the
enumeration is over `q⁴ + q⁵` polynomials and is complete.

| `q` | signatures | strict pairs | ties | `φ`-blind strict | `φ`-violating | interior contacts | **3-cycles** | all-smooth |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 11 | 296 | 43660 | 0 | 7399 | 244 | 87.8% | **132** | 28 |
| 13 | 698 | 243253 | 0 | 39648 | 936 | 85.7% | **1475** | 176 |

Edge patterns of the cycles, against `φ`:

| pattern | `q = 11` | `q = 13` |
|---|---:|---:|
| blind, blind, blind | 73 | 450 |
| blind, consistent, violating | 54 | 1025 |
| consistent, violating, violating | 5 | 0 |
| consistent, consistent, consistent | **0** | **0** |

No cycle is `φ`-consistent throughout, as the lemma of the companion paper
requires; what is new is the first row, where `φ` is *constant* on the triple and
supplies no information at all. Those cycles live entirely inside one `φ`-class,
and inside a class every rate is forced into the interior: both endpoints give
the ratio `1` exactly, so a rate there is either `1` or strictly below it, and
over the 14798 same-class pairs at `q = 11` and 79296 at `q = 13` **not one rate
equals `1`** and not one contact is at an endpoint.

The grid engine used for the census was checked against the package solver on
300 random pairs per field: maximum deviation `1.1·10⁻⁹` at `q = 11` and
`8.5·10⁻¹⁰` at `q = 13`, against margins of `10⁻³`.

### Sampled: mixed families at larger `q`

Hyperelliptic and superelliptic pencils `y^r = P(x) + c` of several degrees,
quadratic twist families `P(x)y²`, additive maps `P(x) + Q(y)` and dense
bivariate polynomials, searched inside each class of equal largest fiber.

| `q` | signatures | `φ`-classes | 3-cycles in the four largest classes | widest margin |
|---:|---:|---:|---|---:|
| 31 | 5960 | 49 | 80, 166, 147, 58 | `3.111·10⁻⁴` |
| 101 | 5946 | 75 | 221, 106, 175, 411 | `1.525·10⁻⁴` |
| 211 | 6250 | 117 | 309, 248, 267, 320 | `1.259·10⁻⁴` |

(200 signatures analysed per class, so these counts are lower bounds.) The
widest cycle at `q = 211`
mixes three genera — pencils of degree 5, 8 and 9 — and the one at `q = 31` mixes
degrees 5, 5 and 7.

Cycles persist at every `q` tried. The widest margins fall `3.1 → 1.5 →
1.3·10⁻⁴` while `1/(√q log q)` falls `5.2 → 2.2 → 1.3·10⁻²`; the two are
consistent in order but the comparison is loose, because a maximum over a fixed
200-signature sample is not a scale-invariant statistic. The median-based table
in the next section is the clean measurement of the exponent.

---

## Why the regime permits it, quantitatively

This is the question the brief asked to answer if no cycle were found. It has an
answer anyway, and the answer explains the abundance.

Write `a_c = q − N_c`, `α_c = −a_c/√q` and

```
Λ_f(β) = log( (1/q) Σ_c (N_c/q)^β ),      log Z_f(β) = (1+β) log q + Λ_f(β)
```

exactly. Substituting `β = τ√q` turns `Λ` into the cumulant generating function
of the normalised traces, and with `Ψ_f(τ) = Λ_f(τ)/τ`,

```
C(u→v) = 1 + inf_β (Λ_u − Λ_v)/((1+β) log q) + …
       = 1 + inf_τ (Ψ_u − Ψ_v)/(√q log q)    + …
```

`Ψ_f` increases from `Ψ(0) = 0` — the first moment `Σ_c a_c` vanishes
identically — to `Ψ(∞) = α_max = (max_c N_c − q)/√q`, which is the only value
`φ` reads. **Two scales decide everything:**

* the endpoint gap between neighbouring `φ`-classes is
  `log((N+1)/N) ≈ 1/q`, so in `u`-units it is `≈ ΔN/(q log q)`;
* the interior excursion of `Ψ_u − Ψ_v` is `Θ(1)`, so in `u`-units it is
  `Θ(1/(√q log q))`.

The interior beats the endpoint by a factor `√q/ΔN`. Measured on mixed pools:

| `q` | median contact `β*` | `β*/√q` | median `\|1−C\|` | `× √q log q` | midrange law, rel. error |
|---:|---:|---:|---:|---:|---:|
| 31 | 7.53 | 1.353 | 1.54·10⁻³ | 0.0295 | 3.95% |
| 61 | 12.30 | 1.575 | 5.75·10⁻⁴ | 0.0185 | 2.76% |
| 101 | 16.51 | 1.643 | 6.55·10⁻⁴ | 0.0304 | 2.28% |
| 211 | 29.38 | 2.023 | 3.23·10⁻⁴ | 0.0251 | 1.58% |
| 401 | 45.34 | 2.264 | 2.19·10⁻⁴ | 0.0263 | 1.14% |

Over a thirteenfold range of `q`: `|1−C|` falls sevenfold while
`|1−C|·√q log q` stays in `0.019–0.030` with no trend, and `β*` grows sixfold
while `β*/√q` stays of order one, drifting only from 1.35 to 2.26. The contact
sits at `τ = β/√q` of order one and the interior correction is
`Θ(1/(√q log q))`. The last column is the median relative error of the leading
term against the computed margins; it falls like `1/√q`, with
`error·√q = 0.220, 0.216, 0.229, 0.230, 0.228` across the five fields.

**So the brief's exponent was off by `√q`.** The `O(1/(q log q))` correction it
quotes is the `β = O(1)` bottleneck of T2.2, where the second moment enters at
`β* = √2 − 1`. That regime is real but it is not where the comparison is decided:
the operative structure is at `β ≍ √q`, and it is `√q` larger.

The consequence is that `φ`-violations are generic rather than marginal. An
interior tangency can overturn a largest-fiber gap of `ΔN ≈ √q`:

| `q` | signatures | `φ`-violating ordered pairs | largest `ΔN` overturned | `√q` |
|---:|---:|---:|---:|---:|
| 31 | 333 | 617 | 6 | 5.6 |
| 61 | 505 | 2052 | 13 | 7.8 |
| 101 | 528 | 1778 | 8 | 10.0 |
| 211 | 732 | 3017 | 11 | 14.5 |
| 401 | 798 | 3886 | 19 | 20.0 |

(The last column is a maximum over a sample, so it is noisy; the point is the
order of magnitude, which is `√q` and not `1`.)

---

## Addendum 1, checked

`research/session_briefs/B_cycles_among_curve_families_add_1.md` was written in a
later session and asks (§6) to be verified before it steers anything. It is,
by `addendum.py`. Its diagnosis of *where* to look is right and its conclusion
about *what would be found there* is wrong.

**§4(c) is the finding, in its own words.** The addendum names the crossover
`β ~ √q` as the one region neither of its expansions covers, and says: histogram
the argmin `β` of every rate, and "if the histogram has mass anywhere other than
`{0} ∪ [0,1] ∪ {∞}`, that is the finding". Inside a class of equal largest
fiber — the only pairs that can carry a cycle — the histogram is

| `q` | rates | `β = 0` | `0 < β ≤ 1` | `1 < β < √q/3` | `√q/3 ≤ β ≤ 3√q` | `β > 3√q` | `β = ∞` | median `β/√q` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 31 | 14280 | 0 | 4738 | 462 | **6316** | 2764 | 0 | 1.376 |
| 101 | 14280 | 0 | 4213 | 329 | **5584** | 4154 | 0 | 1.671 |
| 211 | 14280 | 0 | 3546 | 378 | **5254** | 5102 | 0 | 2.055 |
| 401 | 14280 | 0 | 3813 | 367 | **4401** | 5699 | 0 | 2.292 |

Not one contact at either endpoint, and the mass at `β ≍ √q`.

**§3 is false.** The addendum derives
`φ̃(f) = M_f − ((3−2√2)/2)·m₂(f)`, `M = max_c(−a_c)`, and concludes "this is
still a total order, so it still forbids cycles". On the certified cycle
`φ̃(A) = 6.8284 < φ̃(B) = 6.8596 < φ̃(C) = 7.9220`, a chain, and the edge
`C ≺ A` contradicts it — the same edge `φ` gets wrong, and for the same reason.

**And φ̃ is anti-predictive where it adds anything.** `φ̃` differs from `φ` only
on pairs of equal largest fiber, where it reduces to *the larger `m₂` precedes*.
On exactly those pairs:

| `q` | strict pairs | `φ` wrong | `φ` blind | `φ̃` wrong | pairs `φ̃` decides that `φ` cannot | `φ̃` wrong there |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | 43660 | 244 | 7399 | 6006 | 6371 | **90.4%** |
| 13 | 243253 | 936 | 39648 | 31399 | 34184 | **89.1%** |
| 101 | 139128 | 1778 | 2475 | 3727 | 2414 | **80.7%** |

Independently confirmed with the package solver on 300 random such pairs at
`q = 11`: the larger `m₂` precedes in 14% of them. So `φ̃` is a *worse* order
than `φ` — 13.8% of strict pairs wrong against `φ`'s 0.56% at `q = 11` — because
it converts `φ`'s honest ties into confident errors. The reason is visible in the
`Ψ` picture: a signature whose `Ψ` rises faster at small `τ` (larger `m₂`) must
flatten before the common endpoint `α_max`, and the later, negative excursion is
usually the deeper one.

**No subpopulation rescues it.** Stratifying the pairs `φ̃` decides that `φ`
cannot, three ways at `q = 11` and `q = 13`: by largest fiber, every class runs
between 2.4% and 14.6% correct; by multiplicity of the maximum, 7.0% when the
multiplicities agree and 15.8% when they differ. And by `|Δm₂|` the accuracy
falls *monotonically as the kept term grows*:

| `\|Δm₂\|` quintile (`q = 13`) | pairs | `φ̃` correct |
|---|---:|---:|
| `[0.000, 0.154)` | 6504 | 25.0% |
| `[0.154, 0.308)` | 5222 | 20.4% |
| `[0.308, 0.462)` | 7163 | 11.0% |
| `[0.462, 0.615)` | 7146 | 3.2% |
| `[0.615, 3.154)` | 8149 | **0.2%** |

That is the shape of a reversed sign, not of a missing correction: `φ̃` is most
confidently wrong exactly where the term it keeps is largest. The correct
leading rule on these pairs is the opposite one, *smaller `m₂` precedes*.

**And no rule of that shape can work either — provably, not statistically.**
There are three-cycles inside a *single* largest-fiber class whose three `m₂`
are distinct: 48 of them at `q = 11`, 332 at `q = 13`. The widest at `q = 11`:

```
A:  M = 5,  m2 = 1.454545,  sigma = {16,15,14,14,12,12,10,10,9,8,1}
B:  M = 5,  m2 = 1.090909,  sigma = {16,15,14,14,13,10,10,10,7,6,6}
C:  M = 5,  m2 = 0.909091,  sigma = {16,16,13,12,11,11,11,10,9,6,6}
```

with `A ≺ B ≺ C ≺ A` at margins `9.69·10⁻⁴`, `2.24·10⁻³`, `8.78·10⁻⁴`
(package solver). On such a triple `M` is constant, so **every** scalar
`F(M, m₂)` is a function of `m₂` alone, hence a total order, hence not `≺` — for
any coefficient and either sign. The addendum's programme of adding one more
moment cannot be repaired by choosing a better constant.

**§6, first check.** `D = (C(f→g) − 1)·q log q` against
`min(−0.0858 Δm₂, ΔM)`, over pairs of 40 genus-two pencils with `Δm₂ > 0`:

| `q` | population | pairs | median `\|D\|` | median rel. error | `1/√q` | median `β*/√q` |
|---:|---|---:|---:|---:|---:|---:|
| 211 | largest fibers differ | 712 | 0.0198 | 12.4% | 0.0688 | 0.029 |
| 211 | largest fibers tie | 52 | 0.0239 | 15.4% | 0.0688 | 0.029 |
| 1009 | largest fibers differ | 738 | 0.3119 | 10.0% | 0.0315 | 0.015 |
| 1009 | largest fibers tie | 30 | 0.1183 | **94.6%** | 0.0315 | **3.306** |

Where the largest fibers differ the formula is a fair description and its error
is roughly flat in `q`. Where they tie the error *grows* with `q` and the contact
walks out from `β*` to `3.3√q`. The predicted `O(q^{−1/2})` relative accuracy
holds nowhere and fails in the wrong direction exactly where cycles live,
because the dropped crossover term grows like `√q` against the terms kept.

**§4(b) confirmed and placed.** With the largest fiber attained twice instead of
once, `(R(β) − 1)·β·log(max)` converges to `log 2 = 0.693147` — 0.693139 already
at `β = 100`. It decays like `1/β`, so at `β ≍ √q` it contributes `≈ log 2/√q`:
the same order as the crossover itself and `√q` above the moment terms. It is
one contribution to the crossover, not a separate regime.

**§7 is moot, and brief C is safe.** The fallback statement was to be "no scalar
invariant *beyond the `φ̃` truncation*". Cycles exist, `φ̃` is refuted directly,
and it is refuted on nine tenths of the pairs it was introduced to decide — so
the unqualified statement stands for curve families.

Reproduce: `python research/curve_family_cycles/addendum.py`.

---

## Traps, learned here

* **The package can report a spurious interior contact.** For a pair whose
  infimum is the `β = ∞` endpoint the ratio can be constant to double precision
  from `β ≈ 200` on, and the golden-section refinement then returns an interior
  minimiser — `C(A→C)` above was first reported at `β = 258.9` when the infimum
  is the endpoint `log 18/log 19`, confirmed at 50 digits. Test the rate against
  `log(max_g)/log(max_f)` before publishing a contact temperature.
* **The `β` horizon must scale with `q`.** These signatures have adjacent
  entries differing by one, so the largest fiber is not isolated until
  `β ≳ 36q`; the grids here run to `360q`. This is the same trap as T1.2's
  `β_max ≳ 500`, in the form it takes for curve families.
* **Interval arithmetic needs the `q^β` factored out**, see above.
* **Genus one is useless and genus two is plentiful.** The structural collision
  of T2.2 — 400 random elliptic fibrations give `gcd(4,q−1)+1` signatures —
  leaves nothing to compare. Genus two already gives 296 signatures at `q = 11`
  and 698 at `q = 13`, exhaustively.
* **The leading-order midrange law is only asymptotic.** At `q = 11` it gets the
  sign of the edge `C ≺ A` wrong — midrange `+4.2·10⁻⁴` in the exact-`Λ` form
  and `+5.7·10⁻³` in the `α`-scaled one against an exact `−1.07·10⁻³`, so both
  predict `A ≺ C` — because the term dropped from the denominator,
  `Λ_v/((1+β) log q) → log(max_v/q)/log q`, is nominally `1/(√q log q) = 0.13`
  at that field but reaches `0.23`, and
  because `log(1 + α/√q) ≠ α/√q` when `α/√q = 0.72`. The law directs the search;
  the exact `mid Δ` decides, and the interval certificate proves.

---

## What this settles, and what it does not

Settled: the exchange comparison **does** cycle among families of curves, no
scalar invariant reproduces it — not `φ`, and not the next-order `φ̃` of
addendum 1, which is wrong on nine tenths of the pairs it was introduced to
decide — the
smallest witness is over `F_11`, the phenomenon is not rare (132 cycles in a
complete enumeration of one genus at one small field), and the margins are wide:
the narrowest edge of the certified cycle is `1.7·10⁻³`, five times the
`3·10⁻⁴` of the integer three-cycle, and its widest is `4.7·10⁻³`.

Not settled: whether cycles persist for a *fixed* genus as `q → ∞`. Everything
here scales as `1/(√q log q)`, so the margins shrink; the counts stay in the
hundreds in every sample taken, but at `q ≥ 31` the search is a sample and not a
census, so that is evidence and not a theorem.

## Files

| file | what |
|---|---|
| `common.py` | signature pool and the vectorised rate engine |
| `search.py` | the census; writes `cycles.csv`, `headline_cycle.json` |
| `certify.py` | interval-arithmetic proof of the `F_11` cycle |
| `regime.py` | the two scales; writes `regime.csv` |
| `addendum.py` | the checks addendum 1 asks for; writes `addendum.csv` |
| `plot_oscillations.py` | `curve_family_cycle_f11.svg`, the picture of `mid Δ` |
