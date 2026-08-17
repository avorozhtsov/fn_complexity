# T2.3 — Symmetry type, and how far the exchange matrix detects arithmetic

Scripts: `t2_3_common.py`, `t2_3_cm_families.py`, `t2_3_symmetry.py`,
`t2_3_detector_search.py`. Data: `t2_3_*.csv`.

Everything below is for `f : A² → A¹` over `F_q`, `q` an odd prime, signature
`{N_c}_{c∈F_q}` with `N_c = q − a_c` and `Σ_c N_c = q²`. `L = (q,…,q)` is the
flat signature, `X = (2q−1, (q−1)^{q−1})` the split conic, `X⁻ = ((q+1)^{q−1},1)`
the anisotropic conic.

---

## 0. Summary

1. **The established `q ≡ 2 mod 3` example is a special case of a general
   theorem.** For `f = y² − P(x)` the signature is flat — hence
   `C(f→L) = C(L→f) = 1` **exactly** — *if and only if `P` is a permutation
   polynomial of `F_q`*. Verified with no exceptions on 2093 (P, q) pairs over
   all primes `q ≤ 500`. `P = x^d` gives the known criterion
   `gcd(d, q−1) = 1`; Dickson `P = D_n` gives `gcd(n, q²−1) = 1`, a genuinely
   different congruence type (it involves `q+1`). The superelliptic version
   `f = y^r − x^d` is flat exactly when `min(gcd(r,q−1), gcd(d,q−1)) = 1`
   (4368 triples, 0 violations). This is the clean detector T2.3 was asked for, and
   it has an infinite supply of instances.

2. **The CM-by-`Z[i]` analogue behaves differently from CM-by-`Z[ζ₃]`, and the
   naive prediction is false.** The quartic-twist family (correctly set up as
   `f = x y² − x²`, or equivalently `f = y² − x⁴`) does **not** become flat at
   `q ≡ 3 mod 4`. It becomes the **split conic**, exactly. So the detector is
   `C(f→X) = C(X→f) = 1  ⟺  q ≡ 3 mod 4`, and `C(f→L) ≠ 1` there.

3. **The signature of a monomial fibration reads a gcd exactly.** For
   `f = y² − x^d` and `e = gcd(d, q−1)`:
   `Σ_c a_c = 0` and `Σ_c a_c² = (e−1)·q(q−1)`, *exactly*, for every `q` and `d`.
   The signature depends on `d` only through `e`; `e = 1 ⟺ L`, `e = 2 ⟺ X`.

4. **Symmetry type is not visible, and there is a structural reason.**
   `Σ_c a_c = 0` holds identically for *every* fibration of the affine plane.
   The first moment is exactly the Katz–Sarnak statistic that separates
   orthogonal from symplectic, so it is unavailable by construction. What the
   rate does see beyond the genus is the **second moment**, and `m₂ = 1` for
   symplectic *and* orthogonal large-monodromy families alike.

5. **What *is* visible beyond genus is the monodromy rank.** For `f = y² − P(x)`,
   `m₂ = ν(P)/q − 1` exactly, where `ν(P) = #{(x,x′) : P(x) = P(x′)}`, and
   `ν(P)/q` converges to the number of orbits of the monodromy group of `P` on
   ordered pairs. Cyclic branch maps give `m₂ = e−1`, dihedral (Dickson) ones
   `m₂ ≈ ⌊n/2⌋`, 2-transitive ones `m₂ = 1` at *every* genus. So the rate
   separates "CM / small monodromy" from "large monodromy" at fixed genus, and
   does so with an exact integer, which is a much sharper reading than the genus
   itself (see §5.3 and T2.1).

6. **Non-congruence arithmetic is detectable too.** The quadratic-twist family
   `f = P(x) y²` of a fixed curve `E : y² = P(x)` has an *exactly three-valued*
   signature and coincides with the split conic exactly when `a_E = 0`. For a
   non-CM `E` that is the set of supersingular primes — for `y² = x³+x+1` it is
   `{17, 179, 227}` below 500, which is *not* a union of residue classes mod
   ≤ 60. So exchange-matrix entries detect conditions strictly beyond congruences.

7. **Exact values other than 1 essentially do not occur.** A scan of 69 426
   ordered rates over all pairs from a 58-map pool and all primes `q ≤ 97`,
   at tolerance `1e−12` against every `p/r` (`r ≤ 12`) and every
   `log a / log b` (`a,b ≤ 24`), returns only the values `1` (signature
   coincidence), `1/2` (`C(f→const)`, again equivalent to flatness) and `0`.
   Arithmetic enters the exchange matrix through *coincidences of signatures*,
   not through exotic exact rate values.

---

## 1. What a rate is allowed to see

Two exact identities and one consequence.

**(E1) The first moment vanishes identically.**
`Σ_c a_c = q·q − Σ_c N_c = q² − q² = 0` for every `f : A² → A¹` with full image.

**(E2) The second moment is the fiber-square count.**
`Σ_c a_c² = Σ_c N_c² − q³ = Z_f(2) − q³`, so
`m₂ := q^{−2} Σ_c a_c² = Z_f(2)/q³ − 1`. By Lang–Weil, `m₂ + 1` is the number of
`F_q`-rational irreducible components of `X ×_Y X`, up to `O(q^{−1/2})`.

**(E3) For `f = y² − P(x)`**, `a_c = −Σ_x χ(P(x)+c)`, and orthogonality of `χ` gives

```
Σ_c a_c²  =  q · ( ν(P) − q ),      ν(P) = #{(x,x′) ∈ F_q² : P(x) = P(x′)}.
```

Checked exactly (integer arithmetic, zero mismatches) on 63 random `P` of
degrees 2…9 at `q = 101, 211, 401`. (T2.2 found the same identity independently
and writes `ν(P)` as `K_P`.)

The two directions of the `L`-comparison then read two different things:

| entry | closed form | reads |
|---|---|---|
| `C(L→f)` | `log q / log(max_c N_c)`, attained at `β = ∞` (T2.1) | the **largest fiber**, i.e. `max_c(−a_c)` |
| `C(f→L)` | `1 − κ·m₂/(2 q log q) + O(q^{−3/2})`, `κ = 3−2√2`, at `β* = √2−1` | the **second moment** |

Both equal `1` exactly iff the signature is flat: `log Z_f − log Z_L` is convex,
vanishes at `β = 1` and is `≤ 0` at `β = 0`, so it is `≤ 0` on `[0,1]` and `≥ 0`
on `[1,∞)`, with equality throughout only for a flat full-image signature.

---

## 2. The flatness theorem: permutation polynomials

> **Theorem.** Let `q` be an odd prime and `P ∈ F_q[x]`. For
> `f(x,y) = y² − P(x)`, the signature of `f` is flat — equivalently
> `C(f→L) = C(L→f) = 1` exactly — **iff `P` is a permutation polynomial of `F_q`**.

*Proof.* Flat ⟺ `Σ_c a_c² = 0` ⟺ (by E3) `ν(P) = q`. Writing `n_u` for the
number of preimages of `u`, `Σ_u n_u = q` and `Σ_u n_u² = q` with `n_u ≥ 0`
integers force `n_u ∈ {0,1}` and hence `n_u ≡ 1`. ∎

The exchange matrix entry `C(y² − P(x) → L)` being exactly `1` is therefore
*equivalent* to a classical, well-studied arithmetic condition on `(P, q)`.
Verified with **0 violations on 2093 `(P,q)` pairs** (23 branch polynomials ×
all 91 primes `11 ≤ q ≤ 500`).

Instances, all exact:

| branch map `P` | flat (`C = 1`) exactly when | note |
|---|---|---|
| `x^d` | `gcd(d, q−1) = 1` | `d = 3` recovers the known `q ≡ 2 mod 3` |
| `x^5` | `q ≢ 1 mod 5` | genus-2 fibers, CM by `Z[ζ₅]` |
| `x^7` | `q ≢ 1 mod 7` | genus-3 fibers |
| `D_n(x)` (Dickson) | `gcd(n, q²−1) = 1` | involves `q+1` as well as `q−1`; 0 violations for `n ≤ 9`, all 91 primes |
| `x³ + a x`, `a = 1, 2` | never, for any `q ≤ 500` | negative control: the signature is never flat |

The Dickson row is the qualitatively new one: `D_5` is flat exactly for
`q ≡ ±2 mod 5`, i.e. when `q` has order 4 in `(Z/5)^×` — a condition on
`q mod 5` that the monomial families cannot produce, because they only ever
see `q−1`.

**Superelliptic extension.** The same argument in the `y`-variable gives, for
`f = y^r − x^d`,

> the signature is flat ⟺ `min( gcd(r, q−1), gcd(d, q−1) ) = 1`,

because if either monomial permutes `F_q` the fiber count is `q` for every `c`.
Verified on **4368 triples `(r ≤ 7, d ≤ 9, q ≤ 500)` with 0 violations.**
So `y³ − x⁵` is flat exactly for `q ≢ 1 mod 15`, and the corresponding
exchange-matrix entry detects the simultaneous splitting of `q` in `Z[ζ₃]` and
`Z[ζ₅]`.

---

## 3. CM families

### 3.1 The three-line summary

| family | model as `A² → A¹` | exact coincidence | condition |
|---|---|---|---|
| sextic twists, CM `Z[ζ₃]` | `y² − x³` | signature `= L` | `q ≡ 2 mod 3` |
| quartic twists, CM `Z[i]` | `x y² − x²`  (= `y² − x⁴`) | signature `= X` (split conic) | `q ≡ 3 mod 4` |
| quintic, CM `Z[ζ₅]`, `g=2` | `y² − x⁵` | signature `= L` | `q ≢ 1 mod 5` |
| septic, `g=3` | `y² − x⁷` | signature `= L` | `q ≢ 1 mod 7` |

### 3.2 The quartic twists in detail — the asked-for check, and why the naive answer is wrong

The quartic twists `E_c : Y² = x³ + c x` are not the fibers of any polynomial
`y² − P(x)`. The correct model is `f(x,y) = x y² − x²`: over `c ≠ 0` the
substitution `Y = xy` is a bijection of `{x ≠ 0}` and turns the fiber into
`{Y² = x³ + cx, x ≠ 0}`, i.e. `E_c` minus its rational 2-torsion point. This map
has the *same signature* as the quartic model `y² − x⁴` at every `q` tested
(all six probe rates agree to `0.000e+00`), which is the consistency check that
the model is right.

At `q ≡ 3 mod 4` all these twists are supersingular, so `N_c = q − 1` for
`c ≠ 0`; but the degenerate fiber over `c = 0` is a *pair of lines*, of size
`2q − 1`, not a cusp of size `q`. Hence the signature is the split conic, not
the flat one:

```
q ≡ 3 (mod 4):   sig(x y² − x²) = (2q−1, (q−1)^{q−1}) = sig(xy)     exactly.
```

| `q` | `q mod 4` | signature | `m₂` | `C(f→L)` | `C(L→f)` | `C(f→X)` | `C(X→f)` |
|---:|---:|---|---:|---:|---:|---:|---:|
| 101 | 1 | 5-valued | 2.9703 | 0.999500191602 | 0.870234805815 | 0.999623820678 | 0.994680237238 |
| 211 | 3 | **X** | 0.9953 | 0.999948757006 | 0.885683158448 | **1.000000000000** | **1.000000000000** |
| 401 | 1 | 5-valued | 2.9925 | 0.999904207560 | 0.896513025736 | 0.999928229957 | 0.998387630102 |
| 421 | 1 | 5-valued | 2.9929 | 0.999909618812 | 0.897253039052 | 0.999932317198 | 0.998468345505 |
| 461 | 1 | 5-valued | 2.9935 | 0.999918671657 | 0.898605809779 | 0.999939094580 | 0.998583595706 |
| 491 | 3 | **X** | 0.9980 | 0.999980971257 | 0.899525143151 | **1.000000000000** | **1.000000000000** |
| 601 | 1 | 5-valued | 2.9950 | 0.999940212748 | 0.902365859260 | 0.999955230607 | 0.998879007147 |
| 661 | 1 | 5-valued | 2.9955 | 0.999946450304 | 0.903649240075 | 0.999959905391 | 0.998971530866 |

So the requested statement "`C = 1` against `L` exactly when `q ≡ 3 mod 4`" is
**false**; the true statement is "`C = 1` against the *split conic* exactly when
`q ≡ 3 mod 4`". The difference is entirely the shape of the degenerate fiber:
`y² = x³` is a cuspidal cubic with exactly `q` affine points, `y² = x⁴` is two
parabolas with `2q − 1`.

At `q ≡ 1 mod 4` the signature is exactly 5-valued: the degenerate fiber plus
the four quartic-twist traces `±a, ±b`, each with multiplicity `(q−1)/4`, and
`a² + b² = 4q`. The signature therefore records the splitting of `q` in `Z[i]`
itself, not merely the congruence — e.g. at `q = 61` the fiber sizes are
`(121; 72,70,50,48)`, giving `a = 12, b = 10`, `144 + 100 = 244 = 4·61`.
The sextic family does the same for `4q = a² + 3b²`.

### 3.3 The exact gcd law

For `f = y² − x^d` with `e = gcd(d, q−1)`, `(E3)` and `ν(x^d) = e(q−1) + 1` give

```
Σ_c a_c = 0,     Σ_c a_c² = (e − 1) · q (q − 1)      exactly, for all q, d.
```

Verified for `d = 1…14` at `q = 101, 211, 401`: 0 mismatches. Equivalently
`m₂ = (e−1)(1 − 1/q)` — an integer up to the exact factor `1 − 1/q`. Observed:

| `d` | 2 | 3 | 4 | 5 | 6 | 7 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|
| `m₂·q/(q−1)` at `q = 421` | 1 | 2 | 3 | 4 | 5 | 6 | 9 | 11 |
| `= e − 1` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

and the whole signature depends on `d` only through `e` (237 `(q,e)` classes
over 48 primes, 0 violations of "one signature per `e`"). In particular the
genus of the fibers is *not* what `m₂` measures: `y² − x⁵` (genus 2) and
`y² − x⁶` (genus 2) have `m₂ = 4` and `5`.

### 3.4 Other CM discriminants

There is no analogue for `D = −7, −8, −11, …`. The mechanism used at `D = −3`
and `D = −4` is that `j = 0` and `j = 1728` carry extra automorphisms and hence
a *one-parameter twist family* that is itself a fibration of `A²`. Curves with
other CM discriminants are isolated points of the `j`-line; their twist families
are only quadratic, which is the `f = P(x) y²` construction of §4 — and that one
detects `a_E = 0`, not the CM order.

---

## 4. Quadratic-twist families: exact three-valued signatures, non-congruence detection

For a fixed `E : y² = P(x)` put `f(x,y) = P(x) y²`, `z = #{x : P(x) = 0}`,
`a_E = q − #{(x,y) : y² = P(x)}`. Then **exactly**

```
N_0 = q + z(q − 1),        N_c = q − z − χ(c)·a_E   for c ≠ 0,
```

so the signature has at most three distinct values, with multiplicities
`1, (q−1)/2, (q−1)/2`. Verified on 93 primes `5 ≤ q ≤ 500` for three curves:
0 violations.

Consequences:

* the signature reads `|a_E|` exactly, as half the gap between the two large
  multiplicity classes;
* it equals the split conic `X` exactly iff `z = 1` and `a_E = 0`;
* `z = 0` forces `#E(F_q)` odd, hence `a_E` odd, hence `a_E ≠ 0`: a
  quadratic-twist family is never flat.

| curve | `q` with `a_E = 0`, `q ≤ 500` | signature there | congruence |
|---|---|---|---|
| `y² = x³ + x` (CM `Z[i]`) | 7, 11, 19, …, 499 (49 primes) | `X` (`z = 1`) | `q ≡ 3 mod 4` |
| `y² = x³ + 1` (CM `Z[ζ₃]`) | 5, 11, 17, …, 491 (48 primes) | `X` (`z = 1`) | `q ≡ 2 mod 3` |
| `y² = x³ + x + 1` (non-CM) | **17, 179, 227** | `X` at 17, 179 (`z=1`); `(4q−3,(q−3)^{q−1})` at 227 (`z=3`) | **none up to mod 60** |

The last row is the point: a single entry of the exchange matrix equalling `1`
is equivalent to `q` being a **supersingular prime** for a fixed non-CM elliptic
curve — a condition of density zero (Elkies: infinitely many; Serre: density 0)
that is provably not a congruence. Arithmetic detection by the exchange matrix
is therefore not confined to splitting conditions in abelian extensions.

---

## 5. Symmetry type

### 5.1 The obstruction

`Σ_c a_c = 0` identically (§1, E1) — checked at `q = 601` across generic
hyperelliptic families of genus 1–3, both CM families, the quadratic-twist
(orthogonal) family, the Fermat pencil `x³+y³` and the split conic: `sum a_c = 0`
in every case.

In Katz–Sarnak the statistic that separates the symmetry types at the level of
moments is the *first* one: `∫ tr = 0` for `U`, `USp` and `SO(even)`, but
`∫ tr = +1` for `SO(odd)` and `−1` on the `O⁻` coset; whereas `∫ tr² = 1` for
symplectic **and** orthogonal alike. Since
`Σ_c a_c = 0` is forced for every fibration of `A²` by `Σ_c N_c = q²`, the
exchange rate never sees a non-zero first moment, and the leading term it does
see, `m₂`, does not separate `Sp` from `O`.

> **Verdict.** The exchange matrix does not distinguish orthogonal from
> symplectic symmetry, and the obstruction is structural, not numerical: it is
> the same identity `Σ_c N_c = q²` that makes the whole framework work.

### 5.2 What it does see: the monodromy rank

`m₂ = ν(P)/q − 1` exactly, and `ν(P)/q` → number of orbits of the monodromy
group `Gal(P(x) − t / F̄_q(t))` on ordered pairs (its *rank* as a permutation
group). At `q = 601`:

| `P` | monodromy | rank | `ν(P)/q` | `m₂` |
|---|---|---:|---:|---:|
| `x²` | cyclic `Z/2` | 2 | 1.9983 | 0.9983 |
| `x³` | cyclic `Z/3` | 3 | 2.9967 | 1.9967 |
| `x⁴` | cyclic `Z/4` | 4 | 3.9950 | 2.9950 |
| `x⁵` | cyclic `Z/5` | 5 | 4.9933 | 3.9933 |
| `x⁶` | cyclic `Z/6` | 6 | 5.9917 | 4.9917 |
| `D₃` | dihedral | 2 | 1.9950 | 0.9950 |
| `D₄` | dihedral | 3 | 2.9900 | 1.9900 |
| `D₅` | dihedral | 3 | 2.9834 | 1.9834 |
| `D₆` | dihedral | 4 | 3.9750 | 2.9750 |
| random deg 4 | `S₄` | 2 | 1.9651 | 0.9651 |
| random deg 5 | `S₅` | 2 | 1.9617 | 0.9617 |
| random deg 6 | `S₆` | 2 | 1.9750 | 0.9750 |

`rank = 2` ⟺ 2-transitive monodromy ⟺ `m₂ = 1`, the semicircle value, **at every
genus**. So `m₂` is a monodromy invariant that is completely independent of the
genus, and this is the "beyond the genus" content of the rate.

### 5.3 Genus and monodromy are read by *different* entries, and they disagree

At `q = 421` (chosen so that `q − 1 = 420` is divisible by 3,4,5,6,7):

| family | `g` | `m₂` | `max N_c − q` | `C(f→L)` | `C(L→f)` | `m₂` from `C(f→L)` |
|---|---:|---:|---:|---:|---:|---:|
| generic `y²=P₃+c` | 1 | 0.998 | 38 | 0.999966216299 | 0.985900364099 | 1.002 |
| generic `y²=P₅+c` | 2 | 1.055 | 59 | 0.999964319819 | 0.978756427901 | 1.058 |
| generic `y²=P₇+c` | 3 | 1.045 | 74 | 0.999964834531 | 0.973902261034 | 1.043 |
| generic `y²=P₉+c` | 4 | 0.898 | 53 | 0.999969659615 | 0.980754666769 | 0.900 |
| CM `y²=x³+c` | 1 | 1.995 | 41 | 0.999932532578 | 0.984853545020 | 2.001 |
| CM `y²=x⁴+c` | 1 | 2.993 | 420 | 0.999909618812 | 0.897253039052 | 2.680 |
| CM `y²=x⁵+c` | 2 | 3.991 | 64 | 0.999864394298 | 0.977116324372 | 4.021 |
| CM `y²=x⁷+c` | 3 | 5.986 | 113 | 0.999793991057 | 0.962142035092 | 6.109 |
| Dickson `y²=D₅+c` | 2 | 1.976 | 73 | 0.999932999975 | 0.974219787527 | 1.987 |
| twist `(x³+x+1)y²` | 1 | 0.021 | 3 | 0.999999280805 | 0.998826291835 | 0.021 |

Two readings, in opposite directions of the same reference pair, and they order
the families differently: **10 of 13 positions disagree at `q = 421`, 12 of 13
at `q = 601`**. Neither direction alone is a total order on families, and the
pair `(C(L→f), C(f→L))` is genuinely two-dimensional data.

Note the honest limitation, consistent with T2.1: `C(L→f)` reads `max_c N_c`
*exactly*, but its convergence to the Weil value `2g` is extreme-value
statistics and at `q ≤ 10³` it does not resolve genus at all (`g = 2,3,4` all
give `(1−C)√q log q ≈ 2.4–3.2` against targets 4, 6, 8). The `m₂` reading, by
contrast, is exact to 3–4 digits at `q = 421`. **At reachable `q` the monodromy
reading is far sharper than the genus reading.**

### 5.4 Matched pairs

At `q = 601`:

* Two generic genus-1 families with the **same** `ν(P)` (hence identical `m₂`
  and `m₄`) and opposite `m₃` (`±0.0326`): `C(f→L)` differs by `3.1e−08`, which
  is the `m₃` term at `O(q^{−3/2})` and is five orders above the `1e−13`
  numerical floor. So `m₃` is visible, just faint.
* Genus 1, `m₂ = 1` (generic, symplectic) vs genus 1, `m₂ = 2` (CM): `C(f→L)`
  differs by `2.2e−05`. **Monodromy at fixed genus is visible.**
* Quartic model `y² − x⁴` vs quartic-twist family `x y² − x²`: all six probe
  rates agree to `0.0e+00` — the same signature, as expected.

### 5.5 The one place symmetry type does show up: fluctuation, not value

For a quadratic-twist family with irreducible branch cubic, `m₂ = (1−1/q)a_E²/q`
exactly, and `a_E/√q = 2cos θ` equidistributes by Sato–Tate. Over the 100 primes
`q < 2000` with `z = 0` for `y² = x³+x+1`:

| bin for `m₂` | [0,.5) | [.5,1) | [1,1.5) | [1.5,2) | [2,2.5) | [2.5,3) | [3,3.5) | [3.5,4) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| observed | 54 | 12 | 11 | 9 | 7 | 4 | 2 | 1 |
| Sato–Tate | 44.1 | 16.8 | 11.9 | 9.1 | 7.0 | 5.4 | 3.8 | 2.0 |

mean `0.844` (ST prediction 1), sd `0.901` (ST prediction 1); the excess in the
first bin is the bias from conditioning on `z = 0`. Generic symplectic families
instead give `m₂ = 0.9969 ± 0.0035` — deterministic.

> So a single `q` cannot tell a twist family from a big-monodromy family: both
> have `m₂` near 1 on average. The *distribution of `1 − C(f→L)` over `q`* can:
> deterministic for symplectic, Sato–Tate-spread over `[0,4]` for the twist
> family. Symmetry type is a statement about a family of families, and it shows
> up in the exchange rate only in that form.

---

## 6. Systematic search for exact entries, `q ≤ 500`

Pool: 58 structured maps (`L`, the constant map, both conics, `y²−x^d` for
`d ≤ 12`, `y²−D_n(x)` for `n ≤ 8`, `y²−(x^d+ax)`, `y^r−x^d` for `r ≤ 5`,
`x^a+y^b`, monomial pushforwards, the quartic-twist family, three
quadratic-twist families). All 91 primes `11 ≤ q ≤ 500`.

### 6.1 Which maps have a one- or two-valued signature, and when

Every condition found is a congruence, and every one is a gcd condition
(full table in `t2_3_degenerate.csv`):

| map | flat (`= L`) | two-valued |
|---|---|---|
| `y² − x³` | `q ≡ 2 mod 3` | never |
| `y² − x⁴` | never | `q ≡ 3 mod 4` (`= X`) |
| `y² − x⁵` | `q ≢ 1 mod 5` | never |
| `y² − x⁶` | never | `q ≡ 2 mod 3` (`= X`) |
| `y² − x⁷` | `q ≢ 1 mod 7` | never |
| `y² − x¹¹` | `q ≢ 1 mod 11` | never |
| `y² − x¹²` | never | `q ≡ 11 mod 12` (`= X`) |
| `y² − D₅(x)` | `q ≡ ±2 mod 5` | never |
| `y² − D₇(x)` | `q ≡ 2,3,4,5 mod 7` | never |
| `y³ − x⁵` | `q ≢ 1 mod 15` | never |
| `x² + y²` | never | always (`X` or `X⁻`) |
| `x y² − x²` | never | `q ≡ 3 mod 4` (`= X`) |
| `(x³+x) y²` | never | `q ≡ 3 mod 4` (`= X`) |
| `(x³+1) y²` | never | `q ≡ 2 mod 3` (`= X`) |
| `(x³+x+1) y²` | never | **not a congruence** (`q = 17, 179, 227`) |

The Dickson rows are the ones a monomial family cannot produce: `D₅` is flat
exactly for `q ≡ ±2 mod 5` (i.e. `gcd(5, q²−1) = 1`), which mixes `q−1` and
`q+1`.

### 6.2 Coincidences with the named references

Reading `C(f→R) = C(R→f) = 1` off the table (`t2_3_collisions.csv`):

| map | `= L` | `= X` (split conic) | `= X⁻` (anisotropic) |
|---|---|---|---|
| `y² − x³` | `q ≡ 2 mod 3` | never | never |
| `y² − x⁴`, `y² − x⁸` | never | `q ≡ 3 mod 4` | never |
| `y² − x⁶` | never | `q ≡ 2 mod 3` | never |
| `y² − x¹²` | never | `q ≡ 11 mod 12` | never |
| `y² − D₅` | `q ≡ ±2 mod 5` | never | never |
| `x² + y²` | never | `q ≡ 1 mod 4` | `q ≡ 3 mod 4` |
| `x² + y⁴`, `x⁴ + y⁴` | never | never | `q ≡ 3 mod 4` |
| `x y² − x²` | never | `q ≡ 3 mod 4` | never |
| `(x³+x+1) y²` | never | **not a congruence** | never |

`x² + y²` is the cleanest small example after the monomial ones: one map, two
different exact coincidences, alternating with `q mod 4` — the exchange matrix
reading the splitting of `q` in `Z[i]` off a binary quadratic form.

### 6.3 Coincidence classes that are not a named reference

30 further exact-coincidence classes appear. Of the 25 largest, 20 are
congruences and 5 are not. Examples:

```
q ≡ 3, 5 mod 8    y²−D₄(x)  ==  y²−D₈(x)
q ≡ 1   mod 3     x³+y³     ==  y³−x³
q ≡ 7   mod 12    y²−x¹²    ==  y²−x⁶     ==  y⁴−x⁶
q ≡ 1   mod 24    x y²−x²   ==  x²+y⁴     ==  y²−x⁴   ==  y⁴−x²
q ≡ 7, 31 mod 36  x²+y³ == x³+y⁴ == y²−x³ == y²−x⁹ == y³−x² == y³−x⁴ == y⁴−x³
q ≡ 11  mod 60    x²+y⁵ == x⁴+y⁵ == y²−x⁵ == y⁴−x⁵ == y⁵−x² == y⁵−x⁴ == y⁵−x⁶
```

All of the congruence classes are explained by the same mechanism: the
signature of `y^r − x^d` depends on `(r, d)` only through the pair
`(gcd(r,q−1), gcd(d,q−1))`, and the conditions are the congruences that
equalise those gcds. Four of the five non-congruence classes involve the branch
maps `x^d + a x` or `D₃`, whose value distribution is not a gcd condition; the
fifth is `(x³+x)y² == y⁴−x⁴` at 22 of the 43 primes `q ≡ 1 mod 4`, an
elliptic-curve coincidence rather than a congruence.

### 6.4 Exact rate values other than 1

Scanning every ordered pair from the pool at every prime `11 ≤ q ≤ 97` —
**69 426 rates** — against all `p/r` with `r ≤ 12` and all `log a / log b` with
`a, b ≤ 24`, at tolerance `1e−12`:

* 1006 hits with value `1`: signature coincidences (§6.2–6.3);
* 3084 hits with another value, but **2999 of them occur at a single small `q`**
  and are artefacts of the target list — at `q = 11`, `C(L→X⁻) = log 11/log 12`
  is a "ratio of logarithms of small integers" only because `q` itself is small;
* the persistent hits, present at every eligible `q`, are exactly two families,
  both involving the degenerate constant map `f ≡ 0` with signature `(q²)`:

```
C(f → const) = 1/2   exactly  ⟺  max_c N_c = q  ⟺  f is flat
C(const → f) = 0     exactly  for every f with at least two fibers
```

The first of these is a **third equivalent form of the flatness detector**: for
`f = y² − P(x)`, `C(f → const) = 1/2` exactly iff `P` is a permutation
polynomial. Nothing else survives.

> **Conclusion for (3).** Over `F_q` the exchange matrix takes exact recognisable
> values only at `0`, `1/2` and `1`, and each of them is a statement about the
> *shape* of a signature (empty comparison, flatness, coincidence). The
> arithmetic is carried entirely by which signatures coincide; there is no
> exchange-matrix entry equal to `log 2/log 3`-type constants for arithmetic
> reasons. The reason is structural: an exact value needs the infimum at
> `β = 0` or `β = ∞` with commensurable logarithms, and within maps `A² → A¹`
> the numbers of fibers and the largest fibers are all `q + O(√q)` or small
> multiples of `q`, which are commensurable only in the degenerate cases.

## 7. Ledger

| arithmetic property | detected? | how |
|---|---|---|
| `P` is a permutation polynomial of `F_q` | **exactly** | `C(y²−P → L) = 1` |
| `gcd(d, q−1)` for `y² − x^d` | **exactly** | signature class; `Σ a_c² = (e−1)q(q−1)` |
| `q ≡ 2 mod 3` (splitting in `Z[ζ₃]`) | **exactly** | `y²−x³` flat |
| `q ≡ 3 mod 4` (splitting in `Z[i]`) | **exactly** | quartic-twist family `= X`, **not** `= L` |
| `gcd(n, q²−1) = 1` (Dickson) | **exactly** | `y²−D_n` flat |
| `q ≡ 1 mod 4` vs `3 mod 4` via `x²+y²` | **exactly** | `= X` vs `= X⁻` |
| `q` supersingular for a fixed non-CM `E` | **exactly**, and not a congruence | `P(x)y² = X` |
| the ideal factorisation `4q = a²+b²`, `a²+3b²` | **exactly** | the 5- resp. 7-valued signature itself |
| monodromy rank of the branch map | exactly (`m₂ = ν(P)/q − 1`) | `C(f→L)` at `β* = √2−1` |
| `m₃`, `m₄` | yes, at `O(q^{−3/2})`, `O(q^{−2})` | same entry, higher order |
| genus `g` | in principle (`2g` limit), **not at reachable `q`** | `C(L→f)`; see T2.1 |
| CM order / discriminant of an isolated curve | **no** | no one-parameter family to fiber |
| Katz–Sarnak symmetry type at a single `q` | **no** | `Σ a_c = 0` is forced; `m₂ = 1` for `Sp` and `O` alike |
| symmetry type from the ensemble over `q` | yes | fluctuation law of `1 − C(f→L)`: deterministic vs Sato–Tate |

---

## Reproduce

```
python research/m_and_e_and_a_c/t2_3_cm_families.py
python research/m_and_e_and_a_c/t2_3_symmetry.py
python research/m_and_e_and_a_c/t2_3_detector_search.py 500 100
```
