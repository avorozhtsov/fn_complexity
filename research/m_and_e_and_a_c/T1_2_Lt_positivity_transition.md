# T1.2 — The L^t positivity transition between the tropical and the Hilbert end

Scripts: `t1_2_common.py`, `t1_2_part1_psd_in_t.py`,
`t1_2_part2_lt_interpolation.py`, `t1_2_part2b_window_width.py`,
`t1_2_part3_structure.py`.
Data: `t1_2_part1_lambda_min.csv`, `t1_2_part1_thresholds.csv`,
`t1_2_distance_matrices.json`, `t1_2_part2_price_limit.csv`,
`t1_2_part2_l2_spectrum.csv`, `t1_2_part2_phase.csv`,
`t1_2_part2b_window_width.csv`, `t1_2_part3_bipartite.csv`,
`t1_2_part3_disconnection_search.csv`.

Conventions: `Z_a(β) = Σ a_i^β`, `φ_a = log Z_a > 0`, `u_a = log φ_a`,
`C(a→b) = inf_β φ_a/φ_b`, `d(a,b) = −log(C(a→b)C(b→a)) = osc_β(u_a − u_b)`,
`K_t = exp(−t·d)` entrywise. Eigenvalues below `1e-10` in modulus are read as
zero. Because `λ_min` of a centred distance matrix is 1-homogeneous in `d`, all
cross-family comparisons use the scale-free **defect**
`δ(d) = −λ_min(−½JdJ)/λ_max(−½JdJ) ≥ 0`, which vanishes exactly on
negative-type matrices.

---

## Headline

* **The PSD set is always a closed ray.** `S(d) = {t>0 : K_t ⪰ 0}` came out as
  `[t*, ∞)` in every one of the 26 signature families examined and in every one
  of 6805 randomly generated non-negative-type metrics: exactly one sign change
  of `λ_min(K_t)`, never more. `S` is provably a *closed additive semigroup*
  containing a ray (§3.2), which rules out most shapes but not, by itself, a
  gap; the ray shape is **observed, not proved**.
* **`t*` has an exact closed form on a solvable family.** For the complete
  bipartite graph metric `K_{m,n}` (distance 1 across the parts, 2 inside),
  `t*(K_{m,n}) = ½·log((m−1)(n−1))` — **proved** in §3.1 and verified
  numerically on all 55 non-degenerate pairs with `2 ≤ m ≤ 8`, `m ≤ n ≤ 12`,
  to ≤ 1.01·10⁻⁹ relative error. In particular `t*` is
  *unbounded at fixed diameter*, so no bound `t* ≤ f(diam)` exists.
* **`t*` for the exchange families is predicted to a few per cent by a
  second-order formula**, `t_2 = max_{x ⟂ 1, |x|=1} 2(xᵀdx)/(xᵀd^{∘2}x)`:
  `t_2/t* = 0.98, 0.95, 0.91` for the three weakly-violating families and
  `0.54–0.64` for the strongly-violating ones (§1.3).
* **The two ends are connected by a single continuous curve `s*(t)`, but the
  hardest point is in the middle, not at the tropical end.** On the `L^t`
  bridge `D_t` (§2), `s*(t)` is `0` for tiny `t` (Hilbert end), rises to a
  maximum near `t ≈ 15` (cert13) or `t ≈ 60` (greedy25) — a factor 4 to 200
  above `t*(d)`, depending on the family — then falls back to `t*(d)` as
  `t → ∞`. Positivity is *not* monotone in the interpolation parameter.
* **The Hilbert end is a singular limit.** The `t → 0` limit of `D_t/t` is a
  squared `L²(μ)` distance, hence of negative type by construction — but the
  margin is the smallest eigenvalue of the `L²` Gram of `{u_a}`, which decays
  geometrically (factor ≈ 9 per extra signature) and hits double precision by
  `N ≈ 11`. The critical `t_c` below which negative type genuinely survives is
  estimated at `t_c ≲ 5·10⁻⁶` for the 13-signature certificate and shrinks like
  `≈ 3^{−N}`. So the "unconditional positivity" region is real but
  exponentially thin.

---

## 1. Part 1 — for which `t` is `exp(−t·d)` PSD?

`λ_min(K_t)` was scanned on 601 log-spaced `t ∈ [10⁻³, 50]` for 26 families:
the published 13-signature certificate, 14 pseudo-random families of sizes
5..40 drawn from the pool of decreasing signatures with lengths 2..6 and
entries 1..12, and 11 families grown by hill-climbing the defect over a cached
400-signature pool.

### 1.1 Thresholds

`t*` was located by scanning from above for the largest `t` with a resolvable
negative eigenvalue and then bisecting (bisecting from below cannot work:
`K_t → J` as `t → 0`, so `λ_min(K_t) → 0` and the failure drops under any
absolute tolerance). Full table in `t1_2_part1_thresholds.csv`; the violating
families:

| family | N | defect δ(d) | t* | sign changes | diam d | t*·diam | t_2 (2nd order) | t_2/t* |
|---|---|---|---|---|---|---|---|---|
| cert13    | 13 | 1.228e-03 | 0.124065 | 1 | 0.3659 | 0.0454 | 0.121775 | 0.982 |
| rand30_61 | 30 | 6.002e-04 | 0.198118 | 1 | 0.5011 | 0.0993 | 0.189071 | 0.954 |
| rand40_71 | 40 | 4.865e-04 | 0.274464 | 1 | 0.5736 | 0.1574 | 0.250582 | 0.913 |
| greedy20  | 20 | 1.247e-01 | 12.37874 | 1 | 0.0951 | 1.1767 | 7.975810 | 0.644 |
| greedy25  | 25 | 1.703e-01 | 14.74243 | 1 | 0.1097 | 1.6172 | 8.695661 | 0.590 |
| greedy30  | 30 | 1.686e-01 | 15.86237 | 1 | 0.1009 | 1.6001 | 9.004583 | 0.568 |
| greedy40  | 40 | 1.632e-01 | 16.71336 | 1 | 0.1043 | 1.7430 | 9.099169 | 0.544 |

The other 19 families have `δ(d) ≤ 5·10⁻¹⁶` (numerically zero) and no sign
change anywhere on the grid: they are of negative type and `K_t ⪰ 0` for all
`t > 0` by Schoenberg. Random draws from the pool almost never violate: 2 of the
14 random families did, both at `N ≥ 30`. Directed
search finds violations from `N = 20` upward in this pool but not at `N ≤ 16`,
consistent with the repo's earlier observation that no subfamily of size ≤ 12
violates.

### 1.2 The PSD set is a single ray

For every violating family:

* the 601-point global scan on `[10⁻³, 50]` shows exactly **one** sign change;
* a refined 4001-point scan on `[t*/2, 2t*]` also shows exactly one;
* the Schur/semigroup consistency check (200 random pairs `s, u ≥ t*` per
  family) never produced a non-PSD `K_{s+u}` — smallest observed
  `λ_min(K_{s+u})` was `+4.5·10⁻⁴` (cert13) up to `+2.4·10⁻¹` (greedy30).

So `S(d) = [t*, ∞)` in all cases. **Observed**, not proved (see §3.2 for what
*is* proved and §3.3 for the wider search).

### 1.3 What `t*` depends on

`t*` is exactly `1/c`-homogeneous: `t*(c·d) = t*(d)/c`. The scale-free product
`t*·diam` ranges over 0.045 … 1.74 across the seven violating families and is
*not* a function of `N` — it tracks the defect instead
(`δ = 1.2e-3 → t*·diam = 0.045`; `δ = 0.17 → t*·diam = 1.62`). The bipartite
family of §3.1 shows `t*·diam` is unbounded, so no clean formula in `N` and
`diam` alone can exist.

What does work is the second-order expansion. For centred `x` (`1ᵀx = 0`),
`xᵀJx = 0`, so

    xᵀ K_t x = Σ_{k≥1} (−t)^k xᵀ d^{∘k} x / k!,

and truncating after two terms shows failure for `t < 2(xᵀdx)/(xᵀd^{∘2}x)`.
Maximising over the eigenvectors of `−½JdJ` gives

    t_2 = max_{x ⟂ 1, |x| = 1}  2 (xᵀ d x) / (xᵀ d^{∘2} x),

which reproduces `t*` to 2 % when `t*·diam` is small and degrades to 46 % when
`t*·diam ≈ 1.7` — precisely as the truncation error predicts. `t_2` is a
heuristic scale, **not** a bound: the null vector at `t*` is only *nearly*
centred for the weak violators (`|⟨v,1⟩|/√N ≈ 10⁻⁴`) and clearly non-centred
for the strong ones (`≈ 10⁻²`), so the centred truncation eventually misses the
critical direction. Two bounds that *are* proved bracket `t*` but very loosely:
the Gershgorin diagonal-dominance point (29 … 147 across the table) above, and
the crude certificate bound `t·q > N(e^{t·diam} − 1 − t·diam)`
(2.1e-4 … 0.38) below.

### 1.4 Verification

Every `d` was recomputed independently as `osc_β(u_a − u_b)` on 2·10⁶ points of
`[0, 60]` together with the `β = ∞` endpoint. Agreement with the solver:
`max|d_solver − d_grid| ≤ 4.2·10⁻⁷` (mostly `< 10⁻⁹`), and the recomputed `t*`
agrees to 9–11 significant digits (e.g. cert13: `0.124064840` vs
`0.124064836`; greedy40: `16.713364198` vs `16.713364234`).

---

## 2. Part 2 — the `L^t` bridge

### 2.1 Definitions and the two exact limits

For a probability measure `μ` on the spectrum,

    G_t(i,j) = ∫ (φ_i φ_j)^{t/2} dμ                Gram of φ^{t/2}: PSD for every t
    P_t(i,j) = ( ∫ (φ_i/φ_j)^t dμ )^{1/t}          normalised price
    D_t(i,j) = log P_t(i,j) + log P_t(j,i)         soft irreversibility

`D_t(i,i) = 0` and `D_t ≥ 0` by Jensen. With `g = u_i − u_j`,
`log P_t(g) = Σ_k t^{k−1} κ_k(g)/k!` (cumulants under `μ`), and
`κ_k(−g) = (−1)^k κ_k(g)`, so **all odd cumulants cancel** in `D_t`:

    D_t(i,j) = t·κ₂(g) + t³·κ₄(g)/12 + t⁵·κ₆(g)/360 + …

Hence two exact limits (**proved**):

    D_t / t  →  Var_μ(u_i − u_j)        as t → 0    (a squared L²(μ) distance)
    D_t      →  osc_supp(u_i − u_j)     as t → ∞    (the exchange metric d)

The first is a squared Hilbert-space distance, so it is of negative type and
`exp(−s·D_t)` is PSD for every `s > 0`; the second is not of negative type. The
phase boundary is `s*(t) = inf{s : exp(−s D_t) ⪰ 0}`, with `s*(∞) = t*(d)`.

The `t³` coefficient prediction is confirmed to four digits: `max|D_t/t − E|`
divided by `t²` is constant at `1.196·10⁻⁵` over `t ∈ [0.02, 0.32]`
(`1.1973, 1.1962, 1.1962, 1.1961, 1.1960` ×10⁻⁵).

`G_t` was spot-checked at nine `t` per family: `λ_min` of the unit-diagonal
rescaling never exceeded `10⁻⁹` in modulus (worst `−4.8·10⁻¹⁰` at `t = 3·10⁴`,
pure roundoff on a matrix whose entries are within `10⁻¹⁵` of 1). PSD
unconditionally, as it must be.

`D_t` is **not** a metric for finite `t`: `max_{ijk}(D_ij − D_ik − D_kj) > 0`
throughout, decaying to 0 only as `t → ∞` where `D_t → d`. For cert13 the
triangle defect peaks at `6.8·10⁻²` near `t ≈ 21` and is `3.9·10⁻⁴` at
`t = 3·10⁴`.

### 2.2 The Laplace / Varadhan limit — verified

`P_t(a,b) → 1/C(b→a) = sup_β φ_a/φ_b`. With `μ` the push-forward of the uniform
law on `(0,1)` under `β = x/(1−x)` (so the whole spectrum `[0,∞]` is charged),
4·10⁵ nodes:

| a | b | P_100 | P_1e4 | P_1e6 | 1/C (solver) | rel. err at t = 10⁶ |
|---|---|---|---|---|---|---|
| {5,5,5,1} | {6,3,2}   | 1.2241389623 | 1.2608907533 | 1.2618439891 | 1.2618595071 | 1.23e-05 |
| {7,7,6,1} | {5,4,4}   | 1.2305325633 | 1.2608477791 | 1.2618435437 | 1.2618595071 | 1.27e-05 |
| {6,6,2,2} | {6,5,5,5} | 0.9948779323 | 1.0186102920 | 1.0190976227 | 1.0191049177 | 7.16e-06 |
| {9,1}     | {3,2,2,2} | 1.9543129999 | 1.9994011165 | 1.9999931731 | 2.0000000000 | 3.41e-06 |
| {2,2}     | {3,1}     | 1.0207995799 | 1.0354660928 | 1.0358504674 | 1.0358567600 | 6.07e-06 |

Convergence is the expected `O(log t / t)` (Laplace with a smooth interior
maximum), so `t = 10⁶` still leaves `10⁻⁵`; the *direction* and *rate* both
match the claim. The residual at `t = 10⁶` is dominated by the quadrature, not
by the principle.

### 2.3 The phase curve `s*(t)` (cert13, `μ` = full spectrum)

Reference: `t*(d) = 0.124065`; the grid's own `t → ∞` limit
`t*(osc) = 0.122537` (10⁴ quadrature nodes, `max|osc − d| = 5.7·10⁻⁵`).

| t | δ(D_t) | s*(t) | s*·max D_t | triangle defect |
|---|---|---|---|---|
| 0.030 | 5.36e-10 | — (unresolvable) | — | 1.76e-04 |
| 0.0846 | 1.25e-07 | — | — | 4.96e-04 |
| 0.169 | 8.04e-07 | — | — | 9.90e-04 |
| 0.238 | 1.81e-06 | 0.121007 | 0.000347 | 1.40e-03 |
| 0.475 | 1.05e-05 | 0.478373 | 0.002735 | 2.79e-03 |
| 0.949 | 6.87e-05 | 0.941544 | 0.010734 | 5.56e-03 |
| 1.893 | 3.21e-04 | 1.990585 | 0.045160 | 1.10e-02 |
| 3.777 | 1.25e-03 | 4.460862 | 0.199816 | 2.16e-02 |
| 7.536 | 3.76e-03 | 7.972234 | 0.684434 | 4.01e-02 |
| **15.04** | 4.50e-03 | **10.178902** | 1.525015 | 6.22e-02 |
| 30.0 | 5.70e-03 | 8.638756 | **1.895903** | 6.72e-02 |
| 59.9 | 1.01e-02 | 4.408169 | 1.207750 | 5.37e-02 |
| 119.4 | 9.99e-03 | 2.350693 | 0.730202 | 3.78e-02 |
| 336.6 | 8.27e-03 | 1.229642 | 0.419898 | 1.80e-02 |
| 948.7 | 5.20e-03 | 0.617208 | 0.219457 | 8.09e-03 |
| 2674 | 3.04e-03 | 0.329384 | 0.119118 | 3.33e-03 |
| 7536 | 2.01e-03 | 0.208556 | 0.075950 | 1.34e-03 |
| 30000 | 1.47e-03 | 0.149367 | 0.054576 | 3.88e-04 |
| ∞ | 1.228e-03 | 0.124065 | 0.045394 | 0 |

Same shape for the other two families scanned on the full spectrum
(`t1_2_part2_phase.csv`):

| family | t*(d) | first resolvable s* | peak s* | at t ≈ | s* at t = 3·10⁴ |
|---|---|---|---|---|---|
| cert13    | 0.124065 | 0.1210 at t = 0.238 | 10.18 | 15 | 0.1494 |
| rand30_61 | 0.198118 | 0.0490 at t = 0.169 | 39.37 | 30 | 0.2436 |
| greedy25  | 14.74243 | 0.4224 at t = 0.337 | 58.35 | 60 | 14.9997 |

**So the answer to "does the symmetrised object built from `P_t` lose
positivity at some finite `t`, and does that `t` relate to `t*`?" is:** the
positivity that is lost is measured by `s*(t)`, which is `0` at the `L²` end
and `t*(d)` at the tropical end; it is *finite everywhere in between* and
attains a maximum `4×` (greedy25), `82×` (cert13) or `199×` (rand30_61)
larger than `t*(d)`, at an intermediate `t` of order 15–60. The two ends are joined by a **single continuous but
non-monotone curve** — one connected transition, with the worst positivity
strictly inside the interpolation.

### 2.4 The `L²` end is singular

`s*(t) → 0` as `t → 0` is guaranteed by the exact `t → 0` limit, but the margin
is tiny. `−½J·E·J` is the doubly-centred Gram of `{u_a}` in `L²(μ)`, so its
smallest *nonzero* eigenvalue is the negative-type margin. Its relative
spectrum for nested subfamilies of cert13 (`t1_2_part2_l2_spectrum.csv`):

| N | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| λ_min/λ_max | 3.1e-05 | 6.2e-06 | 2.4e-07 | 1.9e-07 | 1.7e-08 | 6.7e-09 | 2.8e-10 | 2.2e-10 | 2.6e-13 | 1.7e-13 | 4.7e-14 | 8.0e-16 |

Geometric decay, roughly a factor 9 per added signature: the functions
`u_a = log log Z_a` are an analytic family and their `L²` Gram has
Kolmogorov-width-type spectrum. By `N ≈ 11` the margin is at double precision.

Combining with §2.1: `−½J(D_t/t)J = −½JEJ − (t²/24)·JK₄J + O(t⁴)`, so negative
type survives while `t² · ‖JK₄J‖ / 24 < λ₂(−½JEJ)`, i.e.

    t_c  ≈  sqrt( λ₂(−½JEJ) / c ),      c = δ(D_t)/t² in the asymptotic regime.

For cert13, `c = 3.19·10⁻⁵` (from `δ = 1.81·10⁻⁶` at `t = 0.238`) and
`λ₂/λ_max = 8.0·10⁻¹⁶` (an over-estimate — it sits on the machine floor), so

    t_c(cert13)  ≲  5·10⁻⁶,

and since `λ₂` decays like `9^{−N}`, `t_c ≈ 3^{−N}`. **This is an estimate, not
a measurement**: the resolvable boundary in the table is `t ≈ 0.2`, and below
`t ≈ 0.03` the defect is indistinguishable from roundoff. What is *proved* is
that `t_c > 0` (the `t → 0` limit is exactly negative type and `λ₂ > 0` when the
`u_a` are linearly independent modulo constants) and that it is far smaller
than any `t` at which the phenomenon is numerically visible.

### 2.5 Windows matter, and how much

The task asked for uniform `μ` on `[0.05, 20]` and other windows. **Every
bounded window tried destroyed the cert13 violation outright**: with
`μ` uniform on `[0.05,20]`, `[0.01,5]`, `[0.5,60]` or `[0.001,200]`, the `t→∞`
limit `osc_window` has defect `≤ 5·10⁻¹⁷`, i.e. it is of negative type, and
`s*(t) = 0` for every `t`. The reason is quantitative, not structural: the
truncation error `max|osc_window − d|` is `3.5·10⁻³` even for `[0.001,200]`,
which is larger than the violation itself (`λ_min(−½JdJ) = −5.45·10⁻⁴`).

The sweep in `β_max` (`t1_2_part2b_window_width.csv`, `β_min = 10⁻⁴`,
2·10⁵ log-spaced nodes) locates the crossover:

| β_max | 5 | 10 | 20 | 50 | 100 | 200 | 500 | 1000 | 2000 | 5000 | 20000 | 100000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **cert13** max\|osc−d\| | 1.24e-01 | 6.60e-02 | 3.36e-02 | 1.36e-02 | 6.85e-03 | 3.45e-03 | 1.41e-03 | 7.28e-04 | 3.87e-04 | 1.82e-04 | 8.00e-05 | 7.67e-05 |
| **cert13** t*(osc) | neg.type | neg.type | neg.type | neg.type | neg.type | neg.type | 0.042344 | 0.083788 | 0.104476 | 0.116877 | 0.122972 | 0.123758 |
| **rand30_61** t*(osc) | neg.type | neg.type | neg.type | neg.type | neg.type | 0.073835 | 0.141035 | 0.167458 | 0.181543 | 0.190314 | 0.195377 | 0.196697 |
| **greedy25** t*(osc) | neg.type | neg.type | 5.250106 | 11.994962 | 13.562174 | 14.173083 | 14.530385 | 14.645764 | 14.698406 | 14.729780 | 14.743006 | 14.740446 |

So the violation reappears exactly when the truncation error drops below the
size of the violation: `β_max ≳ 500` for cert13, `≳ 200` for rand30_61, `≳ 20`
for the strong greedy25 violator, and `t*(osc) ↗ t*(d)` from below thereafter.
**Practical consequence: any Weil-style or spectral-window computation on these
objects must reach `β` of order 10³ before it can see the failure of negative
type at all.** The full-spectrum measure `β = x/(1−x)`, `x ~ U(0,1)` is the
right default and is what §2.2–2.3 use.

---

## 3. Part 3 — structure

### 3.1 An exactly solvable family (proved)

Let `d` be the graph metric of the complete bipartite graph `K_{m,n}`:
`d = 1` across the parts, `d = 2` inside a part, `d = 0` on the diagonal.

**Claim.** `exp(−t·d) ⪰ 0 ⟺ e^{−2t}(m−1)(n−1) ≤ 1`, i.e.

    S(K_{m,n}) = [ ½ log((m−1)(n−1)) , ∞ ),        t* = ½ log((m−1)(n−1)).

*Proof.* Put `p = e^{−2t}`, `q = e^{−t}`, `A` = the `m`-part, `B` = the
`n`-part. Write `x = a·1_A + u₀ ⊕ b·1_B + v₀` with `1ᵀu₀ = 1ᵀv₀ = 0`. Since
`Σ_{i≠j∈A} x_i x_j = (Σ_A x)² − |x_A|²`,

    xᵀK_t x = (1−p)(|x_A|² + |x_B|²) + p(m²a² + n²b²) + 2q(ma)(nb).

The centred parts enter only through `|x_A|² = ma² + |u₀|²` with `1 − p > 0`, so
the minimum is at `u₀ = v₀ = 0`. Setting `A' = ma²`, `B' = nb²`, the form
becomes `(1−p+pm)A' + (1−p+pn)B' − 2q√(mn)·√(A'B')`, which is `≥ 0` for all
`a, b` iff `(1−p+pm)(1−p+pn) ≥ q²mn = p·mn`. With `r = 1−p`,

    (r+pm)(r+pn) − pmn = r² + rp(m+n) − mnpr = r[ r − p(mn−m−n) ] ≥ 0
    ⟺ 1 − p ≥ p(mn−m−n) ⟺ p·(m−1)(n−1) ≤ 1. ∎

Numerics (`t1_2_part3_bipartite.csv`, all `2 ≤ m ≤ 8`, `m ≤ n ≤ 12`, 56 pairs,
55 of them non-degenerate): relative error `5.21·10⁻¹¹ … 1.01·10⁻⁹`, and
**exactly one sign change** of `λ_min(K_t)` on a 1201-point log grid over
`[10⁻⁴, 3·10⁴]` in every non-degenerate case (`K_{2,2}`: zero sign changes).
`K_{2,2}` (the 4-cycle) is of negative type, as it must be — every metric on ≤ 4
points is `ℓ₁`-embeddable — and the formula returns `t* = ½ log 1 = 0`.

Two corollaries. (i) `t*(K_{n,n}) = log(n−1) → ∞` at fixed diameter 2, so
**`t*` is unbounded on diameter-normalised metrics** and no bound of the form
`t* ≤ f(diam)` can exist; any sharp bound must involve `N` as well, and `log N`
is the right order here since `K_{n,n}` has `N = 2n`.
(ii) For `m ≠ n` the critical direction is *not* centred (`a·m + b·n ≠ 0` at the
optimum), confirming that PSD-ness of `K_t` at a fixed `t` is strictly stronger
than negative type of any linearisation. This is visible numerically:
`|⟨v,1⟩|/√N` at `t*` is `0.029` for `K_{2,3}` and `0.143` for `K_{3,7}`, but
exactly `0` for `K_{3,3}` and `K_{5,5}`.

### 3.2 What is proved about the shape of `S(d)`

Let `S(d) = {t > 0 : K_t ⪰ 0}`.

1. **`S` is closed.** `t ↦ K_t` is continuous and the PSD cone is closed.
2. **`S` is an additive semigroup.** `exp(−s·d) ∘ exp(−u·d) = exp(−(s+u)·d)`
   entrywise, and the Schur product of two PSD matrices is PSD. So
   `s, u ∈ S ⟹ s + u ∈ S`.
3. **`S` contains a ray.** `K_t` is strictly diagonally dominant as soon as
   `max_i Σ_{j≠i} e^{−t·d_ij} ≤ 1`; that holds for all large `t` (bisected value
   in `t1_2_part1_thresholds.csv`, 10.67 … 146.83 across the families).
4. **If `d` is not of negative type then `t* = inf S > 0`,** and `t* ∈ S`.
   Failure at small `t` follows from the certificate: with `x` centred,
   `xᵀK_t x = −t·xᵀdx + Σ_{k≥2}(−t)^k xᵀd^{∘k}x/k!`, and
   `|xᵀ d^{∘k} x| ≤ N·diam^k`, so `K_t` is not PSD whenever
   `t·(xᵀdx) > N(e^{t·diam} − 1 − t·diam)`.
5. **Consequences of 1–4 for gaps.** `t* ∈ S ⟹ k·t* ∈ S` for all `k ∈ ℕ`.
   If `S ⊇ [α, β]` with `β > α` then `S ⊇ [nα, nβ]` for all `n`, and these
   overlap once `n ≥ α/(β−α)`, so `S ⊇ [Nα, ∞)`: **`S` can have only finitely
   many gaps, all below `⌈α/(β−α)⌉·t*`.**

Points 1–5 do **not** force `S = [t*, ∞)`: `{1} ∪ [2, ∞)` is a closed additive
semigroup containing a ray. So "the PSD set is an interval" remains a
conjecture, supported by §3.3.

### 3.3 Search for a disconnected `S` — none found

`λ_min(K_t)` on 1201 log-spaced `t ∈ [10⁻⁴, 3·10⁴]` (diameter-normalised), with
sign changes counted (`t1_2_part3_disconnection_search.csv`):

| ensemble | trials | non-negative-type | with >1 sign change | max sign changes |
|---|---|---|---|---|
| random metric, N = 5 (shortest paths of a random weighted `K_5`) | 4000 | 14 | 0 | 1 |
| random metric, N = 6 | 4000 | 134 | 0 | 1 |
| random metric, N = 8 | 3000 | 803 | 0 | 1 |
| random metric, N = 12 | 2000 | 1865 | 0 | 1 |
| two `K_{m,n}` violators at separated scales, glued at a large gap | 4000 | 3989 | 0 | 1 |

The last ensemble is the adversarial one: it superposes two thresholds that
differ by up to `50×`, which is where a gap would be expected if one exists.
6805 non-negative-type metrics, every one with `S` a single ray.

### 3.4 Conjecture

> For every finite pseudometric `d`, `{t > 0 : exp(−t·d) ⪰ 0}` is a closed ray
> `[t*, ∞)`; equivalently `λ_min(exp(−t·d))` changes sign exactly once.

Proved for `K_{m,n}` (§3.1) and for negative-type `d` (Schoenberg, `t* = 0`);
observed in 26 exchange families and 6805 random metrics; constrained but not
implied by the semigroup structure of §3.2. A proof would presumably come from
showing that the semigroup `S` is *divisible-below-the-ray* — that `t ∈ S`
forces `t/2 ∈ S` or `[t, ∞) ⊆ S` — which is exactly what the Schur argument
does *not* give.

---

## 4. Answers to the three questions as posed

1. **Is the PSD set always an interval `[t*, ∞)`? Is it ever disconnected?**
   Always an interval in everything tested (26 signature families, 55 bipartite
   metrics, 6805 random metrics); never disconnected. Proved to be a closed
   additive semigroup containing a ray, which permits at most finitely many
   gaps but does not exclude them. `t*` values: `0.124` (cert13), `0.198`
   (rand30_61), `0.274` (rand40_71), `12.4–16.7` (the greedy families);
   19 of the 26 families are of negative type and have `S = (0, ∞)`.
   `t*` depends on the family through the *defect*, not through `N`; it scales
   as `1/c` under `d ↦ c·d`; it is unbounded at fixed diameter (`K_{n,n}`);
   and it is captured to a few per cent by the second-order formula `t_2`
   whenever `t*·diam ≪ 1`.

2. **The `L^t` family.** `P_t → 1/C` verified to `10⁻⁵` relative at `t = 10⁶`
   with the expected `O(log t/t)` rate. `G_t` is PSD for every `t`
   (`|λ_min| ≤ 10⁻⁹`, roundoff). The symmetrised object `D_t = log P_t(i,j) +
   log P_t(j,i)` loses positivity in the sense that `s*(t)` becomes positive;
   the *proved* structure is `D_t = t·κ₂ + t³·κ₄/12 + …`, so `D_t/t → ` the
   squared `L²` distance (negative type, `s* = 0`) and `D_t → d` (`s* = t*`).
   The explicit boundary is the curve `s*(t)` tabulated in §2.3: it is `0` only
   for `t ≲ t_c ≈ 5·10⁻⁶` (estimated), positive and rising for
   `t ≳ 0.2`, maximal near `t ≈ 15–60` where it exceeds `t*(d)` by a factor
   `4` (greedy25), `82` (cert13) or `199` (rand30_61), and decreasing to
   `t*(d)` thereafter.

3. **Clean structural findings.** (a) The exact formula
   `t*(K_{m,n}) = ½ log((m−1)(n−1))`, with proof, plus the corollary that `t*`
   is unbounded at fixed diameter. (b) `S(d)` is a closed additive semigroup
   containing a ray, hence has finitely many gaps — the strongest general
   statement I could prove towards "the PSD set is an interval". (c) The
   cumulant identity `D_t = Σ_{k even} 2 t^{k−1} κ_k/k!`, which is why the
   bridge starts at a *squared* `L²` distance rather than at an `L¹` one.
   (d) The `L²` Gram of `{u_a}` has geometrically decaying spectrum, so the
   Hilbert end of the bridge is a singular limit and its negative-type margin
   is exponentially small in `N`.

## 5. Proved vs observed

**Proved.** The two limits of `D_t` and the cumulant expansion; negative type
of the `t → 0` limit; unconditional PSD-ness of `G_t`; the semigroup structure,
closedness and ray-containment of `S`; the Gershgorin upper bound and the
certificate lower bound on `t*`; `t*(K_{m,n}) = ½ log((m−1)(n−1))`;
`t*(c·d) = t*(d)/c`.

**Observed (numerics, tolerance `10⁻¹⁰`, distances certified against an
independent 2·10⁶-point grid to `≤ 4.2·10⁻⁷`).** That `S` is always a ray;
the `t*` values and their relation to the defect; the accuracy of `t_2`; the
shape and the peak of `s*(t)`; the geometric decay of the `L²` Gram spectrum;
the `β_max` thresholds at which windows recover the violation.

**Estimated, below numerical resolution.** `t_c`, the true upper end of the
unconditional-positivity region in `t` (`≲ 5·10⁻⁶` for cert13, `≈ 3^{−N}`
in general). The estimate rests on the measured `t²` law for the defect and on
an `L²` Gram eigenvalue that sits on the double-precision floor; higher
precision would be needed to confirm it.
