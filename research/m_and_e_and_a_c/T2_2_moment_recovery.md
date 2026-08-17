# T2.2 — What the exchange rates know about the trace data `{a_c}`

**Question.** For `f : A² → A¹` over `F_q` with fiber counts `N_c` and traces
`a_c = q − N_c`, how much of `{a_c}` is recoverable from the exchange rates of
the signature `{N_c}` against a small fixed reference family?

**Answer in one line.** The rates read the trace data through exactly two kinds
of channel: the two *endpoints*, which give `|image(f)|` and `max_c N_c` (hence
`min_c a_c`) **exactly, as integers**; and the *interior tangencies*, each of
which is one exact evaluation of `log Z_f(β) = log Σ_c N_c^β` at a pinned `β`.
The interior is not blind past `m₂` — it sees the whole moment ladder — but each
successive moment is damped by a factor `≈ 0.6/√q`. So `m₂` is recoverable to
relative accuracy `O(q^{−1/2})` from a single rate, `m₃` only through a
`q^{−1/2}`-suppressed term, and `max_c a_c` (the *smallest* fiber) is not
recoverable at all beyond what the moments leak, because the rate curve has no
endpoint at `β < 0`.

Scripts: `t2_2_common.py` (rate solver + moments), `ffmaps.py` (map families),
`t2_2_survey.py`, `t2_2_separation.py`, `t2_2_collisions.py`, `t2_2_interior.py`,
`t2_2_vector.py`, `t2_2_validate.py`.

---

## 0. Setup, conventions, and solver certification

`Z_f(β) = Σ_c N_c^β`, `C(g → f) = inf_{β ∈ [0,∞]} log Z_g(β)/log Z_f(β)`.
Normalised moments `m_k = q^{−1−k/2} Σ_c a_c^k`, so `m_k = O(1)` under
Hasse–Weil. `Σ_c N_c = q²` forces `m₁ = 0` whenever every fiber is non-empty.

Reference family (each is a genuine map `A² → A¹` with `q²` points):

| name | map | signature |
|---|---|---|
| `L` | `f = x` | `(q)^q` |
| `Xsplit` | `f = xy` | `(2q−1, (q−1)^{q−1})` |
| `Xaniso` | `f = x² − d y²`, `d` a non-residue | `((q+1)^{q−1}, 1)` |
| `Sq` | `f = x²` | `((2q)^{(q−1)/2}, q)` |

Signatures have length `~q`, so bulk work uses a vectorised
(value, multiplicity) solver. `t2_2_validate.py` certifies it against
`fn_complexity.exchange_rate_result` on 720 rate computations:

```
max |C_fast − C_exact| = 4.441e-16     (9.7x faster)
```

Every headline number below is re-certified with the repo solver.

One trap worth recording: past the saturation horizon `β ≳ 36/min gap` both
partition functions are affine and the ratio is flat *to machine precision*, so a
naive local-minimum scan reports thousands of spurious "interior" minima whose
value equals the `β = ∞` endpoint. `t2_2_common.classify` filters these; without
it, `C(f → Xsplit)` looks interior at `β ≈ 600` when it is really the endpoint.

### An exact handle on `m₂`

For `f_P(x, y) = y² − P(x)` the fiber over `c` is `y² = P(x) + c` and
`a_c = −Σ_x χ(P(x)+c)`. Orthogonality of the quadratic character gives, exactly,

```
Σ_c a_c²  =  q·K_P − q²,        K_P = #{(x, x') : P(x) = P(x')},
```

so **`m₂ = K_P/q − 1` is an exact rational invariant of `P`**. Two consequences
used throughout: `m₂` can be matched *exactly* between different maps (same
`K_P`), and the known limits fall out — `P` with big monodromy has `K_P ≈ 2q`,
so `m₂ → 1`; `P = x³` with `q ≡ 1 (3)` has `K_P = 3q − 2`, so `m₂ → 2`; with
`q ≡ 2 (3)` cubing is a bijection, `K_P = q`, `m₂ = 0`.

---

## 1. Which probes are endpoints, which are interior

Measured over the "Weil subpool" (all fibers non-empty and `max_c|a_c| ≤ 6√q`;
≈ 570 maps per `q` drawn from hyperelliptic, superelliptic, additive and dense
random bivariate families):

| probe | q = 101 | q = 211 | q = 503 | tangency `β` range (q = 211) |
|---|---|---|---|---|
| `C(L → f)` | endpoint `∞` | endpoint `∞` | endpoint `∞` | — |
| `C(f → Xsplit)` | endpoint `∞` | endpoint `∞` | endpoint `∞` | — |
| `C(Xaniso → f)` | endpoint `∞` | endpoint `∞` | endpoint `∞` | — |
| `C(Sq → f)` | endpoint `0` | endpoint `0` | endpoint `0` | — |
| `C(f → L)` | **interior** 91 % | **interior** 94 % | **interior** 90 % | `[0.4097, 0.4148]` |
| `C(Xsplit → f)` | **interior** 100 % | **interior** 100 % | **interior** 100 % | `[0.4424, 5.4446]` |
| `C(f → Xaniso)` | **interior** 91 % | **interior** 94 % | **interior** 90 % | `[0.5299, 1.6990]` |
| `C(f → Sq)` | **interior** 91 % | **interior** 94 % | **interior** 90 % | `[9.63, 41.81]` |

The percentages are over the Weil subpool; the complement in every row is
exactly the supersingular maps (`m₂ = 0`, signature literally equal to `L`),
which is why `C(Xaniso → f)` is reported interior on precisely the 9 % / 6 % /
10 % where the other rows are not.  For every genuine curve fibration the four
top rows are pure endpoints and the four bottom rows are genuine tangencies.

**Endpoint channel.** `C(L → f) = log q / log N_max` inverts to an integer:
`N_max` was recovered exactly for **1718/1718** Weil-regime maps across
`q ∈ {101, 211, 503}`. Likewise `C(Sq → f)` gives `|image(f)|`.

**Interior channel — one rate is one evaluation of `log Z_f`.** If the infimum is
attained at an interior `β₀` then `log Z_f(β₀) = C · log Z_g(β₀)` identically.
`t2_2_interior.py` confirms this to `0.0e+00` (q = 211, `f : y² = x³ + x + c`):

```
probe             C                    beta0    log Z_f(b0) from C        direct              err
C(f->L)           0.999924851379929    0.413844   7.566121768776055   7.566121768776055   0.0e+00
C(Xsplit->f)      0.999984859216071    1.531817  13.551829279030120  13.551829279030120   0.0e+00
C(f->Xaniso)      0.999955514649888    1.460264  13.168555535299504  13.168555535299504   0.0e+00
C(f->Sq)          0.900043135660118   34.887118 194.001392892496540 194.001392892496540   0.0e+00
C(f->Cb)          0.845544579741347   57.514227 317.284104050078156 317.284104050078156   0.0e+00
```

**In principle everything is recoverable — with an unbounded reference family.**
`C(f → flat(n, m)) = inf_β log Z_f(β)/(log n + β log m)` is the support function
of the epigraph of the *convex* function `log Z_f`. Sweeping `(n, m)` over
integers sweeps the lines `A + Bβ` with `log n / log m` dense in `(0, ∞)`, so the
family of flat rates reconstructs `log Z_f` on `[0, ∞)`, hence `Z_f(β)` for all
`β`, hence the multiset `{N_c}`. `t2_2_interior.py` exhibits tangencies at
`β = 0.414, 18.77, 35.07, 347.89` from four flat probes. **A small reference
family gives only finitely many tangent lines**, and that is where the losses
come from.

---

## 2. The moment ladder: what a single interior rate sees

Expanding `N_c = q(1 − a_c/q)` and using `m₁ = 0`,

```
log Z_f(β) = (β+1) log q + Σ_{k≥2} (−1)^k C(β,k) m_k q^{−k/2} + O(q^{−2}),
```

so `log Z_f/log Z_L = 1 + Σ_k (−1)^k C(β,k) m_k / ((β+1) q^{k/2} log q)`.
Minimising the leading term over `β` means minimising `C(β,2)/(β+1)`, giving

> **`β* = √2 − 1 = 0.414213562…`, and `1 − C(f → L) = κ m₂ / (2q log q) + …`
> with `κ = 3 − 2√2 = 0.171572875…`.**

Measured `argmin β` of `C(f → L)` over the Weil subpool:

| q | range of `argmin β` | distance to `β*` |
|---|---|---|
| 101 | `[0.409232, 0.415568]` | `< 5.0e-3` |
| 211 | `[0.409741, 0.414776]` | `< 4.5e-3` |
| 503 | `[0.413498, 0.414537]` | `< 7.2e-4` |

**The tangency is pinned, not free.** That is the crux: the interior contributes
essentially *one* fixed linear functional, not a curve's worth of them. Writing
`1 − C(f → L) = Σ_{k≥2} c_k m_k + O(q^{−2})` with
`c_k = (−1)^{k+1} C(β*,k) / ((β*+1) q^{k/2} log q)`, the successive ratio is

```
c_{k+1}/c_k = −(β* − k)/((k+1)√q):     0.5286/√q,  0.6464/√q,  0.7172/√q, …
```

Numerical verification (rms of `1 − C(f→L)` after subtracting the predicted
contribution of moments `2..K`; ~720 hyperelliptic maps of genus 1–6 per `q`):

| K | q = 101 | q = 211 | q = 503 |
|---|---|---|---|
| raw | 1.848e-04 | 7.681e-05 | 2.748e-05 |
| 2 | 2.993e-06 | 6.199e-07 | 8.517e-08 |
| 3 | 1.957e-06 | 3.922e-07 | 5.707e-08 |
| 4 | 2.005e-07 | 3.063e-08 | 3.672e-09 |
| 5 | 1.658e-07 | 2.779e-08 | 3.627e-09 |
| 6 | 1.192e-07 | 2.324e-08 | 3.357e-09 |
| 8 | 1.167e-07 | 2.311e-08 | 3.354e-09 |

Moments 2, 3 and 4 are all genuinely present with the predicted coefficients.
The plateau from `K = 5` on is exactly the dropped `−u²/2` term of `log(1+u)`:
predicted `C(β*,2)² m₂² / (2q²(β*+1) log q) = 3.30e-09` at `q = 503` (with
`m₂ ≈ 1`), observed plateau `3.36e-09`.

Direct least-squares fit of `1 − C(f → L)` on `(m₂, m₃)` over ~600 Weil-regime
maps (`t2_2_collisions.py`):

| q | fitted / predicted `c₂` | fitted / predicted `c₃` | rms resid, `m₂` only | rms resid, `m₂` + `m₃` |
|---|---|---|---|---|
| 101 | 1.010413 | 1.031781 | 2.90e-06 | 1.97e-06 |
| 211 | 1.004915 | 1.015203 | 5.69e-07 | 3.78e-07 |
| 503 | 1.002032 | 1.007455 | 8.30e-08 | 5.62e-08 |

Both coefficients converge to the predicted values like `O(q^{−1/2})`, and `m₃`
accounts for 56 % of the variance of the `m₂`-only residual at `q = 211`.

---

## 3. (1) Inversion: what the 8-rate vector predicts

Rate vector `(C(L→f), C(f→L), C(Xsplit→f), C(f→Xsplit), C(Xaniso→f),
C(f→Xaniso), C(Sq→f), C(f→Sq))`. Linear regression and leave-one-out 1-NN over
the Weil subpool (`t2_2_survey.py`):

**q = 211** (572 maps)

| target | spread | linear R² | max resid | 1-NN median err |
|---|---|---|---|---|
| `m₂` | 7.962 | 0.9999488 | 3.34e-02 | 9.48e-03 |
| `m₃` | 1.784 | 0.7323367 | 4.86e-01 | 5.96e-02 |
| `m₄` | 191.93 | 0.9917175 | 8.83e+00 | 1.47e-01 |
| `a_min = q − N_max` | 77 | **0.9999999** | 2.89e-02 | **0** |
| `a_max = q − N_min` | 76 | 0.8004081 | 2.50e+01 | 3 |
| `max_c abs(a_c)` | 77 | 0.9039857 | 2.42e+01 | 0 |

**q = 503** (576 maps)

| target | spread | linear R² | max resid | 1-NN median err |
|---|---|---|---|---|
| `m₂` | 1.245 | 0.9999988 | 3.95e-03 | 2.70e-03 |
| `m₃` | 0.625 | 0.9788029 | 1.54e-01 | 4.14e-02 |
| `m₄` | 4.745 | 0.9934640 | 4.17e-01 | 8.55e-02 |
| `a_min = q − N_max` | 111 | **1.0000000** | 3.34e-02 | **0** |
| `a_max = q − N_min` | 96 | 0.9298594 | 2.52e+01 | 4 |
| `max_c abs(a_c)` | 111 | 0.9608790 | 2.50e+01 | 1 |

Single-rate model inversion `m₂_hat = 2q log q (1 − C(f→L))/κ`:

| q | median rel. err | max rel. err | `q^{−1/2}` |
|---|---|---|---|
| 101 | 9.94e-03 | 6.06e-02 | 9.95e-02 |
| 211 | 4.18e-03 | 3.83e-02 | 6.88e-02 |
| 503 | 1.79e-03 | 8.40e-03 | 4.46e-02 |

Two-probe inversion for `(m₂, m₃)` (two most separated interior tangencies plus
the `1/q` expansion):

| q | median err in `m₂` | median err in `m₃` | true `m₃` spread |
|---|---|---|---|
| 101 | 4.61e-03 | 2.65e-01 | 1.441 |
| 211 | 2.25e-03 | 1.94e-01 | 1.784 |
| 503 | 9.26e-04 | 1.24e-01 | 0.625 |

**Verdict on (1).**

* `|image(f)|` and `N_max` (equivalently `min_c a_c`): recovered **exactly**, as
  integers, from the two endpoints. Nothing else is needed.
* `m₂`: recovered to relative accuracy `O(q^{−1/2})` from a *single* rate
  (`1e-3`–`1e-2` in practice), and to `~1e-3` absolute from the vector.
* `m₃`: present but damped by `0.5286/√q`. Recoverable in principle (the solver's
  `1e-13` precision leaves ample headroom) but with a conditioning penalty of
  `√q`; in practice the 2-probe inversion lands within 10–20 % of the `m₃` spread
  of the sample, and the 8-rate regression is markedly worse for `m₃` than for
  `m₂` at every `q`.
* `max_c a_c` (i.e. `min_c N_c`): **not recovered**. `R² ≈ 0.80–0.93`, 1-NN
  median error 3–4 on a spread of 76–96, while `min_c a_c` is exact. See §5.

---

## 4. (2) Separation: matched `m₂`, different `m₃`

Search over 2400 hyperelliptic maps `y² = P(x) + c`, `deg P ∈ {3,5,7}`, grouped
by the exact invariant `K_P` (`t2_2_separation.py`). `P` is listed low-to-high.

### Pair A — `q = 211`, identical `m₂` (exactly), different `m₃`

| | `P` | `K_P` | `m₂` | `m₃` | `m₄` | `N_max` | image |
|---|---|---|---|---|---|---|---|
| `f` | `(2, 24, 166, 4, 13, 1)` | 431 | `1.042654028436` | `−0.313218604` | 2.777291 | 246 | 211 |
| `g` | `(45, 182, 50, 17, 24, 70, 194, 1)` | 431 | `1.042654028436` | `+0.464205558` | 3.606235 | 246 | 211 |

| probe | `f` | `g` | difference | `β₀(f)` |
|---|---|---|---|---|
| `C(f→L)` | 0.999921302100922 | 0.999918988496203 | **2.314e-06** | 0.414564 |
| `C(Xsplit→f)` | 0.999976127589515 | 0.999974976043013 | **1.152e-06** | 1.655805 |
| `C(f→Xaniso)` | 0.999961569475610 | 0.999962735818382 | **1.166e-06** | 1.415300 |
| `C(L→f)` | 0.972122768364702 | 0.972122768364702 | 0.000e+00 | `∞` |
| `C(f→Xsplit)` | 0.911081590998706 | 0.911081590998706 | 1.110e-16 | `∞` |
| `C(Xaniso→f)` | 0.972981597876619 | 0.972981597876619 | 0.000e+00 | `∞` |

Repo-solver certification: `C(f→L) = 0.999921302100922` at `β = 0.414563800`;
`C(g→L) = 0.999918988496203` at `β = 0.412252812`.

Predicted from the ladder: `c₃ Δm₃ = 2.7645e-06 × 0.777424 = 2.149e-06` against
observed `2.314e-06` (the balance is `c₄ Δm₄`). **The interior sees `m₃`; every
endpoint probe is exactly blind to it.**

### Pair B — `q = 503`, same structure

`f: P = (178, 452, 157, 289, 70, 1)`, `g: P = (284, 488, 427, 184, 123, 269, 77, 1)`,
`K_P = 1035`, `m₂ = 1.057654075547` for both, `m₃ = −0.145730292` vs `+0.341243115`,
`N_max = 568` and image `= 503` for both.

`ΔC(f→L) = 3.272e-07`, `ΔC(Xsplit→f) = 1.611e-07`, `ΔC(f→Xaniso) = 1.661e-07`;
all three endpoint probes agree to `0.0e+00`.

### Converse — the blind direction

A single rate is a single functional of `(m₂, m₃, …)`, so it must be blind to the
direction `δm₂ = −(0.5286/√q) δm₃`. Adjacent pairs in the `C(f→L)`-ordered list
of 1797 distinct signatures at `q = 211` (`r = 0.036390`):

| difference in `C(f→L)` | `m₂(f)` | `m₂(g)` | `Δm₂` | `Δm₃` | `Δm₂ + rΔm₃` |
|---|---|---|---|---|---|
| 3.285e-11 | 1.061611374 | 1.042654028 | −1.896e-02 | +5.102e-01 | −3.903e-04 |
| 3.314e-11 | 0.995260664 | 0.985781991 | −9.479e-03 | +2.278e-01 | −1.189e-03 |
| 3.904e-11 | 0.890995261 | 0.900473934 | +9.479e-03 | −2.486e-01 | +4.315e-04 |

The combination `Δm₂ + rΔm₃` is 20–50× smaller than `Δm₂` alone on every
near-collision: the collisions do lie along the predicted blind line.
Certified example (`q = 211`):

```
f: P=(50,193,0,196,158,151,84,1)  m2=1.061611374408  m3=-0.388851248  N_max=260
g: P=(0,169,103,196,78,1)         m2=1.042654028436  m3=+0.121372209  N_max=249
exact C(f->L) = 0.999919993210414
exact C(g->L) = 0.999919993243267      difference = 3.285e-11
```

11 such pairs at `q = 211` and 53 at `q = 503` fall below the `1e-10` resolution
threshold. **But the vector separates them**: because `N_max` differs (260 vs
249), the endpoint probes differ by `7.5e-03`. A single rate is far from
injective on `(m₂, m₃)`; the vector is much better.

---

## 5. (3) The sharpest question: does the interior add anything to the endpoints?

**Yes — and the cleanest witness also exposes what the rates cannot see.**

Search over 1600 full-image hyperelliptic maps at `q = 211`, grouped by the two
endpoint readouts `(image size, N_max)`; 80 734 admissible pairs; report the pair
minimising the sup-norm distance between the *full 8-probe rate vectors*
(`t2_2_vector.py`).

| | `P` | `m₂` | `m₃` | `m₄` | `N_max` | `N_min` | image |
|---|---|---|---|---|---|---|---|
| `f` | `(151,33,27,120,193,35,61,108,199,1)` | `0.966824644550` | `+0.029280744` | 2.351921 | 248 | **173** | 211 |
| `g` | `(200,202,15,180,196,178,98,39,4,1)` | `0.966824644550` | `+0.094949028` | 2.764025 | 248 | **167** | 211 |

Identical image size, identical largest fiber, identical `m₂` (exactly), but the
*smallest* fiber differs by 6, i.e. `max_c a_c` differs by 6 out of ≈ 38.

| probe | `f` | `g` | difference |
|---|---|---|---|
| `C(f→L)` | 0.999926158704570 | 0.999925916828638 | 2.419e-07 |
| `C(Xsplit→f)` | 0.999987373634344 | 0.999987255783059 | 1.179e-07 |
| `C(f→Xaniso)` | 0.999953240785310 | 0.999953358172294 | 1.174e-07 |
| `C(f→Sq)` | 0.900685396874103 | 0.900685413125429 | 1.625e-08 |
| `C(L→f)` | 0.970695075582268 | 0.970695075582268 | **0.000e+00** |
| `C(f→Xsplit)` | 0.912421604607920 | 0.912421604607920 | **0.000e+00** |
| `C(Xaniso→f)` | 0.971552643787758 | 0.971552643787758 | **0.000e+00** |
| `C(Sq→f)` | 0.871368219748219 | 0.871368219748219 | **0.000e+00** |

Certified: `C(f→L) = 0.999926158704570` at `β = 0.413664037`,
`C(g→L) = 0.999925916828638` at `β = 0.413373844`.

Reading of this example:

1. **The interior is not redundant.** With both endpoints and `m₂` matched
   exactly, four interior probes still separate, by `2.4e-07` — four orders of
   magnitude above the `1e-10` resolution floor.
2. **What separates them is `m₃` and `m₄`, not the extreme trace.** Predicted
   `c₃Δm₃ + c₄Δm₄ = 1.816e-07 + 0.507e-07 = 2.32e-07`; observed `2.419e-07`
   (4 % agreement). The 6-unit change in the *smallest* fiber enters only
   through the moments.
3. **The rate curve is structurally asymmetric between the two extreme traces.**
   `β = ∞` isolates `max_c N_c`, so `min_c a_c` is read off as an exact integer.
   There is no `β < 0`, so nothing isolates `min_c N_c`: `max_c a_c` is only ever
   seen as a `q^{−k/2}`-damped moment contribution. This is why `a_min` regresses
   at `R² = 1.0000000` while `a_max` regresses at `R² = 0.80`.

Sensitivity check confirming the asymmetry is structural, not numerical:
`∂ log Z_f(β)/∂N_c ∝ N_c^{β−1}`, so at `β = β*` a small fiber is weighted *more*
per unit change than a large one (ratio `0.850` at `q = 211`) — but that is a
smooth `O(1)` weight shared with every other fiber, whereas `β → ∞` puts weight
`1` on the largest fiber and `0` on all the rest.

---

## 6. How many distinct signatures? (a collision census)

Distinct signatures among 400 independent random maps per family
(`t2_2_collisions.py`):

| family | q = 101 | q = 211 | q = 503 |
|---|---|---|---|
| `y² = P₃(x) + c` (elliptic) | **5** | **3** | **3** |
| `y² = P₅(x) + c` (genus 2) | 394 | 398 | 399 |
| `y² = P₇(x) + c` (genus 3) | 400 | 400 | 400 |
| `a + bx + cy + dxy` (bilinear) | **2** | 2 | 2 |
| dense bidegree (2,2) | 400 | 400 | 400 |
| dense bidegree (3,3) | 400 | 400 | 400 |
| dense bidegree (5,5) | 400 | 400 | 400 |
| `P₄(x) + Q₄(y)` | 396 | 400 | 400 |

The elliptic collapse is **not** an accident and is exactly predicted: the
substitution `x ↦ u²x`, `y ↦ u³y`, `t ↦ u⁶t` carries the fibration
`y² = x³ + ax + t` to `y² = x³ + au^{−4}x + t'` and permutes the base, so the
signature depends only on `a mod (F_q^*)⁴`. That gives `gcd(4, q−1)` classes plus
the `a = 0` (sextic-twist / CM) class:

* `q = 101 ≡ 1 (4)`: `gcd = 4`, predicted `4 + 1 = 5` — observed 5.
* `q = 211, 503 ≡ 3 (4)`: `gcd = 2`, predicted `2 + 1 = 3` — observed 3.

So the exchange rate cannot possibly distinguish two members of the same quartic
twist class of the elliptic fibration — but that is because their *signatures*
coincide, not a failure of the rate. From genus 2 up, random maps have distinct
signatures with probability `≈ 1 − O(1/N)`; the only large collision classes in
the whole study are the elliptic twist classes and the degenerate bilinear maps.

---

## 7. Summary

| statistic | recoverable from a small rate vector? | channel | accuracy |
|---|---|---|---|
| image size | **yes, exactly** | `β = 0` endpoint | exact integer |
| `max_c N_c = q − min_c a_c` | **yes, exactly** | `β = ∞` endpoint | exact integer, 1718/1718 |
| `m₂` | **yes** | interior at `β* = √2 − 1` | rel. `O(q^{−1/2})`; 1.8e-3 median at `q = 503` |
| `m₃` | **yes, weakly** | same tangency, damped `0.5286/√q` | 10–20 % of the sample spread |
| `m₄` | yes, very weakly | damped a further `0.6464/√q` | visible in the ladder, not invertible |
| `min_c N_c = q − max_c a_c` | **no** | none — no endpoint at `β < 0` | 1-NN error 3–4 on a spread of 76–96 |
| the multiset `{N_c}` | yes **only** with an unbounded flat reference family | convex-envelope reconstruction of `log Z_f` | exact in principle |

**Negative results worth recording.**

* The interior tangency of `C(f → L)` is *pinned* at `√2 − 1` for every
  curve-like `f`. Enlarging the reference family with more `q²`-point,
  all-fibers-`≈ q` signatures does **not** buy new tangency points at leading
  order; it re-samples nearly the same functional. The probes that do move
  (`Xaniso`, `Sq`) move because they are *not* Hasse–Weil normalised.
* A single rate is one real functional and is blind to the one-parameter
  direction `δm₂ = −(0.5286/√q) δm₃`; explicit near-collisions at `3.3e-11`
  (`q = 211`) and `4.7e-12` (`q = 503`) with `Δm₂` of `1.9e-02` and `4.0e-03`.
* The positive extreme trace `max_c a_c` is invisible to the endpoints. Two maps
  can share image size, largest fiber and `m₂` exactly while their smallest
  fibers differ by 6, with the entire 8-probe rate vector agreeing to `2.4e-07`.

**Not** the case, contrary to the framing of the question: the rates do *not*
see only the extremes. They see the whole moment ladder, with a clean and
completely explicit damping constant `−(β* − k)/((k+1)√q)` per moment order.
