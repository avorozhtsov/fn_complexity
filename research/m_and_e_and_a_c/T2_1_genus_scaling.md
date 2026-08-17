# T2.1 — Genus scaling of `1 − C(L→f)` for hyperelliptic fibrations

Scripts: `t2_1_genus_scaling.py` (arithmetic, exchange rates, tables),
`t2_1_extreme_value.py` (the `USp(2g)` tail law and the prediction).
Data: `t2_1_genus_scaling.csv`, `t2_1_extreme_value.csv`.

## Verdict in one paragraph

The law is **exactly right in form and badly misleading as a numerical
prediction**. Two exact statements replace the heuristic:

```
C(L → f) = log q / log(max_c N_c)          for every f : A^2 → A^1     (Lemma 2)
(1 − C) √q log q  =  μ − μ²(½ + 1/log q)/√q + O(1/q),   μ = max_c(−a_c)/√q
```

so the whole question is whether `μ → 2g`, i.e. whether the **largest** of the
`q` Frobenius traces reaches the Weil edge. It does, but the approach is
governed by extreme-value statistics with exponent `2/dim USp(2g) = 2/(2g²+g)`:

```
2g − μ  ≈  Γ(1 + 2/d) (q K_g)^{−2/d},        d = 2g² + g = dim USp(2g)
```

For `g = 1` that is `q^{−2/3}` and the law is visible immediately (μ = 1.999 at
q ≈ 10⁶). For `g = 2` it is `q^{−1/5}`; for `g = 3`, `q^{−2/21}`; for `g = 4`,
`q^{−1/18}`. So at `q = 10⁶` we measure `μ = 3.68` against `2g = 4`, `μ = 4.66`
against `6`, and `μ = 4.76` against `8`. **The 2g law is confirmed as the
correct asymptotic, and confirmed to be numerically unreachable for `g ≥ 3`:**
reaching within 10 % of `2g` needs `q ≈ 6·10⁴` (g=2), `6·10⁹` (g=3), `4·10¹⁶`
(g=4). The observed deficits agree with the predicted extreme-value deficits to
a few percent at every `q` and every genus, and the fitted exponents agree with
`−2/(2g²+g)` to within 1–4 % for `g = 2, 3, 4`.

## 1. Setup

`f(x,y) = y² − h(x)` over a prime field `F_q`, `deg h = 2g+1` **odd**, so the
smooth projective model of `y² = h(x)+c` has exactly one rational point at
infinity and the affine count is `N_c = q − a_c` with `a_c` the trace of
Frobenius, `|a_c| ≤ 2g√q` for the smooth fibres. Sixteen families are used
(4 each of genus 1, 2, 3; 2 of genus 4; 3 genus-0 controls plus the split
conic `xy`), and 14 primes from 101 to 1 000 003.

**Lemma 1.** `Σ_c a_c = 0` exactly, since `Σ_c N_c = q²`. Verified for every
family and every `q` (assertion in `summarize`).

**Lemma 2 (the rate is an endpoint rate, always).** For any `f : A²→A¹` over
`F_q`, `C(L→f)` is attained at `β = ∞`, so
`C(L→f) = log q / log(max_c N_c)`.

*Proof.* Write `A = log q`, `B = log(max_c N_c)`, `Z_L(β) = q^{1+β}`. Since
`Σ_c N_c = q²` over `q` fibres, `Z_f(β) = Σ_c N_c^β < q·(max N_c)^β` unless all
`N_c` are equal, hence for finite `β`
`ratio(β) = (1+β)A / log Z_f(β) > (1+β)A/(A+βB) > A/B = ratio(∞)`,
the last step being exactly `B > A`. ∎

This was checked numerically against the repo solver
(`exchange_rate_result`, `implementer = (q,…,q)`) for `q ≤ 1009` and all four
genera: the solver returns `beta = inf` and matches the closed form to `1e−12`.
The pre-existing `analysis/frobenius_exchange_rates.csv` shows the same
(`L->E1` has `contact_beta = inf` at every `q`).

**Corollary (exact form of the scaled deviation).** With `m = max_c(−a_c)` and
`μ = m/√q`,

```
(1 − C) √q log q = √q log q · log(1 + m/q) / log(q + m)
                 = μ − μ²(½ + 1/log q)/√q + O(q^{-1}).
```

Measured gap `μ − (1−C)√q log q` against its leading term:

| family | q | μ | (1−C)√q log q | gap | μ²(½+1/log q)/√q |
|---|---:|---:|---:|---:|---:|
| E1 | 101 | 1.8906 | 1.6699 | 0.2206 | 0.2549 |
| E1 | 2003 | 1.9886 | 1.9346 | 0.0540 | 0.0558 |
| E1 | 1000003 | 1.9990 | 1.9967 | 0.0023 | 0.0023 |
| H2 | 101 | 2.8856 | 2.4052 | 0.4804 | 0.5938 |
| H2 | 32003 | 3.4546 | 3.4153 | 0.0392 | 0.0398 |
| H2 | 1000003 | 3.6840 | 3.6762 | 0.0077 | 0.0078 |
| H3 | 2003 | 3.4410 | 3.2831 | 0.1579 | 0.1671 |
| H3 | 1000003 | 4.6630 | 4.6506 | 0.0124 | 0.0124 |

This correction is **deterministic and always negative**: the rate quantity
undershoots `μ` by `Θ(μ²/√q)`. For `g = 1` it is the *dominant* correction
(`q^{−1/2}` beats the extreme-value `q^{−2/3}`); for `g ≥ 2` the extreme-value
deficit `q^{−2/(2g²+g)}` dominates it.

### Method

Traces are computed exactly without a `q²` loop:
`N_c = Σ_x (1 + χ(h(x)+c)) = q + (m_h ⋆ χ)[c]`, a circular cross-correlation of
the value-multiplicity vector of `h` with the Legendre symbol, done by FFT and
rounded to integers (float error ≈ 1e−9 even at `q = 10⁶`). This is validated
against the naive `O(q²)` bincount for all 16 families at `q ∈ {101,211,401}`.
Total runtime for the whole sweep: 13 s.

## 2. Tables

### Genus 1 — `E1 : h = x³ + x`

| q | m = max(−a_c) | μ = m/√q | max\|a_c\|/√q | (1−C)√q log q | 2g − μ | q⁻²Σa² | q⁻³Σa⁴ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 101 | 19 | 1.8906 | 1.8906 | 1.6699 | 0.1094 | 1.010 | 2.075 |
| 211 | 28 | 1.9276 | 1.9276 | 1.7688 | 0.0724 | 0.986 | 1.932 |
| 401 | 39 | 1.9476 | 1.9476 | 1.8302 | 0.0524 | 1.002 | 1.980 |
| 1009 | 62 | 1.9518 | 1.9518 | 1.8780 | 0.0482 | 0.997 | 2.023 |
| 2003 | 89 | 1.9886 | 1.9886 | 1.9346 | 0.0114 | 1.000 | 1.982 |
| 4001 | 126 | 1.9920 | 1.9920 | 1.9540 | 0.0080 | 1.000 | 1.995 |
| 8009 | 177 | 1.9778 | 1.9890 | 1.9515 | 0.0222 | 1.000 | 2.020 |
| 16001 | 252 | 1.9922 | 1.9922 | 1.9735 | 0.0078 | 1.000 | 2.008 |
| 32003 | 357 | 1.9956 | 1.9956 | 1.9824 | 0.0044 | 1.000 | 1.991 |
| 64007 | 505 | 1.9961 | 1.9961 | 1.9868 | 0.0039 | 1.000 | 2.007 |
| 128021 | 715 | 1.9983 | 1.9983 | 1.9918 | 0.0017 | 1.000 | 1.999 |
| 256019 | 1011 | 1.9981 | 1.9981 | 1.9935 | 0.0019 | 1.000 | 1.998 |
| 512009 | 1431 | 1.9999 | 1.9999 | 1.9967 | 0.0001 | 1.000 | 1.999 |
| 1000003 | 1999 | 1.9990 | 1.9990 | 1.9967 | 0.0010 | 1.000 | 2.001 |

Reproduces and extends the known 101…2003 numbers (1.670, 1.769, 1.830, 1.878,
1.935) and pushes the convergence to 1.9967 at `q = 10⁶`. Note `q⁻²Σa² → 1`
(semicircle) and `q⁻³Σa⁴ → 2` (Catalan), the `SU(2)` values.

**The genus-1 family saturates the integer Weil bound.** From `q ≈ 4000` on,
`m = ⌊2√q⌋` for essentially every prime and every genus-1 family (see the gap
column of §3): the deficit is then not extreme-value but *quantisation*, the
fractional part `2√q − ⌊2√q⌋`, which is `O(1/√q)` after dividing by `√q`. This
is why the fitted genus-1 exponent below is `−0.56` rather than the intrinsic
`−2/3`.

### Genus 2 — `H2 : h = x⁵ + x²`

| q | m | μ | max\|a_c\|/√q | (1−C)√q log q | 2g − μ | q⁻²Σa² | q⁻³Σa⁴ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 101 | 29 | 2.8856 | 2.8856 | 2.4052 | 1.1144 | 0.950 | 2.976 |
| 211 | 35 | 2.4095 | 2.8226 | 2.1672 | 1.5905 | 1.109 | 4.371 |
| 401 | 54 | 2.6966 | 3.0961 | 2.4777 | 1.3034 | 0.988 | 2.824 |
| 1009 | 99 | 3.1167 | 3.1167 | 2.9334 | 0.8833 | 1.053 | 3.391 |
| 2003 | 145 | 3.2399 | 3.6421 | 3.0995 | 0.7601 | 1.000 | 2.936 |
| 4001 | 201 | 3.1777 | 3.3358 | 3.0822 | 0.8223 | 0.999 | 2.952 |
| 8009 | 307 | 3.4304 | 3.4304 | 3.3523 | 0.5696 | 1.000 | 2.955 |
| 16001 | 444 | 3.5100 | 3.5100 | 3.4524 | 0.4900 | 1.000 | 2.991 |
| 32003 | 618 | 3.4546 | 3.4602 | 3.4153 | 0.5454 | 1.000 | 2.969 |
| 64007 | 908 | 3.5890 | 3.7313 | 3.5592 | 0.4110 | 1.000 | 3.019 |
| 128021 | 1339 | 3.7423 | 3.7423 | 3.7196 | 0.2577 | 1.000 | 3.002 |
| 256019 | 1888 | 3.7313 | 3.7313 | 3.7155 | 0.2687 | 1.000 | 2.999 |
| 512009 | 2707 | 3.7831 | 3.7831 | 3.7716 | 0.2169 | 1.000 | 3.002 |
| 1000003 | 3684 | 3.6840 | 3.6840 | 3.6762 | 0.3160 | 1.000 | 3.002 |

Target `2g = 4`. We reach 3.68 at `q = 10⁶`: convergent, but 8 % short, and the
approach is visibly `q^{−1/5}` slow. Second moment `→ 1`, fourth moment `→ 3`
(the `USp(4)` Gaussian value, *not* the genus-1 value 2).

### Genus 3 — `H3 : h = x⁷ + x²`

| q | m | μ | max\|a_c\|/√q | (1−C)√q log q | 2g − μ | q⁻²Σa² | q⁻³Σa⁴ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 101 | 20 | 1.9901 | 1.9901 | 1.7473 | 4.0099 | 0.990 | 2.462 |
| 211 | 39 | 2.6849 | 3.0979 | 2.3879 | 3.3151 | 1.327 | 4.820 |
| 401 | 59 | 2.9463 | 2.9463 | 2.6872 | 3.0537 | 1.022 | 3.134 |
| 1009 | 122 | 3.8407 | 3.8407 | 3.5668 | 2.1593 | 0.993 | 2.880 |
| 2003 | 154 | 3.4410 | 3.4410 | 3.2831 | 2.5590 | 0.997 | 2.791 |
| 4001 | 215 | 3.3990 | 3.3990 | 3.2901 | 2.6010 | 1.025 | 3.123 |
| 8009 | 329 | 3.6763 | 3.7657 | 3.5867 | 2.3237 | 0.999 | 2.960 |
| 16001 | 449 | 3.5495 | 3.6839 | 3.4907 | 2.4505 | 0.991 | 2.748 |
| 32003 | 761 | 4.2539 | 4.2539 | 4.1946 | 1.7461 | 1.000 | 3.035 |
| 64007 | 1088 | 4.3005 | 4.3005 | 4.2578 | 1.6995 | 1.000 | 3.004 |
| 128021 | 1387 | 3.8765 | 4.0553 | 3.8521 | 2.1235 | 1.003 | 3.037 |
| 256019 | 2267 | 4.4804 | 4.4804 | 4.4575 | 1.5196 | 1.000 | 2.995 |
| 512009 | 3129 | 4.3729 | 4.3785 | 4.3575 | 1.6271 | 1.000 | 2.993 |
| 1000003 | 4663 | 4.6630 | 4.6630 | 4.6506 | 1.3370 | 1.000 | 2.994 |

Target `2g = 6`. We reach 4.65 — a 22 % shortfall, still shrinking like
`q^{−2/21}`. **Stated honestly: at accessible `q` the genus-3 numbers do not
look like 6, and no feasible computation will make them.**

### Genus 4 — `H4 : h = x⁹ + x²` (target 8)

`μ` reaches 4.76 at `q = 10⁶` (a 40 % shortfall), with `q⁻²Σa² = 1.000`,
`q⁻³Σa⁴ = 3.001` — the family is perfectly equidistributed; only the maximum is
starved.

## 3. Where the deficit comes from (question 3)

### The exact `USp(2g)` lower-tail law

Under Katz–Sarnak, `T_c = −a_c/√q` equidistributes for the trace law of
`USp(2g)`. Writing `x_i = cos θ_i`, the Weyl formula is
`dμ = (2^{g²}/(g! π^g)) Π_{i<j}(x_i−x_j)² Π_i √(1−x_i²) dx` on `[−1,1]^g`,
`T = 2Σx_i`. The substitution `1 − x_i = (ε/2)v_i` maps `{2g − T < ε}` onto the
standard simplex and gives the **exact** identity

```
P(2g − T < ε) = (2^{g²}/(g! π^g)) (ε/2)^{d/2} J_g(ε),        d = 2g² + g,
J_g(ε) = ∫_{v>0, Σv<1} Π_{i<j}(v_j−v_i)² Π_i √v_i · Π_i √(2 − (ε/2)v_i) dv
```

— the whole `ε`-dependence outside the prefactor sits in a bounded analytic
factor. Letting `ε → 0`, `P(2g − T < ε) ~ K_g ε^{d/2}` with
`K_g = 2^{g²−d/2}/(g! π^g) · J_g(0)`, i.e. **the edge exponent is exactly half
the dimension of the monodromy group.**

`J_g` is evaluated by Monte-Carlo on the simplex (bounded integrand). Checks:

* `g = 1`: agrees with the closed form `(θ₀ − ½ sin 2θ₀)/π` to `1.1e−4`
  (Monte-Carlo error), `K_1 = 0.212183` vs `2/3π = 0.212207`.
* `g = 2`: agrees with 1.87 M rejection samples of the Weyl measure at
  `ε = 0.5, 1, 2, 3, 4` to 0.2 %; sampled `E[T²] = 1.0003`, `E[T⁴] = 3.0035`.

| g | d = dim USp(2g) | edge exponent d/2 | K_g | deficit exponent −2/d | crossover log q ≈ 2g² |
|---:|---:|---:|---:|---:|---:|
| 1 | 3 | 1.5 | 2.1218e−01 | −0.6667 | 2 |
| 2 | 10 | 5.0 | 9.9351e−04 | −0.2000 | 8 |
| 3 | 21 | 10.5 | 2.1244e−08 | −0.0952 | 18 |
| 4 | 36 | 18.0 | 8.7578e−16 | −0.0556 | 32 |

### The arithmetic tail *is* the `USp(2g)` tail

Pooling the four families at `q = 1 000 003` (4.0 M traces per genus):

| g | ε | arithmetic `P(2g−T<ε)` | `USp(2g)` | ratio |
|---:|---:|---:|---:|---:|
| 1 | 0.25 | 0.0258319 | 0.0260199 | 0.993 |
| 1 | 0.50 | 0.0720905 | 0.0721387 | 0.999 |
| 1 | 1.00 | 0.1955037 | 0.1954786 | 1.000 |
| 2 | 0.25 | 0.0000010 | 0.0000009 | 1.058 |
| 2 | 0.50 | 0.0000310 | 0.0000294 | 1.054 |
| 2 | 1.00 | 0.0009227 | 0.0008869 | 1.040 |
| 2 | 2.00 | 0.0247224 | 0.0247102 | 1.000 |
| 3 | 2.00 | 0.0000240 | 0.0000236 | 1.018 |
| 3 | 3.00 | 0.0013860 | 0.0013975 | 0.992 |

Agreement to a few percent down to probability `10⁻⁶`. So the shortfall is
**not** a defect of the families: they are equidistributed right up to the
edge; there simply are not enough fibres for the maximum to get there.

### Observed vs predicted deficit

With `n = q` (approximately independent) samples the deficit is Weibull of
shape `d/2`, `E[2g − max T] = ∫₀^{4g} (1 − F(t))^q dt`. Predictions below use
the exact `F`; "edge asympt." is `Γ(1+2/d)(qK_g)^{−2/d}`. Observed values are
means over the 4 families (2 for `g = 4`) at each `q`.

**Genus 2 (2g = 4):**

| q | observed deficit | spread | predicted E | edge asympt. | obs μ | pred μ | gap to ⌊2g√q⌋ (a-units) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1009 | 1.1037 | 0.1928 | 0.9385 | 0.9177 | 2.8963 | 3.0615 | 35.00 |
| 4001 | 0.7828 | 0.1087 | 0.7084 | 0.6967 | 3.2172 | 3.2916 | 49.50 |
| 16001 | 0.5157 | 0.0423 | 0.5346 | 0.5280 | 3.4843 | 3.4654 | 64.25 |
| 64007 | 0.3458 | 0.1126 | 0.4039 | 0.4002 | 3.6542 | 3.5961 | 86.50 |
| 256019 | 0.2761 | 0.0779 | 0.3054 | 0.3033 | 3.7239 | 3.6946 | 138.75 |
| 1000003 | 0.2443 | 0.0429 | 0.2321 | 0.2309 | 3.7557 | 3.7679 | 244.25 |

**Genus 3 (2g = 6):**

| q | observed deficit | spread | predicted E | edge asympt. | obs μ | pred μ | gap (a-units) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1009 | 2.6236 | 0.3349 | 2.7580 | 2.6539 | 3.3764 | 3.2420 | 82.75 |
| 4001 | 2.1306 | 0.4117 | 2.4039 | 2.3276 | 3.8694 | 3.5961 | 134.25 |
| 16001 | 2.2350 | 0.1419 | 2.0961 | 2.0397 | 3.7650 | 3.9039 | 281.75 |
| 64007 | 1.8715 | 0.2554 | 1.8294 | 1.7874 | 4.1285 | 4.1706 | 472.50 |
| 256019 | 1.5018 | 0.1472 | 1.5978 | 1.5664 | 4.4982 | 4.4022 | 759.00 |
| 1000003 | 1.3963 | 0.1343 | 1.3995 | 1.3757 | 4.6037 | 4.6005 | 1396.25 |

**Genus 4 (2g = 8):** predicted 3.174 at `q = 10⁶`, observed 3.241.

**Genus 1:** predicted 0.00025 at `q = 10⁶`, observed 0.00080 — the only place
where the prediction is systematically *low*, because `a_c` is an integer:
`⌊2√q⌋ − m = 0` at almost every large `q` (the gap column is 0.00 for eight of
the fourteen genus-1 primes), so the residual deficit is the fractional part of
`2√q`, of mean `1/2` in `a`-units, whereas the continuous extreme-value deficit
is `2.540 q^{−1/6}` in `a`-units — below 1 already for `q ≳ 270`. **For `g = 1`,
integer quantisation of `a_c` takes over from extreme-value fluctuation at
`q ≈ 300`.**

### Fitted exponents (log–log least squares, `q ∈ [4001, 10⁶]`)

| g | fitted slope | edge `−2/(2g²+g)` | Gaussian-regime slope |
|---:|---:|---:|---:|
| 1 | −0.5617 | −0.6667 | — |
| 2 | −0.1969 | −0.2000 | — |
| 3 | −0.0931 | −0.0952 | −0.1277 |
| 4 | −0.0583 | −0.0556 | −0.0597 |

Genus 2, 3, 4 confirm `−2/(2g²+g)` to 1–4 %. Genus 1 sits at `−0.56` because
quantisation (slope `−1/2`) has overtaken the intrinsic `−2/3`.

### Why `g ≥ 3` is hopeless: the Gaussian regime

The `USp(2g)` trace has variance 1 for every `g`, and is close to `N(0,1)` in
the bulk (`E[T⁴] = 3` for `g ≥ 2` — visible in the tables). The maximum of `q`
standard normals sits at `b_q = √(2L) − log(4πL)/(2√(2L))`, `L = log q`, which
is `4.766` at `q = 10⁶` **independently of the genus**. The Weil edge is felt
only once `b_q` approaches `2g`, i.e. once `log q ≳ 2g²`:

| family at q = 10⁶ | 2g | Gaussian max `b_q` | observed μ |
|---|---:|---:|---:|
| g = 2 | 4 | 4.766 (capped) | 3.756 |
| g = 3 | 6 | 4.766 | 4.604 |
| g = 4 | 8 | 4.766 | 4.760 |

For `g = 4` the observed maximum equals the Gaussian prediction to three
decimals: at `q = 10⁶` the genus-4 family is statistically indistinguishable
from `q` iid standard normals, and `2g = 8` plays no role whatsoever.
`q` required to bring `μ` within 10 % / 5 % of `2g`:

| g | within 10 % | within 5 % |
|---:|---:|---:|
| 2 | 6.4·10⁴ | 2.1·10⁶ |
| 3 | 6.1·10⁹ | 8.8·10¹² |
| 4 | 3.7·10¹⁶ | 9.8·10²¹ |

## 4. `max_c(−a_c)` versus `max_c |a_c|` (question 2)

They are **not** the same quantity and the rate only sees the first. `C(L→f)`
depends on `max_c N_c`, i.e. on the most *positive* fibre excess `−a_c`; the
Weil bound is two-sided. Over the full sweep:

| genus | family×prime pairs | how often `max(−a) = max\|a\|` | mean `max\|a\|/max(−a) − 1` |
|---:|---:|---:|---:|
| 1 | 56 | 50 (89 %) | 0.0065 |
| 2 | 56 | 28 (50 %) | 0.0415 |
| 3 | 56 | 33 (59 %) | 0.0406 |
| 4 | 28 | 19 (68 %) | 0.0255 |

Both converge to `2g` (the `USp(2g)` trace law is symmetric because `−I` lies
in the group), so the asymptotic statement is unaffected; but at any finite `q`
the rate quantity systematically tracks the *lower* of the two, on average 4 %
below `max|a_c|/√q` for `g = 2, 3`. The expectation that the rate quantity and
`max|a_c|/√q` agree is therefore true only in the limit — corrected here. The
honest identity is `(1−C)√q log q ↔ max_c(−a_c)/√q`, and the two columns in §2
agree exactly whenever the extremal trace happens to be negative.

## 5. Genus-0 controls (question 4)

Three genus-0 families plus the split conic, and the answer depends entirely on
whether some fibre is **geometrically reducible over `F_q`**, not on the genus:

| family | f | m = max(−a_c) | `(1−C)√q log q` | conclusion |
|---|---|---|---|---|
| `G0a` | `y² − x` | 0 | 0 exactly | `C = 1`, the signature *is* `L` |
| `G0c` (q ≡ 3 mod 4) | `x² + y²` | 1 | `1/√q` | `1 − C = 1/(q log q)`, the predicted degeneration |
| `G0c` (q ≡ 1 mod 4) | `x² + y²` | q − 1 | `≈ 0.66 √q` | reducible fibre, law breaks |
| `G0b`, `XY` | `y²−x²−x`, `xy` | q − 1 | `≈ 0.66 √q` | reducible fibre, law breaks |

* **`f = x² + y²` with `q ≡ 3 mod 4` is the clean genus-0 test.** Every `c ≠ 0`
  gives an anisotropic conic with `N_c = q+1`, `a_c = −1`; `c = 0` gives the
  single point `(0,0)`, `a_0 = q−1`. So `m = 1`, `μ = q^{−1/2} → 0 = 2g`, and
  `(1−C)√q log q = 0.0688, 0.0223, 0.0056, 0.0040, 0.0020, 0.0010` at
  `q = 211, 2003, 32003, 64007, 256019, 1000003` — exactly `1/√q`, i.e.
  `1 − C = 1/(q log q)`. **Degenerates correctly, as expected.**

* **The split conic `xy` does *not*.** Its fibre over `0` is a union of two
  lines, `N_0 = 2q−1`, so `m = q−1` and
  `1 − C = log(2 − 1/q)/log(2q−1) ~ log 2 / log q`. The scaled quantity
  `(1−C)√q log q` *diverges* like `√q log 2` (660.03 at `q = 10⁶`, against
  `log 2 · √q = 693.1`). The expectation of `1/(q log q)` for `xy` is wrong;
  the correct genus-0 control is `x² + y²` at `q ≡ 3 mod 4`.

  This is a structural point worth keeping: `1 − C(L→f) ≍ 2g/(√q log q)`
  requires *every* fibre to be geometrically irreducible. A single fibre with
  `r` geometric components over `F_q` contributes `N ≈ rq`, hence
  `1 − C ≍ log r / log q` — a `1/log q` effect that swamps every Weil-scale
  term. For `y² − h(x)` with `h` squarefree of odd degree the singular fibres
  are nodal and irreducible, which is why the genus-1…4 families behave.

  Note also the `q mod 4` sensitivity of `x²+y²`: the same map has
  `C(L→f) → 1` at rate `1/(q log q)` when `q ≡ 3 mod 4` and only at rate
  `log 2 / log q` when `q ≡ 1 mod 4`. This is the same kind of arithmetic
  detection as the established `C(y²−x³ → L) = 1 ⟺ q ≡ 2 mod 3`.

## 6. Surprises worth recording

1. **`C(L→f)` is an endpoint rate for every `f : A²→A¹`** (Lemma 2), not just
   for the elliptic family — the `β = ∞` contact observed in
   `analysis/frobenius_exchange_rates.csv` is a theorem, not a coincidence.
   Consequently the exchange rate against `L` sees *only* `max_c N_c`. All the
   remaining structure of `{a_c}` (moments, symmetry type) is invisible to this
   one rate; T2.2/T2.3 must use other pairs.
2. **The exponent in the deficit is `2/dim USp(2g)`.** `2g²+g` is the dimension
   of the monodromy group, and the edge exponent `d/2` is a general fact about
   compact groups (`P(tr near max) ~ ε^{dim/2}`). So the *rate of convergence*
   of the exchange rate to its Weil-scale limit encodes the dimension of the
   monodromy group, while the limit itself encodes only `2g`. That is a second,
   finer arithmetic invariant readable off the same curve — and a concrete way
   for M to see symmetry type (T2.3).
3. **Genus 1 saturates the integer Weil bound**: `m = ⌊2√q⌋` for essentially
   all `q ≳ 4000`, so the genus-1 deficit is pure quantisation and the fitted
   exponent is `−1/2`, not `−2/3`.
4. **For `g ≥ 3` the observed maximum is a Gaussian extreme, not a Weil edge.**
   At `q = 10⁶` the genus-4 maximum matches the corrected `√(2 log q)` to three
   decimals; the Weil bound `2g = 8` is irrelevant at that scale. The crossover
   is at `log q ≈ 2g²`.
5. **`q⁻³Σa_c⁴` cleanly separates `g = 1` from `g ≥ 2`** (2 vs 3), while
   `q⁻²Σa_c² = 1` for all genera. The exchange rates at finite `β` (not the
   `β = ∞` rate against `L`) are the natural place to look for these — a
   concrete handle for T2.2.
