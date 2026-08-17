# Continued-fraction hyperbola analysis

This analysis concerns

```text
C({2,2} | {3,1}) = log(2) / log(3)
```

and the exact points `(n, k_max(n))` for `1 <= n <= 100`.

## Why the curves are hyperbolas

For a fixed rational slope `e`, define

```text
d_n = e*n - k_max(n).
```

All points with the same value of `d_n` lie on the line

```text
k_max(n) = e*n - d
```

in the `(n, k_max)` plane. Dividing by `n` gives the exact hyperbola

```text
c_n = e - d/n.
```

A single point does not determine `e`: it lies on one such curve for every
candidate slope. The analysis therefore restricts `e` to continued-fraction
convergents of the limiting exchange rate.

## Candidate slopes

The continued fraction begins

```text
log(2)/log(3) = [0; 1, 1, 1, 2, 2, 3, 1, ...].
```

The positive convergents with denominator at most 100 are

```text
1, 1/2, 2/3, 5/8, 12/19, 41/65, 53/84.
```

For each slope, equal values of `d_n` partition all 100 points. Consequently,
every point belongs to exactly seven candidate hyperbolas, one for each slope.
Across all slopes there are 259 distinct curves, including singleton curves.

## Greedy versus minimum cover

The literal first-uncovered-point algorithm is still ambiguous until a curve
selection rule is specified. Using the following deterministic rule:

1. take the smallest uncovered `n`;
2. among its candidate curves, maximize the number of uncovered points;
3. break ties by total curve length, then by smaller `|d|`;

produces 15 curves.

A global coverage-first greedy pass finds 12 curves. An exact branch-and-bound
set-cover search proves that no cover of 11 or fewer candidate curves exists.
The minimum cover has a particularly simple form: all 12 curves use the single
convergent

```text
e = 5/8.
```

These 12 curves are disjoint and partition the data, so the color assignment is
unambiguous. Ordered by their first point, they are:

| Order | Hyperbola | Number of points |
| ---: | --- | ---: |
| 1 | `c_n = 5/8 - 5/(8n)` | 8 |
| 2 | `c_n = 5/8 - 1/(4n)` | 13 |
| 3 | `c_n = 5/8 - 7/(8n)` | 3 |
| 4 | `c_n = 5/8 - 1/(2n)` | 11 |
| 5 | `c_n = 5/8 - 1/(8n)` | 12 |
| 6 | `c_n = 5/8 - 3/(4n)` | 5 |
| 7 | `c_n = 5/8 - 3/(8n)` | 12 |
| 8 | `c_n = 5/8` | 12 |
| 9 | `c_n = 5/8 + 1/(8n)` | 10 |
| 10 | `c_n = 5/8 + 1/(4n)` | 7 |
| 11 | `c_n = 5/8 + 3/(8n)` | 5 |
| 12 | `c_n = 5/8 + 1/(2n)` | 2 |

The structural reason `5/8` wins is the tradeoff between denominator and
approximation error. Its denominator is only eight, while it is already close
to `log(2)/log(3)`. Thus the residue-class branches persist for many steps
before a phase slip changes `d`. Later convergents are closer but have too many
residue classes; earlier convergents have smaller denominators but drift too
quickly.

## Generated files

- `hyperbolas_2-2_over_3-1.csv`: all 259 candidate hyperbolas.
- `point_hyperbola_map_2-2_over_3-1.csv`: every point and its seven candidate
  curves, plus its minimum-cover curve and plot color.
- `minimum_cover_2-2_over_3-1.csv`: the 12-curve optimal partition.
- `first_uncovered_greedy_2-2_over_3-1.csv`: the 15-step literal greedy result.

Regenerate all tables from the project root:

```bash
./analysis/hyperbola_cover.py
```

The curated PDF-only plots and their regeneration commands are documented in
[`images/README.md`](../images/README.md).

## Relativistic radiation energy--entropy diagram

Generate the white-background vector figure comparing photon and
ultrarelativistic electron--positron sectors with:

```bash
python analysis/plot_relativistic_species_energy_entropy.py
```

The script writes
`paper_exchange_rate/figures/relativistic-radiation-energy-entropy.pdf`.  It uses
\(g_\gamma=2\) and \(g_{e^\pm}=(7/8)4=7/2\), so the dotted \(7/4\)
homothety of the photon curve coincides exactly with the electron--positron
curve.

## Opposite direction through `n < 100000`

For `C({3,1} | {2,2})`, feasibility has the exact binomial-tail form

```text
sum(C(n,j), j >= ceil(k log(2)/log(3))) >= 2^k.
```

Generate the best branch for each of the first eight convergents and the
complete `e=1` partition with:

```bash
./analysis/reverse_hyperbola_cover.py
```

This writes:

- `best_convergent_hyperbolas_3-1_over_2-2_n-1-99999.csv`;
- `e1_partition_3-1_over_2-2_n-1-99999.csv`.

The latter uses both sign conventions: `c_n=1-m/n`, or equivalently
`c_n=1+d/n` with `d=-m`.

### First 30 complete `e=1` branches

The total numbers of points on `c_n=1-d/n`, with no upper cutoff on `n`, are:

| `d` | Points | `d` | Points | `d` | Points |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 6 | 11 | 28 | 21 | 29 |
| 2 | 16 | 12 | 28 | 22 | 27 |
| 3 | 21 | 13 | 27 | 23 | 30 |
| 4 | 23 | 14 | 27 | 24 | 27 |
| 5 | 26 | 15 | 28 | 25 | 28 |
| 6 | 25 | 16 | 28 | 26 | 28 |
| 7 | 27 | 17 | 28 | 27 | 29 |
| 8 | 27 | 18 | 28 | 28 | 27 |
| 9 | 26 | 19 | 28 | 29 | 30 |
| 10 | 27 | 20 | 28 | 30 | 28 |

Their last occurrence is `n=789`.  Completeness follows independently of the
`n < 100000` data cutoff: for every `d <= 30`, the binomial Chernoff bound is
already below `2^-30` at `n=900` and decreases thereafter, so `k=n-d` is
infeasible for all `n >= 900`.

The continued fraction of the asymptotic branch spacing is

```text
1 / (1-C) = 28.8887439648...
          = [28; 1, 7, 1, 84, 3, 2, ...].
```

It predicts the limiting average of about `28.8887` points per branch.  It
does not determine each exact count on its own.  Here `k_max` is defined by a
binomial-tail inequality, rather than by a single floor of an irrational
multiple.  Its second-order behavior contains a logarithmic correction and a
lattice phase.  Exact counts require the binomial-tail test in
`reverse_hyperbola_cover.py`; continued fractions organize only the leading
Beatty-like spacing.

## Small exchange-rate matrix

Generate the matrix for all decreasing length-3 signatures with entries below
4, together with the length-4 signature whose entries are below 2:

```bash
./analysis/exchange_matrix_small.py
```

Rows are implementers `g`, columns are implemented signatures `f`, and cells
are `C(g -> f)`. Outputs:

- `exchange_matrix_small.csv` with 15-digit values;
- `exchange_matrix_small.md` with a readable six-decimal table.

The ordering uses the directed comparison requested for the project:

```text
b -> a  when  C(a -> b) > C(b -> a).
```

A deterministic topological sort displays lower-complexity signatures first.
Equal rates add no edge and can produce multi-signature layers. For this
particular set, all 55 pairs compare strictly, the graph is acyclic, and every
topological layer contains one signature.

## Extended exchange-rate matrix

The extended matrix contains:

- length 1 with entries at most 6;
- length 2 with entries at most 6;
- length 3 with entries at most 5;
- length 4 with entries at most 4.

Regenerate its 97 signatures, all 9,409 directed rates, and the cycle-aware
ordering with:

```bash
./analysis/exchange_matrix_extended.py
```

Outputs:

- `exchange_matrix_extended.csv`, with 15-digit values, ranks, layers, and
  strongly connected component metadata;
- `exchange_matrix_extended.md`, with cyclic components, ordered
  condensation-DAG layers, and the six-decimal matrix.

The edge orientation remains `b -> a` when `C(a -> b) > C(b -> a)`. Equal rates
create no edge. The special zero-complexity signature `{1}` produces infinite
rates in its implemented column, matching the operational definition. Rates
are read through `src/fn_complexity/exchange_rates_cache.json`; only missing
entries are computed, and the updated cache is saved atomically. Existing
values from a previous matrix CSV are imported once before that CSV is
replaced.

## Appendix B: first 69 non-special signatures

Generate the paper's 69-signature table and its backward-comparison exception
table with:

```bash
PYTHONPATH=src python analysis/appendix_b_signatures.py
```

Singletons and all-ones signatures are omitted. The remaining signatures are
exhausted by the natural budgets

```text
S_B = {a : len(a) >= 2, a_1 > 1, a_1 + 2 len(a) <= B}.
```

The strict comparison graph is contracted into strongly connected components
and ordered through its condensation DAG. Power-sum screening gives the same
first 69 signatures for `B=18` and `B=19`; the generator stores that stabilized
cutoff and recomputes all displayed rates and exception relations with the
high-accuracy persistent cache before writing
`paper_exchange_rate/appendix_b_signatures.tex`.

Re-run the large stabilization check (approximately 400 MB peak memory) with:

```bash
PYTHONPATH=src:analysis python analysis/stabilize_signature_order.py
```

For a numbered signature `x`, an exception is an earlier signature `a` for
which `x ≺ a`. Each exception prints both `C(x -> a)` and `C(a -> x)`. All
comparisons use and extend the persistent exchange-rate cache.

Appendix B also formulates the global strongly connected component problem for
arrows `a -> b` when `C(a -> b) >= C(b -> a)`. Finite induced subgraphs do not
determine this global cluster. Nested-shell screening currently supports, but
does not prove, the conjecture that the cluster of `{3,1,1}` contains every
non-special signature except `{2,1}` and `{2,2}`.

## Exchange positivity and Weil numbers

`analysis/exchange_positivity.py` runs the three computations behind
`paper_finite_fields_maps/docs/exchange_positivity_and_weil.md`:

```bash
python analysis/exchange_positivity.py
```

1. The isometry identity `d(a,b) = osc_beta(u_a - u_b)` for
   `u_a = log log Z_a`, where `d(a,b) = -log(C(a -> b) C(b -> a))`.
2. An explicit thirteen-signature family whose exchange metric violates
   negative type (`x^T D x = +1.09e-3` with `sum(x) = 0`), so `exp(-t d)` is not
   positive semidefinite for small `t` and the exchange geometry admits no
   Hilbert-space embedding. Weights are written to
   `exchange_negative_type_certificate.csv`.
3. Split and anisotropic conics over `F_q`: the exact endpoint rate
   `C(anisotropic -> split) = log(q+1)/log(2q-1)` and the expansion
   `C(split -> anisotropic) = 1 - kappa/(q log q) + O(q^-2)` with
   `kappa = max_beta (2 beta - 2^beta)/(beta + 1) = 0.068755890904...`.
   Rates are written to `finite_field_conic_rates.csv`.

## Finite-field exchange matrix

`analysis/finite_field_exchange_matrix.py` computes the exchange matrix of the
quadratic-map classes over `F_q`:

```bash
python analysis/finite_field_exchange_matrix.py
```

It tabulates `C(g -> f)` for the four non-degenerate signature classes `S`
(pure square, or `x^2+x` in even characteristic), `L` (linear and parabolic),
`A` (anisotropic) and `X` (split), checks the eight closed forms attained at an
endpoint, reports the two interior constants

```text
lambda = max_beta (beta + 1 - 2^beta)/(beta + 1) = 0.057915307318...   (X -> L)
kappa  = max_beta (2 beta - 2^beta)/(beta + 1)   = 0.068755890904...   (X -> A)
```

and verifies the comparison order `S < L < A < X`, which holds for every
tabulated `q >= 4`; `q = 3` is the single exception, with `L < S`. Output goes
to `finite_field_exchange_matrix.csv`.

## Frobenius traces and exchange rates

`analysis/frobenius_exchange_rates.py` computes the fiber signatures of two
elliptic fibrations over `F_q` and their exchange rates:

```bash
python analysis/frobenius_exchange_rates.py
```

It checks that `sum_c a_c = 0` exactly, tabulates the normalised second moment
`q^-2 sum_c a_c^2` (about 1 for the large-monodromy family `y^2 = x^3 + x + c`,
about 2 for the CM family `y^2 = x^3 + c` when `q = 1 mod 3`, and 0 when
`q = 2 mod 3` where every fiber is supersingular), verifies that `Z_f(k)` counts
the points of the k-fold fiber power, and confirms

```text
(1 - C(linear -> f)) * sqrt(q) * log q  ->  2g
```

for genus-g fibers. It also reproduces the observation that
`C(y^2 - x^3 -> linear) = 1` in both directions exactly when `q = 2 mod 3`.
Output goes to `frobenius_exchange_rates.csv`.

## Finite-field figures

`analysis/plot_finite_field_exchange.py` draws two figures per field order: the
classical degeneration poset beside the total order the exchange rate puts on
the fiber signatures, and one exchange rate as a containment of energy-entropy
regions.

```bash
python analysis/plot_finite_field_exchange.py
```

By default it emits `q = 2, 4, 8, 9, 16, 25`, matching the poset gallery;
pass orders explicitly to choose others.

```bash
python analysis/plot_finite_field_exchange.py 16 25
```

Files are written into `paper_finite_fields_maps/images/` as
`quadratic-map-exchange-order-q{q}.svg` and
`quadratic-map-gibbs-regions-q{q}.svg`, hand-authored in the same style as the
poset gallery, with no plotting dependency. Three layouts are used. Odd `q` has six classes,
two of which share the flat signature `L`. Even `q` has seven, three of which
share it, because Frobenius is a bijection and `x^2` becomes flat. `q = 2` has
only three classes, since `x^2 = x` as a function, and its poset is already a
chain; the contact in its region figure sits at `T = 0` rather than `T = \u221e`,
which the annotation detects from the attained endpoint rather than assuming.

## Cycles in the tensor families over F_3

`analysis/tensor_cycles_f3.py` searches the homogeneous tensor families over
`F_3` for strict three-cycles of the exchange comparison:

```bash
python analysis/tensor_cycles_f3.py
```

It computes the fiber signature of every orbit representative, restricts to
non-special signatures, and reports cycles together with pairs on which the
endpoint index `phi = log(#fibers) * log(max fiber)` fails. Five of the families
have neither. Quadratic homogeneous maps `F_3^3 -> F_3^3` -- the fifty-orbit
case -- have seven distinct three-cycles and twenty-one `phi`-violating pairs,
and every cycle is closed by exactly one `phi`-violating edge whose rate is
attained at an interior temperature. Output goes to `tensor_cycles_f3.csv`.

## Cycle search in the two large cubic families

`analysis/cubic_cycles_search.py` covers the families with no computed orbit
list, by sampling the map space and keeping one witness map per fiber signature:

```bash
python analysis/cubic_cycles_search.py
```

Cubic maps over `F_8` yield 35 distinct signatures -- a count that does not move
between 150k, 600k and 2.4M samples -- with four `phi`-violating pairs and no
cycle. Cubic homogeneous maps `F_3^3 -> F_3^3` yield 586 strict three-cycles,
the widest with a minimum margin of `3.9e-2`; 578 of them are closed by exactly
one `phi`-violating edge and the other eight by two, which is the most a cycle
can have. Output goes to `cubic_cycles_search.csv`.

## The completed zeta function as a resource

`analysis/xi_gibbs_curve.py` verifies Riemann's kernel
`xi(1/2+b) = int Phi(u) e^{bu} du` with `Phi > 0` even, checks that `log Z` is
convex for every real `b`, and renders the energy-entropy curve
(`riemann-xi-gibbs-curve.svg`). `analysis/xi_versus_euler_factors.py` computes
exchange rates between `xi` and truncated Euler factors `{1, p, ..., p^K}`:

```bash
python analysis/xi_gibbs_curve.py
python analysis/xi_versus_euler_factors.py
```

Unrestricted, both rates vanish -- `xi` has an unbounded spectrum with a
different growth exponent. On a temperature window they are finite, and as
`K -> infinity` the product of the two rates converges to a limit independent of
`p` and `K`, equal to the ratio of the extreme chord slopes of `log Z_xi`. That
limit is the irreversibility floor: `d = 2.2837, 2.2364, 2.1143` on the windows
`[0.5,5]`, `[1,10]`, `[2,20]`.

`analysis/xi_versus_igusa_profiles.py` runs the same comparison against p-adic
Igusa profiles `Phi = -log integral |f|_p^s`, after reproducing the local-currency
note's constant `C_zeta(x^2 -> x^2-y^2) = 0.939702787545916` as a check. Igusa
profiles are linear at `s = 0` with slope `log p * E[v_p(f)]` and *saturate* at
`-log measure{|f|_p = 1}`, so both unrestricted rates vanish again and the
windowed irreversibility is markedly worse than for the Euler factors:
`d = 3.62, 4.00, 4.20` against `2.28, 2.24, 2.11` on the same windows. Copies of
a map cancel from the product, so the number measures shape alone.

## Load-bearing verifiers for the paper claims

Two self-contained scripts verify the claims the second paper makes that came
out of the research in `research/m_and_e_and_a_c/`. The exploratory searches stay
in that directory; these only verify.

```bash
python analysis/negative_type_certificate.py
python analysis/frobenius_bottleneck.py
```

The first checks the five-signature witness that the exchange metric is not of
negative type (`x^T D x = +9.81e-4`), recomputes the distance matrix on an
independent `2e6`-point beta-grid out to `beta = 600` (the violation is
invisible to grids truncated below ~500), confirms the triangle inequality, the
four-point minimality bound `MET_4 = CUT_4`, and the pentagonal violation that
places the metric in `MET \ HYP`.

The second checks the two directions of the rate against a linear map over
`F_q`: the forward rate is always the endpoint `log q / log(max_c N_c)`; the
reverse has the universal bottleneck `beta* = sqrt(2) - 1` with
`1 - C(f->L) = (3 - 2 sqrt 2) m_2 / (2 q log q)`, independent of family and
genus; and a signature is flat exactly when `P` is a permutation polynomial.
