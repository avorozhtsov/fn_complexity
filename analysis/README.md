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
are `C(g | f)`. Outputs:

- `exchange_matrix_small.csv` with 15-digit values;
- `exchange_matrix_small.md` with a readable six-decimal table.

The ordering uses the directed comparison requested for the project:

```text
b -> a  when  C(a|b) > C(b|a).
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

The edge orientation remains `b -> a` when `C(a|b) > C(b|a)`. Equal rates
create no edge. The special zero-complexity signature `{1}` produces infinite
rates in its implemented column, matching the operational definition. Rates
are read through `src/fn_complexity/exchange_rates_cache.json`; only missing
entries are computed, and the updated cache is saved atomically. Existing
values from a previous matrix CSV are imported once before that CSV is
replaced.
