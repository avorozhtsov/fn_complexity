# fn-complexity

For an onto map between finite sets, its **signature** is the decreasing
multiset of its non-empty fiber sizes:

```text
(n_1, ..., n_m),  n_i = |f^{-1}(y_i)|.
```

If the input and output relabellings are injective, a signature `a` is
implemented by a signature `b` exactly when the sorted fibers of `a` fit,
one-by-one, into distinct sorted fibers of `b`.

For Cartesian powers, fiber sizes multiply.  The operational definition is

```text
C(g -> f) = lim  max{k : f^k is implemented by g^n} / n.
            n→∞
```

The entropy method proves the equivalent computational formula

```text
                            log(sum_j g_j^beta)
C(g -> f) = inf
          0 <= beta < ∞   log(sum_i f_i^beta),
```

together with the endpoint at infinity,
`log(max(g)) / log(max(f))`.  Thus the infimum is a theorem/algorithm for
computing the operationally defined limit, rather than a second definition.
It is the supporting-line form of the Gibbs entropy-energy construction with
energy levels `-log(n_i)`.

The compiled paper is available as
[exchange_rates_finite_map_signatures.pdf](paper/exchange_rates_finite_map_signatures.pdf).

In map notation the direct product is
`f^k : X_f^k -> Y_f^k`.  For the equation
`f = h_out ◦ g ◦ h_in` to type-check, `h_out` maps the used subset of `Y_g`
bijectively to `Y_f`.

## Requested examples

```python
from fn_complexity import exchange_rate, exchange_rate_result

exchange_rate_result((2, 2), (3, 1))
# ExchangeRateResult(rate=0.630929...,
#                    beta=inf)

exchange_rate_result((3, 1), (2, 2))
# ExchangeRateResult(rate=0.965384...,
#                    beta=0.403680...)
```

The first displayed result follows directly: every fiber of `(2, 2)^n` has
size `2^n`, so `(3, 1)^k` fits iff both `2^k <= 2^n` and
`3^k <= 2^n`.

## Exact finite-power check

```python
from fn_complexity import k_max

k_max((2, 2), (3, 1), 20)
# 12
```

`k_max` uses a compressed `fiber_size -> multiplicity` representation and is
intended as an independent verifier for modest powers.

## CLI and plots

`kmax_cli` uses the argument order `C(g -> f)`: implementer first, implemented
map second.

```bash
./cli/kmax_cli 3,1 2,2 --n-max 120
```

It prints exact `k_max(n)` values and creates a dependency-free SVG convergence
plot of `c_n = k_max(n) / n` in `images/` by default. The asymptotic exchange
rate is shown as a horizontal reference line. The CLI enumerates slopes from
the continued-fraction convergents of the rate and draws one maximum-coverage
curve `c_n = e - d/n` for every convergent, using a different color for each
value of `e`. Exact data points use neutral diamond markers because selected
curves can overlap and need not cover every point. When
the largest-fiber endpoint determines the rate, the legend also shows its analytic
`log(max(g)) / log(max(f))` value. Use `--output result.svg` to select a file,
`--output some-directory/` to select a directory, and `--quiet` to suppress
the table. For large ranges, `--max-convergents 7` restricts the overlay to the
first seven convergents. The exact polyline always contains every requested
point; visible markers are automatically sampled above 1,200 values. Use
`--log-x` when plotting several orders of magnitude in `n`.

The entropy-energy curve itself is available as:

```python
from fn_complexity import gibbs_point

point = gibbs_point((3, 1), temperature=2.0)
point.energy, point.entropy
```

Here `p_i = n_i ** (1 / T) / sum_j n_j ** (1 / T)`, energy is
`sum_i p_i * (-log(n_i))`, and entropy is `-sum_i p_i * log(p_i)`.
Throughout the project, `T` always denotes physical temperature and
`beta = 1 / T` denotes inverse temperature.

## Cached exchange matrix

`analysis/exchange_matrix_extended.py` computes the requested 97-signature
matrix through `ExchangeRateCache`. Its 9,409 directed rates are stored in
`src/fn_complexity/exchange_rates_cache.json` and are reused on later runs.
The same cache also stores the additional comparisons used to construct the
69-signature table in Appendix B, so the cache can contain more entries than
the matrix itself.
The cache records the numerical algorithm and precision; a mismatch is
rejected rather than silently mixing values computed under different
settings.

## Comparison-cluster CLI

`cluster_cli` searches the strongly connected component of a non-special
signature for arrows
`a -> b` when `C(a -> b) >= C(b -> a)`, so arrows point toward the
less-complex signature. It uses the exhaustive finite shells
`S_B = {a : len(a) >= 2, a_1 > 1, a_1 + 2 len(a) <= B}` and the same
vectorized partition-function screening used by the Appendix B analysis.

```bash
./cli/cluster_cli 3,1,1 --n-max 100 --max-b 18
```

With `--n-max`, shells are searched in increasing order and the command stops
when that many distinct members have been found, or at `--max-b`. Without
`--n-max`, it prints the complete component in `S_MAX_B`. Near-equal numerical
comparisons receive arrows in both directions; tune this with `--tolerance`.
Pass `--strict` to replace the default relation
`C(a -> b) >= C(b -> a)` by `C(a -> b) > C(b -> a)`. In strict mode, differences
no larger than `--tolerance` produce no arrow.
The special CLI requires the analysis dependencies:

```bash
python -m pip install -e '.[analysis]'
```

## Development

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

Signatures must be non-empty collections of positive integers.  The signature
`(1,)` is the zero-complexity identity: its exchange rate from any non-empty
target is infinite, and consequently its finite `k_max` is unbounded.
