# Attempts at affine exchange rates over \(\mathbb Q_2\)

Let

\[
q_d(x,y)=x^2+d y^2
\]

and define

\[
k_{d\to e}(n)=\max\{k:q_e^k\text{ is implemented by }q_d^n
\text{ using affine processors over }\mathbb Q_2\}.
\]

The affine exchange rate is the operational limit

\[
C_{\rm aff}(q_d\to q_e)=\lim_{n\to\infty}\frac{k_{d\to e}(n)}n.
\]

Block-diagonal composition makes \(k_{d\to e}\) superadditive, so the limit
exists. The experiments below concern the pair

\[
q_1=x^2+y^2,
\qquad
q_2=x^2+2y^2.
\]

They have different determinant square classes and therefore are not
one-copy equivalent.

## Attempt 1: exact elementary blocks

Two copies of either form implement one copy of the other: restrict the two
source inputs to \((x,0)\) and \((y,0)\), obtaining \((x^2,y^2)\), and take
the required output combination. Consequently

\[
k_{1\to2}(2r)\ge r,
\qquad
k_{2\to1}(2r)\ge r,
\]

and both affine rates are at least \(1/2\).

For distinct nondegenerate classes, an \(n\)-copy to \(n\)-copy conversion is
impossible. Indeed, equality of dimensions forces the input and output linear
maps to be invertible. The determinant of the generic target pencil is a
constant times \(\prod_i t_i^2\). Equality with the transformed source pencil
forces the output transformation to permute and rescale the coordinate
hyperplanes. Each target block would then be similar to one source block,
contradicting the distinct determinant square classes. Hence

\[
k_{1\to2}(n),k_{2\to1}(n)\le n-1.
\]

This excludes every finite rate-one conversion, but does not exclude an
asymptotic rate of one.

## Attempt 2: the first cross-copy case modulo \(4\)

The first possible improvement is

\[
q_1^3\longrightarrow q_2^2,
\]

which would prove a lower bound of \(2/3\). We exhaustively tested the
quadratic part of this identity over \(\mathbb Z/4\mathbb Z\), with standard
integral lattices. For each source copy the search enumerates every matrix

\[
A_i\in\operatorname{Mat}_{2\times4}(\mathbb Z/4\mathbb Z)
\]

and every pulled-back quadratic form \(q_d(A_i z)\), then every pair of output
coefficients. A meet-in-the-middle three-sum tests the two target forms
simultaneously.

| direction | distinct pullbacks | contribution terms | pair sums | solution |
|---|---:|---:|---:|---:|
| \(q_1^3\to q_2^2\) | 136 | 1,486 | 777,226 | none |
| \(q_2^3\to q_1^2\) | 256 | 1,486 | 193,216 | none |

This is an exact obstruction for processors integral in the displayed
lattices. It is not yet an obstruction over all of \(\mathbb Q_2\): a
rational processor may use denominators, equivalently changing the source or
target lattices before reduction.

## Attempt 3: one shared source summand

A natural \(3\)-to-\(2\) construction tries to share one source form \(H\):

\[
F_1=R_1+aH,
\qquad
F_2=R_2+bH,
\]

where \(F_1,F_2\) are the two \(q_2\) target blocks and \(H,R_1,R_2\) are
rank-two forms of source type \(q_1\). Rank at most two for both residuals
forces the restrictions of \(H\) to the two target blocks to have rank one.
Writing \(\Delta\) for the determinant square class of \(H\), a direct
Schur-complement calculation gives

\[
[\det R_1]=[\det R_2]=[-2/\Delta].
\]

For a source-type shared form \([\Delta]=[1]\), both residuals therefore have
class \([-2]\), not \([1]\). Thus this whole shared-summand family fails. It
does not cover a general \(2\times3\) output-mixing matrix, so it is a rigorous
negative result only for this ansatz.

## Attempt 4: finite residue signatures

Reducing the two maps modulo \(2^m\) produces finite maps. Their ordinary
fiber-signature exchange rates are straightforward to compute, but depend on
the chosen integral lattice and allow arbitrary finite-set processors rather
than affine processors.

| \(m\) | \(C_m(q_1\to q_2)\) | \(C_m(q_2\to q_1)\) |
|---:|---:|---:|
| 2 | 0.792481250361 | 0.666666666667 |
| 4 | 0.954242509439 | 0.991176112003 |
| 6 | 0.991534341727 | 0.998646536327 |
| 8 | 0.998413561081 | 0.999748319425 |
| 10 | 0.999688024570 | 0.999949423599 |
| 12 | 0.999936014975 | 0.999989399891 |
| 14 | 0.999986455906 | 0.999997716974 |

The first direction has the exact finite-level formula

\[
C_m(q_1\to q_2)
=\frac{\log(2^{m-1}+1)}{\log(2^{m-1}+2)}.
\]

Both residue rates converge rapidly to one. This agrees with the fact that
the normalized local-zeta data are too coarse to distinguish the pair; it is
not evidence for an affine construction with rate one.

## Current conclusion

The computations improve the evidence, but not the global bound:

\[
\boxed{\frac12\le C_{\rm aff}(q_1\to q_2)\le1},
\qquad
\boxed{\frac12\le C_{\rm aff}(q_2\to q_1)\le1}.
\]

No \(3\)-to-\(2\) construction was found. The exhaustive modulo-four result
and the shared-summand calculation explain why the most immediate attempts
fail, while still leaving general cross-copy mixing over \(\mathbb Q_2\)
open.

Regenerate the numerical tables and the exact modulo-four searches with

```bash
python analysis/p_adic_exchange_attempts.py
```

The generated files are the complete
[\(m=2,\ldots,14\) convergence sequence](../../analysis/p_adic_residue_rates_q1_q2.csv),
the [\(7\times7\) residue-rate matrix](../../analysis/p_adic_residue_rate_matrix_m10.csv)
for all anisotropic square classes at \(m=10\), and the
[modulo-four search counts](../../analysis/p_adic_integral_3_to_2_mod4.csv).

For a distinct zeta-profile diagnostic, including the explicit value
\(C_\zeta(x^2\to x^2-y^2)=0.939702787545916\ldots\) and its relation to the
local Euler factors of the Riemann zeta function, see
[Local zeta profiles and Riemann-zeta currency](p_adic_zeta_exchange_currency.md).
