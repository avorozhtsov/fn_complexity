# Riemann hypothesis as positivity of exchange matrices

This note records an exact reformulation of the Riemann hypothesis in an
extended exchange language. It also states why the existing scalar local
quantity \(C_\zeta\) is not sufficient by itself.

## Why the local scalar matrix does not see RH

The local diagnostic

\[
C_\zeta(f\to g)
=\inf_{s>0}\frac{-\log Z_f(s)}{-\log Z_g(s)}
\]

uses only positive real \(s\) and compresses two complete profiles to one real
number. The Riemann hypothesis concerns complex zeros of the global completed
zeta function

\[
\xi(s)
=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s).
\]

It depends on the interaction of every finite prime with the archimedean
factor and on the reflection \(s\leftrightarrow1-s\). Neither that complex
phase information nor the global functional equation is retained by one
fixed-prime, real-axis infimum. No equivalence with RH is currently claimed
for the finite-field or \(p\)-adic quadratic-map rate matrices in this
project.

## Adelic dilation portfolios

Let \(g\in C_c^\infty(\mathbb R_{>0})\) be an admissible test profile and let

\[
\widehat g(s)=\int_0^\infty g(x)x^{s-1}\,dx
\]

be its Mellin transform. One can interpret \(g\) as a virtual portfolio of
dilation processors

\[
D_a:x\longmapsto ax,
\qquad
T_g=\int_0^\infty g(a)D_a\,\frac{da}{a}.
\]

Mellin transformation diagonalizes these multiplicative-scale processors:
the spectral response of \(T_g\) at \(s\) is \(\widehat g(s)\). These are
linear combinations of processors, not individual polynomial maps.

## The Weil exchange form

For an admissible profile \(g\), define the Weil quadratic form

\[
Q_W(g)
=\sum_\rho
\widehat g(\rho)
\overline{\widehat g(1-\overline\rho)},
\]

where \(\rho\) ranges over the nontrivial zeros of \(\zeta\), with
multiplicity and the standard symmetric convergence convention. Let
\(B_W(g,h)\) be the corresponding sesquilinear form obtained by
polarization.

For a finite family of exchange profiles \(g_1,\ldots,g_N\), define

\[
E(g_1,\ldots,g_N)
=\bigl(B_W(g_i,g_j)\bigr)_{i,j=1}^N.
\]

Then Weil's positivity criterion gives the exact matrix equivalence

\[
\boxed{
\mathrm{RH}
\quad\Longleftrightarrow\quad
E(g_1,\ldots,g_N)\succeq0
\text{ for every finite admissible family }g_1,\ldots,g_N.
}
\]

This is a restatement, by polarization, of
[Weil's explicit-formula positivity criterion](https://doi.org/10.1070/IM1972v006n01ABEH001866),
not a new proof of RH.

If RH holds, then every nontrivial zero satisfies
\(1-\overline\rho=\rho\). Therefore

\[
B_W(g_i,g_j)
=\sum_\rho
\widehat g_i(\rho)\overline{\widehat g_j(\rho)},
\]

so \(E\) is the Gram matrix

\[
E=\sum_\rho v_\rho v_\rho^*,
\qquad
v_\rho=
\begin{pmatrix}
\widehat g_1(\rho)\\
\vdots\\
\widehat g_N(\rho)
\end{pmatrix}.
\]

It is consequently positive semidefinite. Conversely, Weil's criterion says
that positivity for every test profile forces all the nontrivial zeros onto
the critical line. If RH is false, some admissible profile has
\(Q_W(g)<0\), and a corresponding finite exchange matrix has a negative
eigenvalue.

The explicit formula makes this criterion non-circular: the same entries can
be evaluated from prime-power terms weighted by the von Mangoldt function,
the archimedean gamma contribution, and the pole and normalization terms,
rather than from an assumed list of zeros.

## Normalized exchange coefficients

When the diagonal entries are positive, define

\[
R_{ij}=\frac{E_{ij}}{\sqrt{E_{ii}E_{jj}}}.
\]

Under RH, every such finite matrix satisfies

\[
R\succeq0,
\qquad
R_{ii}=1,
\qquad
|R_{ij}|\le1.
\]

Thus RH may be phrased as consistency of every finite collection of adelic
exchange portfolios with a Hilbert-space, or Gram-matrix, valuation. The
matrix is Hermitian and behaves more like an exchange covariance matrix than
like the directed scalar matrix \((C_\zeta(f_i\to f_j))\).

## Li's diagonal alternative

A simpler but less structural matrix reformulation follows from Li's
criterion. Define

\[
\lambda_n
=\left.
\frac1{(n-1)!}\frac{d^n}{ds^n}
\left[s^{n-1}\log\xi(s)\right]
\right|_{s=1}.
\]

Li proved that RH is equivalent to
\(\lambda_n\ge0\) for every \(n\ge1\). Equivalently,

\[
L_N=\operatorname{diag}(\lambda_1,\ldots,\lambda_N)\succeq0
\quad\text{for every }N.
\]

See [Li, *The positivity of a sequence of numbers and the Riemann
hypothesis*](https://doi.org/10.1006/jnth.1997.2137). This diagonal encoding
is exact but does not expose pairwise exchange structure as Weil's matrices
do.

## Open bridge to polynomial-map exchange

The exact positivity theorem uses virtual adelic dilation portfolios. To turn
it into a theorem solely about the actual polynomial maps studied elsewhere
in this project, one would need to prove that their global zeta profiles, or
a natural closure of them under allowed processors, generate a sufficiently
rich class of Weil test functions. No such density or representation theorem
is established here.

Accordingly, the defensible conclusion is that RH is exactly positivity of
every finite Weil exchange matrix for adelic dilation portfolios. Whether an
equally exact criterion exists inside the restricted category of polynomial
maps and the original operational exchange rates remains open.

The local Riemann-factor numeraire that motivates this global formulation is
developed in
[Local zeta profiles and Riemann-zeta currency](p_adic_zeta_exchange_currency.md).

The structural comparison between this Hermitian matrix \(E\) and the directed
exchange matrix \(M_{ij}=C(f_i\to f_j)\) --- including the sup-norm
representation of \(M\), the unconditional failure of the corresponding
positivity for finite maps, and the finite-field register in which Weil numbers
appear directly in the exchange rates --- is in
[Two positivities: the exchange matrix and the Weil matrix](exchange_positivity_and_weil.md).
