# Local zeta profiles and Riemann-zeta currency

This note records a zeta-profile diagnostic for maps over \(p\)-adic fields.
It is separate from the operational affine exchange rate
\(C_{\rm aff}\): no theorem below identifies the zeta-profile quantity with
an achievable affine conversion rate.

## Normalized local profile

Fix \(|p|_p=p^{-1}\), the standard lattice \(\mathbb Z_p^n\), and additive
Haar measure normalized by

\[
\mu(\mathbb Z_p^n)=1.
\]

For a primitive integral polynomial map
\(f\colon\mathbb Z_p^n\to\mathbb Z_p\), define

\[
Z_{f,p}(s)
=\int_{\mathbb Z_p^n}|f(x)|_p^s\,dx,
\qquad
\Phi_{f,p}(s)=-\log Z_{f,p}(s),
\qquad s>0.
\]

The zeta-profile exchange diagnostic is

\[
C_\zeta(f\to g)
=\inf_{s>0}\frac{\Phi_{f,p}(s)}{\Phi_{g,p}(s)}.
\]

The base of the logarithm is immaterial. These profiles are invariant under
\(\operatorname{GL}_n(\mathbb Z_p)\) input changes and multiplication of the
output by a \(p\)-adic unit. They are not invariant under unrestricted
\(\operatorname{GL}_n(\mathbb Q_p)\) changes: such a change can move the
chosen lattice. Thus a lattice or canonical primitive integral representative
is part of the definition.

## The two forms over \(\mathbb Q_2\)

Put \(t=2^{-s}\). For the degenerate quadratic form \(R(x,y)=x^2\), the
second coordinate integrates out and

\[
\begin{aligned}
Z_{R,2}(s)
&=\sum_{n\ge0}\Pr(v_2(x)=n)t^{2n}\\
&=\sum_{n\ge0}2^{-n-1}t^{2n}
=\frac1{2-t^2}.
\end{aligned}
\]

For the split form \(S(x,y)=x^2-y^2\), separating the cases in which
\(x,y\) have opposite parity, are both even, or are both odd gives

\[
Z_{S,2}(s)=\frac{t^2-2t+2}{(2-t)^2}.
\]

Consequently

\[
\boxed{
C_\zeta(x^2\to x^2-y^2)
=\inf_{0<t<1}
\frac{\log(2-t^2)}
{\log\!\left((2-t)^2/(t^2-2t+2)\right)}
}.
\]

Numerical minimization gives the interior minimizer

\[
t_*=0.696541929482172327177326187\ldots,
\qquad
s_*=-\log_2t_*=0.521717894056140064532492062\ldots,
\]

and hence

\[
\boxed{
C_\zeta(x^2\to x^2-y^2)
=0.939702787545916334905213210\ldots .
}
\]

This decimal is the numerical minimum of the displayed elementary
transcendental function. It is not a named mathematical constant and is not,
without an additional coding theorem, the operational affine exchange rate.

## The Riemann local factor

The Riemann zeta function has Euler product

\[
\zeta(s)=\prod_p\zeta_p(s),
\qquad
\zeta_p(s)=\frac1{1-p^{-s}}.
\]

The normalized local profile of the identity map \(u(x)=x\) is

\[
U_p(s):=Z_{x,p}(s)
=\int_{\mathbb Z_p}|x|_p^s\,dx
=\frac{1-p^{-1}}{1-p^{-s-1}}
=(1-p^{-1})\zeta_p(s+1).
\]

Equivalently, if multiplicative Haar measure is normalized by
\(\operatorname{vol}(\mathbb Z_p^\times)=1\), then

\[
\zeta_p(s)
=\int_{\mathbb Q_p^\times}
\mathbf 1_{\mathbb Z_p}(x)|x|_p^s\,d^\times x.
\]

The product of these multiplicative local integrals gives the finite-prime
part of \(\zeta(s)\). The real Gaussian integral supplies the gamma factor of
the completed zeta function in the adelic formulation of
[Tate's thesis](https://sites.math.rutgers.edu/~alexk/2023S572/Tate1950.pdf).

The positive local currency cost curve appropriate to the present additive
normalization is

\[
\Phi_{U_p}(s)
=-\log\!\left((1-p^{-1})\zeta_p(s+1)\right).
\]

The raw quantity \(-\log\zeta_p(s)\) is negative for real \(s>0\), so the
normalization and shift are essential for the exchange interpretation.

## Cases in which this is an exact building-block currency

For a monomial in independent coordinates,

\[
Z_{x_1^{a_1}\cdots x_r^{a_r},p}(s)
=\prod_{j=1}^r U_p(a_js),
\qquad
\Phi_{x_1^{a_1}\cdots x_r^{a_r},p}(s)
=\sum_{j=1}^r\Phi_{U_p}(a_js).
\]

Thus local Riemann-zeta factors form an exact common currency for monomial
and normal-crossing profiles.

For every odd prime \(p\), the integral change of variables

\[
u=x-y,\qquad v=x+y
\]

has unit determinant. It follows that

\[
Z_{x^2-y^2,p}(s)=U_p(s)^2
=(1-p^{-1})^2\zeta_p(s+1)^2,
\qquad p\ne2.
\]

At the bad prime \(2\), the determinant is not a unit. Here

\[
U_2(s)=\frac1{2-t}=\frac12\zeta_2(s+1)
\]

and the exact correction is

\[
Z_{x^2-y^2,2}(s)
=(t^2-2t+2)U_2(s)^2
=\frac{t^2-2t+2}{4}\zeta_2(s+1)^2.
\]

Similarly,

\[
Z_{x^2,2}(s)=U_2(2s)=\frac12\zeta_2(2s+1).
\]

The polynomial \(t^2-2t+2\) is therefore the precise bad-prime correction
which prevents the split quadratic form from being two unmodified units of
local Riemann-zeta currency at \(p=2\).

## Profile-valued rather than scalar prices

Using \(U_p\) as a numeraire, define the whole price curve

\[
P_f(s)=\frac{\Phi_{f,p}(s)}{\Phi_{U_p}(s)}.
\]

Then the pairwise diagnostic is recovered exactly by

\[
C_\zeta(f\to g)
=\inf_{s>0}\frac{P_f(s)}{P_g(s)}.
\]

A single scalar price suffices only when the relevant profiles are
proportional. In general the universal quote is the entire function \(P_f\),
analogous to a yield curve. Igusa local zeta functions are rational functions
of \(p^{-s}\), but their numerator corrections and sums of rational pieces
mean that ordinary Riemann factors alone are not a scalar universal currency;
see the [computational account of Igusa rationality](https://arxiv.org/abs/2006.08926).

The global, complex-valued extension of this currency viewpoint and its exact
relation to the Riemann hypothesis are recorded in
[Riemann hypothesis as positivity of exchange matrices](riemann_hypothesis_exchange_matrices.md).
