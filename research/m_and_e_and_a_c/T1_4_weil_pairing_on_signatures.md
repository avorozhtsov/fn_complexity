# T1.4 — The Weil pairing evaluated on signatures

**Status: theorem, derived and numerically confirmed.** This is the sharpest
form of the M ↔ E connection found so far, and it replaces the earlier
"dictionary" (a shape analogy) with a computation on *the same objects*.

## The two matrices live on one family

A signature \(a=(a_1,\ldots,a_r)\) **is** a Weil test measure. The atomic
measure \(\sum_i\delta_{a_i}\) on \(\mathbb R_{>0}\) has Mellin transform

\[
\widehat{\mu_a}(s)=\sum_i a_i^{\,s}=Z_a(s),
\]

which is exactly the partition function of the first paper. So both matrices of
the project are defined on the same finite family of signatures:

\[
M_{ab}=C(a\to b)=\inf_{\beta\in[0,\infty]}\frac{\log Z_a(\beta)}{\log Z_b(\beta)},
\qquad
E_{ab}=\sum_{\rho}Z_a(\rho)\,\overline{Z_b(\rho)} .
\]

Truncated to the first \(N\) zeros, \(E\) is the Gram matrix of the vectors
\(v_a=(Z_a(\rho_n))_{n\le N}\in\mathbb C^N\) and is positive semidefinite by
construction (confirmed: minimum eigenvalue \(-1.7\cdot10^{-11}\), i.e. zero).

## Theorem (Landau evaluation of the pairing)

Write \(\rho=\tfrac12+i\gamma\), so
\(Z_a(\rho)=\sum_i\sqrt{a_i}\,e^{i\gamma\log a_i}\) and

\[
E_{ab}=\sum_{i,j}\sqrt{a_ib_j}\;S\!\left(\frac{a_i}{b_j}\right),
\qquad
S(x)=\sum_{0<\gamma\le T}x^{\,i\gamma}.
\]

Landau's theorem evaluates \(S\): it equals \(N\) at \(x=1\), and
\(-(T/2\pi)\Lambda(y)/\sqrt y+O(\log T)\) otherwise, where
\(y=\max(x,1/x)\) and \(\Lambda\) is von Mangoldt. Since
\(\sqrt{a_ib_j}/\sqrt y=\min(a_i,b_j)\), everything collapses to

\[
\boxed{\;
E_{ab}=N\cdot O(a,b)\;-\;\frac{T}{2\pi}\cdot A(a,b)\;+\;O(rs\log T),}
\]
\[
O(a,b)=\!\!\sum_{i,j:\,a_i=b_j}\!\! a_i,
\qquad
A(a,b)=\!\!\sum_{i,j:\,a_i\ne b_j}\!\!\min(a_i,b_j)\,
\Lambda\!\left(\frac{\max(a_i,b_j)}{\min(a_i,b_j)}\right).
\]

\(O\) is the **multiset overlap**; \(A\) is a **prime-power ratio correlation**,
nonzero only for pairs of entries whose ratio is exactly \(p^k\).

### Numerical confirmation

With the first \(N=1200\) zeros (\(T=1648.27\), \(T/2\pi=262.33\)) on a
sixteen-signature family, over all 136 pairs:

| quantity | value |
|---|---:|
| maximum absolute error | \(49.9\) |
| mean absolute error | \(15.5\) |
| predicted error size \(O(rs\log T)\) | \(\sim67\) |
| diagonal scale of \(E\) | \(\sim24000\) |

Every error sits inside the Landau remainder. Individual checks: for
\((4,2)\) against \((8,4)\) the prediction \(3345.33\) against the actual
\(3344.96\); for \((2,2)\) against \((6,1)\), where the ratios \(6/2=3\) and
\(2/1=2\) are both prime, prediction \(-940.07\) against actual \(-926.59\).

## Consequence: Weil-orthogonality is generic, and unrelated to the exchange metric

If two signatures share no entry **and** no pair of their entries has a
prime-power ratio, then \(O=A=0\) and the pair is **Weil-orthogonal** — while
being perfectly comparable in the exchange sense. Examples, against a scale of
\(9600\):

| pair | \(E_{ab}\) | exchange distance \(d\) |
|---|---:|---:|
| \((4,2)\) vs \((5,3)\) | \(3.14\) | \(0.163\) |
| \((4,2)\) vs \((5,5)\) | \(-3.59\) | \(0.267\) |
| \((4,2)\) vs \((7,5)\) | \(-5.80\) | \(0.364\) |
| \((2,2)\) vs \((7,5)\) | \(12.62\) | \(1.032\) |

The first row is the sharpest: \((4,2)\) and \((5,3)\) are among the *closest*
pairs in the exchange metric and are *orthogonal* in the Weil geometry.

Across the family, the correlation between the exchange distance \(d(a,b)\) and
the Weil angle \(\arccos R_{ab}\) is

\[
\operatorname{corr}\bigl(d,\ \text{angle}\bigr)=+0.19 ,
\]

i.e. **essentially none**. The angles cluster at \(\pi/2\): \((10,5)\) against
\((11,7)\) has \(d=0.065\), the closest pair in the family, and angle
\(1.5712\approx\pi/2\), orthogonal.

## What this settles for the M ↔ E question

> \(M\) reads the fiber sizes through their **logarithms** — a real,
> order-theoretic, metric datum, finite for every pair. \(E\) reads them through
> **multiplicative coincidences** — equality of entries and prime-power ratios —
> and is zero for a generic pair. The two matrices are built from the same data
> and extract disjoint parts of it.

That is why the earlier "dictionary" could not be upgraded to a derivation:
there is nothing to derive, because the two functionals are close to
independent. It also says precisely what a bridge would have to do — relate
the additive structure of \(\{\log a_i\}\) to the multiplicative structure of
\(\{a_i/b_j\}\), which is the explicit formula itself.

The appearance of \(\Lambda\) is not an accident of the atomic test measures:
it is the explicit formula in the smallest possible instance, one delta per
fiber.

## Caveats

* Atomic measures are **not** admissible Weil test functions (no smoothness or
  decay), so \(E\) here is the finite-rank truncation at \(N\) zeros. It is a
  legitimate matrix — PSD, Hermitian, computable — but the \(N\to\infty\)
  behaviour is not addressed, and \(|Z_a(1/2+i\gamma)|\) does not decay, so the
  full sum does not converge. Smoothing the deltas would restore admissibility
  at the cost of blurring the \(\Lambda\) term.
* Landau's theorem is applied for fixed \(x\) as \(T\to\infty\); the family here
  has \(r s\le 9\) pairs, so the uniformity is not at issue, but a family with
  many entries would need care.
* Only \(\gamma>0\) is summed; \(E\) is Hermitian and the reported values are
  real parts.

Reproduce with `research/m_and_e_and_a_c/weil_pairing_on_signatures.py`
(zeros cached in `zeta_zeros_1200.npy`).

## Sequel

T1.5 (`T1_5_multiplicative_design.md`) asks whether the \(+0.19\) is an
artefact of choosing generic signatures. It is: designed families reach
\(+0.86\) and \(-0.98\). But the positive branch is produced by the overlap
term alone — an \(A\equiv0\) control ladder matches the geometric ladders cell
for cell — and the \(\Lambda\) term produces agreement of the opposite sign.
Two exact invariances (\(R\) under \(a\mapsto\lambda a\); \(d\) under
\(a\mapsto a^{\otimes k}\)) show no functional relation between the two
geometries can exist.
