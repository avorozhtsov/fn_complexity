# Frobenius traces, fiber powers, and the exchange matrix

The exchange matrix over \(\mathbb F_q\) is built from fiber point counts, and
point counts over a finite field are Weil numbers. This note makes the
consequence precise: **the partition function of a map is the point count of its
fiber powers, its expansion in \(q\) is a generating function for the moments of
the Frobenius traces, and the exchange rate reads off the extreme Frobenius
trace --- whose approach to the Weil edge encodes \(\dim USp(2g)\) --- and, in one
clean example, a congruence on \(q\)**.

Everything is reproduced by `analysis/frobenius_exchange_rates.py`.

## The partition function counts fiber powers

Let \(f\colon X\to Y\) be a map of finite sets with fibers of sizes \(N_c\).
For an integer \(k\ge1\),

\[
Z_f(k)=\sum_c N_c^{\,k}
=\#\bigl\{(P_1,\ldots,P_k)\in X^k:\ f(P_1)=\cdots=f(P_k)\bigr\}
=\#\bigl(X\times_YX\times_Y\cdots\times_YX\bigr).
\]

So the monotones of the first paper, restricted to integer inverse temperature,
are the point counts of the **fiber powers** of \(f\); the real parameter
\(\beta\) interpolates them. Over \(\mathbb F_q\) those counts are given by
Grothendieck--Lefschetz as traces of Frobenius on the cohomology of the fiber
powers, which is exactly where the monodromy of the family lives. The whole
exchange matrix \(M\) is therefore a functional of the zeta functions of the
fiber powers of the maps involved.

The same fact in generating-function form: the signature is the pole set of

\[
\sum_{k\ge0}Z_f(k)\,T^k=\sum_c\frac1{1-N_cT},
\]

a rational function whose poles are the reciprocals of the fiber sizes. The
largest fiber is the radius of convergence and the number of fibers is the
number of poles --- the two endpoints \(\beta=\infty\) and \(\beta=0\) of the
rate formula, read off the same series. Verified for \(k=2,3\) at \(q=101\) in
the script.

## The expansion in \(q\) is a moment generating function

Take \(f\colon\mathbb A^2\to\mathbb A^1\) whose fibers are curves, write
\(N_c=q-a_c\) for the affine count, and put

\[
m_1=\frac1q\sum_ca_c,\qquad m_2=\frac1{q^2}\sum_ca_c^2 .
\]

Expanding \(Z_f(\beta)=q^\beta\sum_c(1-a_c/q)^\beta\) gives

\[
\boxed{\ \log Z_f(\beta)=(\beta+1)\log q
+\frac1q\Bigl[\tfrac{\beta(\beta-1)}2m_2-\beta m_1\Bigr]+O(q^{-3/2}).\ }
\]

Consequently, for two such maps,

\[
C(f\to g)=1-\frac1{q\log q}\max_{\beta}
\frac{P_g(\beta)-P_f(\beta)}{\beta+1},
\qquad
P(\beta)=\tfrac{\beta(\beta-1)}2m_2-\beta m_1 ,
\]

so **to first order the exchange rate between two maps over \(\mathbb F_q\) is
determined by the first two moments of their Frobenius traces**, and the higher
orders bring in the higher moments in the same way. These are precisely the
quantities Katz--Sarnak equidistribution computes from the monodromy group.

Two of these moments are forced:

* \(\sum_ca_c=0\) **exactly**, for every \(f\colon\mathbb A^2\to\mathbb A^1\),
  because \(\sum_cN_c=q^2\). So \(m_1=0\) always and the leading correction is
  carried by \(m_2\).
* \(m_2\) is the normalised second moment. The script finds

| family | \(q\) | \(q\bmod3\) | \(m_2\) | distinct fiber sizes | \(\min a_c\) vs \(-2\sqrt q\) |
|---|---:|---:|---:|---:|---|
| \(y^2=x^3+x+c\) | 1009 | 1 | \(0.9970\) | 120 | \(-62\) vs \(-63.5\) |
| \(y^2=x^3+x+c\) | 2003 | 2 | \(1.0005\) | 179 | \(-89\) vs \(-89.5\) |
| \(y^2=x^3+c\) | 1009 | 1 | \(1.9980\) | 7 | \(-62\) vs \(-63.5\) |
| \(y^2=x^3+c\) | 2003 | 2 | \(0.0000\) | 1 | \(0\) |

The first family has large monodromy and \(m_2\to1\), the semicircle value. The
second is the family of sextic twists, has complex multiplication, and gives
\(m_2\to2\) when \(q\equiv1\pmod3\) --- twice the generic value, with only seven
distinct fiber sizes instead of a hundred. The exchange rate sees the
difference in the first correction term.

## Three scales, three pieces of geometry

The order of magnitude of \(1-C\) is itself an invariant. Writing \(L\) for the
linear map \(\{q,\ldots,q\}\) and \(X\) for the split conic:

| comparison | \(1-C\) | what it reads |
|---|---|---|
| \(f\) against \(X\) | \(\asymp\log2/\log q\) | a **singular fiber**: the split conic's \(2q-1\) |
| \(L\to f\), *irreducible* fibers of genus \(g\) | \(\asymp 2g/(\sqrt q\,\log q)\) asymptotically | the **extreme trace**; see the caveats below |
| conic against conic | \(\asymp\kappa/(q\log q)\) | the \(\pm1\) of genus-zero counts |

The middle row needs care, and stating it properly turns out to be more
interesting than the loose version.

**Proposition (the rate against a linear map is an endpoint, always).** For
every \(f\colon\mathbb A^2\to\mathbb A^1\) over \(\mathbb F_q\),

\[
C(L\to f)=\frac{\log q}{\log\bigl(\max_cN_c\bigr)}\qquad\text{attained at }\beta=\infty .
\]

*Proof.* \(Z_f(\beta)=\sum_cN_c^\beta\le q\,(\max_cN_c)^\beta\), so with
\(A=\log q\) and \(B=\log\max_cN_c\) we get
\(R(\beta)=(1+\beta)A/\log Z_f\ge(1+\beta)A/(A+\beta B)>A/B\) for every finite
\(\beta\) whenever \(B>A\); and \(R(\infty)=A/B\). \(\square\)

So **this rate sees the largest fiber and nothing else** --- no other trace
statistic enters. Verified to ten decimals against the solver for genus
\(1,2,3\) at \(q=101,211,1009\).

Consequently \((1-C)\sqrt q\log q=\mu-\mu^2(\tfrac12+1/\log q)/\sqrt q+O(1/q)\)
with \(\mu=\max_c(-a_c)/\sqrt q\), and the whole question becomes the
**extreme-value statistics of the traces**. Two things must be said that the
loose statement hides.

*The law needs geometrically irreducible fibers.* A fiber with \(r\) components
has \(N_c\approx rq\), giving \(1-C\asymp\log r/\log q\), which swamps every
Weil-scale term. That is exactly the first row of the table: the split conic's
fiber over \(0\) is two lines. The correct genus-zero control is
\(x^2+y^2\) with \(q\equiv3\pmod4\), where all \(a_c=\pm1\) and
\((1-C)\sqrt q\log q=1/\sqrt q\) exactly; at \(q\equiv1\pmod4\) the form
splits and the law breaks --- a \(q\bmod4\) sensitivity parallel to the
\(q\bmod3\) phenomenon below.

*The limit \(2g\) is asymptotic and, for \(g\ge2\), unreachable in practice.*
The maximum of \(q\) samples approaches the edge of the \(USp(2g)\) trace
distribution only as fast as its lower tail allows,
\(P(2g-T<\varepsilon)\sim K_g\varepsilon^{d/2}\) with
\(d=\dim USp(2g)=2g^2+g\), giving

\[
2g-\mu\;\approx\;\Gamma(1+2/d)\,(qK_g)^{-2/d}.
\]

Fitted against predicted exponents over \(q\in[4\cdot10^3,10^6]\):
\(-0.1969\) vs \(-0.2000\) for \(g=2\), \(-0.0931\) vs \(-0.0952\) for
\(g=3\), \(-0.0583\) vs \(-0.0556\) for \(g=4\). At \(q=10^6\) the observed
\(\mu\) is \(3.684\), \(4.663\), \(4.760\) against targets \(4,6,8\);
reaching within \(10\%\) of \(2g\) would need \(q\approx6\cdot10^4\),
\(6\cdot10^9\), \(4\cdot10^{16}\). For \(g\ge3\) at reachable \(q\) the
observed maximum is a *Gaussian* extreme rather than a Weil edge --- the
crossover sits at \(\log q\approx2g^2\).

So the honest statement is stronger than the original one:

> **The exchange rate against a linear map measures the extreme Frobenius trace,
> normalised by \(\sqrt q\). Its limit is \(2g\), but its *rate of approach*
> encodes \(\dim USp(2g)=2g^2+g\) --- a finer invariant than the genus, read off
> the same curve.**

For \(y^2=x^3+x+c\), where genus \(1\) saturates the integer Weil bound
\(m=\lfloor2\sqrt q\rfloor\) for essentially every prime above \(4000\):

| \(q\) | 101 | 211 | 401 | 1009 | 2003 |
|---|---:|---:|---:|---:|---:|
| \((1-C)\sqrt q\log q\) | \(1.670\) | \(1.769\) | \(1.830\) | \(1.878\) | \(1.935\) |

One more caveat: \(\max_c(-a_c)\) and \(\max_c|a_c|\) are *not* the same at
finite \(q\) --- they agree only \(50\)--\(68\%\) of the time for
\(g=2,3,4\). The rate tracks the former, the largest-fiber side only.

## The reverse rate, and a universal bottleneck temperature

The rate *into* the linear map is the mirror image of the proposition above, and
it is where the trace statistics actually live. Since \(m_1=0\) always,

\[
C(f\to L)=\inf_\beta\frac{\log Z_f}{\log Z_L}
=1+\frac{m_2}{2q\log q}\,\inf_{\beta\ge0}\frac{\beta(\beta-1)}{\beta+1}
+O(q^{-3/2}),
\]

and the infimum is elementary: \(\frac{d}{d\beta}\frac{\beta(\beta-1)}{\beta+1}
=\frac{\beta^2+2\beta-1}{(\beta+1)^2}\), so the minimiser is the positive root of
\(\beta^2+2\beta-1=0\).

**Proposition (universal bottleneck).** For every map \(\mathbb A^2\to\mathbb A^1\)
over \(\mathbb F_q\) with geometrically irreducible fibers,

\[
\boxed{\ \beta_*=\sqrt2-1=0.414213562\ldots,
\qquad
1-C(f\to L)=\frac{(3-2\sqrt2)\,m_2}{2q\log q}+O(q^{-3/2}).\ }
\]

Both constants are independent of the family and of the genus. Measured:

| \(q\) | \(g=1\) | \(g=2\) | \(g=3\) |
|---|---:|---:|---:|
| \(\beta_*\) at \(q=1009\) | \(0.414089\) | \(0.414127\) | \(0.413946\) |
| \((1-C)\,2q\log q/m_2\) at \(q=1009\) | \(0.171802\) | \(0.171671\) | \(0.172101\) |

against \(\sqrt2-1=0.4142136\) and \(3-2\sqrt2=0.1715729\).

**The interior sees the whole moment ladder.** Beyond the leading \(m_2\) term,
\(1-C(f\to L)=\sum_{k\ge2}c_km_k+O(q^{-2})\) with each order damped by a further
factor \(\approx0.6/\sqrt q\). Subtracting moments \(2,\ldots,K\) in turn drops
the residual r.m.s. from \(2.75\cdot10^{-5}\) to \(8.5\cdot10^{-8}\),
\(5.7\cdot10^{-8}\), \(3.7\cdot10^{-9}\), and then plateaus at exactly the
first dropped term (predicted \(3.30\cdot10^{-9}\), observed
\(3.36\cdot10^{-9}\)). So \(m_2\) is recoverable from a *single* rate to
relative \(O(q^{-1/2})\); \(m_3\) only through its damped term; \(m_4\) is
visible but not invertible.

**A sharp separation.** At \(q=211\) there are two hyperelliptic maps with
*identical* image size, *identical* largest fiber and *identical* \(m_2\), whose
smallest fibers differ (\(173\) against \(167\)). All four endpoint probes agree
to \(0.000\cdot10^{0}\); the interior probes separate them by
\(2.42\cdot10^{-7}\), and the separation is predicted in advance by
\(c_3\Delta m_3+c_4\Delta m_4=2.32\cdot10^{-7}\), a \(4\%\) match. **The
interior of the rate curve is not redundant: what it adds is the higher
moments.**

**What is invisible, and why.** The *smallest* fiber --- equivalently
\(\max_ca_c\) --- is not determined: regression \(R^2=0.80\) against
\(1.0000000\) for \(\min_ca_c\). Nothing isolates it, because that would need
\(\beta<0\), and negative \(\beta\) is exactly what the first paper excludes on
the grounds that \(Z_a\) is order-preserving only for \(\beta\ge0\). So the cost
of that exclusion has a precise arithmetic price: **the exchange rate sees the
largest fiber exactly and the smallest one not at all.**

**Signature collisions are structural.** Four hundred random elliptic
fibrations \(y^2=P_3(x)+c\) produce only \(5,3,3\) distinct signatures at
\(q=101,211,503\) --- exactly \(\gcd(4,q-1)+1\), the number of quartic twist
classes. From genus \(2\) upward, \(398\)--\(400\) of \(400\) are distinct.
Relatedly, \(m_2\) has the closed form \(m_2=K_P/q-1\) with
\(K_P=\#\{(x,x'):P(x)=P(x')\}\) --- the fiber-square count again.

## One exchange rate that detects a congruence

The cleanest example the framework produces so far.

**Observation.** For the family \(f(x,y)=y^2-x^3\), whose fibers are the sextic
twists \(y^2=x^3+c\),

\[
C(f\to L)=C(L\to f)=1
\qquad\Longleftrightarrow\qquad
q\equiv2\pmod 3 .
\]

Verified at \(q=101,401,2003\) (both rates exactly \(1\)) against
\(q=211,1009\) (rates \(0.9765\) and \(0.9915\) in one direction). The reason is
classical: for \(q\equiv2\pmod3\) every curve \(y^2=x^3+c\) is supersingular,
\(a_c=0\), so every fiber has exactly \(q\) points and the signature is the flat
\(\{q,\ldots,q\}\) --- the same resource as a linear map, hence asymptotically
reversible with product \(1\). For \(q\equiv1\pmod3\) the traces are nonzero and
reversibility fails.

So a single entry of the exchange matrix being exactly \(1\) is equivalent to a
congruence condition on \(q\), i.e. to the splitting behaviour of \(q\) in
\(\mathbb Z[\zeta_3]\).

The same phenomenon reverses a comparison. Between the two families,

\[
E_1\prec E_0 \ \ (q\equiv1\bmod3),
\qquad
E_0\prec E_1 \ \ (q\equiv2\bmod3),
\]

because at \(q\equiv2\) the CM family degenerates to the flat signature and
becomes the weaker resource, while at \(q\equiv1\) its heavier trace
distribution makes it the stronger one. The direction of one arrow in the
comparison graph is controlled by \(q\bmod3\).

## What this suggests

Three questions the paper can pose with evidence behind them.

1. **Is \(m_2\) --- hence the monodromy group --- recoverable from the exchange
   matrix?** The expansion says the first correction is a function of
   \((m_1,m_2)\) and \(m_1=0\). If the map \((m_2,m_3,\ldots)\mapsto\) rate is
   injective on some natural class, the exchange matrix determines the
   Katz--Sarnak symmetry type.
2. **Does vertical Sato--Tate give a limit law for the rate?** For a family with
   large monodromy the \(a_c/2\sqrt q\) equidistribute towards the semicircle, so
   \((1-C(L\to f))\sqrt q\log q\to2g\) should hold with an explicit error, and
   the fluctuation should be governed by the extreme-value statistics of the
   trace distribution.
3. **Is there a family whose exchange comparison graph has a cycle?** Over
   \(\mathbb F_q\) the quadratic classes are totally ordered; the elliptic
   families already show one arrow reversing with \(q\). Three families whose
   pairwise comparisons rotate would be the arithmetic counterpart of the
   integer three-cycle of the first paper.
