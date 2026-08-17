# Two positivities: the exchange matrix \(M\) and the Weil matrix \(E\)

This note proposes the backbone of the second paper. It answers one question:
what exactly is the relation between Weil's positivity criterion
(\(E\succeq0\iff\mathrm{RH}\)) and the directed exchange matrix
\(M_{ij}=C(f_i\to f_j)\) of the first paper?

The short answer is a dictionary with one honest gap, and it is sharper than
an analogy. Both matrices are *price matrices of the same commodity*
--- a resource evaluated against a family of multiplicative monotones ---
and they differ in exactly one respect:

* \(M\) evaluates the monotones on the **real** part of the spectrum and
  compares by an **infimum of ratios**. That is an \(L^\infty\), order-theoretic
  object: it produces a quasi-metric, a no-arbitrage law, and a comparison
  that runs in circles.
* \(E\) evaluates them on the **complex** spectrum and compares by a
  **Hermitian pairing**. That is an \(L^2\), Hilbertian object: it produces a
  Gram matrix, an angle metric, and positive semidefiniteness.

The Riemann hypothesis is precisely the statement that, for the arithmetic
resource, the complex spectrum that matters is again real. Everything below
makes this sentence precise, and two of the steps are new theorems about the
first paper's own objects.

---

## 1. The exchange matrix is a sup-norm distance

Fix the notation of the first paper: \(Z_a(\beta)=\sum_ia_i^\beta\),
\(\varphi_a=\log Z_a\), and
\(C(a\to b)=\inf_{\beta\in[0,\infty]}\varphi_a(\beta)/\varphi_b(\beta)\).
Since \(\varphi_a>0\) throughout, put

\[
u_a(\beta)=\log\varphi_a(\beta)=\log\log Z_a(\beta),
\qquad \beta\in[0,\infty].
\]

**Theorem 1 (isometric representation).**

\[
-\log C(a\to b)=\sup_{\beta}\bigl(u_b(\beta)-u_a(\beta)\bigr),
\qquad
d(a,b):=-\log\bigl(C(a\to b)\,C(b\to a)\bigr)
=\operatorname*{osc}_\beta\bigl(u_a-u_b\bigr),
\]

where \(\operatorname{osc}f=\max f-\min f\). Consequently
\(a\mapsto u_a\) is an isometry from signatures, modulo Cartesian powers, into
the quotient space \(C([0,\infty])/\mathbb R\) with the sup-norm, because
\(u_{a^{\otimes k}}=u_a+\log k\) is a *constant* shift.

The proof is one line each way: \(C(a\to b)=\exp\inf_\beta(u_a-u_b)\).
The identity is verified numerically in `analysis/exchange_positivity.py`.

Four facts that the first paper states separately, or uses without proof,
collapse into corollaries of Theorem 1.

1. \(C(a\to a)=1\).
2. \(C(a\to b)C(b\to a)\le1\), with equality iff \(u_a-u_b\) is constant, i.e.
   iff the Gibbs regions are homothetic. So \(d\) is a pseudometric whose zero
   set is exactly asymptotic reversibility.
3. **Supermultiplicativity** \(C(a\to c)\ge C(a\to b)\,C(b\to c)\): this is the
   triangle inequality \(\sup(u_c-u_a)\le\sup(u_c-u_b)+\sup(u_b-u_a)\).
   In particular \(L(a,b)=-\log C(a\to b)\) is a quasi-metric.
4. **No arbitrage.** Around any cycle,
   \(C(f_1\to f_2)C(f_2\to f_3)\cdots C(f_k\to f_1)\le C(f_1\to f_1)=1\).
   Trading in a circle never creates resources.

Item 4 deserves emphasis, because it is the fact a reader will look for after
the three-cycle of the first paper: the *comparison* \(\prec\) runs in circles,
but the *market* does not. Those are different statements, and Theorem 1
separates them cleanly. The cycle is the non-transitivity of comparing
\(\sup(u_b-u_a)\) with \(\sup(u_a-u_b)\) --- two one-sided suprema of the same
function --- which is generic behaviour for asymmetric sup-norms, not an
accident of three small signatures.

---

## 2. The exchange metric is not Hilbertian

Weil's matrix is a Gram matrix. Is ours? The precise question is whether the
pseudometric \(d\) of Theorem 1 is of **negative type**, that is whether

\[
\sum_{i,j}x_ix_j\,d(f_i,f_j)\le0
\qquad\text{whenever}\qquad \sum_ix_i=0 .
\]

By Schoenberg's theorem this is equivalent to \(\sqrt d\) embedding
isometrically into a Hilbert space, and to
\(\bigl(e^{-t\,d(f_i,f_j)}\bigr)_{ij}=\bigl((C(f_i\to f_j)C(f_j\to f_i))^{t}\bigr)_{ij}\)
being positive semidefinite for every \(t>0\). That last matrix is the exact
structural counterpart of Weil's \(E\): a symmetric matrix of pairwise prices,
with unit diagonal, whose positivity would say that the whole family of
resources is consistent with one inner product.

**Theorem 2 (failure of positivity, unconditional, and minimal).** The exchange
metric is not of negative type. The smallest possible witness has **five**
signatures, and one is

\[
\begin{aligned}
a_1&=\{12,10,8,8,2,1\}, &a_2&=\{11,9,7,7,4,1\}, &a_3&=\{12,12,6,5,4,4\},\\
a_4&=\{12,10,7,4,3,3\}, &a_5&=\{11,11,7,7,4,3\},
\end{aligned}
\]

with \(x=(302626,-510642,-576418,330027,454407)/10^6\), so \(\sum_ix_i=0\)
exactly, and

\[
\sum_{i,j}x_ix_j\,d(a_i,a_j)=+9.8126948851\cdot10^{-4}>0 .
\]

**Five is the floor, and it is proved rather than searched.** Every cut
semimetric satisfies \(x^{\mathsf T}\delta_Sx=-2(\sum_{i\in S}x_i)^2\le0\), so
\(\mathrm{CUT}_n\subseteq\mathrm{NEG}_n\); and \(\mathrm{MET}_4=\mathrm{CUT}_4\),
so no four-point metric can violate negative type at all. The family above
attains the bound.

The certificate is verified four ways: the package solver, a dense
\(2\cdot10^6\)-point \(\beta\)-grid on \([0,600]\) (agreeing to
\(1.5\cdot10^{-9}\)), the same on \([0,2000]\), and 40-digit `mpmath`
(\(3.8\cdot10^{-16}\)). The margin \(\sim10^{-3}\) is six orders above the
solver tolerance. Triangle slack over all triples is \(+4.87\cdot10^{-5}\), so
it is a genuine metric, and the twenty infima are attained at \(\beta=0\) (four
pairs), at interior \(\beta\) up to \(12.73\), and at \(\beta=\infty\) (three
pairs) --- the violation uses the whole temperature range.

Two further facts place the metric exactly in the Deza--Laurent hierarchy. The
**pentagonal inequality** \(b=(1,1,1,-1,-1)\) is violated, by
\(+5.3933\cdot10^{-2}\), on
\(\{10,6\},\{8,8,1,1,1,1\},\{10,10,6,5,4,4\}\mid\{9,9,2\},\{10,5,5,3,3,1\}\);
and since \(\mathrm{CUT}_5=\mathrm{HYP}_5\), \(\ell_1\)-embeddability also
breaks at exactly five points. The triangle inequality holds throughout
(minimum slack \(-4.4\cdot10^{-16}\) over all \(8\cdot10^9\) triples of a
2000-signature core).

Certificate and reproduction: `analysis/negative_type_certificate.py`, which
verifies the witness, the triangle inequality, the four-point minimality bound
and the pentagonal violation. The search that found the family is the
exploratory record in `research/m_and_e_and_a_c/`; the thirteen-signature family
reported earlier was a search artefact, superseded by this one.

The interpretation is the point of the note. In the finite-map world the
spectrum is real by construction, so the theory is an order theory and the
natural geometry is \(L^\infty\); \(L^\infty\) geometries are not Hilbertian,
and Theorem 2 exhibits the obstruction concretely. **The analogue of Weil
positivity is false for finite maps, and the reason it is false is exactly the
reason exchange rates exist at all**: an infimum of ratios is a Chebyshev
construction, not a quadratic one.

---

## 3. Where the complex spectrum comes from

The first paper's Remark on monotones proves: the \(\otimes\)-multiplicative,
\(\oplus\)-additive, \(\preccurlyeq\)-monotone maps on signatures are exactly
\(a\mapsto Z_a(\beta)\) for \(\beta\ge0\). Read the proof again and note which
hypothesis does which work.

* Multiplicativity and additivity force \(\varphi(\{n\})=t_n\) to be
  completely multiplicative.
* **Monotonicity** --- and only monotonicity --- forces \(t_n=n^\beta\) with
  \(\beta\) *real and non-negative*.

Drop the order and keep the algebra: the completely multiplicative characters
of \(\mathbb N\) are \(n\mapsto n^{s}\) for \(s\in\mathbb C\), together with
their twists by Dirichlet characters. So:

> The asymptotic spectrum of the signature semiring is the *real* half-line
> \([0,\infty]\); its complexification is the *complex* plane of Dirichlet
> characters. The order picks out the real points.

That single sentence is the bridge. \(M\) is what you get from the real
points; \(E\) is what you get from the complex ones.

### The zeta resource

The bridge is not only formal. Let \(P_{p,k}=\{1,p,p^2,\ldots,p^k\}\), the
signature whose entries are the powers of a prime. In the cost convention
\(s=-\beta\),

\[
Z_{P_{p,k}}(s)=\sum_{j=0}^{k}p^{-js}
\xrightarrow[k\to\infty]{}\frac1{1-p^{-s}}=\zeta_p(s),
\]

and because \(\otimes\) multiplies entries, unique factorization reads

\[
\bigotimes_{p}P_{p,\infty}=\{1,2,3,4,\ldots\},
\qquad
Z(s)=\prod_p\zeta_p(s)=\zeta(s).
\]

**The Euler product is a \(\otimes\)-factorization of a resource in the
signature semiring, and \(\zeta\) is its partition function.** The completed
\(\xi\) adds the archimedean factor, which in the local language of
[Local zeta profiles and Riemann-zeta currency](p_adic_zeta_exchange_currency.md)
is the real place's own processor.

Two honest caveats, and both are informative:

* The product converges only for \(\Re s>1\). The critical strip is exactly the
  region where the monotone *diverges*, so RH is a statement about the analytic
  continuation of a monotone, not about the monotone. Any bridge from the
  operational theory to RH has to cross this gap, and no argument here crosses
  it.
* Each finite \(\otimes\)-subproduct is zero-free in the half-plane, so the
  analogue of RH is vacuous for every finite portfolio of local currencies.
  **The zeros of \(\zeta\) are a phenomenon of the infinite tensor product.**
  In resource language: every finite portfolio of local currencies is
  arbitrage-free for trivial reasons; RH asserts that the completed global
  currency remains so in the limit. This locates the difficulty precisely
  where number theory says it is, and it says why Theorem 2 is not a statement
  about RH.

### Exchange rates against the completed zeta function

The completion repairs \(\zeta\) as a partition function --- with
\(\Phi>0\) even and doubly exponentially decaying,
\(\xi(\tfrac12+\beta)=\int\Phi(u)e^{\beta u}du\) is the Laplace transform of a
positive measure, so \(\log Z\) is finite and convex for every real \(\beta\)
and the Gibbs curve is symmetric about the entropy axis. It does **not** repair
comparability. Normalising the measure to mass one,

\[
\log Z_\xi(\beta)\sim\tfrac12(\log Z)''(0)\,\beta^2=0.0231\,\beta^2
\ \ (\beta\to0),
\qquad
\log Z_\xi(\beta)\sim\tfrac\beta2\log\beta\ \ (\beta\to\infty),
\]

against \(\log Z_{P}(0)=\log(K+1)>0\) and \(\log Z_P\sim K\beta\log p\) for a
truncated Euler factor \(P_{p,K}=\{1,p,\ldots,p^K\}\). The ratio therefore tends
to \(0\) at one end and to \(\infty\) at the other, and **both unrestricted rates
vanish**: \(C(\xi\to P)=C(P\to\xi)=0\). This is the first limitation of the
companion paper --- an unbounded spectrum with a different growth exponent gives
nothing --- so the comparison must be restricted to a temperature window, exactly
as for the Standard Model plasma there.

On a window the rates are finite and attained at opposite ends. What is
striking is the limit. As \(K\to\infty\) the profile \(\log Z_P\to K\beta\log p\)
becomes a straight line through the origin, and both \(K\) and \(p\) cancel from
the *product* of the two rates:

\[
C(\xi\to P)\,C(P\to\xi)\ \longrightarrow\
\frac{\min_{\beta\in W}\log Z_\xi(\beta)/\beta}
     {\max_{\beta\in W}\log Z_\xi(\beta)/\beta},
\qquad
d_W(\xi)=\log\frac{\max_W\log Z_\xi/\beta}{\min_W\log Z_\xi/\beta}.
\]

The quantity being compared is the slope of the chord from the origin, so
\(d_W\) is exactly the failure of \(\log Z_\xi\) to be linear through the origin.
Computed values (`analysis/xi_versus_euler_factors.py`):

| window \(W\) | limiting product | \(d_W(\xi)\) |
|---|---:|---:|
| \([0.5,5]\) | \(0.101908\) | \(2.2837\) |
| \([1,10]\) | \(0.106848\) | \(2.2364\) |
| \([2,20]\) | \(0.120723\) | \(2.1143\) |

and the approach is monotone: on \([1,10]\) against \(p=2\) the irreversibility
falls \(2.6968\to2.3538\to2.2441\) at \(K=1,8,128\), and larger primes converge
faster. **So longer ladders and larger primes are better currency, but there is
a floor no Euler factor can cross**, and the floor is an invariant of \(\xi\)
alone --- the curvature of its log-partition profile on the window. In the
language of the local note below, this is the precise sense in which a single
scalar price cannot quote \(\xi\): the yield curve is not a straight line.

### And against the local Igusa profiles

The same question for the \(p\)-adic side has a different answer, and the
difference is instructive. An Igusa profile
\(\Phi_{f,p}(s)=-\log\int_{\mathbb Z_p^n}|f(x)|_p^s\,dx\) is positive and
increasing with \(\Phi(0)=0\), and both of its ends are classical invariants of
the map:

\[
\Phi_{f,p}'(0)=\log p\cdot\mathbb E\bigl[v_p(f)\bigr],
\qquad
\Phi_{f,p}(\infty)=-\log\mu\{|f|_p=1\}.
\]

Both are **finite**, and the second is the obstruction: an Igusa profile
*saturates*. Computed values confirm the closed forms exactly --- for \(f=x\) on
\(\mathbb Z_p\) the slope is \(\log p/(p-1)\) and the ceiling is
\(-\log(1-1/p)\); for the monomial \(x^2y\) on \(\mathbb Z_2\) the slope is
\(3\log2\) and the ceiling \(2\log2\).

So the mismatch with \(\xi\) is a *new* pair of failures: quadratic against
linear at \(s=0\), and divergent against **bounded** at \(s=\infty\). Both
unrestricted rates vanish again. On a window, and with \(r\) independent copies
cancelling from the product exactly as \(K\) did:

| window | best Euler factor \(d\) | best Igusa profile \(d\) | which |
|---|---:|---:|---|
| \([0.5,5]\) | \(2.2837\) | \(3.6166\) | \(x\) on \(\mathbb Z_2\) |
| \([1,10]\) | \(2.2364\) | \(4.0034\) | \(x\) on \(\mathbb Z_2\) |
| \([2,20]\) | \(2.1143\) | \(4.2029\) | \(x\) on \(\mathbb Z_2\) |

**A long \(p\)-adic ladder is markedly better currency for \(\xi\) than any
single \(p\)-adic map**, by a factor of about \(e^{1.77}\approx5.9\) in the
product on \([1,10]\), and the reason is structural rather than numerical: the
truncated Euler factors at least grow linearly and so can track \(\xi\)'s
divergence up to a slope, whereas one Igusa profile has a ceiling. Within the
Igusa family the ranking is by that ceiling --- smaller \(p\) and lower degree
are better, \(\mathbb Z_2\) best --- and \(xy\) over \(\mathbb Z_3\) scores
exactly the same as \(x\) over \(\mathbb Z_3\), because
\(\Phi_{xy}=2\Phi_x\) is two copies and copies cancel.

Reproduced by `analysis/xi_versus_igusa_profiles.py`, which first recovers the
note's constant \(C_\zeta(x^2\to x^2-y^2)=0.939702787545916\) to fifteen digits
as a check on the profile formulas. The moral for the currency programme is
sharp: a single place cannot quote \(\xi\), because a single place has finite
bandwidth. That is the resource-theoretic reason the global object needs *all*
of them.

### The dictionary

| | exchange theory (papers 1--2) | Weil / RH |
|---|---|---|
| resources | finite maps, polynomial maps | admissible test profiles, adelic dilation portfolios |
| parallel composition | \(\otimes\) (Cartesian product) | product of local factors |
| alternative | \(\oplus\) (disjoint union) | addition of profiles |
| monotones | \(Z_a(\beta)=\sum a_i^\beta\) | \(\widehat g(s)\), Mellin transform |
| spectrum | \(\beta\in[0,\infty]\), real | \(s\in\mathbb C\), critical strip |
| local factors | fiber counts, Igusa \(Z_{f,p}\) | Euler factors \(\zeta_p\) |
| global object | \(\bigotimes_pP_{p,\infty}\) | \(\xi(s)\) |
| pairing | \(M_{ij}=C(f_i\to f_j)\), directed | \(E_{ij}=\sum_\rho\widehat g_i(\rho)\overline{\widehat g_j(1-\overline\rho)}\), Hermitian |
| induced geometry | quasi-metric \(-\log M\), sup-norm | angle metric \(\arccos R\), spherical |
| positivity law | cycle products \(\le1\) (Theorem 1.4, unconditional) | \(E\succeq0\) (\(\iff\) RH) |
| Hilbert embedding | **fails** (Theorem 2) | holds \(\iff\) RH |
| what it asserts | the spectrum is real by construction | the spectrum is real as a theorem |

The last two rows are the payload. Under RH, \(1-\overline\rho=\rho\) and
\(E=\sum_\rho v_\rho v_\rho^{*}\) becomes a Gram matrix --- the same shape as
\(G_{ij}=\int\varphi_{f_i}\varphi_{f_j}\,d\mu\), which in our setting is
positive semidefinite *unconditionally*, for any positive measure \(\mu\) on
the spectrum, simply because it is a Gram matrix in \(L^2(\mu)\). And the two
pairings are the two ends of one family: by Laplace's principle,

\[
C(a\to b)^{-1}
=\lim_{t\to\infty}
\left(\int_0^\infty
\Bigl(\frac{\varphi_b(\beta)}{\varphi_a(\beta)}\Bigr)^{t}d\mu(\beta)\right)^{1/t},
\]

for any \(\mu\) of full support. **The exchange rate is the tropical
(\(t\to\infty\)) limit of a family of \(L^t\) pairings whose \(t=2\) member is a
Gram matrix.** Weil positivity lives at \(t=2\); the exchange matrix lives at
\(t=\infty\); the passage between them is Maslov dequantization. This is the
precise sense in which \(M\) and \(E\) are the same object.

---

## 3b. Four structural connections between \(E\) and \(M\)

The dictionary above is an analogy of shape. These four are tighter: each is a
classical theorem in which both matrices appear as objects of the same kind.

### (i) They are two levels of one classical hierarchy of matrix cones

The cones of matrices on \(N\) points nest:

\[
\text{cut cone}\ \subset\
\text{negative-type cone}\ \subset\
\text{metric cone},
\]

with the middle one equal, by Schoenberg, to the \(\ell_2^2\)-embeddable
metrics, i.e. exactly those \(d\) for which \((e^{-td})\succeq0\) for all
\(t>0\). This is the standard hierarchy of Deza and Laurent, *Geometry of Cuts
and Metrics*.

* **Weil's \(E\succeq0\)** places the arithmetic price data in the
  negative-type cone: RH says the family of test profiles is
  \(\ell_2\)-realisable.
* **No arbitrage** places \(-\log M\) in the metric cone, and nothing more.
* **Theorem 2 above** says the containment is strict for our data: \(M\) is in
  the metric cone and not in the negative-type cone.

So the two positivity statements are not analogous, they are *comparable*:
they are the two adjacent levels of one hierarchy, and the exchange theory
occupies the weaker one. This also says exactly what a proof of a "Weil
positivity for maps" would have to supply --- an \(\ell_2\) structure that
Theorem 2 rules out for the plain rate, hence a different pairing.

### (ii) The transpose is the functional equation

Weil's entries pair \(\widehat g_i\) at \(\rho\) against \(\widehat g_j\) at the
*reflected* point \(1-\overline\rho\). Our symmetrised entries pair
\(C(f_i\to f_j)\) against the *reversed* rate \(C(f_j\to f_i)\). In both cases
the matrix is built from a quantity and its image under an involution, and in
both cases positivity is the statement that the involution acts trivially:

| | involution | acts trivially when | then the matrix is |
|---|---|---|---|
| \(E\) | \(s\mapsto1-\overline s\) | RH: \(\rho=1-\overline\rho\) | a Gram matrix \(\sum_\rho v_\rho v_\rho^{*}\) |
| \(G=M\circ M^{\mathsf T}\) | transpose | reversibility: \(C(a\to b)C(b\to a)=1\) | the all-ones matrix, rank one |

Hence the exact correspondent of "a zero off the critical line" is
**irreversibility**:

\[
\delta=\bigl|\Re\rho-\tfrac12\bigr|
\qquad\longleftrightarrow\qquad
d(a,b)=-\log\bigl(C(a\to b)C(b\to a)\bigr).
\]

Both are non-negative, both vanish exactly in the self-dual case, and both are
what the corresponding positivity would forbid. Asked what in the exchange
theory plays the role of the Riemann hypothesis, the answer is: *asymptotic
reversibility of every pair* --- which is false, demonstrably and by a wide
margin, and whose failure is the entire content of the first paper.

**Status of this entry — now settled, and negatively.** It is a dictionary, not
a derivation, and the reason is no longer a matter of taste: a functional
relation between the two geometries is *impossible*. Each is invariant under a
group the other is not.

* **The Weil angle is invariant under \(a\mapsto\lambda a\).** Since
  \(Z_{\lambda a}(\rho)=\lambda^\rho Z_a(\rho)\), one has
  \(E_{\lambda a,\lambda b}=\lambda^{2\Re\rho}E_{a,b}\), and the factor is
  *constant across zeros only because every \(\Re\rho=\tfrac12\)*. Verified to
  \(6\cdot10^{-14}\); moving \(60\) of \(1200\) zeros to \(\Re=0.7\) breaks it
  to \(1.3\cdot10^{-1}\). **The scale-invariance of the Weil geometry is a
  manifestation of the critical line.**
* **The exchange metric is invariant under \(a\mapsto a^{\otimes k}\)**, since
  \(u_{a^{\otimes k}}=u_a+\log k\) is a constant shift and \(d\) is an
  oscillation. Verified to \(2.5\cdot10^{-16}\).

Under a common rescaling \(d\) moves while the angles do not; under Cartesian
powers the angles move while \(d\) does not. Hence neither is a function of the
other, and any correlation between them is an artefact of the family chosen.
Research record: `research/m_and_e_and_a_c/T1_5_multiplicative_design.md`. There is no map carrying \(E\) to \(G\) or back: they are attached
to different semirings, and no property of \(E\) is expressible through \(G\).
Two further cautions belong here.

*\(G\) is not the structural analogue of \(E\).* The analogue of a Weil matrix
in our setting is the \(L^2\) Gram matrix
\(\Gamma_{ij}=\int\varphi_{f_i}\varphi_{f_j}\,d\mu\), which is positive
semidefinite for every positive \(\mu\), unconditionally --- and that
unconditional positivity *is* the content of the dictionary, since it holds for
exactly the reason RH would: the spectrum is real. \(G=M\circ M^{\mathsf T}\) is
a different, tropical object, sitting at the \(t\to\infty\) end of the
interpolation rather than at \(t=2\).

*The two positivities do not propagate into one another.* If \(d\) were of
negative type, positivity at one \(t\) would give positivity at every \(t\) by
Schoenberg. Theorem 2 says it is not, so \(\Gamma\succeq0\) and the failure of
\(e^{-td}\succeq0\) coexist without contradiction. The correct summary is that
positivity holds at \(t=2\) and fails in the tropical limit, and that RH is a
statement about the \(t=2\) end alone.

### (iii) Both are cycle conditions on a price matrix, from opposite sides

For a positive semidefinite matrix with unit diagonal, the \(3\times3\)
determinant condition reads

\[
1+2R_{12}R_{23}R_{31}-R_{12}^2-R_{23}^2-R_{31}^2\ \ge\ 0,
\]

a *lower* bound on the cycle product \(R_{12}R_{23}R_{31}\). No arbitrage is an
*upper* bound on the cycle product of \(M\). So both positivities constrain the
same functional --- the product of prices around a loop --- and they constrain
it from opposite sides. In the \(E\) world a Hilbert structure forbids cycle
products that are too negative; in the \(M\) world the resource order forbids
cycle products that exceed one.

### (iv) A quantitative bridge: distortion \(O(\log N)\)

By Bourgain's embedding theorem every \(N\)-point metric embeds into Hilbert
space with distortion \(O(\log N)\). Applied to \(d\):

> for every finite family of \(N\) resources there is a Gram matrix --- a
> \(E\)-shaped object --- reproducing the exchange metric to within a factor
> \(O(\log N)\).

So \(M\) is not Hilbertian, but it is never far from Hilbertian, and the
obstruction of Theorem 2 is at worst logarithmic. Bourgain's bound is attained
on expanders, which suggests where to look for worse families: the certificate
of Theorem 2 should be compared against the expansion of the comparison graph.
That is a concrete experiment the project can run with the existing cache.

### A fifth, softer one

Strassen's asymptotic spectrum and Bochner's theorem are the same kind of
statement: a comparison is decided by a family of characters --- via an
infimum in the first case, via an integral in the second. The \(L^t\)
interpolation of Section 3 is the passage between them, and \(t=2\) versus
\(t=\infty\) is precisely Bochner versus Strassen.

---

## 4. Weil's theorem, in the register where it is proved

Over finite fields Weil's Riemann hypothesis is a theorem, and there the
relation to \(M\) is not an analogy but a computation: **the entries of the
exchange matrix over \(\mathbb F_q\) are functions of Weil numbers**, because
fiber sizes are point counts.

For binary quadratic forms over \(\mathbb F_q\) (odd \(q\)) the signatures are

\[
\sigma(xy)=\{2q-1,\underbrace{q-1,\ldots,q-1}_{q-1}\},
\qquad
\sigma(x^2-\nu y^2)=\{\underbrace{q+1,\ldots,q+1}_{q-1},1\},
\]

the \(\pm1\) being exactly the trace of Frobenius on the conic. Then:

**Theorem 3 (conic exchange rates).** For every odd prime power \(q\),

\[
C\bigl(\text{anisotropic}\to\text{split}\bigr)=\frac{\log(q+1)}{\log(2q-1)}
\qquad(\text{attained at }\beta=\infty),
\]

and

\[
C\bigl(\text{split}\to\text{anisotropic}\bigr)
=1-\frac{\kappa}{q\log q}+O(q^{-2}),
\qquad
\kappa=\max_{\beta>0}\frac{2\beta-2^{\beta}}{\beta+1}=0.068755890904\ldots,
\]

the maximum being attained at \(\beta_*=1.478296901967\ldots\), which is also
the limit of the exact bottleneck temperatures.

Both statements are checked against the solver for
\(q\le509\) in `analysis/exchange_positivity.py`; the endpoint formula agrees to
ten digits at every \(q\), and \((1-C)q\log q\) rises to \(0.0680\) at \(q=509\)
against the predicted \(0.06876\).

The structure of the expansion is the general phenomenon and should be the
theorem the paper actually proves:

**Proposed Theorem 3\('\) (Weil-controlled exchange rates).** Let
\(f,g:\mathbb A^n\to\mathbb A^1\) be defined over \(\mathbb Z\) with
geometrically irreducible generic fibers. Then, over \(\mathbb F_q\),

\[
M_q(f\to g)=1+\frac{1}{q\log q}\,W(f,g)+O\!\left(\frac{1}{q^{3/2}\log q}\right),
\]

where \(W(f,g)=\max_\beta\bigl(w_g(\beta)-w_f(\beta)\bigr)/(\beta+1)\) and
\(w_f\) is assembled from the Frobenius traces of the fibers of \(f\).
Lang--Weil supplies the error term; Weil's theorem is what makes the
expansion uniform in \(q\).

Read backwards, this says: *to leading order the exchange rate over a large
finite field sees only the dimension of the fibers; the entire arithmetic
content sits in the first correction, and that correction is a Weil number.*
The exchange matrix is a lens that magnifies exactly the term Weil's theorem
bounds.

The parallel closes at the level of proofs, too. Weil's own proof of RH for
curves is a *positivity* statement --- the Castelnuovo--Severi inequality for
the intersection form on correspondences. So on both sides of the dictionary,
"the spectrum is real" is proved, when it is proved, by exhibiting a positive
definite form. That is the reason to expect the \(E\)-side and not the
\(M\)-side to be the arithmetically deep one, and Theorem 2 confirms it from
the other direction: the \(M\)-side geometry cannot support such a form.

---

## 5. Proposed main results of the second paper

The present material is a set of examples showing that one framework
reproduces classical facts. Four theorems would turn it into a paper.

**A. The degeneration order is local--global.** For quadratic maps
\(K^2\to K\), the posets computed here over \(\mathbb F_q\), \(\mathbb R\),
\(\mathbb C\) and \(\mathbb Q_p\) are the local pieces of one global order:
\([R]\preceq[Q]\) over \(\mathbb Q\) iff it holds over every completion. This
follows from Hasse--Minkowski applied to the isotropy of \(Q\perp\langle-b\rangle\),
and it is the statement that makes the collection of diagrams a single object
rather than a gallery. It also predicts the node counts already in the notes:
\(|K^\times/(K^\times)^2|\) anisotropic classes, hence 8 over \(\mathbb Q_2\),
4 over \(\mathbb Q_p\) for odd \(p\), 2 over \(\mathbb R\) and \(\mathbb F_q\),
1 over \(\mathbb C\).

**B. The local monotones are Igusa zeta functions.** For \(p\)-adic polynomial
maps, the \(\otimes\)-multiplicative monotones are the Igusa local zeta
functions \(s\mapsto Z_{f,p}(s)\), so the diagnostic
\(C_\zeta=\inf_s\Phi_f/\Phi_g\) is the exchange rate of the framework, not an
ad-hoc quantity. The bad-prime correction \(t^2-2t+2\) and the constant
\(0.9397027875\ldots\) become the worked example. What is missing today is the
coding theorem identifying \(C_\zeta\) with an achievable affine rate; state
that gap as a conjecture with the evidence from
[the \(\mathbb Q_2\) attempts](p_adic_exchange_rate_attempts.md).

**C. The Euler product is a tensor factorization** (Section 3 above), with the
two caveats stated honestly. Short, quotable, and it is what licences calling
\(\zeta\) a partition function in this framework.

**D. The two positivities** (Sections 1, 2 and 4 above): Theorem 1, Theorem 2
with its certificate, the \(L^t\) interpolation, and Theorem 3/3\('\) over
finite fields. This is the paper's thesis and should be its title material.

A candidate title: *Exchange rates, degeneration posets and two kinds of
positivity*. A candidate one-sentence abstract: *the same resource calculus
that computes exchange rates between finite maps produces, at its real
spectrum, an order theory with a no-arbitrage law and a non-Hilbertian
geometry, and at its complex spectrum the Weil pairing whose positivity is the
Riemann hypothesis; over finite fields the two meet, and the exchange matrix is
a Weil number to first order.*

### What is not claimed

No implication toward RH is asserted in either direction. Theorem 2 says the
finite-map exchange geometry is not Hilbertian; it says nothing about
\(\zeta\), because the zeta resource is an infinite tensor product outside the
convergence region of the monotone. The open bridge stated in
[Riemann hypothesis as positivity of exchange matrices](riemann_hypothesis_exchange_matrices.md)
--- a density theorem placing polynomial-map profiles inside the admissible
Weil test class --- remains the missing step, and Sections 3 and 4 above say
exactly which two obstacles it must clear.
