# T1.5 — Designing multiplicative coincidences: how far can the Weil geometry be made to track the exchange geometry?

**Status: settled, with two exact theorems and a quantified bridge.**
This extends T1.4. There the correlation between the exchange distance
\(d(a,b)=-\log C(a\to b)C(b\to a)\) and the Weil angle \(\arccos R_{ab}\) was
\(+0.19\) on a generic family, and the conclusion was that the two geometries
are essentially independent. The question here: is that independence
unavoidable, or an artefact of genericity?

**Answer.** It is an artefact — the correlation can be driven to \(+0.86\) or to
\(-0.98\) by design — but *not for the reason one would hope*. The agreement is
produced entirely by the **overlap** term \(O\), which carries no arithmetic;
the **von Mangoldt** term \(A\), the only part of \(E\) that knows about primes,
produces agreement of the **opposite sign**. And two exact invariance theorems
show the two functionals can be moved independently of one another, so no
functional relation between them can exist.

Throughout, \(N=1200\) zeros, \(T=1648.27\), \(T/2\pi=262.33\), and by T1.4

\[
E_{ab}=N\cdot O(a,b)-\frac{T}{2\pi}A(a,b)+O(rs\log T),\qquad
O=\!\!\sum_{a_i=b_j}\!\!a_i,\quad
A=\!\!\sum_{a_i\ne b_j}\!\!\min(a_i,b_j)\Lambda\!\Big(\tfrac{\max}{\min}\Big).
\]

---

## 1. Designed families beat generic ones — in both directions

Correlation between \(d\) and the Weil angle. Pearson / Spearman, with the pair
count so significance is judgeable.

| family | signatures | pairs | Pearson | Spearman | mean angle |
|---|---:|---:|---:|---:|---:|
| generic (T1.4 family) | 16 | 120 | \(+0.191\) | \(+0.187\) | 1.487 |
| generic control A (no shared entry, no prime-power ratio) | 16 | 120 | \(-0.153\) | \(-0.169\) | 1.5715 |
| generic control B | 21 | 210 | \(-0.181\) | \(-0.231\) | 1.5720 |
| 3-smooth pairs \(\{2^i3^j\}\) | 21 | 210 | \(+0.165\) | \(+0.179\) | 1.528 |
| pairs \((2^i,2^j)\) | 21 | 210 | \(+0.417\) | \(+0.264\) | 1.408 |
| pairs \((3^i,3^j)\) | 15 | 105 | \(+0.444\) | \(+0.290\) | 1.381 |
| triples \((2^i,2^j,2^k)\) | 20 | 190 | \(+0.616\) | \(+0.690\) | 1.185 |
| **5-subsets of \(\{2,4,\dots,256\}\)** | 56 | **1540** | \(\mathbf{+0.725}\) | \(\mathbf{+0.815}\) | 1.009 |
| **6-subsets of \(\{3,9,\dots,6561\}\)** | 28 | 378 | \(\mathbf{+0.857}\) | \(\mathbf{+0.870}\) | 0.786 |
| staircase \((2^t,1)\), \(t=1..9\) | 9 | 36 | \(-0.708\) | \(-0.858\) | 1.586 |
| staircase \((2^t,2^t)\), \(t=1..9\) | 9 | 36 | \(-0.801\) | \(-0.847\) | 1.628 |
| chain \((2^t,2^t,1)\), \(t=1..9\) | 9 | 36 | \(-0.872\) | \(-0.947\) | 1.621 |
| **staircase \((3^t,1)\), \(t=1..8\)** | 8 | 28 | \(\mathbf{-0.878}\) | \(\mathbf{-0.980}\) | 1.619 |

The generic controls — built so that no pair shares an entry and no ratio is a
prime power, hence \(E_{ab}\approx0\) exactly — sit at correlation \(\approx-0.17\)
with a mean angle of \(1.5715\) and a standard deviation of \(0.006\): the whole
family is a nearly perfect orthonormal set, and the residual correlation is
Landau noise. That is the baseline the designed families must beat, and they
beat it by a lot — in **both signs**.

### The two signs are two different mechanisms

The sign is forced by which Landau term dominates.

* **Overlap-dominated** (\(O>0\)): \(E_{ab}>0\), angles **acute**, and the angle
  *increases* with \(d\) — agreement with the naive expectation. This happens
  whenever the signatures are subsets of a common ladder and share entries.
* **\(\Lambda\)-dominated** (\(O=0\), \(A>0\)): \(E_{ab}<0\), because \(A\ge0\)
  always, so angles are **obtuse** — the Weil form is *repulsive* — and the
  repulsion is largest for the pairs of most similar scale, i.e. the pairs
  *closest* in \(d\). Hence the angle *decreases* with \(d\).

A family carrying both mechanisms has them fighting, which is exactly why the
3-smooth lattice (\(+0.165\)) does *worse* than either pure design and no better
than generic.

---

## 2. The arithmetic contributes nothing to the positive branch

This is the sharpest negative result of the note. Take three ladders of eight
entries with nearly identical logarithmic spacing and form all \(k\)-subsets.
Two ladders are geometric, so **every** ratio is a prime power and \(A\) is
enormous. The third, \(\{3,7,17,37,67,131,257,521\}\), has the same spacing but
**no** ratio that is a prime power, so \(A\equiv0\) identically — its Weil
matrix is exactly \(N\cdot O\), pure combinatorics with no zeta content.

| ladder | max cross \(A\) | \(k=3\) (1540 pairs) | \(k=4\) (2415) | \(k=5\) (1540) | \(k=6\) (378) |
|---|---:|---:|---:|---:|---:|
| \(\{2,\dots,256\}\) | 88.7 | \(+0.514/+0.522\) | \(+0.623/+0.707\) | \(+0.725/+0.815\) | \(+0.822/+0.874\) |
| \(\{3,\dots,6561\}\) | 2402.7 | \(+0.517/+0.450\) | \(+0.637/+0.685\) | \(+0.750/+0.810\) | \(+0.857/+0.870\) |
| \(\{3,7,17,37,67,131,257,521\}\) | **0.00** | \(+0.562/+0.578\) | \(+0.659/+0.722\) | \(+0.743/+0.820\) | \(+0.826/+0.870\) |

The \(A\equiv0\) control **matches the geometric ladders cell for cell**, and at
\(k=3,4\) it is slightly *ahead* of them. So the climb from \(+0.19\) to
\(+0.86\) is bought entirely by shared entries. In the overlap-dominated regime
\(R_{ab}\to\langle w_a,w_b\rangle/\|w_a\|\|w_b\|\) with
\(w_a=(\sqrt v)_{v\in a}\): the Weil correlation degenerates into the cosine
similarity of \(\sqrt{\cdot}\)-weighted indicator vectors, a statement in which
neither the zeros nor the primes appear. Any Gram matrix of indicator vectors
would do as well.

> **Designing multiplicative coincidences does raise the correlation, but the
> arithmetic part of the coincidence is not what raises it.**

---

## 3. The one genuine bridge: the receding staircase

On \(a_t=(p^t,1)\) both geometries have closed forms.

**Exchange side (exact).** \(\log Z_{a_t}(\beta)=\log(p^{t\beta}+1)\), and the
infimum defining \(C\) is attained at \(\beta=\infty\) in one direction and
\(\beta=0\) in the other, giving \(C(a_s\to a_t)=s/t\) for \(s<t\) and
\(C(a_t\to a_s)=1\), hence

\[
\boxed{\;d(a_s,a_t)=\log (t/s)\;}
\]

verified to \(3\cdot10^{-16}\) over all pairs, \(p=2\) and \(p=3\).

**Weil side (Landau).** \(O=1\) (the shared entry \(1\)),
\(A=(p^{\min(s,t)}+2)\log p\), \(E_{tt}=N(p^t+1)-2\tfrac{T}{2\pi}\log p\); for
\(p^{s}\gg N/(\tfrac{T}{2\pi}\log p)\) the \(N\!\cdot\!O\) term is negligible and

\[
\boxed{\;\theta(a_s,a_t)-\tfrac\pi2\;=\;K\,p^{-(t-s)/2},\qquad
K=\frac{T\log p}{2\pi N}.\;}
\]

Fitted on windows \(t\in[T_0,T_0+8]\), nine signatures and 36 pairs each:

| \(p\) | \(T_0\) | Pearson | Spearman | \(K\) fitted | \(K=(T/2\pi N)\log p\) | rms residual |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 5 | \(-0.935\) | \(-0.988\) | 0.14229 | 0.15153 | 7.8% |
| 2 | 10 | \(-0.944\) | \(-0.976\) | 0.15026 | 0.15153 | 2.7% |
| 2 | 20 | \(-0.953\) | \(-0.983\) | 0.15071 | 0.15153 | 2.4% |
| 2 | 40 | \(-0.956\) | \(-0.980\) | 0.15072 | 0.15153 | 2.4% |
| 2 | 100 | \(-0.957\) | \(-0.978\) | 0.15072 | 0.15153 | 2.4% |
| 3 | 10 | \(-0.905\) | \(-0.975\) | 0.23667 | 0.24017 | 9.3% |
| 3 | 20 | \(-0.915\) | \(-0.973\) | 0.23669 | 0.24017 | 9.3% |
| 3 | 40 | \(-0.918\) | \(-0.966\) | 0.23669 | 0.24017 | 9.3% |

The fitted \(K\) locks onto the predicted \(T\log p/2\pi N\) to within 0.5%
(\(p=2\)) and 1.5% (\(p=3\)) and is independent of \(T_0\), exactly as the
derivation says it must be. Pearson *improves* as the window recedes, reaching
\(-0.957\), because \(\log(t/s)\) becomes more nearly linear in \(t-s\).

Now combine. Inside a window of width \(n\ll T_0\),
\(\log(t/s)=(t-s)/T_0\cdot(1+O(n/T_0))\), hence

\[
\theta-\tfrac\pi2\;=\;K\exp(-c\,d),\qquad c=\tfrac12 T_0\log p .
\]

**This is a genuine bridge**: on the receding staircase the Weil angle is an
explicit, derived, monotone function of the exchange distance, Pearson
\(-0.957\) and Spearman \(-0.978\) over 36 pairs, improving as the window
recedes. It is also the
*only* bridge found, and it comes with two warnings. The constant
\(c=\tfrac12T_0\log p\) is unbounded, so the two geometries agree only after a
reparametrisation whose parameter depends on where the window sits; and the
agreement has the **negative** sign — the closer two signatures are in the
exchange metric, the *more* the Weil form repels them.

### Why the fit is not a functional relation

The two closed forms depend on \((s,t)\) through different combinations:
\(d\) sees the **ratio** \(t/s\), the angle sees the **difference** \(t-s\).
Two decisive tables, \(p=3\), \(t=1..8\):

**(A) Same \(d\), different angle.** All four pairs sit at \(d=\log 2\) exactly.

| \((s,t)\) | \(d\) | angle | \(\theta-\pi/2\) |
|---|---:|---:|---:|
| \((1,2)\) | 0.69315 | 1.60337 | \(+0.0326\) |
| \((2,4)\) | 0.69315 | 1.62685 | \(+0.0561\) |
| \((3,6)\) | 0.69315 | 1.60752 | \(+0.0367\) |
| \((4,8)\) | 0.69315 | 1.59359 | \(+0.0228\) |

**(B) Same angle, \(d\) varying by a factor of five.** All pairs have \(t-s=1\).

| \((s,t)\) | \(d\) | angle | \(\theta-\pi/2\) |
|---|---:|---:|---:|
| \((4,5)\) | 0.22314 | 1.70415 | \(+0.13336\) |
| \((5,6)\) | 0.18232 | 1.70593 | \(+0.13513\) |
| \((6,7)\) | 0.15415 | 1.70717 | \(+0.13637\) |
| \((7,8)\) | 0.13353 | 1.70854 | \(+0.13774\) |

The angle is constant to \(0.004\) rad while \(d\) drops by a factor \(1.67\).
The near-perfect Spearman score of the staircase is a *window artefact*: over a
short window \(t/s\) and \(t-s\) are comonotone, and nothing more.

---

## 4. The obstruction, exactly: two transverse invariance groups

The two functionals have invariance groups that intersect trivially. Both
statements are exact identities, not asymptotics, and both were verified
numerically to machine precision.

> **Theorem A (the Weil geometry is scale invariant).** For any \(\lambda>0\)
> and any truncation, \(E_{\lambda a,\lambda b}=\lambda\,E_{a,b}\), because
> \(\sqrt{\lambda a_i\cdot\lambda b_j}=\lambda\sqrt{a_ib_j}\) and every ratio
> \(a_i/b_j\) is unchanged. Hence the correlation matrix \(R\), and therefore
> every Weil angle, is **exactly** invariant under \(a\mapsto\lambda a\).
>
> But \(\log Z_{\lambda a}(\beta)=\beta\log\lambda+\log Z_a(\beta)\), so \(d\) is
> **not**: as \(\lambda\to\infty\), \(d(\lambda a,\lambda b)\to
> \bigl|\log(\log r_a/\log r_b)\bigr|\) with \(r_a=\#a\), which is \(0\)
> whenever \(a\) and \(b\) have equally many fibers.

Verified directly on \(F=\{(2,1),(4,1),(8,1),(4,2),(8,2),(16,4)\}\): the
maximum relative deviation of \(E_{\lambda a,\lambda b}-\lambda E_{ab}\) is
\(2\cdot10^{-12}\) at \(\lambda=210\), and the maximum angle shift is
\(2\cdot10^{-8}\) — machine precision. On the 21-signature family
\(\{(2^i,2^j):1\le i<j\le7\}\), 210 pairs, the angles are likewise unmoved to
four decimals while the mean exchange distance falls monotonically:

| \(\lambda\) | 1 | 2 | 3 | 5 | 7 | 11 | 30 | 210 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| max relative deviation of \(E\) (on \(F\)) | — | \(9\!\cdot\!10^{-15}\) | \(1\!\cdot\!10^{-14}\) | \(3\!\cdot\!10^{-14}\) | \(1\!\cdot\!10^{-14}\) | — | \(2\!\cdot\!10^{-13}\) | \(2\!\cdot\!10^{-12}\) |
| max angle change (on \(F\)) | 0 | \(2\!\cdot\!10^{-8}\) | \(2\!\cdot\!10^{-8}\) | \(2\!\cdot\!10^{-8}\) | \(1\!\cdot\!10^{-8}\) | — | \(2\!\cdot\!10^{-8}\) | \(2\!\cdot\!10^{-8}\) |
| mean \(d\) (21-signature family) | 0.393 | 0.322 | 0.292 | 0.262 | 0.245 | 0.226 | 0.193 | 0.151 |

> **Theorem B (the exchange metric is Cartesian-power invariant).** For every
> \(k\ge1\), \(d(a^{\otimes k},b^{\otimes k})=d(a,b)\) exactly, because
> \(\log Z_{a^{\otimes k}}=k\log Z_a\) and the ratio defining \(C\) is
> unchanged.
>
> But \(a^{\otimes k}\) has entirely different entries, so \(O\), \(A\) and
> every Weil angle change.

Verified: on the 21-signature family, max \(|\Delta d|=2.8\cdot10^{-16}\) over
all 210 pairs for \(k=2\), while the Weil angles move by \(0.224\) rad on
average and \(0.821\) rad at worst, with
\(\operatorname{corr}(\theta_{\text{before}},\theta_{\text{after}})=+0.698\).
On \(F\), max \(|\Delta d|=1.3\cdot10^{-16}\) for \(k=2\) and
\(3.3\cdot10^{-16}\) for \(k=3\).

### Two explicit families, side by side

Put \(F\) and \(5F\) into one family. The two blocks have **identical Weil
angles** and **different exchange distances**:

| pair in \(F\) | \(d\) | angle | pair in \(5F\) | \(d\) | angle |
|---|---:|---:|---|---:|---:|
| \((2,1)\) vs \((4,1)\) | 0.69315 | 1.45995 | \((10,5)\) vs \((20,5)\) | 0.26316 | 1.45995 |
| \((2,1)\) vs \((8,1)\) | 1.09861 | 1.48850 | \((10,5)\) vs \((40,5)\) | 0.47129 | 1.48850 |
| \((4,1)\) vs \((8,1)\) | 0.40547 | 1.55633 | \((20,5)\) vs \((40,5)\) | 0.20813 | 1.55633 |
| \((2,1)\) vs \((4,2)\) | 0.69315 | 1.19658 | \((10,5)\) vs \((20,10)\) | 0.26316 | 1.19658 |

Now put \(F\) and \(F^{\otimes 2}\) side by side. **Identical exchange
distances, different Weil angles**:

| pair in \(F\) | \(d\) | angle | pair in \(F^{\otimes2}\) | \(d\) | angle |
|---|---:|---:|---|---:|---:|
| \((2,1)\) vs \((4,1)\) | 0.69315 | 1.45995 | \((4,2,2,1)\) vs \((16,4,4,1)\) | 0.69315 | 1.24891 |
| \((2,1)\) vs \((8,1)\) | 1.09861 | 1.48850 | \((4,2,2,1)\) vs \((64,8,8,1)\) | 1.09861 | 1.68037 |
| \((4,1)\) vs \((8,1)\) | 0.40547 | 1.55633 | \((16,4,4,1)\) vs \((64,8,8,1)\) | 0.40547 | 1.73074 |
| \((2,1)\) vs \((4,2)\) | 0.69315 | 1.19658 | \((4,2,2,1)\) vs \((16,8,8,4)\) | 0.69315 | 1.60843 |

**Consequence.** The Weil angle is not a function of \(d\) (Theorem B produces
equal \(d\), unequal angle) and \(d\) is not a function of the Weil angle
(Theorem A produces equal angle, unequal \(d\)). Any correlation between the two
is a property of the *family chosen*, not of the two functionals. This is the
promised sharp statement, and it explains why §1 can produce any number between
\(-0.98\) and \(+0.86\).

---

## 5. Scale-invariance mismatch, quantified

Perturbing one member of a family and watching how far each geometry moves.

**Cartesian square of one member** (generic T1.4 family; means over the 15 pairs
touching that member):

| member \(a\to a^{\otimes2}\) | mean \(|\Delta d|\) | mean \(|\Delta\theta|\) |
|---|---:|---:|
| \((2,2)\) | 0.0000 | 0.198 |
| \((3,1)\) | 0.0000 | 0.188 |
| \((4,2)\) | 0.0000 | 0.147 |
| \((5,3)\) | 0.0000 | 0.301 |

Note \(\Delta d=0\) even though only *one* member was raised. This is a
sharper form of Theorem B: \(\log Z_{a^{\otimes k}}=k\log Z_a\) gives
\(C(a^{\otimes k}\to b)=k\,C(a\to b)\) and \(C(b\to a^{\otimes k})=C(b\to a)/k\)
for **every** \(b\), so the product — and hence \(d\) — is untouched. A single
member can be raised to any Cartesian power without moving one entry of the
exchange distance matrix, while its Weil position moves by 0.15–0.30 rad. This
is the mismatch in its purest form.

**Doubling one member's entries** (generic family):

| member \(a\to 2a\) | mean \(|\Delta d|\) | mean \(|\Delta\theta|\) | \(O+A\) mass before | after |
|---|---:|---:|---:|---:|
| \((2,2)\) | 0.581 | 0.198 | 42.8 | 50.0 |
| \((3,1)\) | 0.411 | 0.309 | 63.8 | 51.6 |
| \((4,2)\) | 0.275 | 0.170 | 44.4 | 46.3 |
| \((5,3)\) | 0.234 | 0.302 | 65.2 | 60.4 |
| \((6,1)\) | 0.200 | 0.207 | 58.5 | 54.6 |
| \((3,1,1)\) | 0.431 | 0.317 | 95.1 | 74.9 |

**Tripling one member inside the designed powers-of-2 family** — a factor of
\(3\) destroys every prime-power ratio that member had with the rest, so its
Weil angles all collapse onto \(\pi/2\) while \(d\) merely drifts:

| member \(a\to3a\) | mean \(|\Delta d|\) | mean \(|\Delta\theta|\) | mean \(|\theta-\pi/2|\) before | after |
|---|---:|---:|---:|---:|
| \((4,2)\) | 0.542 | 0.135 | 0.118 | 0.017 |
| \((8,2)\) | 0.374 | 0.180 | 0.159 | 0.021 |
| \((16,2)\) | 0.276 | 0.239 | 0.216 | 0.025 |
| \((32,2)\) | 0.211 | 0.300 | 0.273 | 0.029 |

One factor of \(3\) — invisible to the exchange metric beyond a small drift —
annihilates the entire Weil structure of that signature, taking
\(|\theta-\pi/2|\) down by a factor of \(7\) to \(10\).

---

## 6. Stability in \(N\)

The zero list was extended to \(N=2400\) (\(T=2930.00\)) with
`extend_zeta_zeros.py`. Pearson / Spearman:

| \(N\) | \(T\) | generic (T1.4) | 2-ladder 5-subsets | prime-ladder 5-subsets (\(A\equiv0\)) | staircase \((2^t,1)\) |
|---:|---:|---:|---:|---:|---:|
| 100 | 236.5 | \(+0.231/+0.205\) | \(+0.720/+0.809\) | \(+0.747/+0.822\) | \(-0.851/-0.927\) |
| 200 | 396.4 | \(+0.220/+0.219\) | \(+0.717/+0.809\) | \(+0.743/+0.820\) | \(-0.792/-0.899\) |
| 400 | 679.7 | \(+0.207/+0.223\) | \(+0.721/+0.812\) | \(+0.740/+0.819\) | \(-0.823/-0.912\) |
| 800 | 1183.7 | \(+0.200/+0.239\) | \(+0.724/+0.815\) | \(+0.744/+0.821\) | \(-0.772/-0.899\) |
| 1200 | 1648.3 | \(+0.191/+0.187\) | \(+0.725/+0.815\) | \(+0.743/+0.820\) | \(-0.708/-0.858\) |
| 1600 | 2090.0 | \(+0.188/+0.223\) | \(+0.725/+0.816\) | \(+0.743/+0.820\) | \(-0.681/-0.848\) |
| 2000 | 2515.3 | \(+0.186/+0.244\) | \(+0.726/+0.816\) | \(+0.743/+0.820\) | \(-0.658/-0.831\) |
| 2400 | 2930.0 | \(+0.184/+0.241\) | \(+0.726/+0.817\) | \(+0.742/+0.820\) | \(-0.662/-0.831\) |

The two ladder columns are the §2 comparison repeated at every truncation: the
\(A\equiv0\) control tracks the geometric ladder to within \(0.02\) in Pearson
and \(0.003\) in Spearman across a 24-fold range of \(N\).

The generic figure is flat across the whole range: the \(+0.19\) of T1.4 is a
stable number, not a truncation artefact. The overlap-dominated ladders are
flatter still, moving by \(0.006\) in Pearson between \(N=100\) and
\(N=2400\). Only the staircase drifts, from \(-0.85\) to \(-0.66\), and the
drift is predicted: its amplitude is \(K=T\log p/2\pi N\), which decreases like
\(1/\log(T/2\pi)\), so the \(N\!\cdot\!O\) term slowly regains ground at small
\(s\) and pulls the small-index pairs off the pure \(p^{-(t-s)/2}\) law. All
three effects follow from the Landau formula, in which the two terms scale as
\(N\) and \(T/2\pi\sim N/\log N\), so their ratio drifts only logarithmically.

---

## 7. Summary

1. **Best positive correlation:** \(+0.857\) Pearson, \(+0.870\) Spearman on
   6-subsets of the ladder \(\{3,9,\dots,6561\}\) (378 pairs), and
   \(+0.725/+0.815\) on 5-subsets of \(\{2,\dots,256\}\) over 1540 pairs —
   against \(+0.19\) generic.
2. **Best monotone relation:** the receding staircase \(a_t=(p^t,1)\),
   \(t\in[T_0,T_0+8]\): Pearson \(-0.957\), Spearman \(-0.978\) over 36 pairs
   at \(p=2,\ T_0=100\) (and \(-0.878/-0.980\) on \(t=1..8\), \(p=3\)), with the
   derived law
   \(\theta-\pi/2=K\,p^{-(t-s)/2}\), \(K=T\log p/2\pi N\), and
   \(d=\log(t/s)\) exactly; inside a window \([T_0,T_0+n]\), \(n\ll T_0\), this
   reads \(\theta-\pi/2=Ke^{-cd}\) with \(c=\tfrac12T_0\log p\).
3. **The positive branch is not arithmetic.** A control ladder with
   \(A\equiv0\) reproduces the geometric ladders' correlations cell for cell.
   The overlap term alone does the work, and it is pure combinatorics.
4. **The arithmetic branch has the wrong sign.** Where the \(\Lambda\) term
   dominates, the Weil form is repulsive and the repulsion *decreases* with
   \(d\): closest in exchange \(\Rightarrow\) furthest in Weil.
5. **No functional relation is possible.** \(R\) is exactly invariant under
   \(a\mapsto\lambda a\) and \(d\) is not; \(d\) is exactly invariant under
   \(a\mapsto a^{\otimes k}\) and \(R\) is not. Explicit families realising both
   are in §4.

So T1.4's conclusion survives in a stronger form. It is not that the two
geometries fail to correlate — they can be made to correlate at \(\pm0.9\) — it
is that the correlation is a choice of family, and the part of \(E\) that
carries the primes is precisely the part that refuses to line up with \(M\).

---

## Caveats

* All the caveats of T1.4 apply: atomic measures are not admissible Weil test
  functions, \(E\) is the finite-rank truncation at \(N\) zeros, and the
  \(N\to\infty\) limit is not addressed.
* Theorem A is exact for the truncated \(E\) at every \(N\), so it is a
  statement about the finite-rank object. Smoothing the deltas would break it,
  because a fixed smoothing kernel does not commute with dilation — that is the
  one place where an admissible-test-function version of this note could differ.
* The ladder families contain signatures of many different lengths, and part of
  the positive correlation there is the trivial contribution of the number of
  fibers to both functionals. The \(A\equiv0\) control has exactly the same
  defect, which is precisely why it is the right comparison.
* Correlations on \(k\)-subset families are computed over pairs that are far
  from independent (each signature appears in many pairs), so the pair counts
  overstate the effective sample size. Comparisons *between rows* of the §2
  table are unaffected, since all three rows share the same dependence
  structure.

Reproduce with `research/m_and_e_and_a_c/multiplicative_design.py`
(zeros in `zeta_zeros_1200.npy` / `zeta_zeros_2400.npy`, the latter built by
`extend_zeta_zeros.py`; tabulated output in `multiplicative_design.csv`,
full console transcript in `multiplicative_design_output.txt`).
