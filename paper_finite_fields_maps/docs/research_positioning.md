# Research positioning and novelty ledger

This note separates the paper's new angle from classical ingredients and open
speculation. It should be updated before submission after a dedicated
bibliographic review.

## The new angle

The paper uses one processor equation

\[
f=h_{\mathrm{out}}\circ g\circ h_{\mathrm{in}}
\]

as a narrative route through several structures:

1. invertible processors give left-right equivalence classes;
2. singular processors give exact reachability and a degeneration preorder;
3. quotienting mutual reachability gives a poset and Hasse diagram;
4. Cartesian powers give operational exchange rates;
5. composition gives directed triangle inequalities and no-arbitrage cycles;
6. forgetting to a finite-map fiber signature gives the partition-function
   observer solved in `paper_exchange_rate`.

The point is not that any one item is new in isolation. The candidate new
contribution is using this chain as one operational framework on the same
named polynomial-map classes, over several kinds of fields.

## Classical ingredients that must not be claimed as new

- Left-right equivalence of maps and its orbit language are classical in
  singularity theory. A modern arbitrary-characteristic reference is Dmitry
  Kerner, [*Orbits of the left-right equivalence of maps in arbitrary
  characteristic*](https://arxiv.org/abs/2111.02715).
- Algebraic-group orbit closures and one-parameter degenerations are classical;
  see David Birkes, [*Orbits of linear algebraic
  groups*](https://doi.org/10.2307/1970884). Our exact singular-processor order
  is related but should not be identified with orbit closure without proof.
- Classification of quadratic forms by rank, radical, Witt index,
  discriminant, square class, and the characteristic-two Arf/trace invariant
  is classical.
- Resource preorders, ordered commutative monoids, and asymptotic spectra are
  classical. The no-arbitrage cycle inequality is an elementary consequence
  of composing operational rates.
- Maximum/minimum mean cycles are classical; see Richard Karp,
  [*A characterization of the minimum cycle mean in a
  digraph*](https://doi.org/10.1016/0012-365X(78)90011-0).
- The interpretation of \(\zeta(\beta)\) as a partition function is not new.
  Bost and Connes constructed a quantum statistical system with this partition
  function in [*Hecke algebras, type III factors and phase transitions with
  spontaneous symmetry breaking in number
  theory*](https://doi.org/10.1007/BF01589495).
- The quantum-spectral analogy for zeta zeros is established research and
  remains conjectural; see Berry and Keating,
  [*The Riemann zeros and eigenvalue
  asymptotics*](https://doi.org/10.1137/S0036144598347497).
- Weil positivity and its equivalence to RH for all admissible test functions
  are classical. Sampling known critical-line zeros to make a finite Gram
  matrix is automatically positive semidefinite and is not an RH test.

## Strongest candidate novel results

Subject to a fuller literature search, the strongest candidates are:

1. **A uniform affine-processor atlas for quadratic maps
   \(\mathbb F_q^2\to\mathbb F_q\).** The normal-form ingredients are
   classical, but the exact all-\(q\) Hasse diagrams, class-size formulas at
   every node, consistent arrow convention, and executable SVG generator may
   be a new compilation.
2. **The exhaustive homogeneous and cubic datasets.** In particular, the
   verified \(\mathbb F_3\) class/cover counts and the 110-class generated
   quadratic-input preorder for cubic maps over \(\mathbb F_8\) are plausible
   original computations.
3. **Putting three geometries on the same classes.** The one-shot affine
   Hasse order, the finite-map signature exchange matrix, and its max-times
   cycle losses are computed and compared without conflating their processor
   models.
4. **The explicit novelty boundary.** The paper proves no-arbitrage for the
   operational matrix \(M\), gives a concrete counterexample showing that this
   does not imply PSD, and explains why the Weil matrix \(E\) is currently a
   separate construction.

These are best described as a new framework specialization, atlas, and
computational study. They should not be advertised as a new classification of
quadratic forms or a new result about RH.

## What the cycle spectrum can reveal

For a directed rate matrix, define

\[
r(C)=\left(\prod_{(i,j)\in C}M_{ij}\right)^{1/|C|}.
\]

This is a per-step retention factor. It can:

- quantify the least lossy round trip among a chosen number of distinct
  classes;
- check composition and no-arbitrage numerically;
- detect that a coarse observer has collapsed distinct algebraic classes.

For the seven anisotropic \(\mathbb Q_2\) forms, residue-signature cycles of
mean one show that the modulo-\(2^{10}\) signature observer cannot distinguish
several square classes. They do **not** prove affine reversibility. The exact
block-affine rate matrix is still unknown.

Ordinary eigenvalues of \(M\) answer a different question: powers of \(M\)
sum over alternative walks. A Perron eigenvalue greater than one therefore
does not contradict no-arbitrage, which follows one composed conversion path.

## The disciplined zeta story

For \(\Re s>1\),

\[
\zeta(s)=\prod_p\sum_{k\ge0}e^{-s k\log p}
\]

is exactly the joint partition function of independent prime occupation modes
with energies \(k\log p\). For real \(\beta>1\), this is an ordinary Gibbs
system. For complex \(s=\beta+it\) with \(\beta>1\),

\[
\frac{\zeta(\beta+it)}{\zeta(\beta)}
=\mathbb E_\beta[e^{-it\log N}]
\]

is a characteristic function.

The nontrivial zeros lie outside this normalizable Gibbs half-plane. Their
ordinates are safely described as Mellin frequencies in log scale. Calling
\(1/2+i\gamma\) a complex energy, or a zero a quantum composition of prime
atoms, is presently metaphor rather than a construction.

The completed function adds the archimedean Gaussian Mellin factor and the
Fourier-dual symmetry \(s\leftrightarrow1-s\), as made systematic by Tate's
thesis. Any future bridge to the exchange framework should be adelic and
functorial, not an entrywise transformation of a finite rate matrix.

## Safe one-sentence claim

> A processor framework introduced for comparing finite-map complexity gives
> a uniform, executable language for affine equivalence classes, exact
> degeneration posets, and asymptotic exchange diagnostics of polynomial maps
> over finite, real, complex, and local fields.
