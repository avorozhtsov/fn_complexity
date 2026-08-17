# Session brief E — is the flux arithmetic?

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (work on a branch;
commit or stash first — uncommitted work is invisible to a new worktree).

**Read first, in this order:**
`research/session_briefs/D_quantum_mechanics_of_the_exchange_matrix.md` Part 0
(the gauge decomposition — take it as given, do not re-derive);
`research/curve_family_cycles/FINDINGS.md` (the certified `F_11` cycle, the
census, the `β ≍ √q` regime);
`research/m_and_e_and_a_c/FINDINGS.md` (Notation, T2.1, T2.2, T2.3);
`research/session_briefs/C_no_scalar_invariant.md` (the question this brief
finally has the tools to answer).

## Where the programme stands

Two sessions converged, independently, on the same decomposition. Writing
`L(a,b) = −log C(a→b)` and `u_a = log log Z_a`:

```
L = S + A,   S = ½ osc(u_b − u_a) = d/2   (symmetric, the metric)
             A = mid(u_b − u_a)           (antisymmetric, a lattice 1-form)
a ≺ b  ⟺  A(a,b) > 0
```

Everything published so far — negative type, the five-point certificate, the
`l2`-distortion `1.3375`, the PSD ray — is a statement about `S` alone. **`S` is
comparison-blind.** `A` carries the entire order, the cycles, and by brief A's
accounting the only arithmetic content that is not a monotone relabelling of a
classical statistic. `A` has never been studied.

Brief B then found that `A` has genuine curvature over arithmetic objects: 132
strict three-cycles in a complete enumeration of genus-two pencils
`y² = P(x) + c` over `F_11`, 1475 over `F_13`, one of them certified by interval
arithmetic with margins `1.7·10⁻³` to `4.7·10⁻³`. And brief D showed the
endpoint regime is exactly the flat-connection locus, `A = dψ` with
`ψ = ½ log φ`, verified to `1.1·10⁻¹⁶` on 328 endpoint pairs.

So the object exists, it is non-trivial, and it is defined on curves. Nobody has
asked what it *is*.

## The question

**Is the flux arithmetic?**

Three sub-questions, in order. They are independent enough that a negative
answer to one does not sink the others.

### E1 — Hodge decomposition over an arithmetic pool

On the complete graph of the `F_11` (296 signatures, exhaustive) and `F_13`
(698, exhaustive) genus-two pools, decompose

```
A  =  grad ψ_opt  +  curl-free-residual
```

by least squares (HodgeRank; `ψ_opt` is minus the row-mean of `A`, one line).
Report the energy split `‖grad‖/‖A‖` and `‖residual‖/‖A‖`.

A calibration run on 8/16/24 random *integer* signatures (not curves) gives
`‖grad‖/‖A‖ ≈ 0.996–0.998` and residual `0.065–0.088`, with 0, 3 and 2 strict
3-cycles out of 56, 560 and 2024 triangles. **Reproduce that first as a smoke
test**, then run the arithmetic pools and compare. The interesting outcome is
either "arithmetic pools have far more curl than random ones" (the flux is
arithmetic) or "the same" (the flux is a generic feature of nearly flat
signatures and the arithmetic is only in which signatures occur).

### E2 — what is `ψ_opt`?

The gradient part is the *best possible scalar complexity* for genus-two pencils
over `F_q`. It exists whether or not it is exact. Identify it.

Concretely: regress `ψ_opt` against the statistics that are known to be visible
— `max_c N_c` (equivalently the extreme trace `M = max_c(−a_c)`), `m₂`, `m₃`,
`m₄`, the multiplicity `μ` of the largest fiber, `ν(P) = #{(x,x′):P(x)=P(x′)}`,
genus, and `φ = log q · log max_c N_c`. Report `R²` for nested models.

The honest target is a sentence of the form "the best scalar invariant of a
genus-two pencil over `F_q` is `ψ_opt = α·M + β·m₂ + …` to `R² = …`, and the
residual is the flux." **If `ψ_opt` turns out to be `½ log φ` plus a small
correction, say so** — that would mean the endpoint potential survives as the
best scalar even where the connection is not flat, which is a clean result.

Brief B's addendum predicts `φ̃ = M − ((3−2√2)/2)·m₂` as the two-term scalar,
but derived it in the `β = O(1)` regime; brief B's own measurement puts the
operative scale at `β ≍ √q`, `√q` deeper. **Expect `φ̃` to fail and record how.**

### E3 — brief C's crux, which is now decidable

The signature merges objects the geometry separates. A cycle among *signatures*
is not automatically a cycle among *curve families*. Settle it:

* For the certified `F_11` triple `A, B, C`, compute the actual arithmetic of the
  eleven fibers of each pencil: the Jacobians, their isogeny decomposition where
  computable, the trace vectors `(a_c)`, the monodromy of the branch map.
  **Do the three pencils differ by anything classical?** If two of them have
  isogenous Jacobian fibrations and still compare strictly, the comparison is
  reading something finer than isogeny; if the signature collapses a distinction
  the curves make, say which.
* Then the real question: **is the flux `A(f,g)` an invariant of the pair of
  families, or only of the pair of signatures?** It is manifestly the latter by
  construction. So the content is whether the map `family ↦ signature` is
  injective enough on the pools used. Measure it: how many distinct genus-two
  pencils share each signature at `q = 11, 13`? T2.2 recorded that genus ≥ 2
  gives 398–400 of 400 distinct, so the collapse should be mild — verify, and
  give the fiber-size distribution of the map.
* **The claim to make, and no more.** Brief C proposed "no real-valued invariant
  of such families is compatible with asymptotic conversion". With `A` in hand
  the precise version is: *there is no `φ` with `a ≺ b ⟺ φ(a) < φ(b)`, because
  `A` has non-zero curl; the best scalar approximation is `ψ_opt` and it explains
  `X%` of the comparison; the remaining `Y%` is irreducibly pairwise.* Write
  that, with `X` and `Y` measured. Do not write the unquantified version.

## Why this is the right next step

It is the only question in the programme that survives brief A's objection. The
one-way-flow argument says every *scalar* the rate reports is a monotone
relabelling of a statistic arithmetic already studies. That argument does not
touch `A`, because `A` is not a function of one family — it is a function of a
pair, and no classical invariant of pairs of curve families is on offer. If the
flux correlates with nothing classical, the programme's arithmetic claim is
honestly dead and should be retracted; if it correlates with something, that is
the paper.

## Traps

* `−½JDJ` has the constant vector in its kernel; work in an orthonormal basis of
  `{Σx = 0}` (FINDINGS T1.1).
* The `β` horizon must scale with `q`: these signatures have adjacent entries
  differing by one, so the largest fiber is not isolated until `β ≳ 36q`. Brief B
  ran grids to `360q`. A grid truncated at `β ≈ 500` hides everything.
* The package can report a spurious interior contact when the infimum is the
  `β = ∞` endpoint (brief B saw `β = 258.9` reported for an exact endpoint).
  Test every contact against `log(max_g)/log(max_f)` before believing it.
* `exchange_rate` is accurate to `~1e−13`; treat differences below `1e−10` as
  ties. Margins here are `10⁻³`–`10⁻⁴`, comfortably above.
* Interval arithmetic on these objects needs `q^β` factored out first:
  `log Z = β log q + log S`, `S(β) = Σ_c (N_c/q)^β`. Without it the boxes must be
  `10⁻⁵` wide and branch-and-bound does not close.
* Genus one is useless — 400 random elliptic fibrations give `gcd(4,q−1)+1`
  signatures. Use genus ≥ 2.

## Success criterion

The energy split of `A` on both arithmetic pools and on a matched random-integer
control; an explicit `ψ_opt` with its regression table; and a settled statement
about whether the certified cycle descends to the curves. Negative results
count and must be written as such.

## Reproduce / build on

`research/curve_family_cycles/common.py` (pool + vectorised rate engine),
`search.py`, `certify.py`; `research/m_and_e_and_a_c/gauge_decomposition.py`.
