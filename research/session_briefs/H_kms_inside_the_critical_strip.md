# Session brief H — what replaces the exchange rate inside the critical strip?

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (branch; commit or
stash first).

**Read first:** `research/session_briefs/D_quantum_mechanics_of_the_exchange_matrix.md`
(all of it, and Part 3 item "Bost–Connes" especially — that brief calls this
"the most promising single item" and then defers it to a programme; this brief
is that programme's first session);
`paper_finite_fields_maps/docs/exchange_positivity_and_weil.md` §3;
`paper_finite_fields_maps/docs/riemann_hypothesis_exchange_matrices.md`;
`research/m_and_e_and_a_c/FINDINGS.md` T1.4 (the caveat about admissibility).

## The obstruction, stated exactly

The whole framework rests on

```
C(g→f) = inf_{β ≥ 0} log Z_g(β) / log Z_f(β),    Z_a(β) = Σ_i a_i^β
```

and the project's zeta section observes that the tensor product of all local
prime modes has partition function `ζ(β)`:

```
⨂_p P_{p,∞} = {1, 2, 3, …},     Z(β) = ζ(β)
```

**But `ζ(β)` converges only for `Re β > 1`.** The critical strip — the only place
RH lives — is exactly where the monotone diverges and `C` is undefined. Every
positive result in this project lives at `Re β > 1` or on the real axis; every
statement about RH lives at `Re β = ½`. The two never meet, and no amount of
work on either side closes the gap. This is stated as a caveat in the notes; it
is in fact the central structural obstruction of the entire programme.

**There is a developed theory of exactly this divergence.** The infinite tensor
product above is the *primon gas* / Riemann gas (B. Julia, 1990). Its
non-commutative completion is the **Bost–Connes C\*-dynamical system**
`(A, σ_t)`, whose partition function is `ζ(β)` and which has a **KMS phase
transition at `β = 1`**:

* for `β > 1` there is a unique KMS`_β` state, and it is the Gibbs state — which
  reproduces `Z(β) = ζ(β)` and hence the existing framework;
* at `β ≤ 1` the Gibbs state does not exist; KMS`_β` states still do, they form a
  simplex, and for `β ≤ 1` the symmetry group `Ẑ* ≅ Gal(Q^{ab}/Q)` acts on them —
  spontaneous symmetry breaking, with the class field theory of `Q` as the
  symmetry.

So the divergence is not a defect of the framework. **It is a known physical
phase transition with an arithmetic order parameter, and the object that
replaces the Gibbs state below it is classified.**

## The question

> **Does the KMS condition define a comparison of resources on the critical
> strip, and does it reduce to `C(g→f)` above `β = 1`?**

This is the only route in the whole programme that reaches `Re β = ½` at all.
Whether it reaches anything *useful* there is exactly what this session decides.

## The checkpoint — one session, and stop

Brief D specifies a first checkpoint and it is the right one. Do this and
nothing beyond it before reporting:

1. **Rebuild the primon gas explicitly** as a `C*`-dynamical system, at the
   level of concreteness the repo can compute with: Hamiltonian `H|n⟩ = log n|n⟩`
   on `ℓ²(N)`, `Tr e^{−βH} = ζ(β)`. Confirm the Gibbs state at `β > 1` and its
   partition function.
2. **Verify the framework is the `β > 1` Gibbs restriction.** Compute the
   exchange rate against a truncated Euler factor and check it against
   `analysis/xi_versus_euler_factors.py`, which already exists. If they do not
   agree, that mismatch is the finding and everything else waits.
3. **Identify precisely what breaks at `β = 1` in resource language.** The
   partition function diverges; the Gibbs state ceases to exist; the "resource"
   `⨂_p P_{p,∞}` has infinite entropy at that temperature. Say which of the
   framework's axioms fails first — the infimum, the normalisation, or the
   monotonicity — and at what rate.
4. **State the KMS replacement.** For `β ≤ 1` the extremal KMS`_β` states are
   indexed (Bost–Connes) by embeddings, with a `Ẑ*` action. Write down what
   "`C(g→f)`" would have to mean when the two states are extremal KMS states
   rather than Gibbs states. It will not be an infimum of a ratio of partition
   functions, because there is no partition function. **The candidate to test:
   the relative entropy / Connes cocycle Radon–Nikodym derivative `[Dφ : Dψ]_t`,
   which exists for any pair of states and reduces to the ratio of Gibbs weights
   in the type-I case.** Check that reduction explicitly at `β > 1`; that check
   is the entire scientific content of the checkpoint.
5. **Then decide, in writing, whether to continue.** The failure mode is
   real and likely: the Bost–Connes system is type III at `β ≤ 1`, there is no
   trace, and "exchange rate" may simply have no meaning there. If so, say so and
   stop — a written negative result naming the type-III obstruction is worth
   more than three sessions of analogy.

## What would make this more than a restatement

Two things, and only these two count:

* **A comparison that is defined at `β = ½` and reduces to `C` at `β > 1`.**
  Anything that is only defined where `C` already is, is a renaming.
* **The `Ẑ*` symmetry doing work.** If the extremal KMS states below `β = 1`
  give a *family* of comparisons permuted by `Gal(Q^{ab}/Q)`, then the exchange
  order is Galois-equivariant and that is a genuinely new statement. If the
  Galois action turns out to be invisible to any resource-theoretic quantity,
  that is the honest negative and is also worth writing.

## Prior art a referee will assume you know

Bost–Connes (1995); Connes' trace formula and the adele class space; Connes–Marcolli
on quantum statistical mechanics and number theory; Julia's Riemann gas;
Berry–Keating `xp`; Bender–Brody–Müller (2017). None of it is shortened by
anything in this repo. **Read enough to place the contribution, and if the
placement is "this is Bost–Connes restated in exchange language", write that
sentence and stop.** Brief D already commits the project to claiming only two
new things on the quantum side — the gauge decomposition and the `t = ½`
Szegedy point. Do not add a third unless it survives this scrutiny.

## Standing obstructions that this brief does *not* remove

Both are recorded in FINDINGS and must be repeated in anything written:

* atomic measures `Σ_i δ_{a_i}` are **not** admissible Weil test functions, so
  the matrix `E` is a finite-rank truncation and `|Z_a(½+iγ)|` does not decay;
* the exchange monotone diverges in the critical strip (the subject of this
  brief).

The KMS route addresses the second and not the first. Say so.

## Traps

* Do not build anything on `MᵀM` or `MᵀDM`: PSD for free, informative never
  (brief D, Part 2).
* Two unrelated one-parameter families in the notes are both called `t` — the
  Schoenberg family `e^{−t d}` and the Maslov/idempotent family with `ħ = 1/t`.
  Rename before writing.
* `exchange_rate` is good to `~1e−13`; the Euler-factor comparisons in
  `xi_versus_euler_factors.py` are truncations and their truncation error must be
  reported alongside any agreement claim.
* Numerical work near `β = 1` needs care: `ζ(β) ~ 1/(β−1)`, so ratios of
  logarithms are ill-conditioned there. Use mpmath at 40 digits, as T1.1 did.

## Success criterion

Steps 1–4 completed with numbers, and a written decision on step 5 with its
reason. A one-paragraph honest "this is Bost–Connes restated, the type-III
obstruction is X, do not continue" is a full success for this session. So is a
comparison functional defined at `β = ½`. Anything in between should be labelled
as which.

## Reproduce / build on

`analysis/xi_versus_euler_factors.py`, `analysis/prime_mode_entropy_energy_curves.py`,
`analysis/xi_gibbs_curve.py`, `analysis/zeta_entropy_energy_curve.py`.
