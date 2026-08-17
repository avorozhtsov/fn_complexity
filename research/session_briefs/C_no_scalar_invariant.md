# Session brief C — would a cycle mean no scalar invariant exists?

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (branch; commit or
stash first). **Depends on brief B**: this has content only if arithmetic
families do cycle. If B comes back negative, this becomes the question of
whether the *impossibility* is itself provable for curve families.

**Read first:** `research/m_and_e_and_a_c/FINDINGS.md`, and the cycle sections of
`paper_finite_fields_maps/docs/finite_field_exchange_matrix.md`.

## The claim to make precise

If three families of curves over `F_q` satisfy `f₁ ≺ f₂ ≺ f₃ ≺ f₁`, then **no
real-valued invariant of such families is compatible with asymptotic
conversion**: there is no `φ` with `a ≺ b ⟺ φ(a) < φ(b)`, because `≺` is not
transitive and any such `φ` would make it so.

This is the one statement in the whole programme that is *not* a repackaging of
classical data. Every other result — the extreme trace, the second moment, the
Weil bound — is a monotone function of a statistic that was already studied. A
cycle is different in kind: it says the comparison itself is not a function of
any single number, and that is a statement about arithmetic families which
cannot be phrased in terms of `m₂`, `max_c N_c`, genus, or monodromy separately.

## What the session should establish

1. The exact statement, with the hypotheses spelled out — which conversions are
   allowed, what "invariant" means (a function of the family? of the signature?
   of the isogeny class?), and over what class of families the claim ranges.
2. Whether the obstruction is genuinely arithmetic or an artefact of the
   signature being a coarse invariant. Note that the signature *merges* classes
   the geometry separates (linear and parabolic maps share `(q,…,q)`), so a
   cycle among signatures need not be a cycle among the underlying families.
   This distinction is the crux and must be settled before any claim is made.
3. The relation to the second place the matrix may exceed its scalar shadows:
   the pairwise distance `d(f,g) = −log(C(f→g)C(g→f))` is a joint invariant of
   two families, not reducible to comparing two numbers. Is it arithmetically
   meaningful — does it separate families that classical invariants do not?

## Honest framing to preserve

Non-transitivity is already proved for integer signatures and for `P²` morphisms;
what would be new here is that it happens for *arithmetic* families, where a
scalar complexity is exactly what one would expect to exist. Claim that and
nothing more.
