# Findings — is the exchange rate superadditive under tensor products?

Answer to the question: can `C(a⊗b → c) > C(a→c) + C(b→c)`, for quadratic or
cubic map classes, and for the affine rate `C_aff`?

## Summary

**Yes, always — and that is the uninteresting half.** Superadditivity is a
two-line consequence of the definition and holds for *every* triple. The content
is entirely in **when it is strict**, and there the answer is exact:

> `C(a⊗b→c) = C(a→c) + C(b→c)` **iff the two ratio functions
> `log Z_a/log Z_c` and `log Z_b/log Z_c` share a minimiser.**
> Strict superadditivity ⟺ the two contact temperatures are disjoint.

So the gap is a **temperature-mismatch functional**: it measures the extent to
which `a` and `b` are useful for making `c` at *different* temperatures. Two
resources whose contacts differ are complementary, and combining them beats
using them separately.

For `C_aff` superadditivity also holds, by block-diagonal composition — the same
argument the paper already uses for Fekete's lemma. **Whether it is ever strict
is open, and that is the version of the question worth a session** (§4).

---

## 1. The proposition *(proved)*

The Cartesian product of maps multiplies fiber sizes, so the signature of
`a⊗b` is the multiset `{a_i b_j}` and

```
Z_{a⊗b}(β) = Z_a(β)·Z_b(β),    hence   log Z_{a⊗b} = log Z_a + log Z_b.
```

Therefore, with `F = log Z`,

```
C(a⊗b→c) = inf_β (F_a + F_b)/F_c ≥ inf_β F_a/F_c + inf_β F_b/F_c
         = C(a→c) + C(b→c),
```

because the infimum of a sum is at least the sum of the infima. Equality holds
iff some `β` minimises both `F_a/F_c` and `F_b/F_c` simultaneously. ∎

*Verified:* the minimum gap over every triple of every pool below is
`−8.9·10⁻¹⁶`, i.e. zero to floating point.

---

## 2. Quadratic and cubic classes *(computed, exhaustive)*

Signatures of all maps `F_q^2 → F_q` of the stated degree, non-constant only.

| pool | signatures | triples | strictly superadditive | largest gap |
|---|---|---:|---:|---:|
| `F_3` quadratic | `(6,3) (5,2,2) (4,4,1) (3,3,3)` | 64 | **20** | 0.317394 |
| `F_5` quadratic | `(10,10,5) (9,4,4,4,4) (6,6,6,6,1) (5,5,5,5,5)` | 64 | **20** | 0.267513 |
| `F_3` cubic | the above plus `(7,1,1)` | 125 | **50** | 0.317394 |

The extremal `F_3` example, and it is typical of all of them:

```
a = (6,3)      C(a→c) = 0.630930   contact β = 0
b = (3,3,3)    C(b→c) = 0.682606   contact β = ∞
c = (5,2,2)    C(a⊗b→c) = 1.630930          gap = 0.317394
```

**In every strict case in these three pools the two contacts sit at *opposite
endpoints*, `β = 0` and `β = ∞`.** That is the whole mechanism at this size: `a`
is the better resource at the fiber-count end, `b` at the largest-fiber end, and
the product exploits both. The classes are too small to produce an
interior-driven example.

## 3. Interior contacts *(computed)*

Random integer signatures (2–7 fibers, values ≤ 30), 6000 sampled triples:

* 2287 of 4000 triples strict in the first pool;
* **1238 strict cases have *both* contacts interior and distinct** — e.g.
  `(19,6,2) ⊗ (30,27,24,23,20,7) → (12,7,4)` with `β_a = 0.1744`,
  `β_b = 298.83`, gap `0.193798`, confirmed at 40 digits on a `10⁻⁴…10⁴` grid
  (`0.1937981452492805352100597834161352469248`).

So the phenomenon is not an endpoint artefact; it just needs bigger pools than
the `F_3`/`F_5` classes provide.

**Numerical caveat, and it matters.** 428 sampled triples had reported contacts
differing by as much as `8·10²` yet a gap below `10⁻⁹`. Two causes, both known
to this project: minima of `F/F_c` can be extremely flat, so moving from one
contact to the other costs almost nothing; and the solver is documented to
report *spurious* interior contacts when the ratio is constant to double
precision (the trap recorded in `curve_family_cycles/FINDINGS.md`). **"Contacts
differ" is therefore necessary for strictness but is not a reliable numerical
test of it** — test the gap, not the argmin.

---

## 4. The affine rate `C_aff` *(superadditivity proved; strictness open)*

With `k_{g→f}(r) = max{k : f^{×k} ⪯_aff g^{×r}}` and
`C_aff(g→f) = lim_r k_{g→f}(r)/r`:

**Proposition.** `C_aff(a⊗b→c) ≥ C_aff(a→c) + C_aff(b→c)`.

*Proof.* If `a^{×r}` implements `c^{×k₁}` and `b^{×r}` implements `c^{×k₂}`, put
the two block processors side by side: `(a⊗b)^{×r} = a^{×r} ⊗ b^{×r}` implements
`c^{×(k₁+k₂)}`. So `k_{a⊗b→c}(r) ≥ k_{a→c}(r) + k_{b→c}(r)`; divide by `r` and
take limits. ∎

This is the same block-diagonal argument the paper already uses to get
`k(r+s) ≥ k(r) + k(s)` and hence Fekete's lemma.

**Why the affine version is the interesting one.** `C_aff` has no variational
formula, so there is no "contact temperature" to appeal to, and strictness would
not be a statement about where an infimum sits. It would say that **affine
processors achieve jointly what they cannot achieve separately** — and there is a
concrete mechanism available, because the paper is explicit that *the block
processor is allowed to mix the `r` copies*. A strict gap would have to come from
mixing `a`-copies with `b`-copies, which is exactly the operation unavailable to
either rate alone. That is genuine synergy rather than a temperature artefact.

Note also `C_aff ≤ C`: affine implementation implies signature implementation, so
the signature rate is an upper bound and the results of §2–§3 bound the affine
gap from above but say nothing about whether it is positive.

**Open, and the natural next session:** exhibit quadratic or cubic classes over a
small `F_q` with `C_aff(a⊗b→c) > C_aff(a→c) + C_aff(b→c)`, or prove that affine
mixing never helps. The `F_3` quadratic and cubic classes are small enough that
`k_{g→f}(r)` is computable exactly for small `r`.

---

## 5. Is it interesting? An honest assessment

**The inequality: no.** It is two lines, it is automatic, and it holds for every
triple. It should be stated as a remark and not as a result.

**The strictness criterion: mildly, and as a third instance of one fact.** This
project has now found three separate phenomena that all reduce to *the infimum
being attained at different temperatures*: the cycles (an interior tangency
overturning the endpoint index), the curl of the gauge field, and now the
synergy gap. The gap is the **weakest** of the three as a diagnostic — it is
nonzero whenever two contacts merely differ, which happened in 2287 of 4000
random triples, whereas 3-cycles occur in about 0.07% of random 4-subsets. So it
detects *non-degeneracy*, not *non-scalarity*, and it carries no information
beyond the contact temperature, which is already recorded throughout the project.

**The affine version: yes.** It is the only form of the question whose answer is
not already implied by the variational formula.

**One notion worth keeping.** The gap defines *complementarity of resources
relative to a target*: `a` and `b` are complementary for `c` when their contact
temperatures differ, and then the pair is strictly better than the sum of its
parts. That is a natural resource-theoretic concept, it is exactly computable,
and this project has not named it before.

---

## Files

| file | what |
|---|---|
| `synergy.py` | enumerates the quadratic and cubic classes, computes the gap tables of §2–§3, and the 40-digit confirmation |
