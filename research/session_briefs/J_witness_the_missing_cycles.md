# Session brief J — witness `4·SU(2)` and `SU(2)×3·SU(2)` at genus 4

**Repo:** `fn_complexity`, work on a branch.

**Read first:** `research/sato_tate_limit/REPEATED_FACTOR.md` — all of it,
especially Theorem A and the Open section, which writes out the route this brief
asks you to take. Then `research/sato_tate_limit/FINDINGS.md` (the nine cycles
and the `k`-threshold table).

## The question

Brief F found nine 3-cycles in the `q → ∞` limit; the follow-up session closed
the smallest by exhibiting pencils for all three of its vertices. **Six of the
nine remain unwitnessed**, and they are the six widest. All of them need one or
both of two measures at genus 4:

```
4·SU(2)          Jac ∼ E⁴,  E varying
SU(2)×3·SU(2)    Jac ∼ E₁ × E²... precisely: Jac ∼ J₁ × E³ with E varying
```

Multiplicity two is proved arithmetic (Theorem A: `f` even and palindromic of
degree `2m`, `m` odd ⟹ `Jac ∼ J²`; the pencil
`y² = (x²+1)(x⁴+cx²+1)` has `Jac ∼ E_c²` with `j(E_c) = 256(c+1)³/(c+2)`
non-constant, verified on 32721 fibres over 29 primes with zero mismatches).
Multiplicity three is proved arithmetic at genus 3 (a `(ℤ/2)²`-cover branched at
`{a,1,ζa,ζ,ζ²a,ζ²}` gives `Jac ∼ E³`). **Multiplicity four is unknown in any
construction.**

## The route, already written out

From REPEATED_FACTOR.md's Open section — this is the concrete plan, execute it:

> Take a `(ℤ/2)²`-cover with `|B₁| = |B₂| = 4`, `|B₁ ∩ B₂| = 1`, hence `|B₃| = 6`
> and quotient genera `(1, 1, 2)`. Take `B₃` to be the branch set of the
> `Jac ∼ E²` sextic `(x²+1)(x⁴+cx²+1)`, split it into two triples `S ⊔ T`, and
> put `B₁ = S ∪ {p}`, `B₂ = T ∪ {p}` for a seventh point `p`. Then
> `Jac ∼ J₁ × J₂ × E²` with `J₁, J₂` elliptic and two free parameters `(c, p)`.
> Imposing `J₂ ≅ E` (one condition) gives `SU(2) × 3·SU(2)`; imposing
> `J₁ ≅ J₂ ≅ E` (two) gives `4·SU(2)`.

`SU(2)×3·SU(2)` is therefore a **curve in the `(c,p)`-plane over `ℚ`**. The two
things to determine: is it rational (so that it carries a one-parameter family
over `F_q`), and **does `E` actually vary along it** — an isotrivial answer is
worthless and is the trap that has already cost this project time.

`4·SU(2)` is two conditions on two parameters, so generically a *point*, not a
family. If it comes out zero-dimensional here too, say so: combined with the
recorded fact that the even+palindromic route also gives a zero-dimensional
locus, that is real evidence — though **not a proof** — that `4·SU(2)` may be
unreachable. A genuine non-existence proof would be a strong result in its own
right, since it would permanently cap the six widest cycles.

**Non-hyperelliptic genus-4 curves were not examined at all.** If the
`(ℤ/2)²` route stalls, that is the obvious second front.

## The detector, which is cheap and decisive

With `α_c = −a_c/√q` and `m₂ = (1/q²)Σ_c a_c²`, `m₂ → E[α²]`:

| family | `m₂` |
|---|---:|
| `USp(8)` (generic genus 4) | 1 |
| four independent elliptic factors | 4 |
| `SU(2) × 3·SU(2)` | 10 |
| **`4·SU(2)`** | **16** |

so scan candidate pencils for `m₂ ≈ 16` (resp. 10) and confirm any hit fibre by
fibre against `a_c(C) = 4a_c(E)` (resp. `a_c(J₁) + 3a_c(E)`), **not by the
moment alone**. Two traps recorded by the previous session, both of which cost it
real time:

* a bielliptic split rational only over a quadratic extension gives the *swap*
  measure — `m₂ = 4`, not `8`. This is why the rationality condition
  `(r+1)(r−s) = □` is load-bearing. Check rationality over the base field.
* a symmetry that propagates to the other half gives `SU(2)³ × USp(4)` —
  `m₂ = 4`, not `3`. Verify the factorisation you think you have.

## If you succeed

Close the loop as the previous session did: confirm the limiting measure, plug it
into `research/sato_tate_limit/` and report which of the six cycles are now
witnessed, with midranges at 40 digits and margins. State the theorem precisely,
and **keep the standing caveat explicit**: this is a theorem about the limit.
`sup|Ψ_f − Ψ_μ| = 0.28` at genus 2 even at `q = 4·10⁵` against a `1.2·10⁻²`
margin, so `q₀` is astronomically beyond any census.

## Optional, if the main line finishes early

REPEATED_FACTOR.md lists **higher genus with multiplicity two** as the cheapest
way to a *wider* cycle: the `k ≤ 2` sub-cone was searched only to `α_max = 12`,
and wider margins mean a smaller `q₀`. Cost is `USp(16)`, a few hundred digits
per `τ`. This is lower-risk than the genus-4 hunt and may be the better use of a
second half-session.

## Working standards

Proved vs computed marked on every statement; negative results recorded as
results; isotriviality and reducible/unipotent monodromy checked explicitly
(brief F had to discard 7 of 10 seed measures as unrealisable by Deligne
semisimplicity plus a unipotent argument — apply the same scrutiny). Verify
headline numbers at 40 digits with mpmath.

## Deliverable

`research/sato_tate_limit/GENUS4_WITNESSES.md` in the house style (summary
first, sections marked proved/computed, corrections, open, files table), scripts
and CSVs alongside. **If the harness blocks you from writing that file, put its
complete body in your final message instead and say so** — that happened to an
earlier agent and the coordinator committed the file on its behalf.
Build on `curve_lib.py`, `repeated_factor.py`, `witness_search.py`.
