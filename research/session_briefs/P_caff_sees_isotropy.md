# Session brief P — does `C_aff` see the square class?

**Read first:** `research/affine_rate/FINDINGS.md` in full — especially the
reduction lemma, the support monotone, and Propositions 1–2, which are the
machinery this brief extends. Then
`paper_finite_fields_maps/docs/p_adic_quadratic_map_poset.md` (the seven
anisotropic square classes over `Q_2`) and
`paper_finite_fields_maps/docs/quadratic_map_posets.md`.

## The claim to test

`C_aff` and `C_sig` are incomparable (brief O), but on the `F_3` quadratics the
sign of the difference is decided by one-step implementability, and where
`C_aff` is exactly known off that locus it is always the smaller. The
interesting instance is this pair of rank-two forms over `F_3`:

```
xy      isotropic     C_aff(→ x) = 1       C_sig(→ x) = 0.9840
x²+y²   anisotropic   C_aff(→ x) = 2/3     C_sig(→ x) = 0.9753
```

They differ by **isotropy** — a genuine arithmetic invariant — and `C_aff`
separates them by `1/3` where `C_sig` separates them by `0.0087`. A factor of
38. Moreover the proof of `C_aff(x²+y² → x) ≤ 2/3` *uses* anisotropy: because
`−1` is not a square in `F_3`, `ℓ² + m² = 0` forces `ℓ = m = 0`, so the support
is two-dimensional and costs two atoms. For isotropic `xy` the argument
collapses and the rate is `1`.

> **Hypothesis.** `C_aff(x²+y² → x) = 2/3` when `q ≡ 3 mod 4` and `= 1` when
> `q ≡ 1 mod 4`, jumping at the congruence, while `C_sig` crosses it smoothly.

`C_sig(x²+y² → x)` is already measured and does *not* jump:
`0.975301, 0.993132, 0.988581, 0.992632, 0.998289, 0.998811, 0.995734` at
`q = 3, 5, 7, 11, 13, 17, 19`. Note it does move the right way — isotropic is
easier — but by `0.018` where the hypothesis predicts `0.333`.

**If the hypothesis holds, the headline is: the affine rate detects the square
class of the discriminant, and the signature rate very nearly erases it.** That
is the first arithmetic content in this project that is genuinely invisible to
`C_sig`, and it is worth a section of its own.

## What to do

1. **Settle the congruence.** Compute `C_aff(Q → x)` for the anisotropic binary
   form over `F_q` at `q = 3, 5, 7, 11, 13` (and further if cheap). Note
   `x²+y²` is anisotropic only for `q ≡ 3 mod 4`; the form that is anisotropic
   for **every** `q` is `x² − n y²` with `n` a non-square — use that as the
   uniform family and `x²+y²` as the congruence probe. Signatures are already
   known: `x² − n y²` has `(q+1)` repeated `(q−1)/2` times and a single `1`.
   * Proposition 2 of `affine_rate` should generalise verbatim: the argument only
     needs the form to be anisotropic, so `C_aff(\text{anisotropic} → x) ≤ 2/3`
     over **any** `F_q`. Check that, and check the matching witness `N_2 = 3`
     per `q` — the `F_3` witness is explicit in the findings and may transplant.
   * `x ⪯_aff xy` in one step over any field (`a(x,y) = (x,1)`), so the
     isotropic side is `1` for free. The content is entirely the anisotropic
     side.
2. **Close the six open brackets** of `affine_rate` §2 — `xy → x²+y²`,
   `xy → x²+y`, `x²+y² → xy`, `x²+y² → x²+y`, `x²+y → xy`, `x²+y → x²+y²`. They
   all have `N_1 = 2` and an upper bound of only `1`, and they decide a real
   conjecture: on this pool every exactly-known non-implementable pair has
   `C_aff < C_sig`, and six of the seven brackets top out *above* their `C_sig`.
   So **"`C_aff ≤ C_sig` off the implementable locus" is equivalent to
   `affine_rate`'s open item 2, "can `C_aff = 1` without `f ⪯_aff g`?"** Settle
   that question — it is the crux and it is stated there with the obstruction.
3. **Rationality.** Every `C_aff` computed so far is rational because `N_k` turns
   out eventually *exactly* linear (`N_{2m}(x²+y² → x) = 3m`, proved). Is that
   general? `C_aff = 1/\lim(N_k/k)` is a limit of rationals with denominator `k`
   and need not be rational a priori — this is the affine analogue of the
   rationality of the matrix-multiplication exponent, and saying so honestly is
   worth more than a guess. Report which cases are proved eventually-linear.

## The anisotropic-pair question, and why it is not a finite-field question

Asked directly: what is `C_aff` between two *anisotropic* classes?

Over `F_q` the question is **empty**, and that should be stated rather than
computed around. Binary quadratic forms over `F_q` are classified by their
discriminant modulo squares; the equivalence in use allows scaling by
`α ∈ F_q^×`; so there is exactly **one** anisotropic class, the norm form of
`F_{q²}`, and `C_aff(A → A) = 1` trivially. In three or more variables
Chevalley–Warning makes every form over a finite field isotropic, so there are
no anisotropic classes at all.

**So the question is `p`-adic.** Over `Q_p` anisotropic forms exist in up to four
variables, and `p_adic_quadratic_map_poset.md` already separates **seven**
anisotropic square classes over `Q_2`. That is the right arena, and it is the
natural continuation if item 1 comes back positive:

* Are distinct anisotropic square classes over `Q_p` distinguished by `C_aff`,
  and is the value ever strictly between `0` and `1`?
* The `p`-adic case has an infinite residue tower, so the reduction lemma of
  `affine_rate` §0 does not transplant unchanged — say precisely what breaks
  before attempting numbers. `p_adic_exchange_rate_attempts.md` records that
  both affine rates there are only bracketed in `[1/2, 1]`; that bracket is the
  state of the art and this session should not claim to beat it without a
  genuine argument.

Treat this section as **exploratory**: a clear statement of why the finite-field
version is empty, plus a scoped assessment of what the `p`-adic version would
need, is a full result for it. Do not let it consume the session — items 1 and 2
are the deliverable.

## Traps

* Enumerate **orbits** under the affine group, not maps; prune with the signature
  test as a cheap necessary condition (`affine_rate` traps).
* Fekete gives `C_aff = sup_k k/N_k`, so any computed `k/N_k` is a **rigorous
  lower bound** — report brackets, never point estimates without a matching
  upper bound.
* A claimed `N_k` needs the witness processors or an exhaustive-search argument.
* `C_sig` values here are infima over `β ∈ [0,∞]`; compute at ≥ 40 digits and
  check endpoint values against their closed forms (`log 3/log 6` and friends).
* Do **not** reuse the refuted claim `C_aff ≤ C_sig`; it is false in general and
  the whole point of item 2 is to find the domain where it is true.

## Deliverable

`research/affine_rate/ISOTROPY.md`, house style (summary; every claim marked
proved/computed; corrections; open; files table), scripts and CSVs alongside.
**If the harness blocks the write, paste the complete body in your final message
and say so** — that has happened to several agents here and the coordinator
commits the file on their behalf.
