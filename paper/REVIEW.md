# Working list — *Exchange Rates for Finite Map Signatures*

Status of `paper/main.tex`: 29 pp., builds clean (`latexmk -pdf main.tex`, zero
overfull boxes, zero warnings, all references resolve).

All numerics in the paper were independently recomputed for this review and
agree to the displayed precision: every rate in Table 1, both tangency
temperatures, the pair rates, and the shell cardinalities 129 / 2560 / 6738.

---

## Done

- **Introduction** — opens on `{2,2}` vs `{3,1}` with the products worked by
  hand; three results stated without terminology; abstract rewritten to match.
- **Proposition 1 fully proved** — counting lemma, Chernoff upper bound, method
  of types, uniform gap `δ = inf_β(nφ_g − kφ_f) > 0`. The near-extremal regime
  is handled by the strict `β = ∞` inequality, which keeps the `g`-side
  threshold a fixed distance `η` below `log g₁`; this is stated explicitly as
  the one place the endpoint is used.
- **Machinery moved to the back** — Weihrauch, Strassen, monotones, resource
  theories, thermomajorization now live in "Relations to other work".
- **Remark: power sums are exactly the monotones** — classification via
  singletons, with the honest caveat that `β = ∞` is not a homomorphism
  (`max` is not additive) and enters only as a limit.
- **Physics** — four sections plus a summary table of ten pairs: exact
  homotheties (`c = (g_A/v_A³)/(g_B/v_B³)`, giving `e⁺e⁻` `7/4`, neutrinos
  `21/8`, Debye phonons `3·10¹⁴`), and finite-temperature contact (SM plasma,
  spin-1/2 vs spin-1, qubits, crystal fields, Ising ring). Degenerate cases and
  the temperature-window repair stated honestly.
- **Gauge** — `C` is scale-invariant but not shift-invariant; complexity is
  carried by the pair (curve, marked origin), and the operational model is what
  supplies the origin. Stated in the introduction and in a remark.
- **Convex geometry** — `R_{a⊗b} = R_a ⊕ R_b` (as a remark: immediate from
  Legendre duality), verified numerically on mixed products to `~10⁻⁴`
  (discretization error). Area bound `C(g→f) ≤ √(A_g/A_f)` as a proposition.
- **Negative temperatures** — remark on why the curve stops at `β = 0`:
  `Z(β)` is order-preserving only for `β ≥ 0` (`{1} ⪯ {2}` but
  `Z=1 > Z=0.707` at `β = −½`), `log Z` can go negative, and the horizontal
  closure means an inverted population is never the better resource.
- **Sparaciari–Oppenheim–Fritz + thermomajorization** cited with the novelty
  boundary drawn explicitly.
- **Notation** `C(g|f) → C(g→f)`, 342 replacements across 16 files including
  the generators; `k_max` arrow direction fixed.
- **Housekeeping** — date, email, no affiliation (decided: none), code
  availability section, `paper/README.md` 99 → 69, matplotlib declared in
  `[analysis]`, figures regenerated with fixed legends and `Target` titles.

---

## Open

### 1. Certify the three-cycle --- DONE

Appendix C now proves it.  Two of the three universal lower bounds turned out
to be elementary: `5^b + 3^b >= (6^b+1)^(log5/log6)` from concavity of `x^a`
(the needed `1 - log(5/3)/log6 > 0` is just `6 > 5/3`), and
`6^b + 1 >= (3^b+2)^(log2/log3)` from `c log(t+2) <= log(t+1)`, which holds
with equality at `t = 1`.  The third, `3^b + 2 > (5^b+3^b)^(log2/log3)`,
follows from convexity of `x^(1/c)` plus a monotone split at `s = 12`.

So four of the six rates are exact --- `log2/log3` twice, `log3/log6`,
`log5/log6` --- and only five finite comparisons remain, all with margin at
least `3e-4`: `A(12)`, `B(12)`, two products of logarithms, and the single
evaluation at `beta = 9/20`.  `analysis/certify_cycle.py` rechecks them in
mpmath interval arithmetic.

Still open here: Appendix B's cluster claims elsewhere lean on "tolerance
`2e-6`" without saying what it guarantees, and the length-three cycle
`{6,3,3} -> {7,2,1} -> {6,5,1}` (margins `4-16e-3`) is not certified.

### 2. Appendix B --- DONE

Both longtables are gone.  The 69-signature table cost one page and the
exception table six; the paper went from 30 pages to 23.  The prefix survives
as a compact run of text, the exceptions as a generated summary (22 of the 69
signatures are involved, 112 pairs, every offending `x` of length two, and the
counts thin out as its second entry grows: 7, 5, 4, 3, 2, 1), and both are
written in full to `paper/anc/first_69_signatures.csv` and
`paper/anc/order_exceptions.csv` by the same generator.  The tie-breaking
caveat is now stated up front.

### 3. arXiv metadata

- [x] MSC 2020 in the paper: 94A17 primary; 05A17, 06F25, 16Y60, 82B03
      secondary. 05A17 / 06F25 / 16Y60 verified verbatim against the official
      list; 94A17 and 82B03 from memory — confirm on msc2020.org before
      submitting.
- [x] Keywords in the paper: exchange rate; asymptotic conversion; fiber
      signature; partition function; Gibbs curve; energy--entropy region;
      resource theory.
- [ ] Categories — suggest **cs.IT (= math.IT)** primary, **math.CO** cross-list;
      consider **quant-ph** given the resource-theory audience.
- [ ] License at submission (CC BY 4.0, consistent with the repo LICENSE).
- [ ] arXiv IDs in the bibliography — currently only `cao-et-al-2024` has one.
- [ ] Acknowledgements and a funding / competing-interests statement.
- [ ] Ship `analysis/*.csv` (and the B2 certificates) under `anc/`.

### 4. Open decision: put the regauging table in the paper?

`C` is invariant under `ε → λε` but not under `ε → ε + c`. Numerically, for the
two crystal-field schemes: `0.849525` (top gauge) → `0.883079` (origin `+0.5`)
→ `0.928249` (origin `+2`), and `C → 1` as the origin recedes, since
`Ψ → Ψ + βc`. And the three-cycle itself:

| gauge | `f₁ vs f₂` | `f₂ vs f₃` | `f₃ vs f₁` | cycle? |
|---|---|---|---|---|
| paper's (`ε_i = −log a_i`) | `f₁ ≺ f₂` | `f₂ ≺ f₃` | `f₃ ≺ f₁` | **yes** |
| highest level at 0 | `f₁ ≺ f₂` | `f₂ ≺ f₃` | `f₁ ≺ f₃` | no |
| `E(T=∞) = 0` | `f₁ ≺ f₂` | `f₃ ≺ f₂` | `f₁ ≺ f₃` | no |

No error: for signatures the origin is fixed by the model, so the cycle is a
true theorem about finite maps. The scoping is now stated in the text. What
remains is a choice — include the table (honest, pre-empts a referee, and turns
the point into an argument for the paper's two-part structure) or leave the
prose statement as is. Frame it as "the operational model supplies the origin",
never as "the cycle collapses".

### 5. Content that would raise the level further

1. **Reversibility --- STATED IN THE PAPER, CONJECTURE CORRECTED.**
   `C(a->b)C(b->a) <= 1` with equality iff `R = log Z_a/log Z_b` is constant,
   iff `Z_a = Z_b^c`, iff the regions are homothetic.

   The conjecture as first written here --- that the only reversible pairs
   share a `⊗`-root --- is **false**. An exhaustive screen of the 466,095
   pairs in `S_16` found eight reversible pairs, six of them with no common
   root, and all six are *flat* signatures (all entries equal). A flat
   `{k}^n` has `log Z` affine in `beta`, so its region is the rectangle
   `[-log k, 0] x [0, log n]`, and two rectangles are homothetic exactly when
   `log n/log k = log m/log l`. Hence `{2,2}`, `{3,3,3}`, `{4,4,4,4}`,
   `{5,5,5,5,5}` are mutually reversible, as are `{4,4}` and `{9,9,9}`. A flat
   is never reversible against a non-flat, whose `log Z` is strictly convex.

   Corrected conjecture, now in the conclusion without proof: **among non-flat
   signatures the only reversible pairs are `a = d^p`, `b = d^q`.** Verified
   over `S_16`.

   Toward a proof: substituting `x_j = p_j^beta` for the primes dividing the
   entries makes `Z_a` a polynomial with non-negative integer coefficients,
   and `Z_a^q = Z_b^p` in a UFD gives `Z_a = D^{p'}`, `Z_b = D^{q'}`. The gap
   is showing `D` is itself a signature --- and the flat family is exactly
   where that step fails, so it is a real obstruction, not a technicality.
   Irrational `c` is ruled out separately: `Z_b^c` is then not a finite
   exponential sum.

2. **Which endpoint minimizes.** Three examples are given (`β* = 0`, `∞`,
   interior) with no criterion. The endpoint derivative conditions give a
   checkable one, and it explains why four of six rates in Table 1 are endpoint
   values.
3. **Tangency as a proposition.** The equal-temperature / equal-slope contact
   result is a genuine little theorem still sitting in inline prose.
4. **Finite-size behaviour — already computed, not in the paper.**
   `analysis/README.md` proves that for `{2,2}` over `{3,1}` the exact points
   `(n, k_max(n))` are covered by 12 disjoint hyperbolas all using the single
   convergent `e = 5/8`, and that 11 is impossible. Two figures
   (`images/kmax_*.pdf`) already exist and are unused. A theorem like
   `k_max(n) = ⌊n log2/log3⌋ + O(1)` would add concrete mathematics and justify
   the `k_max` machinery.
5. **Rényi in one display.** With `p_i = a_i/N`, `N = Σ a_i`:
   `log Z_a(β) = β log N + (1−β) H_β(p)`. Turns the Jensen analogy into a change
   of variables.
6. **Easy structural facts.** `C(a→a) = 1`; `C(a→b)C(b→a) ≤ 1` (used without
   proof); degenerate values collected in one remark.
7. **Algorithm and precision.** One short section: how `R(β)` is minimized, is
   it unimodal (if not, how is that handled — this is the soundness question
   behind every table), what precision is certified, what it costs.
8. **Mixed area.** `A_{a⊗b} = A_a + 2V(R_a,R_b) + A_b`. `V` is a joint invariant
   of two signatures the theory does not otherwise see. Noted in the paper as
   not pursued.

### 6. Presentation nits

- Propositions typeset upright: all environments are declared under
  `\theoremstyle{definition}`. Plain style (italic) is conventional.
- Hard-coded "Definition~1" / "Definition~3" cross-references — correct today,
  silently wrong the moment a definition is inserted. Add `cleveref`.
- Figure 2: the tangency is the point of the section and is invisible at print
  size. A zoomed inset at each contact point would fix it.
- Figure 3 uses a blue/purple pair — check grayscale and colour-blind legibility.
- `\sig{}` renders as `{5,3}`, colliding with set-builder braces used two lines
  away (`S_B = {a : …}`). `⟨5,3⟩` would separate them.
- Subtitle "A Simplified Thermodynamic Complexity" undersells and points at the
  weakest material.
- Appendix A's `\boxed` formula — boxes are unusual in math papers.
- Clipped prose: "Its convention is reciprocal to ours after reversing the
  conversion" (§ Relations to other work) still needs a clause of explanation.
  The `{2,2}`/`{3,1}` one is fixed.
- The §6 cycle and the §7 `→ˢ` cycle traverse the same three signatures in
  opposite directions (correctly — `≺` and `→ˢ` are reversed). Half a sentence
  would stop it reading as an inconsistency.
