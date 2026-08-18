# Session brief I — prove `C_4` is unrealisable

**Repo:** `fn_complexity`, work on a branch.

**Read first:** `research/realizability/FINDINGS.md` — all of §1 (the structure
theorem; take it as given, it is proved and verified), then §4.1, §4.2, §4.3.

## The question

Brief G established that the exchange metric is a **line metric plus a bounded
defect**, and then found the one place the framework looks genuinely obstructed:
of nine target metrics, only the 4-point equilateral metric is realisable. `C_4`
— the 4-cycle graph metric, `d = 1` on the sides, `2` on the diagonals — stalls
at distortion `1.255692` under a 24-parameter search at three independent fiber
counts `r = 3, 4, 6`.

**That is a stalled search, not a theorem. Turn it into one.**

> **Prove: no `n = 4` signatures realise `s·C_4` for any `s > 0`.**

A proof converts §4.3 from "we looked and failed" into *the* obstruction theorem
of the programme — the statement that says what the framework cannot do, and
therefore what it is. It is the single most valuable open item in the project.

## Why it should be provable, and the reduction to make first

`C_4` is `ℓ_1` on the 2-cube, isometric to `ℓ_∞` in the plane, so it needs
**exactly two independent directions**. The framework has one — brief G proves
every exchange metric is within `2 log 2` of the line metric `ℓ_ab = |σ_a − σ_b|`.
So the whole of `C_4`'s two-dimensionality has to be carried by the defect, and
the defect is bounded. That is the shape of the argument; the work is making it
exact.

**Reduction.** By §1.3, `d(a,b) = |σ_a − σ_b| + P + Q` with `P, Q ≥ 0`. `C_4` is
not a line metric, so at least part of it must come from `P + Q`, and
`P + Q ≤ 2 log(1 + e^{−|Δσ|})`. Two regimes, and the proof needs both:

* **Large `s`.** The non-line part of `s·C_4` grows linearly in `s` while the
  defect budget is capped at `2 log 2`. Make "non-line part" precise — the
  `ℓ_∞` distance from `s·C_4` to the cone of line metrics on 4 points, which is
  a small explicit LP — and you get an upper bound on `s` immediately. **Do this
  first; it is half a page and it kills the large-scale regime outright.**
* **Small `s`.** Here `σ` may be nearly constant and the metric is essentially
  pure defect. This is the hard half and where the search stalled (the realised
  scale was `0.148`). Specialise §1.2: with all `σ_a` equal, every `w_a` is a
  unimodal bump peaked at the **same** point `s = σ`, of height `h_a ∈ (0, log 2]`,
  1-Lipschitz, vanishing at `±∞`. Then `φ = w_b − w_a` vanishes at both ends, so

  ```
  d(a,b) = osc(w_b − w_a),   with both bumps co-peaked
  ```

  and the question becomes: **which 4-point metrics are oscillations of
  differences of co-peaked unimodal bumps?** That is a clean one-dimensional
  function-space problem with no signatures in it.

**A proved inequality you already have** (§Open of brief G): when the `σ` agree,
`|h_a − h_b| ≤ d(a,b) ≤ h_a + h_b`. Note this alone is *not* enough — `h ≡ s`
satisfies it for `C_4` — so it is a starting point and not the argument. The
missing ingredient is presumably that co-peaked unimodal bumps cannot produce
two *independent* pairs of far-apart points, which is exactly the
two-dimensionality obstruction. Look for a statement of the form "the four
numbers `osc(w_b − w_a)` satisfy a linear relation forced by co-peakedness".

## Also worth settling, same machinery

Both are listed as open in brief G and both are one-parameter-family questions
that the same analysis should reach:

* **Is `ε ≤ log(1+e^{−Δ})` rather than `2 log(1+e^{−Δ})`?** Every search returns
  `0.9889·log 2` at `Δ = 0`, *always with one of `P, Q` exactly zero*. That "one
  of them is always zero" is the clue — prove it and the metric bound halves to
  `d ≤ ℓ + log 2`.
* **Is the maximum triangle curl `(log 2)/2` rather than `3(log 2)/2`?** Three
  independent parametrisations reach `0.9889·(log 2)/2` and no further. A proof
  sharpens brief G §3.4 from *geometric-mean* to *total* asymmetry at most a
  factor 2 around any cycle.
* **Convexity of `U` in `s`.** Verified to `3.6·10⁻⁹` but not proved. A proof
  would turn §4.1's observed "never more than four interior extrema" into a
  theorem bounding them by the number of distinct atom-scales — which may itself
  be the cleanest route to `C_4`.

## If `C_4` turns out to be realisable

Then the search was simply not good enough, and you must produce the witness with
a certified margin (`≥ 1e−6`, verified at 40 digits). That would be an equally
important result — it would mean §4.3's obstruction is illusory and the
framework is metrically universal after all. **Do not assume the negative.**

## Traps

* `exchange_rate` is good to `~1e−13`; treat differences below `1e−10` as ties.
  A "witness" with margin `1e−12` is nothing.
* Work in an orthonormal basis of `{Σx = 0}` for anything eigenvalue-based —
  `−½JDJ` has the constant vector in its kernel and gives no search gradient.
* Windowed `β` computations need `β ~ 10³`; brief G's `common.py` already
  handles the extrema by a certified Lipschitz bracket — reuse it rather than
  gridding.
* Exclude `(1,)` and all-ones signatures; they are degenerate.
* Real-valued atoms are legitimate for a structural question and much easier
  than integers. If you use them, say so, and check whether a witness rounds.

## Success criterion

A proof that `s·C_4` is unrealisable for every `s > 0`, or a certified witness
realising it. Partial credit, in order: the large-`s` bound alone; the co-peaked
reduction stated and solved for `n = 3`; either of the two sharp-constant
questions settled.

## Deliverable

`research/realizability/OBSTRUCTION.md` in the house style (summary first,
sections marked proved/computed, corrections, open, files table), scripts
alongside. Build on `research/realizability/common.py`, `g2_metrics.py`,
`structure.py`, `extremes.py`.
