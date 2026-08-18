# Session brief N — is a processor a Birkhoff contraction?

**Read first:** `research/realizability/OBSTRUCTION.md` — the Hilbert projective
metric identification, which is the whole basis of this brief — then
`research/realizability/FINDINGS.md` §1.6 (Cartesian powers) and
`paper_finite_fields_maps/main.tex` around line 300 (affine implementation).

## The observation

Brief I proved:

> `d(a,b)` is exactly the **Hilbert projective metric** between `F_a = log Z_a`
> and `F_b` on the cone `C = {Φ convex, nondecreasing, Φ ≥ Λ_Φ·β}`, and `C` is
> the projective closure of the achievable set. Cartesian powers are the
> projective rescaling the Hilbert metric quotients out.

The Hilbert metric has one classical theorem attached to it, and it is exactly
about maps:

> **Birkhoff–Hopf.** A linear map `T` preserving a cone contracts its Hilbert
> metric with ratio `tanh(Δ(T)/4)`, where `Δ(T)` is the projective diameter of
> the image. In particular `Δ < ∞` ⟹ strict contraction.

This project's maps are the **processors** — the `h_out ∘ g ∘ h_in` of the
framework, and their affine restriction. The question writes itself.

## What to establish

1. **Which operations on resources act linearly on the cone `C`?** Cartesian
   product is `Φ ↦ Φ + Φ'` (an isometry, by §1.6). Composition with a processor
   is something else. Identify the induced map on `C` for: a fixed input
   processor `h_in`, a fixed output processor `h_out`, and a fixed resource
   tensored in. **Say honestly which of these are linear on the cone and which
   are not** — the theorem applies only to the linear ones and it is easy to
   wave at this step.
2. **Compute `Δ` and the contraction ratio** for whichever maps qualify, on the
   small classes where everything is enumerable (`F_3` quadratics and cubics,
   `research/synergy/synergy.py` builds them). A processor that strictly
   contracts `d` means *conversion loses distinguishability at a definite rate* —
   a quantitative second-law-flavoured statement the project does not have.
3. **The one thing that would be genuinely new.** Birkhoff contraction is the
   classical ancestor of the **contraction coefficients of quantum channels**
   (Hilbert-metric contraction for positive maps on the PSD cone; Reeb–Kastoryano
   and the relative-entropy contraction literature). If processors are Birkhoff
   contractions, the exchange framework and quantum channel contraction are the
   same theorem on two different cones. Establish or refute that, and place it
   against the literature — assume a referee knows Birkhoff, Bushell, and
   Reeb–Kastoryano.
4. **Does contraction explain the `log 2`?** Brief I's sharp defect
   `d ≤ ℓ + log 2` and brief G's `(log 2)/2` comparison bound are constants of
   the cone. Birkhoff-type arguments produce exactly such constants from the
   projective diameter of a sub-cone. **If `log 2` is the diameter of something
   natural, that is the clean proof the project is missing** — and it would also
   settle the open question whether the maximum triangle curl is `(log 2)/2`
   rather than `3(log 2)/2`.

## The likely failure mode, stated up front

Resource conversion here is *not* a linear map on `C` — the exchange rate is
defined by an inequality between profiles, not by applying an operator. If so,
Birkhoff does not apply and the Hilbert-metric identification is a description of
the geometry rather than a dynamical statement. **That is an acceptable and
publishable answer**; write it as a proposition about what the cone structure
does and does not buy, rather than leaving it implicit.

## Deliverable

`research/birkhoff/FINDINGS.md`, house style (summary; proved/computed marking;
corrections; open; files table). **If the harness blocks the write, paste the
full body in your final message.**
