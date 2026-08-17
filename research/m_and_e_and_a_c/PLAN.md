# Research plan: M ↔ E, and what M knows about {a_c}

Two tracks, run in parallel, revised as results arrive.

## Conventions

Repo root is the worktree. In Python:

```python
import sys; sys.path.insert(0, "src")
from fn_complexity import exchange_rate, exchange_rate_result, gibbs_point
```

* `C(g→f) = inf_β log Z_g(β)/log Z_f(β)`, `Z_a(β) = Σ a_i^β`, β ∈ [0,∞].
* Exchange metric `d(a,b) = −log(C(a→b)·C(b→a)) ≥ 0`, a pseudometric.
* Isometry: `u_a = log log Z_a`, `d(a,b) = osc_β(u_a − u_b)`; sup-norm geometry.
* Endpoint index `φ(a) = log(#fibers)·log(max fiber)`; in the endpoint regime
  `a ≺ b ⟺ φ(a) < φ(b)`, so cycles need a φ-violating pair.
* Established: `d` is **not** of negative type (13-signature certificate in
  `analysis/exchange_negative_type_certificate.csv`).

## Track 1 — the M ↔ E connection

The thesis so far: M is an L^∞/tropical object (quasi-metric, no-arbitrage,
cycles), E is an L²/Hilbertian object (Gram, PSD ⟺ RH). Open questions:

1. **T1.1** Minimal negative-type violation; which Deza–Laurent inequality
   (triangle, pentagonal, hypermetric) fails first, and at what size.
2. **T1.2** The L^t interpolation: for which t is `exp(−t·d)` PSD? Locate the
   transition; is there a critical t per family?
3. **T1.3** ℓ₂-distortion of the exchange metric vs Bourgain's `O(log N)`.
4. **T1.4** *(centrepiece, run by the lead)* Each signature `a` is a genuine
   Weil test measure `Σ_i δ_{a_i}`, with Mellin transform `Z_a(s)`. So
   `E_N(a,b) = Σ_{n≤N} Z_a(ρ_n) conj(Z_b(ρ_n))` over the first N zeta zeros is
   a real Weil-shaped Gram matrix **on the same objects as M**. Compare the two
   geometries directly. *(Done: `T1_4_weil_pairing_on_signatures.md`.)*
5. **T1.5** Can multiplicative design make the two geometries agree?
   *(Done: `T1_5_multiplicative_design.md` — no; the invariance groups are
   transverse.)*

## Track 2 — what M knows about {a_c}

Established: `Z_f(k)` counts k-fold fiber powers; `Σ_c a_c = 0` exactly;
`1 − C(L→f) ≍ 2g/(√q log q)`; `C(y²−x³ → L) = 1 ⟺ q ≡ 2 mod 3`.

6. **T2.1** Genus scaling: test `2g/(√q log q)` for g = 2, 3.
7. **T2.2** Injectivity: do the rates determine the trace moments? Families
   with matched m₂ but different m₃.
8. **T2.3** Symmetry type (Katz–Sarnak) — does M distinguish unitary from
   symplectic from orthogonal families?
9. **T2.4** Arithmetic detection beyond q ≡ 2 mod 3: CM, supersingularity.

## Outputs

One markdown note per question in this directory, scripts alongside, plus
`FINDINGS.md` synthesising. Negative results are results; record them.
