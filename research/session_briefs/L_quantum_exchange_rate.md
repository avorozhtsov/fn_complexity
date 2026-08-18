# Session brief L — the quantum exchange rate

**Read first:** `research/realizability/OBSTRUCTION.md` §on the Hilbert
projective metric (the structural identification this brief builds on),
`research/realizability/FINDINGS.md` §1, and
`research/session_briefs/D_quantum_mechanics_of_the_exchange_matrix.md` (what is
already prior art on the quantum side — do not re-derive it).

## The observation that makes this a real question

A signature `a = (a_1,…,a_r)` is a positive diagonal operator `A = diag(a)`, and

```
Z_a(β) = Σ_i a_i^β = Tr A^β
```

so `log Z_a(β)` is a Rényi-type quantity and the exchange rate is

```
C(a→b) = inf_β  log Tr A^β / log Tr B^β .
```

**Written this way the definition never used commutativity.** For positive
operators `A, B` on a finite-dimensional Hilbert space the same formula is
defined verbatim, and reduces to the classical one when `[A,B] = 0`. That is the
quantum exchange rate, and it is a definition rather than an analogy.

## What to establish

1. **Does the structure theorem survive?** Brief G proves `u_a = log log Tr A^β`
   is a kink plus a unimodal bump of height ≤ `log 2`, and brief I proves `d` is
   the Hilbert projective metric on the cone of such profiles. Both proofs use
   only that `F = log Tr A^β` is convex, increasing, and sandwiched between
   `max(log r, β log‖A‖)` and `log r + β log‖A‖`. **Check each step for a
   non-commuting pair.** `Tr A^β` depends only on the spectrum of `A`, so at
   first sight nothing changes — say clearly whether the quantum object is
   therefore *only a function of the two spectra*, because if so the honest
   conclusion is that this quantisation is empty and the interesting definition
   is a different one (see 3).
2. **Cycles.** Do non-commuting triples cycle? Is there a cycle with **no
   classical shadow** — i.e. one that disappears when each operator is replaced
   by its spectrum? If `Tr A^β` sees only spectra, there cannot be, and that is
   the decisive negative to report early.
3. **The definition that is not spectral.** The natural repair is the
   **sandwiched Rényi divergence**
   `D̃_α(A‖B) = (1/(α−1)) log Tr[(B^{(1−α)/2α} A B^{(1−α)/2α})^α]`, which does not
   factor through the spectra. Define `C̃(A→B)` by the same infimum over `α` and
   redo 1–2. **This is where any genuinely quantum phenomenon has to live.**
4. **The prior art that matters, and it is close.** Brandão–Horodecki–Ng–
   Oppenheim–Wehner, *second laws of quantum thermodynamics*: state conversion
   under thermal operations is governed by a **continuum of Rényi α-free
   energies, all of which must decrease** — exactly the β-indexed family of
   monotones this framework infimises over. Read it before claiming novelty. The
   honest question is whether the exchange rate is a known quantity there under
   another name, or a genuinely different functional of the same family (it
   takes a *ratio* of monotones and an infimum, which is not their construction).
   Say which.

## Success criterion

Either a quantum cycle with no classical shadow, with certification; or a proof
that `C` factors through the spectra and the sandwiched version is the only
non-trivial quantisation, with its structure theorem stated. A written placement
against BHNOW either way.

## Traps

`exchange_rate` tolerance is `1e-13`, tie threshold `1e-10`; β horizons must
reach `10³`. Do not build on `MᵀM` or `MᵀDM` (PSD for free, informative never —
brief D Part 2). Verify headline numbers at 40 digits.

## Deliverable

`research/quantum/FINDINGS.md`, house style (summary; sections marked
proved/computed; corrections; open; files table). **If the harness blocks the
write, paste the full body in your final message and say so.**
