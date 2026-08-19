# Session brief O — how far is `C_aff` from `C_sig`?

**Read first:** `paper_finite_fields_maps/main.tex` around lines 295–390 (affine
implementation, `k_{g→f}(r)`, `C_aff`), the remark "No silent identification with
the affine rate" (~line 764), `research/synergy/FINDINGS.md`, and
`research/flux_arithmetic/FINDINGS.md` (the orbit counts, which bound how lossy
the signature is).

## The gap nobody has measured

Everything computed in this project is the **signature rate**
`C = C_sig(g→f) = inf_β log Z_g/log Z_f`. The **operational** object is the
affine rate `C_aff(g→f) = lim_r k_{g→f}(r)/r` with
`k_{g→f}(r) = max{k : f^{×k} ⪯_aff g^{×r}}`.

Affine implementation forces signature implementation, so

```
C_aff ≤ C_sig     always.
```

**Nowhere in this project is the gap measured, and every arithmetic claim is
about `C_sig`.** Brief A's objection ("information flows into the matrix") and
brief E's orbit counts (1744 pencils onto 296 signatures at `q = 11`; every
multiply-realised signature comes from pencils with different isogeny data) both
say the signature is lossy. The affine rate is not.

## What to establish

1. **Compute `k_aff(r)` exactly for small `r`** on the `F_3` quadratic classes
   (4 signatures) and cubic classes (5), where the whole space of affine
   processors is finite and enumerable. `f^{×k} ⪯_aff g^{×r}` is a search over
   affine `a : (F_3^2)^k → (F_3^2)^r` and `b : F_3^r → F_3^k`. Report
   `k(1), k(2), k(3), …` and the Fekete ratios `k(r)/r`, with the
   `C_sig` value alongside. **Even `r ≤ 3` would be the first data on `C_aff`
   this project has.**
2. **Is the gap ever strictly positive?** A single pair with
   `C_aff < C_sig` certified is a result: it shows the signature rate is a strict
   over-estimate and that the paper's separation of the two is necessary rather
   than cautious.
3. **Is `C_aff` strictly superadditive?** `research/synergy/FINDINGS.md` proves
   `C_aff(a⊗b→c) ≥ C_aff(a→c) + C_aff(b→c)` by block-diagonal composition, and
   leaves strictness open. **This is the interesting half**, because `C_aff` has
   no variational formula, so strictness cannot be a contact-temperature artefact
   — it would mean affine processors achieve jointly what they cannot separately,
   using the mixing of copies the paper explicitly permits. Search for it at
   small `r` on the `F_3` classes.
4. **Do the cycles survive?** Every certified cycle in this project is a cycle of
   `C_sig`. If `C_aff` orders the same triples differently, the arithmetic cycle
   results are about the signature shadow and not about the maps — which is
   exactly brief C's crux one level down. Check the `F_3` tensor 3-cycles
   (`analysis/tensor_cycles_f3.py`) against `C_aff` at whatever `r` is reachable.

## Why this is foundational rather than incremental

The programme's headline is that exchange rates give a complexity theory of maps.
`C_sig` is a computable proxy; `C_aff` is the definition. If they differ on the
smallest examples, every arithmetic statement needs the word "signature" in it —
which the paper currently says in a remark and nowhere else. If they agree, that
agreement is itself the theorem that licenses the whole computational programme.

## Traps

`k_aff(r)` grows fast; use symmetry (the affine group acts, so enumerate orbits,
not maps) and prune with the signature test, which is a cheap necessary
condition. Fekete gives `C_aff = sup_r k(r)/r`, so **any computed `k(r)/r` is a
lower bound on `C_aff` and can be reported as such with certainty** — an exact
lower bound plus the `C_sig` upper bound already brackets it.

## Deliverable

`research/affine_rate/FINDINGS.md`, house style. **If the harness blocks the
write, paste the full body in your final message.**
