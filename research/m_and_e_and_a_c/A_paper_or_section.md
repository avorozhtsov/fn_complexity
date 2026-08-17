# A — paper or section? The decision, and where the text went

Answer to session brief `A_is_this_a_paper.md`.

## Verdict: a section, not a paper

Written as `\section{Fibrations of the affine plane over F_q}`
(`\label{sec:frobenius}`) in `paper_finite_fields_maps/main.tex`, placed
immediately before `\section{Prime modes and the zeta partition function}`.
Compiles clean: 42 pages, no overfull boxes, no undefined references.

### Why not a paper

The objection in the brief is fatal to the standalone version and only
survivable inside the framework paper. Every arithmetic statement here is a
change of variables on a statistic that arithmetic already studies:

* `C(L→f) = log q / log max_c N_c` is a strictly monotone function of the
  largest fiber;
* `1 − C(f→L)` is, to leading order, linear in `m₂`;
* the flatness theorem is two lines from the classical character-sum identity
  `Σ a_c² = q(ν(P) − q)`, and permutation polynomials are a well-worked
  classical subject;
* the `USp(2g)` edge exponent `d/2` is a general fact about compact groups.

The rate is a function of the signature and the signature is the list `{N_c}`,
so nothing the rate reports could fail to be computable from `{N_c}` directly.
A number-theory referee sees this on page one, and a standalone paper has
nothing left to be *about*: its subject would have to be the resource theory,
which is not in it.

Inside `One Equation, Many Posets` the same material has a job. That paper
evaluates the signature observer on `F_3` toy classes in §7 and then jumps to
the Euler product in §11; the weakest joint in the paper is exactly that jump.
A section evaluating the same observer on plane-curve fibrations is the missing
rung: arithmetically serious, still elementary, and it ends on
`Z_f(k) = #(X ×_Y ⋯ ×_Y X)`, which hands off directly to the zeta section.

### How the objection is met in the text

Not by writing around it. Subsection `sub:one-way`, "The direction of
information flow", is signposted in the section's opening paragraph ("a reader
who wants the objection before the results should read it first") and states
the reduction in `Remark~\ref{rem:reduction}` in its strongest form, including
the sentence that information flows into the matrix and not out of it, and the
explicit disclaimer that no result in the section is a new theorem about curves.

Four things are then claimed to survive, and nothing more:

1. **The compression is selected, not chosen.** Of all ways to compress
   `{N_c}`, the resource-theoretic comparison against a linear map returns
   exactly the extreme trace and the second moment. A consistency check on the
   framework, not a discovery about curves — stated that way.
2. **`β* = √2 − 1` is a fact about the comparison** and carries no arithmetic;
   its content is that the reverse comparison is decided at one fixed finite
   temperature, and that `Σ a_c = 0` is what puts it there.
3. **The invisibility statements are genuine**, because they are statements
   about the comparison rather than about `f`: the smallest fiber needs `β < 0`
   and is therefore unreachable, and symmetry type is unreachable because
   `Σ a_c = 0` is exactly the Katz–Sarnak separating statistic.
4. **The one escape route is open, and is flagged as open** — a cycle among
   arithmetic families, which would show the comparison is not a function of any
   scalar. Added as `Question~\ref{q:arith-cycle}` with the reason the present
   regime resists it (endpoint-attained pairs are provably acyclic under
   `φ = log(#fibers)·log(max fiber)`, and the interior correction is only
   `O(1/(q log q))`). This is brief B's experiment, stated as the paper's
   question rather than as a result.

The novelty ledger in `sec:novelty` gained a row saying the same thing in the
paper's own accounting format, and the abstract now says plainly that these are
translations of classical statistics.

## What went in

| item | status in the text |
|---|---|
| `C(L→f) = log q/log max N_c`, at `β = ∞` | Proposition, 3-line proof |
| `1 − C(f→L) = (3−2√2)m₂/(2q log q) + O(q^{−3/2})`, `β* = √2−1` | Proposition, derivation |
| moment ladder, damping `−(β*−k)/((k+1)√q)` | displayed, with the residual table |
| `2g` limit, convergence exponent `−2/dim USp(2g)` | Remark, marked numerical, with the honest "unreachable for `g ≥ 3`" |
| flatness ⟺ permutation polynomial | Theorem, with monomial / Dickson / superelliptic instances |
| supersingular-prime entry | Example, with the Dirichlet density argument for "not a congruence" |
| smallest fiber invisible | Proposition + the `q = 211` matched pair (`N_min` 173 vs 167) |
| symmetry type invisible | Proposition, with the Sato–Tate fluctuation escape |
| monodromy rank `m₂ = ν(P)/q − 1` | table at `q = 601` |
| `Z_f(k)` = fiber-power point count | own subsection, bridge to the zeta section |
| geometric irreducibility hypothesis | Remark, with the `xy` counterexample and the `x²+y²` control |

Left out as too fine for a section of this paper: the full genus-scaling tables
(T2.1 §2), the eight-probe inversion regression (T2.2 §3), the coincidence-class
census (T2.3 §6.3). They stay in `FINDINGS.md` and the per-thread notes.

## Reproduce

```
python analysis/frobenius_bottleneck.py
python analysis/frobenius_exchange_rates.py
```

Re-run at the time of writing; `frobenius_bottleneck.py` reproduces `β*` to
`6.9e−10`, the constant to `2.5e−13`, and reports 0 disagreements between
"permutation polynomial" and "flat".
