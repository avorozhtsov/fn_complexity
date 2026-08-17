# Implementing maps by maps

**Working title.** *Implementing maps by maps: degeneration posets and
exchange rates over finite and local fields.*

**Subtitle.** *One resource calculus behind Witt, Hasse--Minkowski and Igusa.*

The paper is expository by intention. Its message is that a single elementary
calculus --- resources, processors, a degeneration order, and multiplicative
monotones --- reaches a string of classical classification results, and that
drawing them all in one language with one set of pictures is worth doing. The
shift from the earlier plan is only in emphasis: the classical material stays,
but the paper now carries one thing the classical theory does not have, so that
the reader is not asked to accept new notation for its own sake.

## The thesis

Three sentences, in order of decreasing familiarity.

1. The relation "\(R\) is implemented by \(Q\)" is classical. Restricted to
   forms it says \(R=Q\circ A\) for a linear map \(A\): the **representation of
   a form by a form**. The two compositions are classical too --- \(\oplus\) is
   the orthogonal sum, and \(\otimes\) of two scalar maps is a pair of forms, so
   the tensor closure of the theory is the theory of pencils of quadrics. Even
   the counting flavour is classical: level, Pythagoras number, \(u\)-invariant.
2. What the classical theory does not do is ask the **asymptotic** question:
   how many copies of \(R\) can \(n\) copies of \(Q\) build, in the limit. That
   limit exists, has a closed form, and produces a real number for every ordered
   pair of classes.
3. Over \(\mathbb F_q\) those numbers are **Weil numbers**, because fiber sizes
   are point counts; and they **strictly order classes that the degeneration
   poset leaves incomparable**. That is the paper's one new result, and it is
   what licences the whole exercise.

## Proposed outline

1. **Implementation.** Processors, the degeneration order, and the dictionary
   to representation of forms. State plainly what is classical.
2. **The posets.** `quadratic_map_posets.md` in full, then `real_complex_…`,
   `p_adic_…`, `q3_adic_…`. **Keep every figure.** The gallery is a deliberate
   part of the paper: the classical sources for this material carry no pictures
   at all, and a student meeting quadratic forms over \(\mathbb F_q\),
   \(\mathbb R\), \(\mathbb C\) and \(\mathbb Q_p\) side by side, drawn to one
   convention, learns something no textbook offers. What the invariants
   \(u(K)\) and \(|K^\times/(K^\times)^2|\) add is a *caption*: they say why the
   diagrams have the shapes they have, so the gallery is a theorem illustrated
   rather than a set of examples. Put the invariant in each caption and the
   local--global principle for field-rational representations at the end.
3. **Tensor closure.** `homogeneous_tensor_posets_q3.md`. This is not a side
   gallery: \(\otimes\) of two scalar-valued maps is vector-valued, so the
   semiring is not closed until pairs and triples of forms are in it. The
   \(\mathbb F_3\) computations are the first cases.
4. **Exchange rates.** `finite_field_exchange_matrix.md` --- the core section.
   The \(4\times4\) matrix, eight exact closed forms, the two constants
   \(\lambda,\kappa\), and the total order \(S\prec L\prec A\prec X\).
5. **Local fields.** `p_adic_zeta_exchange_currency.md` (Igusa zeta functions as
   the local monotones, the bad-prime correction \(t^2-2t+2\), the constant
   \(0.9397\ldots\)) and `p_adic_exchange_rate_attempts.md` (the sharp negative
   results) as one section on where the operational rate stops being computable.
6. **Cubic maps.** `cubic_map_posets_q3.md`, `cubic_map_posets_q8.md` --- an
   appendix or a data note. Include only if the page count allows; they do not
   serve the thesis.
7. **Outlook, two pages.** `riemann_hypothesis_exchange_matrices.md` and
   `exchange_positivity_and_weil.md`, reduced to the dictionary and the two
   named obstructions. Do **not** present the Weil restatement as a section of
   the paper's own results; it is a restatement, it says so, and inflating it
   costs the rest of the paper its credibility.

## Frame, and what to avoid

"One framework reproduces classical results" is a fine message for an arXiv
expository paper, and the illustrative material here is better than what the
classical sources carry. The only thing to avoid is the impression that the
framework is *only* a relabelling. Two sentences prevent it, and both are now
supported: the relation is classical but the asymptotic rate is not, and the
rate decides comparisons the classical order declines to make.

The honest gap, stated once and not repeated: no application of these rates
outside the framework has been found. Section 7 says where the natural
candidate --- Weil positivity --- does and does not connect.
