# What the paper claims, and what is only evidence

`main.tex` now carries the title *One Equation, Many Posets*, so the naming
question is settled. What follows is the part of the earlier planning note that
survives it: the thesis, an outline for the results sections, and a split that
is worth enforcing as the draft grows.

## Three tiers

A dataset is not a result. "We computed the five-, seven-, fifty-, twenty-six-
and nineteen-class diagrams and the 110-class preorder over \(\mathbb F_8\)" is
a list of files; no referee counts it and no reader is persuaded by it. Sort the
material:

* **Results** --- what the paper claims. The endpoint-regime theorem and the
  cycle lemma; the explicit three-cycle of quadratic morphisms
  \(\mathbf P^2\to\mathbf P^2\), which shows no numerical complexity on those
  morphisms is compatible with asymptotic conversion; the endpoint proposition
  \(C(L\to f)=\log q/\log\max_cN_c\) and the universal bottleneck
  \(\beta_*=\sqrt2-1\); flatness as a permutation-polynomial condition; and the
  minimal five-signature failure of negative type.
* **Evidence** --- rate matrices, cycle searches, saturation checks. Tables and
  ancillary CSVs, each attached to the statement it witnesses. A specific
  diagram or matrix earns its place only as the *witness of a general
  statement*: the \(\mathbb F_q\) matrix because eight of its twelve entries are
  exact ratios of logs of point counts; the fifty-class diagram because it is
  the smallest family here that cycles; the \(\mathbb F_8\) search because it
  shows \(\varphi\)-violations occur without producing a cycle. Detach the
  statement and each is a table again.
* **Exposition** --- the poset gallery and the new figures. Classical content,
  drawn better than the sources draw it. That is a legitimate contribution to an
  expository paper, provided it is labelled exposition and not smuggled in as a
  result.

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
