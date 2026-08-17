# Session brief A — is the arithmetic material a paper, or a section?

**Repo:** `/Users/artemvorozhtsov/projects/fn_complexity` (work on a branch; commit
or stash any draft first — uncommitted work is invisible to a new worktree).

**Read first:** `research/m_and_e_and_a_c/FINDINGS.md` (whole file, especially
the Notation section and tracks T2.1, T2.2, T2.3), and
`paper_finite_fields_maps/main.tex` sections "The fiber-signature observer"
(~line 729) and "The Weil matrix E and the rate matrix M" (~line 1147).

## The question

An earlier session established, for `f : A² → A¹` over `F_q` with geometrically
irreducible fibers, `N_c = #f⁻¹(c)`, `a_c = q − N_c`, `m₂ = q⁻²Σa_c²`, and `L`
the flat signature `(q,…,q)`:

* `C(L→f) = log q / log(max_c N_c)`, always attained at `β = ∞` (3-line proof);
* `C(f→L) = 1 − (3−2√2)·m₂/(2q log q) + O(q^{−3/2})`, with the bottleneck at
  `β* = √2 − 1`, independent of family, genus and `q`;
* `Z_f(k) = #(X ×_Y ⋯ ×_Y X)`, the k-fold fiber power, so the partition function
  at integer `β` is a point count;
* a signature is flat iff `P` is a permutation polynomial (subsumes the
  `q ≡ 2 mod 3` result and reaches supersingular primes).

Decide, honestly, whether this is a standalone paper or a section of the
existing one, and then write whichever it is.

## The objection that must be met head-on

`C(L→f)` is a strictly monotone function of `max_c N_c`, and `1 − C(f→L)` is to
leading order an affine function of `m₂`. So "the exchange rate encodes the
extreme Frobenius trace and the second moment" reduces, on inspection, to "a
monotone function of X is a function of X". The rate is a function of the
signature, and the signature is just the list `{N_c}`, so nothing the rate says
about arithmetic could fail to be computable from `{N_c}` directly. Information
flows from the arithmetic into the matrix, not out of it.

A referee in number theory will see this on page one. Do not write around it;
write through it. Either find the place where the matrix carries information its
scalar shadows do not (see brief C), or state the reduction openly and let the
two propositions and `β*` stand on their own modest merits.

## What is defensible either way

The two propositions; `β* = √2 − 1` with `3 − 2√2`; the fiber-power identity;
and the observation that of all ways to compress `{N_c}`, the resource-theoretic
comparison selects exactly the extreme value and the second moment — canonical
rather than arbitrary. Also worth stating: what the rate provably *cannot* see
(the smallest fiber, since that needs `β < 0`) and why symmetry type is invisible
by construction (`Σ a_c = 0` is exactly the Katz–Sarnak separating statistic).

## Reproduce

`analysis/frobenius_bottleneck.py`, `analysis/frobenius_exchange_rates.py`.
