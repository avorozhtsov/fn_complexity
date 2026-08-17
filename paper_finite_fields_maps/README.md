# Paper build

`main.tex` is the narrative draft *One Equation, Many Posets: Polynomial maps
over fields as resources*. The rendered reading copy is
`one_equation_many_posets.pdf`.

Build it with:

```bash
cd paper_finite_fields_maps
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The committed `figures/*.pdf` files are PDF conversions of the source SVGs in
`images/`. The SVGs remain the canonical diagram artifacts. For example:

```bash
rsvg-convert -f pdf \
  -o figures/quadratic-map-poset-q3.pdf \
  images/quadratic-map-poset-q3.svg
```

Regenerate the underlying finite-field and homogeneous SVG diagrams from the
project root with:

```bash
./cli/finite_field_map_poset_cli --no-titles
./cli/homogeneous_tensor_poset_cli --no-titles
./cli/cubic_map_poset_cli --no-titles
./cli/cubic_map_poset_cli --q 8 --case quadratic --no-titles
```

All three diagram CLIs accept `--titles` and `--no-titles`. Titles are on by
default for standalone images. The checked-in figures used by the paper are
generated with `--no-titles`; their captions carry the field, processor
convention, class count, and Hasse-cover count.

Regenerate the signature-exchange matrix and cycle-mean tables used in the
paper with:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
  python3 analysis/polynomial_map_exchange_examples.py
```

The generated CSVs are stored in `anc/`. They deliberately distinguish the
finite-map signature observer from the generally unknown operational affine
exchange rate.

Regenerate the zeta entropy-energy curve and its paper PDF with:

```bash
uv run --extra analysis python analysis/zeta_entropy_energy_curve.py --no-titles
rsvg-convert -f pdf \
  -o paper_finite_fields_maps/figures/zeta-entropy-energy.pdf \
  paper_finite_fields_maps/images/zeta-entropy-energy.svg

uv run --extra analysis python analysis/prime_mode_entropy_energy_curves.py --no-titles
rsvg-convert -f pdf \
  -o paper_finite_fields_maps/figures/prime-mode-entropy-energy-curves.pdf \
  paper_finite_fields_maps/images/prime-mode-entropy-energy-curves.svg

python analysis/xi_gibbs_curve.py
rsvg-convert -f pdf \
  -o paper_finite_fields_maps/figures/riemann-xi-gibbs-curve.pdf \
  paper_finite_fields_maps/images/riemann-xi-gibbs-curve.svg
```

The physical branch exists only for inverse temperature `beta > 1`. The
figure shows the full excluded region `beta <= 1`, including negative
temperatures, rather than treating analytic continuation as a Gibbs state.
The companion `riemann-xi-gibbs-curve.svg` is the same construction for the
*completed* zeta function, where that obstruction disappears: Riemann's kernel
makes `xi(1/2 + beta)` the Laplace transform of a positive even measure, so the
curve exists for every real `beta` and its mirror symmetry is the functional
equation. The two figures are a pair -- the first shows where the Gibbs reading
of `zeta` breaks, the second what repairs it.
The companion prime-mode figure compares the first six spectra
`k log(p)`, for `k = 0, 1, 2, ...`.

The finite-field source is named `quadratic-map-poset-F2.svg`, while the
historical 2-adic source remains `quadratic-map-poset-Q2.svg` (with the
unambiguous alias `quadratic-map-poset-q2-adic.svg`). This prevents one from
overwriting the other on case-insensitive filesystems.
