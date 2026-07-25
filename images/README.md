# Reproducing the figures

The `images/` directory intentionally contains only PDF figures with white
backgrounds. PDF is the most portable vector format for `pdflatex` and arXiv.

## Retained finite-rate plots

The convention is that `C(g | f)` counts the asymptotic number of copies of
`f` implemented per copy of `g`.

| Plot | PDF |
| --- | --- |
| Greedy family `c_n=1-d/n` for `C({3,1}|{2,2})`, `n<100` | `kmax_g-3-1_f-2-2_n-1-99_unit-cover.pdf` |
| First 14 branches `c_n=1/2+d/n` for `C({2,2}|{3,1})`, `n<100` | `kmax_g-2-2_f-3-1_n-1-99_first-14-half-branches.pdf` |
| Transformed ordinate `Y_n=1/(C-c_n)`, through `n=200` | `kmax_g-3-1_f-2-2_n-1-200_unit-cover_inverse-gap.pdf` |

The CLI writes SVG, so regeneration uses a temporary directory and converts
only the final result to PDF:

```bash
tmpdir=$(mktemp -d)

./cli/kmax_cli 3,1 2,2 --n-min 1 --n-max 99 \
  --unit-cover --quiet --output "$tmpdir/"
./cli/kmax_cli 2,2 3,1 --n-min 1 --n-max 99 \
  --half-branches 14 --quiet --output "$tmpdir/"
./cli/kmax_cli 3,1 2,2 --n-min 1 --n-max 200 \
  --unit-cover --inverse-limit-gap --quiet --output "$tmpdir/"

for svg in "$tmpdir"/*.svg; do
  name=$(basename "$svg" .svg)
  rsvg-convert --format pdf --output "images/$name.pdf" "$svg"
done
```

On macOS, install the converter with `brew install librsvg`.

## Gibbs and homothety figures

All plot labels use `T` for physical temperature; inverse temperature is
`beta=1/T`. Generate the white-background PDF figures with:

```bash
PYTHONPATH=src python analysis/plot_energy_entropy.py
PYTHONPATH=src python analysis/plot_pair_homotheties.py
PYTHONPATH=src python analysis/plot_finite_temperature_pair.py
PYTHONPATH=src python analysis/plot_exchange_homotheties.py
```

The solid lines are original Gibbs regions. Dotted lines—including endpoint
closures—belong to scaled homothetic regions; thin solid endpoint closures
belong to the original regions.

Copy the three paper figures and rebuild the paper:

```bash
cp images/exchange-homotheties_2-2_3-1.pdf paper/figures/
cp images/exchange-homotheties_6-5-2-1_6-4-3-2.pdf paper/figures/
cp images/exchange-homotheties_cycle-5-3_3-1-1_6-1.pdf paper/figures/
(cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex)
```

Verify retained PDFs with `pdfinfo images/*.pdf`.
