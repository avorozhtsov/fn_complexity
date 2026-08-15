# Paper build

The directory is self-contained for a pdfLaTeX-based arXiv submission.

Compile with TeX Live:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
cp main.pdf exchange_rates_finite_map_signatures.pdf
```

`main.pdf` remains an ignored build artifact.  The descriptively named copy
`exchange_rates_finite_map_signatures.pdf` is committed with the source.

The four vector figures are regenerated from the project root with:

```bash
PYTHONPATH=src python analysis/plot_pair_homotheties.py
PYTHONPATH=src python analysis/plot_finite_temperature_pair.py
PYTHONPATH=src python analysis/plot_exchange_homotheties.py
python analysis/plot_relativistic_species_energy_entropy.py
```

The first three scripts write PDF versions to `images/`; copy those into
`paper/figures/` before compiling or submitting.  The relativistic-species
script writes its PDF directly to `paper/figures/`.

Recheck the numerical comparisons of Appendix C in interval arithmetic with:

```bash
python analysis/certify_cycle.py
```

Generate the ordered-signature and exception tables in Appendix B with:

```bash
PYTHONPATH=src python analysis/appendix_b_signatures.py
```

The larger NumPy/SciPy screening calculation that verifies the stable first
69 across the \(B=18\) and \(B=19\) candidate shells is run from the project
root with:

```bash
PYTHONPATH=src:analysis python analysis/stabilize_signature_order.py
```
