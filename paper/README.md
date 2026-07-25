# Paper build

The directory is self-contained for a pdfLaTeX-based arXiv submission.

Compile with TeX Live:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

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

Generate the ordered-signature and exception tables in Appendix B with:

```bash
PYTHONPATH=src python analysis/appendix_b_signatures.py
```
