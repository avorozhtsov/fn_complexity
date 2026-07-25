# Paper build

The directory is self-contained for a pdfLaTeX-based arXiv submission.

Compile with TeX Live:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The three vector figures are regenerated from the project root with:

```bash
PYTHONPATH=src python analysis/plot_pair_homotheties.py
PYTHONPATH=src python analysis/plot_finite_temperature_pair.py
PYTHONPATH=src python analysis/plot_exchange_homotheties.py
```

After regeneration, copy their PDF versions from `images/` into
`paper/figures/` before compiling or submitting.
