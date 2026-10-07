# Numerical-range certificate paper

This directory contains the manuscript's canonical LaTeX source, bibliography and compiled PDF. The research abstract contains the research summary; the generative-AI disclosure appears in its dedicated end section.

From the repository root, run `python -B reproduce/build_paper.py` with `pdflatex` and `bibtex` available. The rebuilt PDF is written to `output/pdf/`. Exact requests, certificates and experiment outputs are under `examples/` and `results/`. Programs and an independently implemented ellipse-study oracle are under `experiments/`.

Reference identity/context evidence is recorded in `audit/citation-verification.json`. It is not independent verification of the cited theorems. The current manuscript is a technical note accompanying an alpha prototype, with comparative controls and explicit assurance scope.

Author: Robin Gounder, Vaionex Corporation; robin@gounder.com. Manuscript copyright is separate from the software MIT licence; see `LICENSE` in this directory.
