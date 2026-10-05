# Paper

* `paper.tex` — conference paper draft (article class; swap in the SIAM template when the venue is
  fixed). Open author items are marked `\todo{...}` and listed in `docs/STATUS_<date>.md`.
* `references.bib` — 33 entries, each verified against Crossref, arXiv, zbMATH or the publisher
  (rebuilt 2026-10-01; the August 2026 version had fabricated entries).
* `figures/` — copies of the figures used, from `figures/`.

Build: `pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper` (or `make`).
