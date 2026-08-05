# Paper: Computation of Lagrangian Tori Around NRHOs

This directory contains the LaTeX source for a research paper on computing invariant tori around Near-Rectilinear Halo Orbits.

## Files

- `paper.tex` - Main LaTeX document (SIAM format)
- `references.bib` - BibTeX bibliography
- `Makefile` - Build automation
- `figures/` - Directory for figure files (create as needed)

## Building the Paper

### Requirements

- LaTeX distribution (TeX Live, MiKTeX, etc.)
- SIAM LaTeX class files (`siamart220329.cls`)
- BibTeX or BibLaTeX

### Quick Build

```bash
make
```

This will:
1. Run pdflatex
2. Run bibtex
3. Run pdflatex twice more for cross-references
4. Open the PDF (on Linux with xdg-open)

### Individual Commands

```bash
make pdf      # Build PDF once
make clean    # Remove auxiliary files
make distclean # Remove all generated files including PDF
```

### Manual Build

```bash
pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex
```

## SIAM Document Class

The paper uses the SIAM article class `siamart220329`. If you don't have it installed:

1. Download from: https://www.siam.org/publications/journals/about-siam-journals/information-for-authors
2. Or use CTAN: https://ctan.org/pkg/siam
3. Place `siamart220329.cls` in this directory or your local texmf tree

Alternatively, you can change the document class to a standard class:
```latex
\documentclass[11pt,letterpaper]{article}
```

## Sections to Complete

The paper template includes TODO markers for sections that need completion:

- [ ] Section 4.1: Complete numerical results table
- [ ] Section 4.2: Add torus geometry analysis
- [ ] Section 4.3: Add comparison with Poincaré map iteration
- [ ] Add figures (convergence plots, torus visualizations, cross-sections)
- [ ] Fill in specific numerical values (NRHO period, rotation frequencies, etc.)
- [ ] Generate and include results from the corrected C++ implementation

## Figures

Place figure files in the `figures/` subdirectory. Suggested figures:

1. `nrho_orbit.pdf` - Visualization of the NRHO in configuration space
2. `torus_3d.pdf` - 3D rendering of the computed torus
3. `kam_convergence.pdf` - Plot of invariance error vs iteration
4. `grid_refinement.pdf` - Illustration of adaptive grid refinement
5. `poincare_section.pdf` - Poincaré section showing torus cross-section
6. `rotation_curve.pdf` - Rotation frequencies as function of parameters

## Citation Style

SIAM uses a numbered bibliography style. References are cited as [1], [2], etc.

## Submission

For SIAM conference submission:
- Check specific conference requirements (page limit, format, etc.)
- Ensure all figures are publication quality (vector format preferred)
- Verify math notation is clear and consistent
- Have coauthors review before submission

## Notes

- The paper is currently in "review" mode (line numbers, double spacing)
- To switch to final mode: change `\documentclass[review,...]` to `\documentclass[final,...]`
- AMS subject classifications: 37J40 (Hamiltonian systems), 70F07 (Three-body problems), 70K43 (Quasi-periodic motions), 37M05 (Computational methods)
