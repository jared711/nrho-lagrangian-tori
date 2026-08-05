# Makefile for LaTeX paper compilation

# Main document name (without .tex extension)
MAIN = paper

# LaTeX compiler
LATEX = pdflatex
BIBTEX = bibtex

# Flags
LATEX_FLAGS = -halt-on-error -interaction=nonstopmode

# PDF viewer (Linux)
VIEWER = xdg-open

# All auxiliary files to clean
AUX_FILES = *.aux *.log *.bbl *.blg *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz *.bcf *.run.xml

.PHONY: all pdf clean distclean view help

# Default target
all: pdf

# Build PDF with bibliography
pdf: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex references.bib
	@echo "First LaTeX pass..."
	$(LATEX) $(LATEX_FLAGS) $(MAIN).tex
	@echo "Running BibTeX..."
	$(BIBTEX) $(MAIN)
	@echo "Second LaTeX pass..."
	$(LATEX) $(LATEX_FLAGS) $(MAIN).tex
	@echo "Third LaTeX pass..."
	$(LATEX) $(LATEX_FLAGS) $(MAIN).tex
	@echo "Build complete: $(MAIN).pdf"

# Quick build (single pass, no bibliography)
quick:
	$(LATEX) $(LATEX_FLAGS) $(MAIN).tex

# View the PDF
view: $(MAIN).pdf
	$(VIEWER) $(MAIN).pdf &

# Clean auxiliary files but keep PDF
clean:
	@echo "Cleaning auxiliary files..."
	rm -f $(AUX_FILES)
	@echo "Clean complete."

# Clean everything including PDF
distclean: clean
	@echo "Removing PDF..."
	rm -f $(MAIN).pdf
	@echo "Distclean complete."

# Rebuild from scratch
rebuild: distclean all

# Help
help:
	@echo "LaTeX Paper Makefile"
	@echo "===================="
	@echo ""
	@echo "Targets:"
	@echo "  all (default)  - Build PDF with bibliography (full build)"
	@echo "  pdf           - Same as 'all'"
	@echo "  quick         - Quick build (single LaTeX pass, no bib)"
	@echo "  view          - Open PDF in viewer"
	@echo "  clean         - Remove auxiliary files, keep PDF"
	@echo "  distclean     - Remove all generated files including PDF"
	@echo "  rebuild       - Clean and rebuild from scratch"
	@echo "  help          - Show this help message"
	@echo ""
	@echo "Example usage:"
	@echo "  make          # Full build"
	@echo "  make view     # Build and view"
	@echo "  make clean    # Clean up"
