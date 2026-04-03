#!/usr/bin/env bash
# Prepare a self-contained ApJ submission directory in Paper/ApJ/
# All files are flat (no subdirectories) as required by ApJ submission.
# Run from the repository root: bash Paper/prepare_apj.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
APJ_DIR="$REPO_ROOT/Paper/ApJ"

# Clean and recreate (flat — no subdirectories)
rm -rf "$APJ_DIR"
mkdir -p "$APJ_DIR"

# Copy main manuscript and bibliography
cp "$REPO_ROOT/Paper/adaptive_gridding_for_snec.tex" "$APJ_DIR/"
cp "$REPO_ROOT/Paper/adaptive_gridding_for_snec.bib" "$APJ_DIR/"

# Copy generated tables (flat into ApJ/)
cp "$REPO_ROOT/Paper/generated/"*.tex "$APJ_DIR/"

# Copy figures (flat into ApJ/)
cp "$REPO_ROOT/Analysis/figures/adaptive_runtime_validation.png" "$APJ_DIR/"
cp "$REPO_ROOT/Analysis/figures/lc_comparison_sne.pdf" "$APJ_DIR/"

# Copy built PDF
cp "$REPO_ROOT/Paper/adaptive_gridding_for_snec.pdf" "$APJ_DIR/"

# Bundle AASTeX v7 class and bst (arXiv TeX Live may not have them yet)
cp "$(kpsewhich aastex701.cls)" "$APJ_DIR/"
cp "$(kpsewhich aasjournalv7.bst)" "$APJ_DIR/"

# Rewrite paths in tex: ../Analysis/figures/ -> (flat)
sed -i '' 's|\.\./Analysis/figures/||g' "$APJ_DIR/adaptive_gridding_for_snec.tex"

# Rewrite paths in tex: generated/ -> (flat)
sed -i '' 's|generated/||g' "$APJ_DIR/adaptive_gridding_for_snec.tex"

# Build in ApJ dir to verify and generate .log
cd "$APJ_DIR"
pdflatex -interaction=nonstopmode adaptive_gridding_for_snec.tex > /dev/null 2>&1 || true
bibtex adaptive_gridding_for_snec > /dev/null 2>&1 || true
pdflatex -interaction=nonstopmode adaptive_gridding_for_snec.tex > /dev/null 2>&1 || true
pdflatex -interaction=nonstopmode adaptive_gridding_for_snec.tex > /dev/null 2>&1 || true

# Summary
echo "ApJ submission prepared in Paper/ApJ/ (flat, no subdirectories)"
echo "Contents:"
find "$APJ_DIR" -type f | sort | while read -r f; do
    echo "  ${f#$REPO_ROOT/}"
done
