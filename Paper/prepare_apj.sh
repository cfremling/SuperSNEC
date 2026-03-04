#!/usr/bin/env bash
# Prepare a self-contained ApJ submission directory in Paper/ApJ/
# Run from the repository root: bash Paper/prepare_apj.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
APJ_DIR="$REPO_ROOT/Paper/ApJ"

# Clean and recreate
rm -rf "$APJ_DIR"
mkdir -p "$APJ_DIR/figures" "$APJ_DIR/generated"

# Copy main manuscript and bibliography
cp "$REPO_ROOT/Paper/adaptive_gridding_for_snec.tex" "$APJ_DIR/"
cp "$REPO_ROOT/Paper/adaptive_gridding_for_snec.bib" "$APJ_DIR/"

# Copy generated tables
cp "$REPO_ROOT/Paper/generated/"*.tex "$APJ_DIR/generated/"

# Copy figures
cp "$REPO_ROOT/Analysis/figures/adaptive_runtime_validation.png" "$APJ_DIR/figures/"
cp "$REPO_ROOT/Analysis/figures/lc_comparison_sne.pdf" "$APJ_DIR/figures/"

# Copy built PDF
cp "$REPO_ROOT/Paper/adaptive_gridding_for_snec.pdf" "$APJ_DIR/"

# Rewrite figure paths: ../Analysis/figures/ -> figures/
sed -i '' 's|\.\./Analysis/figures/|figures/|g' "$APJ_DIR/adaptive_gridding_for_snec.tex"

# Summary
echo "ApJ submission prepared in Paper/ApJ/"
echo "Contents:"
find "$APJ_DIR" -type f | sort | while read -r f; do
    echo "  ${f#$REPO_ROOT/}"
done
