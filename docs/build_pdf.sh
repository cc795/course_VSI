#!/usr/bin/env bash
set -euo pipefail
project=$(cd "$(dirname "$0")/.." && pwd)
cd "$project"
# Requires Pandoc, XeLaTeX, PyMuPDF, and Noto Serif/Sans CJK SC fonts.
pandoc docs/section3_implementation.md \
  --pdf-engine=xelatex --resource-path=.:docs -o readme.pdf

python3 docs/style_pdf.py
