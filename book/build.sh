#!/bin/sh
# Generate the LaTeX and compile it to PDF.
#   ./build.sh                     -> catalog.pdf   (2 pages, empty notes)
#   ./build.sh --notes --flow      -> same, notes filled from BSC5
# Any arguments are passed straight through to build_catalog.py.
#
# BasicTeX has no latexmk, and longtable needs two passes to settle its
# column widths, so we just run pdflatex twice.
set -e

cd "$(dirname "$0")"
PY=/Users/dseverin/.venv/bin/python
TEX=/Library/TeX/texbin/pdflatex

OUT=catalog.tex
for a in "$@"; do
    case "$prev" in --out) OUT=$a ;; esac
    prev=$a
done
BASE=$(basename "$OUT" .tex)

"$PY" build_catalog.py "$@"
"$TEX" -interaction=nonstopmode -halt-on-error "$OUT" >/dev/null
"$TEX" -interaction=nonstopmode -halt-on-error "$OUT" >/dev/null
rm -f "$BASE.aux" "$BASE.log" "$BASE.out"

echo "wrote $BASE.pdf"
