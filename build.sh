#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
if [[ "${1:-}" == "--recompute" ]]; then
  python code/run_benchmarks.py --part all
  python code/run_extension.py
  python code/make_assets.py
elif [[ -n "${1:-}" ]]; then
  echo 'Usage: bash build.sh [--recompute]' >&2; exit 2
fi
mkdir -p .build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=.build manuscript.tex >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=.build manuscript.tex >/dev/null
cp .build/manuscript.pdf manuscript.pdf
printf 'Built manuscript.pdf\n'
