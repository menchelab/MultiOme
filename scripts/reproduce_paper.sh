#!/usr/bin/env bash
# Full-scale run on the shipped data (46 layers, 26 Orphanet groups, seed 0).
# Produces the numbers in docs/CHANGES.md and the figures in docs/figures/.
# CV times on a laptop: train ~18 min, paper ~7 min, degree ~67 min.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-results}"
uv run multiome modularity --paper --per-group-plots --out "$OUT/paper"
uv run multiome rank --paper --group Ciliopathy --modularity "$OUT/paper/modularity.tsv" \
  --out "$OUT/paper"
uv run multiome cv --paper --out "$OUT/cv_train"
uv run multiome cv --paper --protocol paper --out "$OUT/cv_paper"
uv run multiome cv --paper --null degree --out "$OUT/cv_degree"
