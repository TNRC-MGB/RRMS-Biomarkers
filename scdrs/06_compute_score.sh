#!/usr/bin/env bash
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham / Harvard Medical School
#
#
# Usage:
#   bash scdrs/06_compute_score.sh <h5ad> <cov> <out-folder> [gs-file]
#
#   bash scdrs/06_compute_score.sh scdrs/work/ms_pbmc.h5ad \
#                                  scdrs/work/ms_pbmc.cov \
#                                  scdrs/output/pbmc
#
#   bash scdrs/06_compute_score.sh scdrs/work/ms_bcells.h5ad \
#                                  scdrs/work/ms_bcells.cov \
#                                  scdrs/output/bcells
#
# R/figure2.R and R/preprocessing_bcells.R read scdrs/output/pbmc/MS.score.gz.
# ==============================================================================

set -euo pipefail

H5AD=${1:?usage: 06_compute_score.sh <h5ad> <cov> <out-folder> [gs-file]}
COV=${2:?usage: 06_compute_score.sh <h5ad> <cov> <out-folder> [gs-file]}
OUTDIR=${3:?usage: 06_compute_score.sh <h5ad> <cov> <out-folder> [gs-file]}

# Defaults to the gene set committed with this repository, so MAGMA can be
# skipped. Pass scdrs/work/MS.gs to use a freshly built one.
GS=${4:-scdrs/MS.gs}

# Number of control gene sets used to build the null e.g. 1000.
N_CTRL=1000

command -v scdrs >/dev/null 2>&1 || {
  echo "scdrs not found on PATH. Activate the environment first, e.g." >&2
  echo "  source /work/scdrs/scdrs/bin/activate" >&2
  exit 1
}

for f in "$H5AD" "$COV" "$GS"; do
  [ -f "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

scdrs compute-score \
  --h5ad-file         "$H5AD" \
  --h5ad-species      human \
  --gs-file           "$GS" \
  --gs-species        human \
  --cov-file          "$COV" \
  --out-folder        "$OUTDIR" \
  --flag-filter-data  True \
  --flag-raw-count    True \
  --n-ctrl            "$N_CTRL"

echo "scDRS scores written -> $OUTDIR/MS.score.gz"
