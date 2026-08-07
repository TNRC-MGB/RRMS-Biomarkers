#!/usr/bin/env bash
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham / Harvard Medical School
#
# Updated: 7-25-2026
#
# scDRS step 04 - build the scDRS gene set from MAGMA gene-level z-scores
#
# Takes the top 1000 genes by MAGMA z-score and weights them by that z-score.
# The resulting MS.gs is a two-line file (header + one MS row) and is included
# in this repository (to skip this step).
#
# Inputs : scdrs/work/MS_magma_symbol_zscore.tsv   (from 03_magma_to_zscores.R)
# Outputs: scdrs/work/MS.gs
# ==============================================================================

set -euo pipefail

WORK=scdrs/work
N_MAX=1000          # number of top genes retained in the gene set

command -v scdrs >/dev/null 2>&1 || {
  echo "scdrs not found on PATH. Activate the environment first, e.g." >&2
  echo "  source /work/scdrs/scdrs/bin/activate" >&2
  exit 1
}

scdrs munge-gs \
  --out-file    "$WORK/MS.gs" \
  --zscore-file "$WORK/MS_magma_symbol_zscore.tsv" \
  --weight      zscore \
  --n-max        "$N_MAX"

echo "gene set written -> $WORK/MS.gs"
echo "(the version used for the manuscript is committed at scdrs/MS.gs)"
