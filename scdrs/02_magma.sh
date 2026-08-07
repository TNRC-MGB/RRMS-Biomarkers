#!/usr/bin/env bash
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham / Harvard Medical School
#
# Updated: 7-25-2026
#
# scDRS step 02 - MAGMA: MS GWAS SNPs -> gene-level association statistics
#
#
# Inputs : scdrs/ref/g1000_eur/g1000_eur.{bed,bim,fam,synonyms}   (LD reference)
#          scdrs/ref/gene_loc/NCBI37.3.gene.loc                   (gene locations)
#          scdrs/work/MS.final.magma.tsv                          (from step 01)
# Outputs: scdrs/work/MS_GWAS.genes.annot
#          scdrs/work/MS_GWAS.genes.out
#          scdrs/work/MS_GWAS.genes.raw
#          scdrs/work/MS_GWAS.log
# ==============================================================================

set -euo pipefail

REF=scdrs/ref
WORK=scdrs/work
MAGMA=${MAGMA:-"$REF/magma/magma"}      # override with MAGMA=/path/to/magma

N_GWAS=41505

# SNP -> gene window, in kb, upstream,downstream
WINDOW="10,10"

mkdir -p "$WORK"

command -v "$MAGMA" >/dev/null 2>&1 || [ -x "$MAGMA" ] || {
  echo "MAGMA binary not found at '$MAGMA'." >&2
  echo "Download v1.10 from https://cncr.nl/research/magma/ and set MAGMA=..." >&2
  exit 1
}

# --- SNP -> gene annotation ---------------------------------------------------
"$MAGMA" \
  --annotate window="$WINDOW" \
  --snp-loc  "$REF/g1000_eur/g1000_eur.bim" \
  --gene-loc "$REF/gene_loc/NCBI37.3.gene.loc" \
  --out      "$WORK/MS_GWAS"

# --- gene analysis (SNPwise-mean, MAGMA's default) ---------------------------
"$MAGMA" \
  --bfile      "$REF/g1000_eur/g1000_eur" \
  --pval       "$WORK/MS.final.magma.tsv" N="$N_GWAS" \
  --gene-annot "$WORK/MS_GWAS.genes.annot" \
  --out        "$WORK/MS_GWAS"

echo "MAGMA gene analysis complete -> $WORK/MS_GWAS.genes.out"
