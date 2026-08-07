
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Rebuilds the EBV factor-derived host gene sets consumed by
# scrna_bcell_nebula_enrichment_analysis.R, extended_data_9.R,
# extended_data_10.R, figure5.R, and bulk/figure5e.R
#
# Inputs : data/genesets/S14 EBV Host Response Correlation Matrix.xlsx
#            (Spearman rho, 22,467 human genes x 79 EBV genes, derived
#            from Arvey et al. 2012, PMID 22901543)
# Outputs: data/genesets/pebv500.RDS   top 500 host genes per EBV factor
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

ENRICH_DIR <- "data/genesets"

# Rank-based selection only: no rho cutoff and no significance filter,
# so every module is the same size
N_TOP <- 500L


# ---------------------------------------------------------------------------- #
# Correlation matrix

ac <- as.data.frame(read_excel(
  file.path(ENRICH_DIR, "S14 EBV Host Response Correlation Matrix.xlsx"),
  sheet = "Correlation_matrix"
))

row.names(ac) <- ac[[1]]
ac <- as.matrix(ac[, -1])

stopifnot(!anyNA(ac), !anyDuplicated(row.names(ac)), nrow(ac) >= N_TOP)


# ---------------------------------------------------------------------------- #
# Top N host genes per EBV factor

pebv500 <- sapply(colnames(ac), function(i) {
  aci <- ac[, i]
  aci <- sort(aci, decreasing = TRUE)
  names(aci[1:N_TOP])
}, simplify = FALSE, USE.NAMES = TRUE)

saveRDS(pebv500, file.path(ENRICH_DIR, "pebv500.RDS"))

message("Wrote pebv500.RDS: ", length(pebv500), " EBV factor modules of ",
        N_TOP, " host genes")
