
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-25-2026
#
# Dependencies
#
# Every analysis script begins with
#     setwd(<repository root>)
#     source('R/packages.R')
#
# ==============================================================================

# Hard requirement
.req <- function(...) {
  for (p in c(...)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Required package '", p, "' is not installed.\n",
           "  See the Requirements section of README.md for install commands.",
           call. = FALSE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

# Soft requirement
.opt <- function(p, needed_by) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
    invisible(TRUE)
  } else {
    warning("Optional package '", p, "' is not installed; ", needed_by,
            " cannot run. Everything else is unaffected.", call. = FALSE)
    invisible(FALSE)
  }
}

.assert_resolves <- function(fn, pkg) {
  if (!exists(fn)) return(invisible(NULL))
  where <- environmentName(environment(get(fn)))
  if (!identical(where, pkg)) {
    stop("Name collision: ", fn, "() resolves to '", where, "' but must resolve ",
         "to '", pkg, "'.\n  The load order in R/packages.R has been disturbed. ",
         "Restart R and source this file again before running any analysis.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# dependencies
.req("Matrix", "data.table")
.opt("mclust", "figure4.R (FlowSOM metacluster consensus)")
.opt("UpSetR", "supplementary_figure1.R (assay UpSet plot)")
.req("SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
     "scater", "scDblFinder", "fgsea")

# (masks dplyr::select - must precede the tidyverse)
.opt("AnnotationDbi", "scdrs/03_magma_to_zscores.R (Entrez/symbol mapping)")
.opt("org.Hs.eg.db",  "scdrs/03_magma_to_zscores.R (Entrez/symbol mapping)")
.opt("msigdbr",       "gene-set build only (see section 9)")

# Seurat and integration
.req("SeuratObject", "Seurat", "harmony")
.opt("Azimuth", "preprocessing_pbmc.R (pbmcref reference mapping)")
.opt("scDist",  "scrna_pbmc_scdist.R, preprocessing_bcells.R, scdist helpers")

# statistics - lme4 MUST be attached before lmerTest
.req("lme4")
.req("lmerTest")
.req("emmeans")
.opt("broom.mixed", "figure4.R (tidy mixed-model output)")

# flow cytometry and qPCR
.opt("FlowSOM",   "figure4.R (flow clustering)")
.opt("Spectre",   "figure4.R (flow clustering)")
.opt("RQdeltaCT", "figure5.R (qPCR delta-Ct analysis)")

# plotting (scales masks purrr::discard - precedes the tidyverse)
.req("scales", "ggplot2", "patchwork", "colorspace", "ggrepel")
.opt("paletteer", "figure2.R (palettes)")
.opt("ggrastr",   "figure2.R (rasterised UMAP layers)")

# parallel / misc
.req("BiocParallel")
.opt("future", "scrna_pbmc_scdist.R (parallel scDist)")

# # --- 9. peak annotation - ONLY for rebuilding the EBNA2 gene sets -------------
# # The builders are not distributed: their products ship in data/genesets/, so
# # nothing here is needed to reproduce a figure. Off unless the option is set.
# if (isTRUE(getOption("rrms.gene_set_build", FALSE))) {
#   .opt("GEOquery",                          "build_ebna2_gmt.R (GEO downloads)")
#   .opt("GenomicRanges",                     "build_ebna2_gmt.R (peak ranges)")
#   .opt("IRanges",                           "build_ebna2_gmt.R (peak ranges)")
#   .opt("rtracklayer",                       "build_ebna2_gmt.R (BED/narrowPeak import)")
#   .opt("ChIPseeker",                        "build_ebna2_gmt.R (peak-to-gene annotation)")
#   .opt("TxDb.Hsapiens.UCSC.hg19.knownGene", "build_ebna2_gmt.R (hg19 annotation)")
# }

# tidyverse LAST, so dplyr and purrr win every remaining collision
.req("tibble", "readr", "readxl", "stringr", "purrr", "tidyr", "dplyr")

if (requireNamespace("lmerTest", quietly = TRUE)) .assert_resolves("lmer",   "lmerTest")
.assert_resolves("select", "dplyr")
.assert_resolves("filter", "dplyr")
.assert_resolves("map",    "purrr")
.assert_resolves("first",  "dplyr")
.assert_resolves("discard","purrr")

message("R/packages.R: dependencies loaded (lmer -> lmerTest, select/filter -> dplyr).")
