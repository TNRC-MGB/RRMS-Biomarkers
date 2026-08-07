
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : scdist_pbmc_relapse_paired.RDS   (scrna_pbmc_scdist.R)
#          scdist_pbmc_prerelapse_paired.RDS   (scrna_pbmc_scdist.R)
#
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

source("R/scdist_functions.R")   # DistPlot2

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

YLIM <- c(0, 12)   # shared x-range so a and b are directly comparable

# --- Inputs -------------------------------------------------------------------
scdist_pbmc_rel_paired <- readRDS(file.path(SCDIST_DIR, "scdist_pbmc_relapse_paired.RDS"))
scdist_pbmc_pre_paired <- readRDS(file.path(SCDIST_DIR, "scdist_pbmc_prerelapse_paired.RDS"))


# ---------------------------------------------------------------------------- #
# a - Relapse vs Remission (paired donors)
p_a <- DistPlot2(scdist_pbmc_rel_paired) +
  ylim(YLIM) +
  labs(title = "Relapse")

ggsave(file.path(out_dir, "extended_data_5a.pdf"), p_a,
       width = 120, height = 100, units = "mm")


# ---------------------------------------------------------------------------- #
# b - PreRelapse vs Remission (paired donors)
p_b <- DistPlot2(scdist_pbmc_pre_paired) +
  ylim(YLIM) +
  labs(title = "Pre-relapse")

ggsave(file.path(out_dir, "extended_data_5b.pdf"), p_b,
       width = 120, height = 100, units = "mm")


# ---------------------------------------------------------------------------- #
# Assembled figure (a | b), legend collected once
p_ab <- p_a + p_b + plot_layout(ncol = 2, guides = "collect")

ggsave(file.path(out_dir, "extended_data_5.pdf"), p_ab,
       width = 240, height = 100, units = "mm")
