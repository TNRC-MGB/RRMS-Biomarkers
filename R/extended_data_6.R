
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : B_first_pass.RDS   (preprocessing_bcells.R, step 3 checkpoint)
#          B_clean.RDS        (preprocessing_bcells.R, cleaned + harmonized)
#
# NOTE: panels a and b are drawn BEFORE contaminant removal so they need the
#       first-pass object (clusters 0-32, contaminants present)
#

# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms") 

source("R/packages.R")
source("R/scrna_preprocessing_functions.R")

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Contaminant / low-quality clusters removed after the first pass. Kept here for
# reference; this must stay in sync with preprocessing_bcells.R
#   10 CD4 T   14 DC   19 RBC   25 T/NK   27 low quality
#   30 monocyte doublet         32 low quality
clusters_remove <- c("10", "14", "19", "25", "27", "30", "32")


# ---------------------------------------------------------------------------- #
# First-pass B-cell object (panels a, b)
B_first <- readRDS(file.path(OBJ_DIR, "B_first_pass.RDS"))

qcB <- plot_cluster_qc_panel(
  B_first,
  cluster_col = "clusters.first_pass",
  reduction   = "umap.first_pass"
)


# ---------------------------------------------------------------------------- #
# a - First-pass B-cell clustering UMAP
p_a <- qcB$umap + coord_fixed()
ggsave(file.path(out_dir, "extended_data_6a.pdf"), p_a,
       width = 500, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# b - Lineage / QC marker dot-plot panel
# Same canonical lineage marker panel as Extended Data Fig. 1b. Orange-highlighted
# rows are the clusters listed in clusters_remove above.
p_b <- qcB$panel
ggsave(file.path(out_dir, "extended_data_6b.pdf"), p_b,
       width = 500, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# Cleaned, harmonized B-cell object (panels c, d, e)
B <- readRDS(file.path(OBJ_DIR, "B_clean.RDS"))


# ---------------------------------------------------------------------------- #
# c - Second-pass clustering UMAP (pre-Harmony), colored by batch
p_c <- DimPlot(B, reduction = "umap.B", group.by = "Batch",
               raster = FALSE, shuffle = TRUE) +
  coord_fixed() +
  labs(title = "Second-pass clustering")

ggsave(file.path(out_dir, "extended_data_6c.pdf"), p_c,
       width = 100, height = 100, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# d - Harmonized UMAP by experimental batch
p_d <- DimPlot(B, reduction = "umap.B.harmony", group.by = "Batch",
               shuffle = TRUE) +
  coord_fixed() +
  labs(title = "Batch")

ggsave(file.path(out_dir, "extended_data_6d.pdf"), p_d,
       width = 100, height = 100, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# e - Final harmonized B-cell cluster identities
p_e <- DimPlot(B, reduction = "umap.B.harmony", group.by = "clusters.B.harmony",
               label = TRUE, shuffle = TRUE) +
  coord_fixed() +
  labs(title = "clusters.B.harmony")

ggsave(file.path(out_dir, "extended_data_6e.pdf"), p_e,
       width = 100, height = 100, units = "mm", dpi = 600)
