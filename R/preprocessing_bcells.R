
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : pbmc_final.RDS            (scrna preprocessing pipeline)
#          scdrs/output/pbmc/MS.score.gz  (scDRS, MS GWAS - see scdrs/README.md)
# Outputs: B_first_pass.RDS         (first-pass B-cell object, pre-contaminant
#                                    removal; Extended Data Fig. 6a/b)
#          B_clean.RDS              (cleaned, annotated B-cell object)
#          B_clean_scdrs.RDS        (B_clean + scDRS norm_score/zscore/ etc)
#          scdist_bcell_prerelapse_v_remission.RDS
#            (scDist, PreRelapse vs Remission)
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

source("R/scrna_preprocessing_functions.R") 
source("R/scdist_functions.R")

OBJ_DIR <- "objects"

# scDRS output for the MS GWAS - produced by the pipeline in scdrs/ and included
# in Zenodo/this repo
SCDRS_SCORE <- "scdrs/output/pbmc/MS.score.gz"


# ---------------------------------------------------------------------------- #
# Load and slim the full PBMC atlas
pbmc <- readRDS(file.path(OBJ_DIR, "pbmc_final.RDS"))

pbmc <- DietSeurat(pbmc, assays = c("RNA", "impADT"), layers = c("counts", "data"))


# ---------------------------------------------------------------------------- #
# Subset B cells from the atlas
B_cells <- c(
  "TCL1A+ Naive B",
  "Naive B",
  "Intermediate B",
  "Switched Memory B",
  "Plasmablast"
)

B <- subset(pbmc, subset = pbmc_annotations %in% B_cells)


# ---------------------------------------------------------------------------- #
# First-pass reclustering of the B-cell compartment
#     run_seurat_standard: NormalizeData -> FindVariableFeatures -> ScaleData ->
#     RunPCA -> FindNeighbors -> FindClusters -> RunUMAP
# ---------------------------------------------------------------------------- #
B <- run_seurat_standard(
  B,
  split_by     = "Batch",
  npcs         = 50,
  dims         = 1:30,
  resolution   = 2,
  pca.name     = "pca.first_pass",
  umap.name    = "umap.first_pass",
  cluster.name = "clusters.first_pass",
  seed         = 1,
  verbose      = TRUE
)

# Checkpoint the first-pass object BEFORE contaminant removal. Extended Data
# Fig. 6a/b are drawn at this stage
saveRDS(B, file.path(OBJ_DIR, "B_first_pass.RDS"))


# ---------------------------------------------------------------------------- #
# Remove contaminant / low-quality clusters (T/NK, monocyte, RBC, DC, platelet)
clusters_remove <- c("10", "14", "19", "25", "27", "30", "32")

B <- subset(B, subset = !(clusters.first_pass %in% clusters_remove))


# ---------------------------------------------------------------------------- #
# Re-run standard preprocessing on the cleaned B-cell object
B <- run_seurat_standard(
  B,
  split_by     = "Batch",
  npcs         = 50,
  dims         = 1:30,
  resolution   = 2,
  pca.name     = "pca.B",
  umap.name    = "umap.B",
  cluster.name = "clusters.B",
  seed         = 1,
  verbose      = TRUE
)


# ---------------------------------------------------------------------------- #
# Harmony integration on Batch + second-pass clustering
set.seed(1)
B <- RunHarmony(
  object         = B,
  group.by.vars  = "Batch",
  reduction.use  = "pca.B",
  dims.use       = 1:30,
  assay.use      = "RNA",
  reduction.save = "B.harmony"
)

B <- FindNeighbors(B, reduction = "B.harmony", dims = 1:30)

set.seed(1)
B <- FindClusters(B, resolution = 1, cluster.name = "clusters.B.harmony")

set.seed(1)
B <- RunUMAP(B, reduction = "B.harmony", dims = 1:30, reduction.name = "umap.B.harmony")


# ---------------------------------------------------------------------------- #
# Annotate harmonized B-cell subclusters (final B_annotations, 14 subtypes)
B_manual_cluster_labels <- c(
  "0"  = "Anergic/Naive-leaning B",
  "1"  = "Activated Memory B",
  "2"  = "TCL1A+ Naive B",
  "3"  = "Activated Switched Memory B",
  "4"  = "Naive B",
  "5"  = "Transitional/Immature-like B",
  "6"  = "Switched Memory B",
  "7"  = "ABC",
  "8"  = "Atypical/ABC-like Memory B",
  "9"  = "CD11c++ Activated Memory B",
  "10" = "Early Activated Memory B",
  "11" = "Plasmablast",
  "12" = "IFN-stimulated Naive B",
  "13" = "Discard_T-cell Contaminant",
  "14" = "Cycling Plasmablast"
)

B_final_levels <- c(
  "Transitional/Immature-like B",
  "TCL1A+ Naive B",
  "Naive B",
  "Anergic/Naive-leaning B",
  "IFN-stimulated Naive B",
  "Switched Memory B",
  "Activated Switched Memory B",
  "Activated Memory B",
  "Early Activated Memory B",
  "Atypical/ABC-like Memory B",
  "CD11c++ Activated Memory B",
  "ABC",
  "Plasmablast",
  "Cycling Plasmablast"
)

B$B_annotations <- unname(B_manual_cluster_labels[as.character(B$clusters.B.harmony)])

# cluster 13 shows a strong T-cell program (BCL11B, CD247, THEMIS, IL7R, CD2, CD6);
# a non-B contaminant, dropped from the final object
B_clean <- subset(B, subset = !(clusters.B.harmony %in% c("13")))
B_clean$B_annotations <- factor(B_clean$B_annotations, levels = B_final_levels)

saveRDS(B_clean, file.path(OBJ_DIR, "B_clean.RDS"))


# ---------------------------------------------------------------------------- #
# add MS-GWAS scDRS scores from the PBMC scDRS run
sc <- fread(SCDRS_SCORE)
setnames(sc, 1, "cell")
stopifnot(identical(sc$cell, colnames(pbmc)))
pbmc <- AddMetaData(pbmc, metadata = as.data.frame(sc), col.name = NULL)

pbmc.b <- subset(pbmc, cells = Cells(B_clean))
stopifnot(identical(Cells(pbmc.b), Cells(B_clean)))

B_clean$norm_score  <- pbmc.b$norm_score
B_clean$zscore      <- pbmc.b$zscore
B_clean$raw_score   <- pbmc.b$raw_score
B_clean$mc_pval     <- pbmc.b$mc_pval
B_clean$pval        <- pbmc.b$pval
B_clean$nlog10_pval <- pbmc.b$nlog10_pval

saveRDS(B_clean, file.path(OBJ_DIR, "B_clean_scdrs.RDS"))


# ---------------------------------------------------------------------------- #
# scDist: PreRelapse vs Remission across B-cell subtypes
pre_all_filt <- prefilter_clusters_for_scdist(
  seu                        = B_clean,
  condition_col              = "Condition",
  patient_col                = "Patient",
  cluster_col                = "B_annotations",
  levels_use                 = c("Remission", "PreRelapse"),
  min_cells_per_condition    = 20,
  min_patients_per_condition = 2
)

scdist_bcell_pre <- run_scdist_standard(
  seu             = pre_all_filt$seu,
  fixed_effect    = "Condition",
  reference_level = "Remission",
  random_effect   = "Patient",
  cluster_col     = "B_annotations",
  assay           = "RNA",
  levels_use      = c("Remission", "PreRelapse"),
  d               = 20
)

saveRDS(scdist_bcell_pre, file.path(SCDIST_DIR, "scdist_bcell_prerelapse_v_remission.RDS"))
