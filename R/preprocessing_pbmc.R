
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 8-4-2026
#
# PBMC scRNA-seq preprocessing pipeline
#   Builds the whole-PBMC atlas used throughout the paper.
#
# Inputs : rrms/zenodo/scrna/counts/<pool>/outs/filtered_feature_bc_matrix
#            (Cell Ranger)
#          rrms/zenodo/scrna/demux/<pool>/vireo/donor_ids.tsv      (Vireo demux)
#          rrms/Supplementary Information/S4 scRNA Libraries.xlsx
#            (library metadata)
# Output : pbmc_final.RDS
#
# ------------------------------------------------------------------------------
# NOTE ON REPRODUCIBILITY
# ------------------------------------------------------------------------------
# Three steps in this pipeline are stochastic and their output depends on the 
# random seed and also package versions
#
#  - scDblFinder, which trains an xgboost classifier on artificial doublets
#  - Harmony integration and the Louvain clustering plus UMAP
#  - Azimuth reference mapping, which also depends on the version of the
#         `pbmcref` reference distributed through SeuratData.
#
# The atlas for the published figures was built under R 4.5.2 in March
# 2026. Re-running the identical code with a later R installation reproduces the
# QC output exactly but diverges at scDblFinder and later steps.
# No value of set.seed() recovers the original result once the package set has 
# changed; see R/published_environment_notes.txt.
#
# This script therefore re-runs the published stochastic outcomes from small
# frozen tables included in this repo. There are several blocks like:
#
#     ## ---- ORIGINAL (as run for the publication) - kept for reference -------
#     # <the original, stochastic code, commented out>
#     ## ---- PUBLISHED VALUES (used when this script is run) ------------------
#     <code that loads the frozen result>
#
# The commented block is the method as it was actually performed; the live block
# reproduces its outcome deterministically.
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")
source("R/scrna_import_functions.R")
source("R/scrna_preprocessing_functions.R")


# ---------------------------------------------------------------------------- #
# Tables for the published result (see NOTE above)
FROZEN_DOUBLETS <- "data/frozen/scDblFinder_doublet_calls.tsv.gz"
FROZEN_META     <- "data/frozen/published_cell_metadata.tsv.gz"

read_published_metadata <- function(path) {
  if (!file.exists(path)) stop("frozen metadata not found: ", path)
  con <- gzfile(path, "rt"); hdr <- readLines(con, n = 1L); close(con)
  levs <- strsplit(sub('^#\\s*annotation_levels\\t', "", hdr), '\\|', fixed = FALSE)[[1]]
  df <- read.delim(gzfile(path), skip = 1L, stringsAsFactors = FALSE,
                   quote = "", check.names = FALSE)
  for (cc in c("cell", "Sample_ID", "Pool_ID", "clusters_harmony", "pbmc_annotations"))
    if (cc %in% names(df)) df[[cc]] <- as.character(df[[cc]])
  attr(df, "annotation_levels") <- levs
  df
}


# ---------------------------------------------------------------------------- #
# Import - Cell Ranger counts + Vireo genotype demux
data_dirs  <- dir("zenodo/scrna/counts", full.names = TRUE)
vireo_dirs <- dir("zenodo/scrna/demux",  full.names = TRUE)
data_dirs  <- file.path(data_dirs,  "outs/filtered_feature_bc_matrix")
vireo_dirs <- file.path(vireo_dirs, "vireo")
names(data_dirs)  <- get_pool_id(data_dirs)
names(vireo_dirs) <- get_pool_id(vireo_dirs)
if (any(is.na(names(data_dirs))))  stop("Could not parse Pool_ID from some data_dirs.")
if (any(is.na(names(vireo_dirs)))) warning("Could not parse Pool_ID from some vireo_dirs; those pools will be treated as no-Vireo if needed.")

# Metadata determines which pools to import.
md <- readxl::read_excel("Supplementary Information/S4 scRNA Libraries.xlsx")
md <- validate_library_meta(md)
pool_ids <- stringr::str_sort(unique(md$Pool_ID), numeric = TRUE)
missing_counts <- setdiff(pool_ids, names(data_dirs))
if (length(missing_counts) > 0) stop("Missing data_dirs for pool(s): ", paste(missing_counts, collapse = ", "))
data_dirs2  <- data_dirs[pool_ids]
vireo_dirs2 <- vireo_dirs[pool_ids]

# Donors expected in each pool (sanity-checks Vireo assignments; empty = single-donor)
expected_ids <- list(W01=c("MS01","MS05","MS06"),
                     W02=c("MS01","MS06","MS09"),
                     W03=c("MS01","MS06","MS09"),
                     W04=c("MS01","MS05","MS06"),
                     W05=c(), # single
                     W06=c(), # single
                     W07=c("HC01","HC02","HC03"),
                     W08=c(), # single
                     W09=c("MS02","MS03","MS04"),
                     W10=c(), # single
                     W11=c(), # single
                     W12=c("MS02","MS03","MS04"),
                     BW1=c("MS12","HC04"),
                     BW2=c("MS13","HC05"),
                     BW3=c("HC06"), # pooled with other study
                     BW4=c("MS07","HC07"),
                     BW5=c("HC08"), # pooled with other study
                     BW6=c("HC09"), # pooled with other study
                     BW7=c("HC10"), # pooled with other study
                     BW8=c("MS14","HC11"),
                     BW9=c("HC12"), # pooled with other study
                     BW10=c("MS08","HC13"),
                     BW11=c("MS15","HC14"),
                     BW12=c("HC15","HC16"),
                     BW13=c("MS10","HC17","HC18"),
                     BW14=c("HC19","HC20"),
                     BW15=c("HC21"), # pooled with other study
                     BW16=c("MS11")) # pooled with other study

# Import raw counts + genotype-based Vireo demux/doublet calls, per pool
seu_list <- setNames(lapply(pool_ids, function(pid) {
  vdir <- vireo_dirs2[[pid]]
  has_vireo <- !is.na(vdir) && nzchar(vdir) && dir.exists(vdir) &&
    file.exists(file.path(vdir, "donor_ids.tsv"))

  donor_unpooled <- if (!has_vireo) infer_unpooled_donor(pid, md) else NULL
  if (!has_vireo && is.na(donor_unpooled)) {
    stop("Pool ", pid, " has no Vireo output and metadata contains >1 Patient_ID. ",
         "Provide donor_id_unpooled for this pool (or fix metadata).")
  }

  import_counts_vireo(
    data_dir10x = data_dirs2[[pid]],
    pool_id = pid,
    library_meta = md,
    vireo_dir = if (has_vireo) vdir else NULL,
    expected_map = expected_ids,        # optional; can be NULL
    donor_id_unpooled = donor_unpooled, # used only when Vireo missing
    allow_missing_vireo = TRUE,
    verbose = TRUE
  )
}), pool_ids)
saveRDS(seu_list, "objects/seu_list.RDS")

# Merge pools into one object
seu <- merge(seu_list[[1]], y = seu_list[-1])
seu <- JoinLayers(seu)
saveRDS(seu, "objects/seu_raw.RDS")


# ---------------------------------------------------------------------------- #
# Minimal QC - per 10x GEM well
# Each Pool_ID is one GEM well (the partitioning reaction where doublets form);
# drop empty droplets / extreme low-feature outliers before doublet detection.
# This step is fully deterministic and reproduces exactly across R versions.
seu <- add_basic_qc(seu)
seu <- qc_filter_per_GEM_well(seu,
                              split_by = "Pool_ID",
                              min_features = 300,
                              min_counts = 500,
                              max_percent_mt = 10,
                              nmads = 2,
                              outlier_type = "lower",
                              log_outliers = TRUE,
                              verbose = TRUE)
saveRDS(seu, "objects/seu_minQC.RDS")


# ---------------------------------------------------------------------------- #
# Expression doublets - scDblFinder (per GEM well)
## ---- ORIGINAL (as run for the publication and kept for reference) -----------
## scDblFinder generates artificial doublets at random and trains an xgboost
## classifier on them. set.seed() alone does not make this reproducible across
## package versions: the classifier itself changed between the March 2026 run
## and any current install, so the same seed yields marginally different calls.
#
# seu <- run_scDblFinder_per_GEM_well(seu,
#                                     split_by = 'Pool_ID',
#                                     assay = 'RNA',
#                                     counts_layer = 'counts',
#                                     dbl_rate = NULL,
#                                     filter = FALSE,
#                                     seed = 1,
#                                     verbose = TRUE)
#
## ---- PUBLISHED VALUES (used when this script is run) ------------------------
message("Applying published scDblFinder calls from ", FROZEN_DOUBLETS)
.frozen <- read.delim(gzfile(FROZEN_DOUBLETS), stringsAsFactors = FALSE, quote = "")
.idx <- match(colnames(seu), .frozen$cell)
if (anyNA(.idx)) {
  stop(sum(is.na(.idx)), " of ", ncol(seu), " cells have no published ",
       "scDblFinder call. The QC step above did not reproduce the published ",
       "cell set. See R/published_environment_notes.txt.")
}
seu$scDblFinder_score <- as.numeric(.frozen$scDblFinder_score[.idx])
seu$scDblFinder_class <- as.character(.frozen$scDblFinder_class[.idx])
seu$scDblFinder_call  <- seu$scDblFinder_class == "doublet"
message("  matched ", ncol(seu), " cells; ", sum(seu$scDblFinder_call),
        " called as doublets")
rm(.frozen, .idx)
saveRDS(seu, "objects/seu_minQC_scDblFinder.RDS")


# ---------------------------------------------------------------------------- #
# Doublet removal - union(Vireo genotype, scDblFinder expression)
seu <- process_doublets(seu,
                        GEM_well = "Pool_ID",
                        filter = TRUE,
                        verbose = TRUE)

# Keep only samples belonging to this study
md <- seu[[]]
keep <- rownames(md)[!is.na(md$Sample_ID)]
seu <- subset(seu, cells = keep)
saveRDS(seu, "objects/seu_minQC_scDblFinder_filt.RDS")


# ---------------------------------------------------------------------------- #
# First-pass clustering - isolate straggler low-quality / doublet-like cells
# Permissive first pass so that ambient-RNA / low-UMI junk, stressed
# or dying cells, residual doublet-like mixtures and RBC/platelet contamination
# collapse into their own clusters for removal.
seu <- run_seurat_standard(seu,
                           split_by = "Batch",
                           npcs = 50,
                           dims = 1:30,
                           resolution = 2,
                           pca.name = "pca.first_pass",
                           umap.name = "umap.first_pass",
                           cluster.name = "clusters.first_pass",
                           seed = 1,
                           verbose = TRUE)


# ---------------------------------------------------------------------------- #
# Second-pass clustering
## ---- ORIGINAL (as run for the publication) - kept for reference -------------
## Clusters 31, 36, 37, 45, 47 and 49 of the first pass were identified by
## marker inspection as low-quality / ambient / residual-doublet populations and
## removed (3.1% of cells). Those cluster NUMBERS are an emergent property of
## the Louvain run above: under a different Seurat/igraph version the same
## numbers point at different populations so they cannot be reapplied blindly.
#
# clusters_to_remove <- c(31, 36, 37, 45, 47, 49)
# seu$first_pass_remove <- ifelse(seu$clusters.first_pass %in% clusters_to_remove, TRUE, FALSE)
#
## ---- PUBLISHED VALUES (used when this script is run) ------------------------
## No cells are removed after this point, so the published cell set itself
## defines the first-pass filter: any cell absent from published_cell_metadata
## is a cell the filter removed. Applying it by barcode is exact and immune to
## renumbering.
.pub <- read_published_metadata(FROZEN_META)
.lost <- setdiff(.pub$cell, colnames(seu))
if (length(.lost) > 0) {
  stop(length(.lost), " cell(s) present in the published atlas are missing from ",
       "this object. Upstream steps did not reproduce - do not proceed.")
}
seu$first_pass_remove <- !(colnames(seu) %in% .pub$cell)
message("First-pass filter (published): removing ", sum(seu$first_pass_remove),
        " of ", ncol(seu), " cells (",
        round(100 * mean(seu$first_pass_remove), 2), "%)")

# Checkpoint. Extended Data 1 (first-pass panels) reads this object. It is
# written here, after `first_pass_remove` has been set, so that the extended
# data script can read the published directly instead of re-deriving it
# from the version-dependent cluster numbers above.
saveRDS(seu, "objects/seu_minQC_scDblFinder_filt_firstPass.RDS")

seu <- subset(seu, subset = !(first_pass_remove))

DefaultAssay(seu) <- "RNA"
DefaultLayer(seu[["RNA"]]) <- "counts"
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL

seu <- run_seurat_standard(seu,
                           split_by = "Batch",
                           npcs = 50,
                           dims = 1:30,
                           resolution = 2,
                           pca.name = "pca.second_pass",
                           umap.name = "umap.second_pass",
                           cluster.name = "clusters.second_pass",
                           seed = 1,
                           verbose = TRUE)
# Checkpoint. Extended Data 1 (second-pass / batch / Harmony panels) reads this object:
saveRDS(seu, "objects/seu_minQC_scDblFinder_filt_secondPass.RDS")


# ---------------------------------------------------------------------------- #
# Harmony for batch integration
# NOTE: the published run used the pre-1.0 harmony in which the input
# reduction was named `reduction` and `assay.use` was required:
#
#   RunHarmony(object, group.by.vars = "Batch", reduction = "pca.second_pass",
#              dims.use = 1:30, assay.use = "RNA", reduction.save = "harmony")
#
# Current (July 2026) harmony renamed that argument to `reduction.use` and dropped
# `assay.use`. The call below is the current equivalent.
set.seed(1)
# seu <- RunHarmony(
#   object         = seu,
#   group.by.vars  = "Batch",
#   reduction.use  = "pca.second_pass",
#   dims.use       = 1:30,
#   reduction.save = "harmony"
# )

seu <- RunHarmony(
  object = seu,
  group.by.vars = "Batch",
  reduction = "pca.second_pass",
  dims.use = 1:30,
  assay.use = "RNA",
  reduction.save = "harmony"
)

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
set.seed(1)
seu <- FindClusters(seu, resolution = 2, cluster.name = "clusters.harmony")
set.seed(1)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
saveRDS(seu, "objects/seu_minQC_scDblFinder_filt_secondPass_harmony.RDS")


# ---------------------------------------------------------------------------- #
# Azimuth - reference mapping to pbmcref (+ imputed ADT)
set.seed(1)
seu <- RunAzimuth.Seurat.v5(
  query = seu,
  reference = "pbmcref",
  umap.name = "azimuth",
  do.adt = TRUE,
  assay = "RNA"
)
saveRDS(seu, "objects/seu_minQC_scDblFinder_filt_secondPass_harmony_azimuth.RDS")


# ---------------------------------------------------------------------------- #
# Manual annotation
## ---- ORIGINAL (as run for the publication) - kept for reference -------------
## Each Harmony cluster was labeled by its majority Azimuth call (>= 0.80),
## and fifteen clusters were then relabeled manually after marker inspection.
## As with the first-pass filter, the keys are Louvain CLUSTER NUMBERS and are
## only meaningful for the specific clustering they were derived from
#
# res <- assign_cluster_labels_from_azimuth(
#   seu = seu,
#   cluster_col = "clusters.harmony",
#   azimuth_col = "predicted.celltype.l2",
#   threshold = 0.80,
#   new_col = "clusters.harmony.annot"
# )
# seu <- res$seu
#
# manual_cluster_labels <- c(
#   "6"  = "Cytotoxic gdT",
#   "11" = "CD8 Naive/TCM",
#   "12" = "CD8 TEM (GZMK+)",
#   "14" = "CD8 Naive-like",
#   "15" = "Activated Treg",
#   "18" = "Switched Memory B",
#   "22" = "FGFBP2+ Cytotoxic gdT",
#   "23" = "Activated cDC2",
#   "27" = "CD4 Naive/TCM",
#   "30" = "CD56bright NK",
#   "37" = "Naive B",
#   "38" = "TCL1A+ Naive B",
#   "39" = "Cycling Cytotoxic Lymphocytes",
#   "40" = "Mature NK",
#   "41" = "CD16 Mono"
# )
#
# label_renames <- c(
#   "B naive" = "Naive B",
#   "B intermediate" = "Intermediate B"
# )
#
# cluster_id <- as.character(seu$clusters.harmony)
# seu$clusters.harmony.annot.manual <- as.character(seu$clusters.harmony.annot)
# idx_cluster <- cluster_id %in% names(manual_cluster_labels)
# seu$clusters.harmony.annot.manual[idx_cluster] <- manual_cluster_labels[cluster_id[idx_cluster]]
# idx_label <- seu$clusters.harmony.annot.manual %in% names(label_renames)
# seu$clusters.harmony.annot.manual[idx_label] <- label_renames[seu$clusters.harmony.annot.manual[idx_label]]
#
# final_levels <- c(
#   "CD4 Naive/TCM", "CD4 TCM", "Activated Treg", "CD8 Naive", "CD8 Naive-like",
#   "CD8 Naive/TCM", "CD8 TEM (GZMK+)", "MAIT", "Cytotoxic gdT",
#   "FGFBP2+ Cytotoxic gdT", "CD56bright NK", "NK", "Mature NK",
#   "Cycling Cytotoxic Lymphocytes", "TCL1A+ Naive B", "Naive B",
#   "Intermediate B", "Switched Memory B", "Plasmablast", "pDC",
#   "Activated cDC2", "CD14 Mono", "CD16 Mono"
# )
# seu$clusters.harmony.annot.manual <- factor(seu$clusters.harmony.annot.manual, levels = final_levels)
# seu$pbmc_annotations <- seu$clusters.harmony.annot.manual
#
.m <- match(colnames(seu), .pub$cell)
if (anyNA(.m)) stop("published metadata does not cover every cell in the object.")

# keep the freshly computed embedding alongside, for comparison
.rerun <- Seurat::Embeddings(seu, "umap.harmony")
colnames(.rerun) <- c("umapharmonyrerun_1", "umapharmonyrerun_2")
seu[["umap.harmony.rerun"]] <- Seurat::CreateDimReducObject(
  embeddings = .rerun, key = "umapharmonyrerun_", assay = DefaultAssay(seu))
seu$clusters.harmony.rerun <- seu$clusters.harmony

.emb <- as.matrix(.pub[.m, c("umap_harmony_1", "umap_harmony_2")])
rownames(.emb) <- colnames(seu)
colnames(.emb) <- c("umapharmony_1", "umapharmony_2")
seu[["umap.harmony"]] <- Seurat::CreateDimReducObject(
  embeddings = .emb, key = "umapharmony_", assay = DefaultAssay(seu))

.cl <- as.character(.pub$clusters_harmony[.m])
.lv <- suppressWarnings(as.numeric(unique(.cl)))
.lv <- if (anyNA(.lv)) sort(unique(.cl)) else as.character(sort(unique(as.numeric(.cl))))
seu$clusters.harmony <- factor(.cl, levels = .lv)

seu$pbmc_annotations <- factor(.pub$pbmc_annotations[.m],
                               levels = attr(.pub, "annotation_levels"))
seu$clusters.harmony.annot.manual <- seu$pbmc_annotations

message("Restored published annotations for ", ncol(seu), " cells across ",
        nlevels(seu$pbmc_annotations), " cell types.")
print(table(seu$pbmc_annotations))
rm(.pub, .m, .emb, .rerun, .lost, .cl, .lv)

Seurat::Idents(seu) <- "pbmc_annotations"


# ---------------------------------------------------------------------------- #
# Save
saveRDS(seu, "objects/pbmc_final.RDS")
