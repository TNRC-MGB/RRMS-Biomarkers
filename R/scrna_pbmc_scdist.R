
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# scDist analysis; whole-PBMC contrasts
#   (scdist_pbmc_pre/scdist_pbmc_rel/scdist_pbmc_ms_hc)
#
# Outputs: scdist_pbmc_prerelapse_v_remission.RDS /
#          scdist_pbmc_prerelapse_paired.RDS
#            PreRelapse vs Remission
#          scdist_pbmc_relapse_v_remission.RDS / scdist_pbmc_relapse_paired.RDS
#            Relapse vs Remission
#          scdist_pbmc_ms_v_healthy.RDS                        MS vs Healthy
#          healthy_split_null.RDS     within-Healthy split null   (ED Fig. 4a,c)
#          remission_split_null.RDS     within-Remission split null (ED Fig. 4b)
#          ms_condition_perm_control.RDS    clinical-label permutation control
# ==============================================================================


# --- environment --------------------------------------------------------------
setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

source("R/scdist_functions.R")

  plan(sequential)
  options(future.globals.maxSize = Inf)





  ### Import PBMC object
  pbmc <- readRDS("objects/pbmc_final.RDS")



  # ### Define background genes (those with detectable expression that could conceivably be detected as distGenes)
  # # Genes detected in at least 1% of cells in at least 20% of patients
  # DefaultAssay(pbmc) <- "RNA"
  # counts <- GetAssayData(pbmc, assay='RNA',layer='counts')
  # thr_cells <- 0.01
  # thr_patients <- 0.20
  # det_by_patient <- sapply(split(seq_len(ncol(counts)), pbmc$Patient), function(ix) {
  #   Matrix::rowMeans(counts[, ix, drop = FALSE] > 0)
  # })
  # frac_patients_passing <- rowMeans(det_by_patient >= thr_cells)
  # pbmc_background_genes <- rownames(counts)[frac_patients_passing >= thr_patients]
  #
  # #saveRDS(pbmc_background_genes,'pbmc_background_genes.RDS')
  # pbmc_background_genes <- readRDS('pbmc_background_genes.RDS')



  ### SCTransform (since scDist author suggests using SCTransform scale.data for scDist)
  #DefaultAssay(pbmc) <- 'RNA'
  #pbmc <- DietSeurat(pbmc,assays='RNA',layers='counts')
  #pbmc <- SCTransform(pbmc, return.only.var.genes = FALSE)


  cluster_map <- pbmc[[]] |>
    dplyr::distinct(clusters.harmony, pbmc_annotations) |>
    dplyr::arrange(clusters.harmony)

  cluster_to_annot <- setNames(
    paste(cluster_map$clusters.harmony,cluster_map$pbmc_annotations,sep="."),
    as.character(cluster_map$clusters.harmony)
  )


  ### 1) PreRelapse vs Remission: all available samples
  pre_all_filt <- prefilter_clusters_for_scdist(
    seu = pbmc,
    condition_col = "Condition",
    patient_col = "Patient",
    cluster_col = "pbmc_annotations",
    levels_use = c("Remission", "PreRelapse"),
    min_cells_per_condition = 50,
    min_patients_per_condition = 2
  )

  scdist_pbmc_pre <- run_scdist_standard(
    seu = pre_all_filt$seu,
    fixed_effect = "Condition",
    reference_level = "Remission",
    random_effect = "Patient",
    cluster_col = "pbmc_annotations",
    assay = "RNA",
    levels_use = c("Remission", "PreRelapse"),
    d = 20
  )

  saveRDS(scdist_pbmc_pre, file.path(SCDIST_DIR, "scdist_pbmc_prerelapse_v_remission.RDS"))


  ### 2) PreRelapse vs Remission: paired-only samples
  paired_patients_pre <- pbmc[[]] |>
    dplyr::filter(Condition %in% c("Remission", "PreRelapse")) |>
    dplyr::count(Patient, Condition) |>
    dplyr::distinct(Patient, Condition) |>
    dplyr::count(Patient) |>
    dplyr::filter(n == 2) |>
    dplyr::pull(Patient)

  paired_pre <- subset(pbmc, subset = Patient %in% paired_patients_pre)

  paired_pre_filt <- prefilter_clusters_for_scdist(
    seu = paired_pre,
    condition_col = "Condition",
    patient_col = "Patient",
    cluster_col = "pbmc_annotations",
    levels_use = c("Remission", "PreRelapse"),
    min_cells_per_condition = 50,
    min_patients_per_condition = 2
  )

  scdist_pbmc_pre_paired <- run_scdist_standard(
    seu = paired_pre_filt$seu,
    fixed_effect = "Condition",
    reference_level = "Remission",
    random_effect = "Patient",
    cluster_col = "pbmc_annotations",
    assay = "RNA",
    levels_use = c("Remission", "PreRelapse"),
    d = 20
  )

  saveRDS(scdist_pbmc_pre_paired, file.path(SCDIST_DIR, "scdist_pbmc_prerelapse_paired.RDS"))







  ### 1) Relapse vs Remission: all available samples
  rel_all_filt <- prefilter_clusters_for_scdist(
    seu = pbmc,
    condition_col = "Condition",
    patient_col = "Patient",
    cluster_col = "pbmc_annotations",
    levels_use = c("Remission", "Relapse"),
    min_cells_per_condition = 50,
    min_patients_per_condition = 2
  )

  scdist_pbmc_rel <- run_scdist_standard(
    seu = rel_all_filt$seu,
    fixed_effect = "Condition",
    reference_level = "Remission",
    random_effect = "Patient",
    cluster_col = "pbmc_annotations",
    assay = "RNA",
    levels_use = c("Remission", "Relapse"),
    d = 20
  )

  saveRDS(scdist_pbmc_rel, file.path(SCDIST_DIR, "scdist_pbmc_relapse_v_remission.RDS"))


  ### 2) Relapse vs Remission: paired-only samples
  paired_patients_rel <- pbmc[[]] %>%
    dplyr::filter(Condition %in% c("Remission", "Relapse")) %>%
    dplyr::count(Patient, Condition) %>%
    dplyr::distinct(Patient, Condition) %>%
    dplyr::count(Patient) %>%
    dplyr::filter(n == 2) %>%
    dplyr::pull(Patient)

  paired_rel <- subset(pbmc, subset = Patient %in% paired_patients_rel)

  paired_rel_filt <- prefilter_clusters_for_scdist(
    seu = paired_rel,
    condition_col = "Condition",
    patient_col = "Patient",
    cluster_col = "pbmc_annotations",
    levels_use = c("Remission", "Relapse"),
    min_cells_per_condition = 50,
    min_patients_per_condition = 2
  )

  scdist_pbmc_rel_paired <- run_scdist_standard(
    seu = paired_rel_filt$seu,
    fixed_effect = "Condition",
    reference_level = "Remission",
    random_effect = "Patient",
    cluster_col = "pbmc_annotations",
    assay = "RNA",
    levels_use = c("Remission", "Relapse"),
    d = 20
  )

  saveRDS(scdist_pbmc_rel_paired, file.path(SCDIST_DIR, "scdist_pbmc_relapse_paired.RDS"))





  ### MS v healthy
  pbmc$Disease <- ifelse(pbmc$Condition=="Healthy","Healthy","MS")
  ms_hc_filt <- prefilter_clusters_for_scdist(
    seu = pbmc,
    condition_col = "Disease",
    patient_col = "Patient",
    cluster_col = "pbmc_annotations",
    levels_use = c("Healthy", "MS"),
    min_cells_per_condition = 50,
    min_patients_per_condition = 2
  )

  scdist_pbmc_ms_hc <- run_scdist_standard(
    seu = ms_hc_filt$seu,
    fixed_effect = "Disease",
    reference_level = "Healthy",
    random_effect = "Patient",
    cluster_col = "pbmc_annotations",
    assay = "RNA",
    levels_use = c("Healthy", "MS"),
    d = 20
  )

  saveRDS(scdist_pbmc_ms_hc,file.path(SCDIST_DIR, "scdist_pbmc_ms_v_healthy.RDS"))




  ### Add scdist scores to seurat object
  pbmc <- add_scdist_to_seurat(pbmc,scdist_pbmc_ms_hc,cluster_col="pbmc_annotations",new_col="scdist.ms_v_healthy")
  pbmc <- add_scdist_to_seurat(pbmc,scdist_pbmc_rel,cluster_col="pbmc_annotations",new_col="scdist.relapse")
  pbmc <- add_scdist_to_seurat(pbmc,scdist_pbmc_pre,cluster_col="pbmc_annotations",new_col="scdist.pre_relapse")

  pbmc <- add_scdist_to_seurat(pbmc,scdist_pbmc_rel_paired,cluster_col="pbmc_annotations",new_col="scdist.relapse_paired")
  pbmc <- add_scdist_to_seurat(pbmc,scdist_pbmc_pre_paired,cluster_col="pbmc_annotations",new_col="scdist.pre_relapse_paired")



  ### Controls (null distributions)
  # healthy v healthy split (not paired within MS patients)
  # remission v remission split (also not paired)
  # patient-label permutation in the MS paired data (paired)
  #
  # Null distributions are generated by repeated random partitioning of healthy
  # or remission samples within condition, and by permutation of clinical-state
  # labels among MS samples, while preserving patient-level aggregation. The
  # observed cell-type scDist distances are then compared against these nulls to
  # test whether relapse-associated signal exceeds within-state heterogeneity
  # and label-randomized expectation. These objects are what Extended Data
  # Fig. 4 is drawn from.
  #
  # WARNING: each call runs scDist 100 times over the full atlas. This takes
  # hours ...

  source("R/scdist_null_controls.R")   # run_scdist_split_null,
                                       # run_scdist_ms_condition_permutation_control


  ### 1) Healthy vs healthy split null
  pbmc_healthy_split_null <- run_scdist_split_null(
    seu = pbmc,
    condition_value = "Healthy",
    n_iter = 100,
    sample_col = "Sample_ID",
    patient_col = "Patient",
    condition_col = "Condition",
    cluster_col = "pbmc_annotations",
    pseudo_col = "pseudo_group",
    levels_use = c("group1", "group2"),
    reference_level = "group1",
    assay = "RNA",
    d = 20,
    min_cells_per_condition = 50,
    min_patients_per_condition = 2,
    seed = 42
  )

  saveRDS(pbmc_healthy_split_null, "objects/healthy_split_null.RDS")


  ### 2) MS remission vs remission split null
  pbmc_remission_split_null <- run_scdist_split_null(
    seu = pbmc,
    condition_value = "Remission",
    n_iter = 100,
    sample_col = "Sample_ID",
    patient_col = "Patient",
    condition_col = "Condition",
    cluster_col = "pbmc_annotations",
    pseudo_col = "pseudo_group",
    levels_use = c("group1", "group2"),
    reference_level = "group1",
    assay = "RNA",
    d = 20,
    min_cells_per_condition = 50,
    min_patients_per_condition = 2,
    seed = 43
  )

  saveRDS(pbmc_remission_split_null, "objects/remission_split_null.RDS")


  ### 3) Clinical-label permutation control within MS samples
  pbmc_ms_condition_perm_control <- run_scdist_ms_condition_permutation_control(
    seu = pbmc,
    ms_levels = c("Remission", "PreRelapse", "Relapse"),
    contrasts = list(
      list(name = "pre", levels_use = c("Remission", "PreRelapse"), reference_level = "Remission"),
      list(name = "rel", levels_use = c("Remission", "Relapse"),    reference_level = "Remission")
    ),
    n_iter = 100,
    sample_col = "Sample_ID",
    patient_col = "Patient",
    condition_col = "Condition",
    perm_col = "Condition_perm",
    cluster_col = "pbmc_annotations",
    assay = "RNA",
    d = 20,
    min_cells_per_condition = 50,
    min_patients_per_condition = 2,
    seed = 44
  )

  saveRDS(pbmc_ms_condition_perm_control, "objects/ms_condition_perm_control.RDS")
