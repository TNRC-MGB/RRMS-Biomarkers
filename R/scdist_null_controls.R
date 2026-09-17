
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Null distributions for the scDist cell-type distances, Extended Data Fig. 4
#
# Two controls, each fitting scDist n_iter times over the full atlas:
#
#   run_scdist_split_null      randomly partitions the samples of ONE condition
#                              into two pseudo-groups and fits the contrast, so
#                              the null is within-state heterogeneity
#   run_scdist_ms_condition_permutation_control
#                              permutes clinical-state labels among MS samples
#                              and refits every contrast
#
# ==============================================================================


make_sample_map <- function(seu,
                            sample_col = "Sample_ID",
                            patient_col = "Patient",
                            condition_col = "Condition") {
  seu[[]] %>%
    dplyr::select(all_of(c(sample_col, patient_col, condition_col))) %>%
    dplyr::distinct()
}

attach_sample_level_metadata <- function(seu,
                                         label_df,
                                         sample_col = "Sample_ID",
                                         new_col) {
  md <- seu[[]]
  stopifnot(sample_col %in% colnames(label_df))
  stopifnot(new_col %in% colnames(label_df))

  md2 <- md %>%
    dplyr::left_join(label_df[, c(sample_col, new_col)], by = sample_col)

  seu[[new_col]] <- md2[[new_col]]
  seu
}


# ---------------------------------------------------------------------------- #
# Permutation 
# Which half receives which label is itself randomized, so the pseudo-groups
# carry no residual ordering from the sample list.
make_sample_split_permutations <- function(samples,
                                           n_iter = 100,
                                           seed = 1,
                                           labels = c("group1", "group2")) {
  set.seed(seed)
  samples <- sort(unique(samples))
  n <- length(samples)
  n1 <- floor(n / 2)

  replicate(n_iter, {
    shuffled <- sample(samples)

    split1 <- shuffled[1:n1]
    split2 <- shuffled[(n1 + 1):n]

    if (sample(c(TRUE, FALSE), 1)) {
      group1 <- split1
      group2 <- split2
    } else {
      group1 <- split2
      group2 <- split1
    }

    list(
      split_table = data.frame(
        Sample_ID = samples,
        pseudo_group = NA_character_,
        stringsAsFactors = FALSE
      ) |>
        dplyr::mutate(
          pseudo_group = dplyr::case_when(
            Sample_ID %in% group1 ~ labels[1],
            Sample_ID %in% group2 ~ labels[2],
            TRUE ~ NA_character_
          )
        ),
      group1 = group1,
      group2 = group2
    )
  }, simplify = FALSE)
}

make_condition_permutations <- function(sample_map,
                                        n_iter = 100,
                                        sample_col = "Sample_ID",
                                        condition_col = "Condition",
                                        seed = 1,
                                        perm_col = "Condition_perm") {
  set.seed(seed)

  samples <- sample_map[[sample_col]]
  labels  <- sample_map[[condition_col]]

  replicate(n_iter, {
    data.frame(
      Sample_ID = samples,
      Condition_perm = sample(labels, length(labels), replace = FALSE),
      stringsAsFactors = FALSE
    )
  }, simplify = FALSE)
}


# ---------------------------------------------------------------------------- #
# One prefilter-and-fit iteration

run_scdist_one <- function(seu,
                           condition_col,
                           patient_col = "Patient",
                           cluster_col = "pbmc_annotations",
                           levels_use,
                           reference_level,
                           assay = "RNA",
                           d = 20,
                           min_cells_per_condition = 50,
                           min_patients_per_condition = 2) {
  filt <- prefilter_clusters_for_scdist(
    seu = seu,
    condition_col = condition_col,
    patient_col = patient_col,
    cluster_col = cluster_col,
    levels_use = levels_use,
    min_cells_per_condition = min_cells_per_condition,
    min_patients_per_condition = min_patients_per_condition
  )

  fit <- run_scdist_standard(
    seu = filt$seu,
    fixed_effect = condition_col,
    reference_level = reference_level,
    random_effect = patient_col,
    cluster_col = cluster_col,
    assay = assay,
    levels_use = levels_use,
    d = d
  )

  rm(filt)
  gc()

  fit
}


# ---------------------------------------------------------------------------- #
# Within-condition split null

# A failed iteration is captured rather than thrown, so one degenerate split
# does not discard the other 99 fits. Can inspect $results[[i]]$scdist for
# condition inherits "error".
run_scdist_split_null <- function(seu,
                                  condition_value,
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
                                  seed = 1) {
  md <- seu[[]]
  cells_keep <- rownames(md)[md[[condition_col]] == condition_value]
  seu_sub <- subset(seu, cells = cells_keep)

  sample_map <- make_sample_map(
    seu_sub,
    sample_col = sample_col,
    patient_col = patient_col,
    condition_col = condition_col
  )

  permutations <- make_sample_split_permutations(
    samples = sample_map[[sample_col]],
    n_iter = n_iter,
    seed = seed,
    labels = levels_use
  )

  results <- lapply(seq_along(permutations), function(i) {
    perm_i <- permutations[[i]]

    seu_i <- attach_sample_level_metadata(
      seu = seu_sub,
      label_df = perm_i$split_table,
      sample_col = sample_col,
      new_col = pseudo_col
    )

    fit_i <- tryCatch(
      run_scdist_one(
        seu = seu_i,
        condition_col = pseudo_col,
        patient_col = patient_col,
        cluster_col = cluster_col,
        levels_use = levels_use,
        reference_level = reference_level,
        assay = assay,
        d = d,
        min_cells_per_condition = min_cells_per_condition,
        min_patients_per_condition = min_patients_per_condition
      ),
      error = function(e) e
    )

    rm(seu_i)
    gc()

    list(
      run_id = i,
      null_type = paste0(tolower(condition_value), "_split"),
      condition_value = condition_value,
      permutation = perm_i,
      scdist = fit_i
    )
  })

  rm(seu_sub)
  gc()

  list(
    sample_map = sample_map,
    permutations = permutations,
    results = results
  )
}


# ---------------------------------------------------------------------------- #
# MS clinical-state permutation control

run_scdist_ms_condition_permutation_control <- function(
    seu,
    ms_levels = c("Remission", "PreRelapse", "Relapse"),
    contrasts = list(
      list(name = "pre", levels_use = c("Remission", "PreRelapse"),
           reference_level = "Remission"),
      list(name = "rel", levels_use = c("Remission", "Relapse"),
           reference_level = "Remission")
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
    seed = 1) {
  md <- seu[[]]
  cells_keep <- rownames(md)[md[[condition_col]] %in% ms_levels]
  seu_ms <- subset(seu, cells = cells_keep)

  sample_map <- make_sample_map(
    seu_ms,
    sample_col = sample_col,
    patient_col = patient_col,
    condition_col = condition_col
  ) %>%
    dplyr::filter(.data[[condition_col]] %in% ms_levels)

  permutations <- make_condition_permutations(
    sample_map = sample_map,
    n_iter = n_iter,
    sample_col = sample_col,
    condition_col = condition_col,
    seed = seed,
    perm_col = perm_col
  )

  results <- lapply(seq_along(permutations), function(i) {
    perm_i <- permutations[[i]]

    seu_i <- attach_sample_level_metadata(
      seu = seu_ms,
      label_df = perm_i,
      sample_col = sample_col,
      new_col = perm_col
    )

    contrast_results <- lapply(contrasts, function(ct) {
      fit_ct <- tryCatch(
        run_scdist_one(
          seu = seu_i,
          condition_col = perm_col,
          patient_col = patient_col,
          cluster_col = cluster_col,
          levels_use = ct$levels_use,
          reference_level = ct$reference_level,
          assay = assay,
          d = d,
          min_cells_per_condition = min_cells_per_condition,
          min_patients_per_condition = min_patients_per_condition
        ),
        error = function(e) e
      )

      list(
        contrast_name = ct$name,
        levels_use = ct$levels_use,
        reference_level = ct$reference_level,
        scdist = fit_ct
      )
    })

    rm(seu_i)
    gc()

    list(
      run_id = i,
      null_type = "ms_condition_permutation",
      permutation = perm_i,
      contrast_results = contrast_results
    )
  })

  rm(seu_ms)
  gc()

  list(
    sample_map = sample_map,
    permutations = permutations,
    results = results
  )
}
