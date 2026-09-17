
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# scDist helpers: cluster prefiltering, the standard fit, and the distance plot
#
#   prefilter_clusters_for_scdist   drop clusters without enough support
#   run_scdist_standard             fit scDist on a Seurat object
#   add_scdist_to_seurat            map per-cluster distances onto cell metadata
#   DistPlot2                       distance and CI plot, colored by adjusted P
#
# Sourced by figure2.R, figure3.R, extended_data_4.R, extended_data_5.R,
# preprocessing_bcells.R and scrna_pbmc_scdist.R. The permutation and
# split-null controls that build on these live in R/scdist_null_controls.R.
#
# ==============================================================================


# ---------------------------------------------------------------------------- #
# Cluster prefiltering
# Requires every retained cluster to be above cell and patient thresholds in
# EVERY requested condition, so a cluster present in only one cond. is dropped
# rather than silently fitted on unbalanced support.
prefilter_clusters_for_scdist <- function(
    seu,
    condition_col = "Condition",
    patient_col = "Patient",
    cluster_col = "pbmc_annotations",
    levels_use = NULL,
    min_cells_per_condition = 50,
    min_patients_per_condition = 2,
    verbose = TRUE
) {
  stopifnot(inherits(seu, "Seurat"))

  req <- c(condition_col, patient_col, cluster_col)
  if (!all(req %in% colnames(seu[[]]))) {
    missing <- setdiff(req, colnames(seu[[]]))
    stop("Missing required metadata columns: ", paste(missing, collapse = ", "))
  }

  md <- seu[[]]
  md$.cell <- rownames(md)

  if (is.null(levels_use)) {
    levels_use <- unique(as.character(md[[condition_col]]))
    levels_use <- levels_use[!is.na(levels_use)]
  }

  md <- md[!is.na(md[[condition_col]]) &
             md[[condition_col]] %in% levels_use, , drop = FALSE]

  if (nrow(md) == 0) {
    stop("No cells remain after filtering to levels_use.")
  }

  cluster_condition_summary <- md %>%
    dplyr::count(
      .data[[cluster_col]],
      .data[[condition_col]],
      .data[[patient_col]],
      name = "n_cells_patient"
    ) %>%
    dplyr::group_by(.data[[cluster_col]], .data[[condition_col]]) %>%
    dplyr::summarize(
      n_cells = sum(n_cells_patient),
      n_patients = dplyr::n_distinct(.data[[patient_col]]),
      .groups = "drop"
    )

  keep_clusters <- cluster_condition_summary %>%
    dplyr::mutate(
      pass = n_cells >= min_cells_per_condition &
        n_patients >= min_patients_per_condition
    ) %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::summarize(
      n_conditions_present = dplyr::n_distinct(.data[[condition_col]]),
      all_conditions_pass = all(pass),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      n_conditions_present == length(levels_use),
      all_conditions_pass
    ) %>%
    dplyr::pull(.data[[cluster_col]]) %>%
    as.character()

  dropped_clusters <- setdiff(unique(as.character(md[[cluster_col]])),
                              keep_clusters)

  if (length(keep_clusters) == 0) {
    stop(
      "No clusters passed filtering with min_cells_per_condition = ",
      min_cells_per_condition,
      " and min_patients_per_condition = ",
      min_patients_per_condition,
      "."
    )
  }

  cells_keep <- md$.cell[as.character(md[[cluster_col]]) %in% keep_clusters]
  seu_filt <- subset(seu, cells = cells_keep)

  cluster_filter_table <- cluster_condition_summary %>%
    dplyr::mutate(
      keep_cluster = as.character(.data[[cluster_col]]) %in% keep_clusters
    ) %>%
    dplyr::arrange(.data[[cluster_col]], .data[[condition_col]])

  if (verbose) {
    message("Conditions retained: ", paste(levels_use, collapse = ", "))
    message("Clusters retained: ", length(keep_clusters))
    message("Clusters dropped: ", length(dropped_clusters))
    if (length(dropped_clusters) > 0) {
      message("Dropped cluster IDs: ", paste(dropped_clusters, collapse = ", "))
    }
    message("Cells retained: ", ncol(seu_filt), " / ", ncol(seu))
  }

  list(
    seu = seu_filt,
    keep_clusters = keep_clusters,
    dropped_clusters = dropped_clusters,
    summary = cluster_filter_table
  )
}


# ---------------------------------------------------------------------------- #
run_scdist_standard <- function(
    seu,
    fixed_effect,
    reference_level,
    random_effect,
    cluster_col,
    assay = "RNA",
    levels_use = NULL,
    genes = NULL,
    seed = 1,
    verbose = TRUE,
    ...
) {
  req <- unique(c(fixed_effect, random_effect, cluster_col))
  if (!all(req %in% colnames(seu[[]]))) {
    missing <- setdiff(req, colnames(seu[[]]))
    stop("Missing required metadata columns: ", paste(missing, collapse = ", "))
  }

  if (!assay %in% Seurat::Assays(seu)) {
    stop("Assay '", assay, "' not found in object.")
  }

  set.seed(seed)

  md0 <- seu[[]]

  if (is.null(levels_use)) {
    levels_use <- unique(as.character(md0[[fixed_effect]]))
    levels_use <- levels_use[!is.na(levels_use)]
  }

  keep <- md0[[fixed_effect]] %in% levels_use
  keep[is.na(keep)] <- FALSE
  cells_keep <- rownames(md0)[keep]

  if (length(cells_keep) == 0) {
    stop(
      "No cells remain after filtering ", fixed_effect, " to levels: ",
      paste(levels_use, collapse = ", ")
    )
  }

  seu_sub <- subset(seu, cells = cells_keep)
  Seurat::DefaultAssay(seu_sub) <- assay

  data_mat <- try(
    Seurat::GetAssayData(seu_sub, assay = assay, layer = "data"),
    silent = TRUE
  )
  if (inherits(data_mat, "try-error")) {
    stop(
      "Failed to retrieve layer='data' from assay '", assay, "'. ",
      "Make sure NormalizeData() has been run and the data layer exists."
    )
  }

  if (is.null(genes)) {
    genes_use <- rownames(data_mat)
  } else {
    genes_use <- intersect(genes, rownames(data_mat))
    if (length(genes_use) == 0) {
      stop(
        "After intersecting with assay='", assay,
        "', layer='data', genes length is 0."
      )
    }
  }

  Y <- as.matrix(data_mat[genes_use, , drop = FALSE])

  meta_cols <- unique(c(fixed_effect, random_effect, cluster_col))
  md <- seu_sub[[meta_cols]]

  fe_vals <- as.character(md[[fixed_effect]])
  fe_levels_present <- unique(fe_vals[!is.na(fe_vals)])

  if (!reference_level %in% fe_levels_present) {
    stop(
      "reference_level '", reference_level,
      "' is not present in the subsetted data for fixed_effect '",
      fixed_effect, "'."
    )
  }

  other_levels <- setdiff(fe_levels_present, reference_level)
  md[[fixed_effect]] <- factor(
    fe_vals,
    levels = c(reference_level, other_levels)
  )

  md[[random_effect]] <- factor(md[[random_effect]])
  md[[cluster_col]] <- factor(md[[cluster_col]])

  stopifnot(identical(colnames(Y), rownames(md)))

  if (verbose) {
    message("Assay used: ", assay)
    message("Layer used: data")
    message("Fixed effect: ", fixed_effect)
    message("Reference level: ", reference_level)
    message("Levels used: ", paste(levels(md[[fixed_effect]]), collapse = ", "))
    message("Random effect: ", random_effect)
    message("Cluster column: ", cluster_col)
    message("Group counts:")
    print(table(md[[fixed_effect]]))
    message("N random-effect levels: ", length(unique(md[[random_effect]])))
    message("N clusters: ", length(unique(md[[cluster_col]])))
    message("N genes used: ", nrow(Y))
  }

  scDist::scDist(
    normalized_counts = Y,
    meta.data = md,
    fixed.effects = fixed_effect,
    random.effects = random_effect,
    clusters = cluster_col,
    ...
  )
}


# ---------------------------------------------------------------------------- #
# Per-cluster distance as a cell-level metadata column

add_scdist_to_seurat <- function(
    seu,
    scdist_out,
    cluster_col = "clusters.harmony",
    value_col = "Dist.",
    new_col = "scdist"
) {
  stopifnot(inherits(seu, "Seurat"))

  if (!cluster_col %in% colnames(seu[[]])) {
    stop("cluster_col not found in seu metadata: ", cluster_col)
  }

  if (is.null(scdist_out$results)) {
    stop("scdist_out does not contain a $results element.")
  }

  res <- scdist_out$results

  if (!value_col %in% colnames(res)) {
    stop("value_col '", value_col, "' not found in scdist_out$results.")
  }

  if (is.null(rownames(res))) {
    stop("scdist_out$results must have row names corresponding to cluster IDs.")
  }

  dist_map <- setNames(res[[value_col]], rownames(res))

  seu[[new_col]] <- unname(dist_map[as.character(seu[[cluster_col]][, 1])])

  seu
}


# ---------------------------------------------------------------------------- #
# Distance plot
# Adapted from the scDist source code
DistPlot2 <- function(
    scd.object,
    color_limits = c(0, 3),
    point_size = 3,
    stroke_size = 2,
    p.adjust_method = "BH",
    annotation_map = NULL,
    factor_levels = NULL
) {
  results <- scd.object$results

  dist_col <- dplyr::coalesce(
    if ("Dist." %in% colnames(results)) "Dist." else NA_character_,
    if ("D.post.med" %in% colnames(results)) "D.post.med" else NA_character_,
    if ("D.hat" %in% colnames(results)) "D.hat" else NA_character_
  )
  if (is.na(dist_col)) {
    stop("Could not find a distance column (tried Dist., D.post.med, D.hat).")
  }

  ci_up_col <- dplyr::coalesce(
    if ("95% CI (upper)" %in% colnames(results)) "95% CI (upper)" else NA_character_,
    if ("D.post.ub" %in% colnames(results)) "D.post.ub" else NA_character_
  )
  ci_lo_col <- dplyr::coalesce(
    if ("95% CI (low)" %in% colnames(results)) "95% CI (low)" else NA_character_,
    if ("D.post.lb" %in% colnames(results)) "D.post.lb" else NA_character_
  )
  if (is.na(ci_up_col) || is.na(ci_lo_col)) {
    stop("Could not find CI columns (tried `95% CI (upper)/(low)` and ",
         "`D.post.ub/lb`).")
  }

  p_col <- dplyr::coalesce(
    if ("p.val" %in% colnames(results)) "p.val" else NA_character_,
    if ("p.sum" %in% colnames(results)) "p.sum" else NA_character_
  )
  if (is.na(p_col)) {
    stop("Could not find p-value column (tried p.val, p.sum).")
  }

  df <- data.frame(
    cell_type = rownames(results),
    dist      = results[[dist_col]],
    se_up     = results[[ci_up_col]],
    se_down   = results[[ci_lo_col]],
    p_value   = results[[p_col]],
    stringsAsFactors = FALSE
  )

  df$adj_p_val <- p.adjust(df$p_value, method = p.adjust_method)

  if (is.null(factor_levels)) {
    df$cell_type <- factor(df$cell_type, levels = df$cell_type)
  } else {
    df$cell_type <- factor(df$cell_type, levels = factor_levels)
  }

  if (!is.null(annotation_map)) {
    df$cell_type <- annotation_map[as.character(df$cell_type)]
  }

  ggplot2::ggplot(df, ggplot2::aes(x = reorder(cell_type, dist), y = dist)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = se_down, ymax = se_up),
      width = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = -log10(adj_p_val)),
      shape = 21,
      size = point_size,
      stroke = stroke_size,
      fill = "lightgray"
    ) +
    ggplot2::xlab("Cell type") +
    ggplot2::ylab(dist_col) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    paletteer::scale_color_paletteer_c(
      "ggthemes::Red-Black Diverging",
      direction = -1,
      limits = color_limits,
      oob = scales::squish
    )
}
