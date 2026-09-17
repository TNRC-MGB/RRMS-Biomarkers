
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# scRNA preprocessing + plotting helper functions
# ==============================================================================

.n <- function(x) message(paste(x, collapse = "\n"))




# (could expand these qc features in future)
add_basic_qc <- function(seu, mt_pattern = "^MT-") {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)
  seu
}




qc_filter_per_GEM_well <- function(seu,
                                   split_by = "Pool_ID",
                                   min_features = 300,
                                   min_counts = 500,
                                   max_percent_mt = 10,
                                   nmads = 2,
                                   outlier_type = c("lower", "both"),
                                   log_outliers = TRUE,
                                   verbose = TRUE) {
  outlier_type <- match.arg(outlier_type)

  md <- seu[[]]
  if (!(split_by %in% colnames(md))) {
    stop("split_by column not found in meta.data: ", split_by)
  }
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% colnames(md))) {
    stop("Expected nCount_RNA and nFeature_RNA in meta.data. Run Seurat QC metrics first.")
  }
  if (!("percent.mt" %in% colnames(md))) {
    stop("Expected percent.mt in meta.data. Compute it with PercentageFeatureSet(...).")
  }
  if (!requireNamespace("scater", quietly = TRUE)) {
    stop("Package 'scater' not installed. Install via Bioconductor: BiocManager::install('scater')")
  }

  obj_list <- SplitObject(seu, split.by = split_by)

  obj_list <- lapply(names(obj_list), function(k) {
    if (verbose) message(split_by, "=", k)
    i <- obj_list[[k]]

    # Hard thresholds (per lane)
    i <- subset(i, subset = percent.mt < max_percent_mt)
    i <- subset(i, subset = nFeature_RNA > min_features)
    i <- subset(i, subset = nCount_RNA > min_counts)

    md_i <- i[[]]
    idx1 <- scater::isOutlier(md_i[["nCount_RNA"]],   nmads = nmads, log = log_outliers, type = outlier_type)
    idx2 <- scater::isOutlier(md_i[["nFeature_RNA"]], nmads = nmads, log = log_outliers, type = outlier_type)
    keep <- !(idx1 | idx2)

    if (verbose) message("  kept ", sum(keep), " / ", length(keep))
    i <- i[, keep]   # <- no drop= for Seurat objects
    i
  })

  # Merge back together
  seu2 <- merge(obj_list[[1]], y = obj_list[-1])
  seu2 <- JoinLayers(seu2)
  seu2
}





run_scDblFinder_per_GEM_well <- function(seu,
                                         split_by = "Pool_ID",
                                         assay = "RNA",
                                         counts_layer = "counts",     # Seurat v5 layer name
                                         fallback_layer = "data",
                                         dbl_rate = NULL,             # optional expected doublet rate (0-1); if NULL, scDblFinder estimates
                                         score_col = "scDblFinder_score",
                                         call_col  = "scDblFinder_call",
                                         class_col = "scDblFinder_class",
                                         filter = FALSE,
                                         seed = 1,
                                         verbose = TRUE) {
  if (!(split_by %in% colnames(seu[[]]))) {
    stop("split_by column not found in meta.data: ", split_by)
  }
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop("Package 'scDblFinder' not installed. Install via Bioconductor: BiocManager::install('scDblFinder')")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' not installed. Install via Bioconductor.")
  }

  set.seed(seed)

  obj_list <- SplitObject(seu, split.by = split_by)

  obj_list <- lapply(names(obj_list), function(k) {
    if (verbose) message(split_by, "=", k)
    i <- obj_list[[k]]

    # counts matrix (prefer counts layer)
    mat <- tryCatch(
      GetAssayData(i, assay = assay, layer = counts_layer),
      error = function(e) NULL
    )
    if (is.null(mat)) {
      if (verbose) message("  counts layer not found; using fallback layer='", fallback_layer, "'")
      mat <- GetAssayData(i, assay = assay, layer = fallback_layer)
    }

    # Build SCE
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))

    # Run scDblFinder (per lane)
    if (is.null(dbl_rate)) {
      sce <- scDblFinder::scDblFinder(sce)
    } else {
      sce <- scDblFinder::scDblFinder(sce, dbr = dbl_rate)
    }

    # Pull results
    cd <- SingleCellExperiment::colData(sce)
    if (!all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(cd))) {
      stop("scDblFinder did not return expected columns. Found: ",
           paste(colnames(cd), collapse = ", "))
    }

    i[[score_col]] <- as.numeric(cd$scDblFinder.score)
    i[[class_col]] <- as.character(cd$scDblFinder.class)
    i[[call_col]]  <- i@meta.data[[class_col]] == "doublet"

    if (verbose) {
      n_call <- sum(i@meta.data[[call_col]], na.rm = TRUE)
      message("  called doublets=", n_call, " / ", ncol(i))
    }

    if (filter) {
      keep <- !i@meta.data[[call_col]]
      i <- i[, keep]   # <- no drop= in Seurat
    }

    i
  })

  seu2 <- merge(obj_list[[1]], y = obj_list[-1])
  seu2 <- JoinLayers(seu2)
  seu2
}



process_doublets <- function(seu,
                             GEM_well = "Pool_ID",
                             vireo_donor_col = "donor_id",
                             vireo_prob_col  = "prob_doublet",
                             vireo_prob_thresh = 0.5,
                             scdbl_class_col = "scDblFinder_class",
                             scdbl_call_col  = "scDblFinder_call",
                             scdbl_doublet_label = "doublet",
                             na_pool_label = "UNKNOWN",
                             filter = TRUE,
                             verbose = TRUE) {

  stopifnot(inherits(seu, "Seurat"))
  md <- seu[[]]

  # ---- column checks ----
  required <- c(GEM_well, vireo_donor_col, vireo_prob_col, scdbl_class_col, scdbl_call_col)
  missing <- required[!required %in% colnames(md)]
  if (length(missing) > 0) {
    stop("Missing required metadata columns: ", paste(missing, collapse = ", "))
  }

  # ---- pool bookkeeping ----
  md[[GEM_well]] <- as.character(md[[GEM_well]])
  md[[GEM_well]][is.na(md[[GEM_well]]) | md[[GEM_well]] == ""] <- na_pool_label

  # ---- vireo doublet call ----
  v_donor <- md[[vireo_donor_col]]
  v_prob  <- md[[vireo_prob_col]]
  md$vireo_doublet <- (v_donor %in% "doublet") |
    (!is.na(v_prob) & (v_prob >= vireo_prob_thresh))

  # ---- scDblFinder doublet call ----
  s_class <- md[[scdbl_class_col]]
  s_call  <- md[[scdbl_call_col]]

  s_call_logical <- rep(FALSE, nrow(md))
  if (is.logical(s_call)) {
    s_call_logical <- (!is.na(s_call) & s_call)
  } else if (is.numeric(s_call)) {
    s_call_logical <- (!is.na(s_call) & s_call == 1)
  } else {
    s_call_chr <- as.character(s_call)
    s_call_logical <- (!is.na(s_call_chr) & tolower(s_call_chr) %in% c("true","t","1","yes","y"))
  }

  md$scDbl_doublet <- s_call_logical |
    (!is.na(s_class) & as.character(s_class) == scdbl_doublet_label)

  # ---- union + reason ----
  md$doublet_union <- md$vireo_doublet | md$scDbl_doublet
  md$doublet_reason <- dplyr::case_when(
    md$vireo_doublet & md$scDbl_doublet ~ "both",
    md$vireo_doublet & !md$scDbl_doublet ~ "vireo_only",
    !md$vireo_doublet & md$scDbl_doublet ~ "scDblFinder_only",
    TRUE ~ "neither"
  )

  # write metadata back (use SeuratObject to avoid export issues)
  seu <- SeuratObject::AddMetaData(
    seu,
    md[, c(GEM_well, "vireo_doublet", "scDbl_doublet", "doublet_union", "doublet_reason"), drop = FALSE]
  )

  # ---- optional filtering ----
  seu2 <- seu
  if (filter) {
    n_drop <- sum(seu[[]]$doublet_union, na.rm = TRUE)
    if (verbose) message("Dropping union doublets: ", n_drop, " cells.")
    seu2 <- subset(seu, subset = doublet_union == FALSE)
  } else {
    if (verbose) message("filter=FALSE; returning unfiltered Seurat object (metadata added).")
  }

  seu2
}





run_seurat_standard <- function(seu,
                                split_by = "Batch",
                                npcs = 50,
                                dims = 1:30,
                                resolution = 2,
                                pca.name = "pca.first_pass",
                                umap.name = "umap.first_pass",
                                cluster.name = "clusters.first_pass",
                                seed = 1,
                                verbose = FALSE) {

  stopifnot(inherits(seu, "Seurat"))
  md <- seu[[]]
  if (!split_by %in% colnames(md)) stop("split_by column not in metadata: ", split_by)

  # split RNA assay into layers by technical batch
  seu[["RNA"]] <- split(seu[["RNA"]], f = md[[split_by]])

  set.seed(seed)
  seu <- NormalizeData(seu, verbose = verbose)
  seu <- FindVariableFeatures(seu, verbose = verbose)
  seu <- ScaleData(seu, verbose = verbose)

  # join before dimensionality reduction / graph construction
  seu <- JoinLayers(seu)

  seu <- RunPCA(seu, npcs = npcs, reduction.name = pca.name, verbose = verbose)
  seu <- FindNeighbors(seu, reduction = pca.name, dims = dims, verbose = verbose)

  set.seed(seed)
  seu <- FindClusters(seu, resolution = resolution, cluster.name = cluster.name, verbose = verbose)

  set.seed(seed)
  seu <- RunUMAP(seu, reduction = pca.name, dims = dims, reduction.name = umap.name, verbose = verbose)

  return(seu)
}



#' Cluster QC panel for PBMC scRNA-seq (DotPlot + cluster size + median QC)
#'
#' @param seu        Seurat object
#' @param cluster_col Metadata column with cluster IDs (default "clusters.first_pass")
#' @param reduction  Reduction for DimPlot (default "umap.first_pass")
#' @param font_sizeS Numeric; passed to FontSize() if you use that helper
#' @param features   Character vector of genes/metadata features for DotPlot
#' @param show_umap  Logical; include a labeled DimPlot in output
#' @param label_size Numeric; label size for DimPlot
#'
#' @return list(panel=patchwork, dotplot=ggplot, sizeplot=ggplot, qcplot=ggplot,
#'              umap=ggplot_or_NULL, clust_qc=data.frame)
plot_cluster_qc_panel <- function(
    seu,
    cluster_col = "clusters.first_pass",
    reduction   = "umap.first_pass",
    font_sizeS  = 10,
    features = c(
      # T
      "TRAC","CD3D","CD3E","IL7R","CCR7",
      # NK/cytotoxic
      "NKG7","GNLY","PRF1","GZMB","FCER1G",
      # B
      "CD19","MS4A1","CD79A","BANK1",
      # Plasma
      "MZB1","JCHAIN","XBP1","IGKC",
      # Mono
      "CD14","S100A8","S100A9","LYZ","FCGR3A",
      # RBC/PLT
      "HBB","HBA1","HBA2","PPBP","PF4","ITGA2B",
      # Endothelial/vascular contamination
      "COL14A1","ANGPT1","EGFL7","IGFBP7",
      # Possible junk and/or dying
      "MALAT1","JUN",
      # Doublet score (metadata feature is OK in DotPlot)
      "scDblFinder_score"
    ),
    show_umap  = TRUE,
    rasterize_dotplot = FALSE,
    label_size = 4
) {
  stopifnot(inherits(seu, "Seurat"))
  if (!(cluster_col %in% colnames(seu[[]]))) {
    stop("cluster_col not found in meta.data: ", cluster_col)
  }

  # set identities
  Seurat::Idents(seu) <- cluster_col

  # UMAP (optional)
  p_umap <- NULL
  if (isTRUE(show_umap)) {
    p_umap <- Seurat::DimPlot(seu, reduction = reduction, label = TRUE, label.size = label_size) +
      Seurat::NoLegend() +
      ggplot2::coord_fixed()
    # If you have FontSize() helper in your environment, use it; otherwise ignore.
    if (exists("FontSize", mode = "function")) {
      p_umap <- p_umap + FontSize(x.text = font_sizeS, y.text = font_sizeS,
                                  x.title = font_sizeS, y.title = font_sizeS)
    }
  }

  # meta.data with explicit cluster column
  md <- seu[[]]
  md$cluster <- md[[cluster_col]]

  # 1) Cluster sizes + median QC
  clust_qc <- dplyr::as_tibble(md) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarize(
      n_cells      = dplyr::n(),
      med_nFeature = stats::median(.data$nFeature_RNA, na.rm = TRUE),
      med_nCount   = stats::median(.data$nCount_RNA,   na.rm = TRUE),
      med_pct_mt   = stats::median(.data$percent.mt,   na.rm = TRUE),
      .groups = "drop"
    )

  # 2) Lineage/QC marker dot plot
  if(rasterize_dotplot) {
    p_dot <- rasterize(Seurat::DotPlot(seu, features = features) +
      Seurat::RotatedAxis() +
      ggplot2::scale_color_viridis_c() +
      ggplot2::theme(legend.position = "left"),dpi=600)
  } else {
    p_dot <- Seurat::DotPlot(seu, features = features) +
      Seurat::RotatedAxis() +
      ggplot2::scale_color_viridis_c() +
      ggplot2::theme(legend.position = "left")
  }

  if (exists("FontSize", mode = "function")) {
    p_dot <- p_dot + FontSize(x.text = font_sizeS, y.text = font_sizeS, x.title = 0, y.title = 0)
  }

  # 3) Cluster size bar plot
  p_size <- ggplot2::ggplot(clust_qc, ggplot2::aes(x = stats::reorder(.data$cluster, -.data$n_cells),
                                                   y = .data$n_cells)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title = ggplot2::element_blank())
  if (exists("FontSize", mode = "function")) {
    p_size <- p_size + FontSize(x.text = font_sizeS, y.text = font_sizeS, x.title = 0, y.title = 0)
  }

  # 4) Median QC scatter (faceted)
  p_qc <- clust_qc %>%
    tidyr::pivot_longer(cols = c("med_nFeature", "med_nCount", "med_pct_mt"),
                        names_to = "metric", values_to = "value") %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$cluster, y = .data$value)) +
    ggplot2::geom_point() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.title = ggplot2::element_blank())
  if (exists("FontSize", mode = "function")) {
    p_qc <- p_qc + FontSize(x.text = font_sizeS, y.text = font_sizeS, x.title = 0, y.title = 0)
  }

  # Combined panel
  panel <- p_dot + p_size + p_qc + patchwork::plot_layout(widths = c(3, 0.5, 1))

  # return everything for flexible reuse
  list(
    panel    = panel,
    dotplot  = p_dot,
    sizeplot = p_size,
    qcplot   = p_qc,
    umap     = p_umap,
    clust_qc = clust_qc
  )
}

# Example:
# res <- plot_cluster_qc_panel(seu, cluster_col="clusters.first_pass", reduction="umap.first_pass", font_sizeS=10)
# res$panel
# res$umap
# res$clust_qc








# RunAzimuth.Seurat <- function (query, reference, query.modality = "RNA", annotation.levels = NULL,
#           umap.name = "ref.umap", do.adt = FALSE, verbose = TRUE, assay = NULL,
#           k.weight = 50, n.trees = 20, mapping.score.k = 100, ...)
# {
#   CheckDots(...)
#   assay <- assay %||% DefaultAssay(query)
#   if (query.modality == "ATAC") {
#     query <- RunAzimuthATAC(query = query, reference = reference,
#                             annotation.levels = annotation.levels, umap.name = umap.name,
#                             verbose = verbose, assay = assay, k.weight = k.weight,
#                             n.trees = n.trees, mapping.score.k = mapping.score.k,
#                             ...)
#   }
#   else {
#     if (dir.exists(reference)) {
#       reference <- LoadReference(reference)$map
#     }
#     else {
#       reference <- tolower(reference)
#       if (reference %in% InstalledData()$Dataset) {
#         reference <- LoadData(reference, type = "azimuth")$map
#       }
#       else if (reference %in% AvailableData()$Dataset) {
#         InstallData(reference)
#         reference <- LoadData(reference, type = "azimuth")$map
#       }
#       else {
#         possible.references <- AvailableData()$Dataset[grepl("*ref",
#                                                              AvailableData()$Dataset)]
#         print("Choose one of:")
#         print(possible.references)
#         stop(paste("Could not find a reference for",
#                    reference))
#       }
#       if (!"num_precomputed_nns" %in% names(Misc(reference[["refUMAP"]])$model)) {
#         Misc(reference[["refUMAP"]], slot = "model")$num_precomputed_nns <- 1
#       }
#       key.pattern = "^[^_]*_"
#       new.colnames <- gsub(pattern = key.pattern, replacement = Key(reference[["refDR"]]),
#                            x = colnames(Loadings(object = reference[["refDR"]],
#                                                  projected = FALSE)))
#       colnames(Loadings(object = reference[["refDR"]],
#                         projected = FALSE)) <- new.colnames
#     }
#     dims <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
#     if (isTRUE(do.adt) && !("ADT" %in% Assays(reference))) {
#       warning("Cannot impute an ADT assay because the reference does not have antibody data")
#       do.adt = FALSE
#     }
#     reference.version <- ReferenceVersion(reference)
#     azimuth.version <- as.character(packageVersion(pkg = "Azimuth"))
#     seurat.version <- as.character(packageVersion(pkg = "Seurat"))
#     meta.data <- names(slot(reference, "meta.data"))
#     if (is.null(annotation.levels)) {
#       annotation.levels <- names(slot(object = reference,
#                                       name = "meta.data"))
#       annotation.levels <- annotation.levels[!grepl(pattern = "^nCount",
#                                                     x = annotation.levels)]
#       annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature",
#                                                     x = annotation.levels)]
#       annotation.levels <- annotation.levels[!grepl(pattern = "^ori",
#                                                     x = annotation.levels)]
#     }
#     query <- ConvertGeneNames(object = query, reference.names = rownames(x = reference),
#                               homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds")
#     if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
#       calcn <- as.data.frame(x = Seurat:::CalcN(object = query[[assay]]))
#       colnames(x = calcn) <- paste(colnames(x = calcn),
#                                    assay, sep = "_")
#       query <- AddMetaData(object = query, metadata = calcn)
#       rm(calcn)
#     }
#     if (any(grepl(pattern = "^MT-", x = rownames(x = query)))) {
#       query <- PercentageFeatureSet(object = query, pattern = "^MT-",
#                                     col.name = "percent.mt", assay = assay)
#     }
#     anchors <- FindTransferAnchors(reference = reference,
#                                    query = query, k.filter = NA, reference.neighbors = "refdr.annoy.neighbors",
#                                    reference.assay = "refAssay", query.assay = assay,
#                                    reference.reduction = "refDR", normalization.method = "SCT",
#                                    features = rownames(Loadings(reference[["refDR"]])),
#                                    dims = 1:dims, n.trees = n.trees, mapping.score.k = mapping.score.k,
#                                    verbose = verbose)
#     refdata <- lapply(X = annotation.levels, function(x) {
#       reference[[x, drop = TRUE]]
#     })
#     names(x = refdata) <- annotation.levels
#     if (isTRUE(do.adt)) {
#       refdata[["impADT"]] <- GetAssayData(object = reference[["ADT"]],
#                                           slot = "data")
#     }
#     query <- TransferData(reference = reference, query = query,
#                           query.assay = assay, dims = 1:dims, anchorset = anchors,
#                           refdata = refdata, n.trees = 20, store.weights = TRUE,
#                           k.weight = k.weight, verbose = verbose)
#     query <- IntegrateEmbeddings(anchorset = anchors, reference = reference,
#                                  query = query, query.assay = assay, reductions = "pcaproject",
#                                  reuse.weights.matrix = TRUE, verbose = verbose)
#     query[["query_ref.nn"]] <- FindNeighbors(object = Embeddings(reference[["refDR"]]),
#                                              query = Embeddings(query[["integrated_dr"]]), return.neighbor = TRUE,
#                                              l2.norm = TRUE, verbose = verbose)
#     query <- NNTransform(object = query, meta.data = reference[[]])
#     query[[umap.name]] <- RunUMAP(object = query[["query_ref.nn"]],
#                                   reduction.model = reference[["refUMAP"]], reduction.key = "UMAP_",
#                                   verbose = verbose)
#     query <- AddMetaData(object = query, metadata = MappingScore(anchors = anchors,
#                                                                  ndim = dims), col.name = "mapping.score")
#   }
#   return(query)
# }





RunAzimuth.Seurat.v5 <- function(
    query,
    reference,
    query.modality = "RNA",
    annotation.levels = NULL,
    umap.name = "ref.umap",
    do.adt = FALSE,
    verbose = TRUE,
    assay = NULL,
    k.weight = 50,
    n.trees = 20,
    mapping.score.k = 100,
    ...
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Package 'Azimuth' is required.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  assay <- assay %||% Seurat::DefaultAssay(query)

  if (!inherits(query, "Seurat")) {
    stop("'query' must be a Seurat object")
  }
  if (!assay %in% Seurat::Assays(query)) {
    stop("Assay '", assay, "' not found in query")
  }

  if (identical(query.modality, "ATAC")) {
    query <- Azimuth:::RunAzimuthATAC(
      query = query,
      reference = reference,
      annotation.levels = annotation.levels,
      umap.name = umap.name,
      verbose = verbose,
      assay = assay,
      k.weight = k.weight,
      n.trees = n.trees,
      mapping.score.k = mapping.score.k,
      ...
    )
    return(query)
  }

  # ---------------------------------------------------------------------------
  # Load / resolve reference
  # ---------------------------------------------------------------------------
  if (dir.exists(reference)) {
    reference <- Azimuth:::LoadReference(reference)$map
  } else {
    if (!requireNamespace("SeuratData", quietly = TRUE)) {
      stop(
        "Package 'SeuratData' is required when 'reference' is not a local directory. ",
        "Install it with install.packages('SeuratData') or remotes::install_github('satijalab/seurat-data')."
      )
    }

    reference <- tolower(reference)

    installed.tbl <- SeuratData::InstalledData()
    available.tbl <- SeuratData::AvailableData()

    if (reference %in% installed.tbl$Dataset) {
      reference <- SeuratData::LoadData(reference, type = "azimuth")$map
    } else if (reference %in% available.tbl$Dataset) {
      SeuratData::InstallData(reference)
      reference <- SeuratData::LoadData(reference, type = "azimuth")$map
    } else {
      possible.references <- available.tbl$Dataset[
        grepl("ref", available.tbl$Dataset)
      ]
      print("Choose one of:")
      print(possible.references)
      stop("Could not find a reference for ", reference)
    }
  }

  if (!inherits(reference, "Seurat")) {
    stop("Resolved reference is not a Seurat object")
  }

  # ---------------------------------------------------------------------------
  # Validate required reference components
  # ---------------------------------------------------------------------------
  required.reductions <- c("refDR", "refUMAP")
  missing.reductions <- setdiff(required.reductions, names(reference@reductions))
  if (length(missing.reductions) > 0) {
    stop(
      "Reference is missing required reductions: ",
      paste(missing.reductions, collapse = ", ")
    )
  }

  if (!"refdr.annoy.neighbors" %in% names(reference@neighbors)) {
    stop("Reference is missing neighbor object 'refdr.annoy.neighbors'")
  }

  if (!"refAssay" %in% Seurat::Assays(reference)) {
    stop("Reference is missing assay 'refAssay'")
  }

  # ---------------------------------------------------------------------------
  # Ensure stored UMAP model has expected metadata
  # ---------------------------------------------------------------------------
  ref.umap.misc <- Seurat::Misc(reference[["refUMAP"]])
  if (is.null(ref.umap.misc$model$num_precomputed_nns)) {
    ref.umap.misc$model$num_precomputed_nns <- 1
    Seurat::Misc(reference[["refUMAP"]]) <- ref.umap.misc
  }

  # ---------------------------------------------------------------------------
  # Harmonize refDR loadings names with reduction key
  # ---------------------------------------------------------------------------
  refdr.key <- Seurat::Key(reference[["refDR"]])
  ref.loadings <- Seurat::Loadings(reference[["refDR"]], projected = FALSE)
  colnames(ref.loadings) <- gsub(
    pattern = "^[^_]*_",
    replacement = refdr.key,
    x = colnames(ref.loadings)
  )
  reference[["refDR"]]@feature.loadings <- ref.loadings

  # ---------------------------------------------------------------------------
  # Dimensionality from reference NN index
  # ---------------------------------------------------------------------------
  dims <- as.integer(reference@neighbors$refdr.annoy.neighbors@alg.info$ndim)
  if (length(dims) != 1 || is.na(dims) || dims < 1) {
    stop("Could not determine dimensionality from reference neighbor object")
  }

  # ---------------------------------------------------------------------------
  # Optional ADT transfer
  # ---------------------------------------------------------------------------
  if (isTRUE(do.adt) && !("ADT" %in% Seurat::Assays(reference))) {
    warning("Cannot impute an ADT assay because the reference does not have antibody data")
    do.adt <- FALSE
  }

  # ---------------------------------------------------------------------------
  # Version capture
  # ---------------------------------------------------------------------------
  reference.version <- Azimuth:::ReferenceVersion(reference)
  azimuth.version <- as.character(utils::packageVersion("Azimuth"))
  seurat.version <- as.character(utils::packageVersion("Seurat"))

  # ---------------------------------------------------------------------------
  # Annotation levels
  # ---------------------------------------------------------------------------
  if (is.null(annotation.levels)) {
    annotation.levels <- colnames(reference[[]])
    annotation.levels <- annotation.levels[!grepl("^nCount", annotation.levels)]
    annotation.levels <- annotation.levels[!grepl("^nFeature", annotation.levels)]
    annotation.levels <- annotation.levels[!grepl("^ori", annotation.levels)]
  }

  missing.annotation.levels <- setdiff(annotation.levels, colnames(reference[[]]))
  if (length(missing.annotation.levels) > 0) {
    stop(
      "Requested annotation.levels not found in reference metadata: ",
      paste(missing.annotation.levels, collapse = ", ")
    )
  }

  # ---------------------------------------------------------------------------
  # Match query gene names to reference
  # ---------------------------------------------------------------------------
  query <- Azimuth:::ConvertGeneNames(
    object = query,
    reference.names = rownames(reference),
    homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds"
  )

  # ---------------------------------------------------------------------------
  # Compute nCount / nFeature from counts layer if missing
  # ---------------------------------------------------------------------------
  needed.meta <- c(paste0("nCount_", assay), paste0("nFeature_", assay))
  if (!all(needed.meta %in% colnames(query[[]]))) {
    counts.mat <- Seurat::GetAssayData(
      object = query,
      assay = assay,
      layer = "counts"
    )

    calcn <- data.frame(
      row.names = colnames(counts.mat),
      nCount = Matrix::colSums(counts.mat),
      nFeature = Matrix::colSums(counts.mat > 0)
    )
    colnames(calcn) <- paste0(colnames(calcn), "_", assay)

    query <- Seurat::AddMetaData(
      object = query,
      metadata = calcn
    )
  }

  # ---------------------------------------------------------------------------
  # Mito percentage
  # ---------------------------------------------------------------------------
  if (any(grepl("^MT-", rownames(query))) && !"percent.mt" %in% colnames(query[[]])) {
    query <- Seurat::PercentageFeatureSet(
      object = query,
      pattern = "^MT-",
      col.name = "percent.mt",
      assay = assay
    )
  }

  # ---------------------------------------------------------------------------
  # Features used for anchor mapping
  # ---------------------------------------------------------------------------
  ref.features <- rownames(
    Seurat::Loadings(reference[["refDR"]], projected = FALSE)
  )

  if (length(ref.features) == 0) {
    stop("No features found in reference refDR loadings")
  }

  # ---------------------------------------------------------------------------
  # Transfer anchors
  # ---------------------------------------------------------------------------
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = query,
    reference.assay = "refAssay",
    query.assay = assay,
    reference.reduction = "refDR",
    reference.neighbors = "refdr.annoy.neighbors",
    normalization.method = "SCT",
    features = ref.features,
    dims = seq_len(dims),
    k.filter = NA,
    n.trees = n.trees,
    mapping.score.k = mapping.score.k,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Metadata / ADT transfer targets
  # ---------------------------------------------------------------------------
  refdata <- stats::setNames(
    lapply(annotation.levels, function(x) reference[[x, drop = TRUE]]),
    annotation.levels
  )

  if (isTRUE(do.adt)) {
    refdata[["impADT"]] <- Seurat::GetAssayData(
      object = reference,
      assay = "ADT",
      layer = "data"
    )
  }

  # ---------------------------------------------------------------------------
  # Transfer labels / values
  # ---------------------------------------------------------------------------
  query <- Seurat::TransferData(
    reference = reference,
    query = query,
    anchorset = anchors,
    refdata = refdata,
    query.assay = assay,
    dims = seq_len(dims),
    n.trees = n.trees,
    k.weight = k.weight,
    store.weights = TRUE,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Integrate projected embeddings
  # ---------------------------------------------------------------------------
  query <- Seurat::IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference,
    query = query,
    query.assay = assay,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    verbose = verbose
  )

  if (!"integrated_dr" %in% names(query@reductions)) {
    stop("Expected reduction 'integrated_dr' was not created")
  }

  # ---------------------------------------------------------------------------
  # Query-to-reference nearest neighbors
  # ---------------------------------------------------------------------------
  query[["query_ref.nn"]] <- Seurat::FindNeighbors(
    object = Seurat::Embeddings(reference[["refDR"]]),
    query = Seurat::Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE,
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Transform neighbor-derived metadata
  # ---------------------------------------------------------------------------
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = reference[[]]
  )

  # ---------------------------------------------------------------------------
  # Project into reference UMAP
  # ---------------------------------------------------------------------------
  query[[umap.name]] <- Seurat::RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = reference[["refUMAP"]],
    reduction.key = "UMAP_",
    verbose = verbose
  )

  # ---------------------------------------------------------------------------
  # Mapping score
  # ---------------------------------------------------------------------------
  warning("Skipping mapping.score because this Azimuth build does not provide MappingScore()")

  # ---------------------------------------------------------------------------
  # run information
  # ---------------------------------------------------------------------------
  run.info <- list(
    reference.version = reference.version,
    azimuth.version = azimuth.version,
    seurat.version = seurat.version,
    assay = assay,
    dims = dims,
    annotation.levels = annotation.levels,
    do.adt = do.adt,
    k.weight = k.weight,
    n.trees = n.trees,
    mapping.score.k = mapping.score.k
  )
  Seurat::Misc(query, slot = "azimuth_run") <- run.info

  return(query)
}






assign_cluster_labels_from_azimuth <- function(
    seu,
    cluster_col = "clusters.harmony",
    azimuth_col = "predicted.celltype.l2",
    threshold = 0.5,
    new_col = "clusters.harmony.annot"
) {
  stopifnot(inherits(seu, "Seurat"))

  md <- seu[[]]

  if (!cluster_col %in% colnames(md)) {
    stop("cluster_col not found in seu metadata: ", cluster_col)
  }
  if (!azimuth_col %in% colnames(md)) {
    stop("azimuth_col not found in seu metadata: ", azimuth_col)
  }

  df <- data.frame(
    cluster = md[[cluster_col]],
    azimuth = md[[azimuth_col]],
    stringsAsFactors = FALSE
  )

  df <- df[!is.na(df$cluster) & !is.na(df$azimuth), , drop = FALSE]

  # counts: rows = clusters, cols = azimuth labels
  tab <- table(df$cluster, df$azimuth)

  # fractions within each cluster
  frac_tab <- prop.table(tab, margin = 1)

  # dominant label and fraction per cluster
  cluster_assignments <- lapply(seq_len(nrow(frac_tab)), function(i) {
    clust <- rownames(frac_tab)[i]
    vals <- frac_tab[i, ]
    max_idx <- which.max(vals)
    top_label <- colnames(frac_tab)[max_idx]
    top_frac <- as.numeric(vals[max_idx])

    assigned_label <- if (top_frac >= threshold) top_label else as.character(clust)

    data.frame(
      cluster = as.character(clust),
      top_azimuth_label = top_label,
      top_fraction = top_frac,
      assigned_label = assigned_label,
      stringsAsFactors = FALSE
    )
  })

  cluster_assignments <- do.call(rbind, cluster_assignments)

  # map assigned cluster label back to every cell
  label_map <- setNames(
    cluster_assignments$assigned_label,
    cluster_assignments$cluster
  )

  seu[[new_col]] <- unname(label_map[as.character(seu[[cluster_col]][, 1])])

  return(list(
    seu = seu,
    summary = cluster_assignments,
    count_table = as.data.frame.matrix(tab),
    fraction_table = as.data.frame.matrix(frac_tab)
  ))
}




make_impadt_heatmap_df <- function(seu, markers, group_col = "pbmc_annotations",
                                   assay = "impADT", layer = "data") {
  stopifnot(all(markers %in% rownames(seu[[assay]])))
  stopifnot(group_col %in% colnames(seu[[]]))

  avg_list <- Seurat::AverageExpression(
    object = seu,
    assays = assay,
    features = markers,
    group.by = group_col,
    layer = layer,
    return.seurat = FALSE,
    verbose = FALSE
  )

  # matrix is features x groups
  mat <- avg_list[[assay]]

  # reorder groups to match factor levels present in the object
  group_levels <- levels(seu[[group_col]][, 1])
  group_levels <- intersect(group_levels, colnames(mat))
  mat <- mat[, group_levels, drop = FALSE]

  # transpose: rows = cell types, cols = markers
  mat_plot <- t(mat)

  # z-score each marker across cell types
  mat_plot_z <- scale(mat_plot)
  mat_plot_z[is.na(mat_plot_z)] <- 0

  df_heat <- as.data.frame(mat_plot_z) %>%
    tibble::rownames_to_column(var = "pbmc_annotations") %>%
    tidyr::pivot_longer(
      cols = -pbmc_annotations,
      names_to = "marker",
      values_to = "z"
    )

  df_heat$pbmc_annotations <- factor(df_heat$pbmc_annotations, levels = group_levels)
  df_heat$marker <- factor(df_heat$marker, levels = markers)

  return(df_heat)
}




plot_impadt_heatmap <- function(df_heat, title = NULL, option="magma") {
  ggplot(df_heat, aes(x = marker, y = pbmc_annotations, fill = z)) +
    geom_tile(color = NA) +
    # scale_fill_gradient2(
    #   low = "#3B4CC0",
    #   mid = "white",
    #   high = "#B40426",
    #   midpoint = 0,
    #   name = "Scaled\nexpression"
    # ) +
    scale_fill_viridis_c(option = option) +
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = "plain"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}
