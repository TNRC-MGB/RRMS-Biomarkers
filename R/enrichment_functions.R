
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
# ==============================================================================


calc_fora_from_ranks <- function(
    ranks,
    pathways,
    universe,
    clusters = NULL,
    cluster_name = NULL,
    score_threshold = NULL,
    direction = c("up", "down", "abs"),
    top_n = NULL,
    top_prop = NULL,
    min_genes = 10,
    pathway_min_size = 1,
    pathway_max_size = Inf,
    p_adjust_scope = c("global", "by_cluster")
) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required.")
  }

  direction <- match.arg(direction)
  p_adjust_scope <- match.arg(p_adjust_scope)

  # exactly one selection rule
  n_rules <- sum(!vapply(
    list(score_threshold, top_n, top_prop),
    is.null,
    logical(1)
  ))
  if (n_rules != 1) {
    stop("Provide exactly one of: score_threshold, top_n, or top_prop.")
  }

  # -----------------------------
  # normalize ranks input
  # -----------------------------
  if (is.numeric(ranks)) {
    if (is.null(names(ranks))) {
      stop("If 'ranks' is a numeric vector, it must be named by gene.")
    }
    if (is.null(cluster_name)) {
      cluster_name <- "cluster1"
    }
    rank_list <- setNames(list(ranks), cluster_name)

  } else if (is.list(ranks)) {
    rank_list <- ranks

    if (is.null(names(rank_list))) {
      if (!is.null(clusters) && length(clusters) == length(rank_list)) {
        names(rank_list) <- as.character(clusters)
      } else {
        stop("If 'ranks' is a list, it must be named, or 'clusters' must be provided with matching length.")
      }
    }

  } else {
    stop("'ranks' must be either a named numeric vector or a named list of named numeric vectors.")
  }

  if (!is.null(clusters)) {
    clusters <- as.character(clusters)
    missing_clusters <- setdiff(clusters, names(rank_list))
    if (length(missing_clusters) > 0) {
      stop("These clusters are missing from 'ranks': ", paste(missing_clusters, collapse = ", "))
    }
    rank_list <- rank_list[clusters]
  }

  # -----------------------------
  # universe can be:
  #   - one character vector used for all clusters
  #   - a named list of character vectors, one per cluster
  # -----------------------------
  get_universe_for_cluster <- function(cl) {
    if (is.list(universe)) {
      if (is.null(names(universe))) {
        stop("If 'universe' is a list, it must be named by cluster.")
      }
      u <- universe[[cl]]
      if (is.null(u)) {
        stop("No universe provided for cluster: ", cl)
      }
    } else {
      u <- universe
    }
    unique(as.character(u))
  }

  # helper: deduplicate genes by keeping the value with largest absolute magnitude
  dedupe_rank_vector <- function(x) {
    if (is.null(names(x))) stop("Rank vector must have gene names.")
    sp <- split(x, names(x))
    out <- vapply(sp, function(v) v[which.max(abs(v))], numeric(1))
    out
  }

  res_list <- lapply(names(rank_list), function(cl) {
    rnk <- rank_list[[cl]]

    if (!is.numeric(rnk) || is.null(names(rnk))) {
      stop("Rank vector for cluster '", cl, "' must be a named numeric vector.")
    }

    rnk <- rnk[is.finite(rnk)]
    rnk <- rnk[!is.na(names(rnk)) & names(rnk) != ""]
    if (length(rnk) == 0) return(NULL)

    if (anyDuplicated(names(rnk)) > 0) {
      rnk <- dedupe_rank_vector(rnk)
    }

    u <- get_universe_for_cluster(cl)
    rnk <- rnk[names(rnk) %in% u]
    if (length(rnk) == 0) return(NULL)

    score <- switch(
      direction,
      up   = rnk,
      down = -rnk,
      abs  = abs(rnk)
    )

    # select genes
    if (!is.null(top_n)) {
      top_n_use <- min(top_n, length(score))
      genes <- names(sort(score, decreasing = TRUE))[seq_len(top_n_use)]

    } else if (!is.null(top_prop)) {
      stopifnot(top_prop > 0 && top_prop < 1)
      k <- max(1, floor(length(score) * top_prop))
      genes <- names(sort(score, decreasing = TRUE))[seq_len(k)]

    } else if (!is.null(score_threshold)) {
      genes <- names(score[score > score_threshold])

    } else {
      stop("Provide one of: score_threshold, top_n, or top_prop.")
    }

    genes <- unique(intersect(genes, u))
    if (length(genes) < min_genes) return(NULL)

    out <- fgsea::fora(
      pathways = pathways,
      genes = genes,
      universe = u,
      minSize = pathway_min_size,
      maxSize = pathway_max_size
    )

    if (is.null(out) || nrow(out) == 0) return(NULL)

    out$cluster <- cl
    out$n_ranked <- length(rnk)
    out$n_genes_selected <- length(genes)
    out$direction <- direction

    out
  })

  res <- do.call(rbind, res_list)
  if (is.null(res) || nrow(res) == 0) return(res)

  # preserve fgsea's per-run padj and add our own
  res$padj_fora <- res$padj

  if (p_adjust_scope == "global") {
    res$padj_custom <- p.adjust(res$pval, method = "BH")
  } else {
    res$padj_custom <- ave(
      res$pval,
      res$cluster,
      FUN = function(p) p.adjust(p, method = "BH")
    )
  }

  res <- res[order(res$padj_custom, res$pval), , drop = FALSE]
  rownames(res) <- NULL
  res
}



calc_fgsea_from_ranks <- function(
    ranks,
    pathways,
    universe = NULL,
    clusters = NULL,
    cluster_name = NULL,
    direction = c("std", "up", "down"),
    pathway_min_size = 10,
    pathway_max_size = 500,
    gsea_param = 1,
    eps = 1e-50,
    sample_size = 101,
    nperm_simple = 1000,
    p_adjust_scope = c("native", "global", "by_cluster")
) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required.")
  }

  direction <- match.arg(direction)
  p_adjust_scope <- match.arg(p_adjust_scope)

  # -----------------------------
  # normalize ranks input
  # -----------------------------
  if (is.numeric(ranks)) {
    if (is.null(names(ranks))) {
      stop("If 'ranks' is a numeric vector, it must be named by gene.")
    }
    if (is.null(cluster_name)) {
      cluster_name <- "cluster1"
    }
    rank_list <- setNames(list(ranks), cluster_name)

  } else if (is.list(ranks)) {
    rank_list <- ranks

    if (is.null(names(rank_list))) {
      if (!is.null(clusters) && length(clusters) == length(rank_list)) {
        names(rank_list) <- as.character(clusters)
      } else {
        stop("If 'ranks' is a list, it must be named, or 'clusters' must be provided with matching length.")
      }
    }

  } else {
    stop("'ranks' must be either a named numeric vector or a named list of named numeric vectors.")
  }

  if (!is.null(clusters)) {
    clusters <- as.character(clusters)
    missing_clusters <- setdiff(clusters, names(rank_list))
    if (length(missing_clusters) > 0) {
      stop("These clusters are missing from 'ranks': ", paste(missing_clusters, collapse = ", "))
    }
    rank_list <- rank_list[clusters]
  }

  # -----------------------------
  # universe can be:
  #   - NULL
  #   - one character vector used for all clusters
  #   - a named list of character vectors, one per cluster
  # -----------------------------
  get_universe_for_cluster <- function(cl) {
    if (is.null(universe)) return(NULL)

    if (is.list(universe)) {
      if (is.null(names(universe))) {
        stop("If 'universe' is a list, it must be named by cluster.")
      }
      u <- universe[[cl]]
      if (is.null(u)) {
        stop("No universe provided for cluster: ", cl)
      }
    } else {
      u <- universe
    }

    unique(as.character(u))
  }

  # keep one value per duplicated gene name
  dedupe_rank_vector <- function(x) {
    sp <- split(x, names(x))
    out <- vapply(sp, function(v) v[which.max(abs(v))], numeric(1))
    out
  }

  res_list <- lapply(names(rank_list), function(cl) {
    rnk <- rank_list[[cl]]

    if (!is.numeric(rnk) || is.null(names(rnk))) {
      stop("Rank vector for cluster '", cl, "' must be a named numeric vector.")
    }

    rnk <- rnk[is.finite(rnk)]
    rnk <- rnk[!is.na(names(rnk)) & names(rnk) != ""]
    if (length(rnk) == 0) return(NULL)

    if (anyDuplicated(names(rnk)) > 0) {
      rnk <- dedupe_rank_vector(rnk)
    }

    # apply cluster-specific tested-gene universe if supplied
    u <- get_universe_for_cluster(cl)
    if (!is.null(u)) {
      rnk <- rnk[names(rnk) %in% u]
    }
    if (length(rnk) == 0) return(NULL)

    # choose score direction
    # std  = standard signed enrichment
    # up   = one-tailed positive enrichment
    # down = one-tailed negative enrichment via sign flip + pos test
    if (direction == "std") {
      stats_use <- rnk
      score_type <- "std"
      direction_label <- "two_sided_signed"
    } else if (direction == "up") {
      stats_use <- rnk
      score_type <- "pos"
      direction_label <- "higher_at_top_of_rank"
    } else if (direction == "down") {
      stats_use <- -rnk
      score_type <- "pos"
      direction_label <- "lower_at_top_of_original_rank"
    }

    # fgsea expects decreasing order
    stats_use <- sort(stats_use, decreasing = TRUE)

    # restrict pathways to genes present in stats
    pathways_use <- lapply(pathways, function(gs) intersect(gs, names(stats_use)))
    pathways_use <- pathways_use[lengths(pathways_use) >= pathway_min_size]
    if (length(pathways_use) == 0) return(NULL)

    fg <- fgsea::fgseaMultilevel(
      pathways = pathways_use,
      stats = stats_use,
      minSize = pathway_min_size,
      maxSize = pathway_max_size,
      eps = eps,
      scoreType = score_type,
      gseaParam = gsea_param,
      sampleSize = sample_size,
      nPermSimple = nperm_simple
    )

    if (is.null(fg) || nrow(fg) == 0) return(NULL)

    fg <- as.data.frame(fg)
    fg$cluster <- cl
    fg$direction <- direction_label
    fg$n_ranked <- length(stats_use)
    fg
  })

  res <- do.call(rbind, res_list)
  if (is.null(res) || nrow(res) == 0) return(res)

  # preserve fgsea native padj
  if (!"padj" %in% colnames(res)) {
    res$padj <- p.adjust(res$pval, method = "BH")
  }
  res$padj_native <- res$padj

  # optional extra adjustment across all clusters/pathways
  if (p_adjust_scope == "native") {
    res$padj_custom <- res$padj_native
  } else if (p_adjust_scope == "global") {
    res$padj_custom <- p.adjust(res$pval, method = "BH")
  } else {
    res$padj_custom <- ave(
      res$pval,
      res$cluster,
      FUN = function(p) p.adjust(p, method = "BH")
    )
  }

  res <- res[order(res$padj_custom, res$pval), , drop = FALSE]
  rownames(res) <- NULL
  res
}


# ---------------------------------------------------------------------------- #
# NEBULA output; ranked gene vectors / gene universe, per B-cell subset
# Genes are ranked by the NEBULA Wald statistic (logFC / SE), descending, which
# is the ranking used for every EBV enrichment analysis in the manuscript.

build_nebula_ranks <- function(res, cluster_col = "B_annotation") {
  split(res, res[[cluster_col]]) |>
    lapply(function(df) {
      x <- df %>%
        dplyr::filter(!is.na(gene), gene != "",
                      is.finite(logFC), is.finite(SE), SE > 0) %>%
        dplyr::mutate(rank_stat = logFC / SE) %>%
        dplyr::distinct(gene, .keep_all = TRUE) %>%
        dplyr::arrange(dplyr::desc(rank_stat))
      stats::setNames(x$rank_stat, x$gene)
    })
}

build_nebula_universe <- function(res, cluster_col = "B_annotation") {
  split(res, res[[cluster_col]]) |>
    lapply(function(df) {
      df %>%
        dplyr::filter(!is.na(gene), gene != "",
                      is.finite(logFC), is.finite(SE), is.finite(PValue)) %>%
        dplyr::pull(gene) %>%
        unique()
    })
}


# ---------------------------------------------------------------------------- #
# EBV lifecycle stages
# Maps each EBV factor to its stage in the viral lifecycle. Used to color and
# order the gene sets in Figure 5d, Extended Data Fig. 9 and Fig. 10.

latent <- c(
  "EBNA-3A", "EBNA-3B/EBNA-3C", "LMP-2B", "A73", "RPMS1",
  "LMP-1", "LMP-2A", "EBNA-1", "EBNA-2", "EBNA-LP"
)

early <- c(
  "BALF5", "BHRF1", "BRLF1", "BRRF1", "BDLF4", "BXLF1", "BGLF3",
  "BGLF3.5", "BNLF2a", "BNLF2b", "BZLF1", "BSRF1", "BKRF3", "BKRF4", "BGLF4",
  "BGLF5", "BBLF2/BBLF3", "BBLF4", "BSLF1", "BLLF3", "BaRF1", "BMRF1",
  "BALF2", "BFRF2", "BORF2", "BALF1", "BARF1", "BLLF2", "BMRF2",
  "BVLF1", "BVRF1", "BHLF1"
)

leaky_late <- c(
  "BXLF2", "BKRF2", "BRRF2", "BLLF1", "BBRF2", "BBRF3", "BDLF3",
  "BLRF1", "BLRF2", "BdRF1", "BBRF1", "BcRF1", "BALF3", "BALF4", "BBLF1"
)

true_late <- c(
  "BTRF1", "BZLF2", "BNRF1", "BILF2", "BGLF2", "BcLF1",
  "BDLF1", "BDLF2", "BVRF2", "BXRF1", "BOLF1", "BPLF1", "BFRF1",
  "BFRF3", "BFRF1A", "BGLF1"
)

other <- c(
  "BFLF1", "BFLF2", "BGRF1/BDRF1", "BSLF2/BMLF1", "EBER1", "EBER2"
)

custom_pathway_order <- c(
  other,
  latent,
  early,
  leaky_late,
  true_late
)


ebv_lifecycle_stage <- function(x) {
  dplyr::case_when(
    x %in% latent     ~ "latent",
    x %in% early      ~ "early",
    x %in% leaky_late ~ "leaky late",
    x %in% true_late  ~ "true late",
    x %in% other      ~ "other",
    TRUE              ~ NA_character_
  )
}

ebv_stage_colors <- c(
  "latent"     = "#E41A1C",
  "early"      = "#377EB8",
  "leaky late" = "#4DAF4A",
  "true late"  = "#984EA3",
  "other"      = "#FF7F00"
)

ebv_stage_labels <- c(
  "latent"     = "Latent",
  "early"      = "Early",
  "leaky late" = "Leaky late",
  "true late"  = "True late",
  "other"      = "Unknown"
)


# ---------------------------------------------------------------------------- #
# scDist-based over-representation, as used for Figure 5d
#
# `dat` is an scDist result object: dat$gene.names, and
# dat$vals[[cluster]]$beta.hat.
calc_fora <- function(dat,
                      clusters,
                      pathways,
                      universe,
                      dist.threshold = NULL,
                      direction = c("up", "down", "abs"),
                      top_n = NULL,
                      top_prop = NULL,
                      min_genes = 10,
                      beta_dim = 1,
                      p_adjust_scope = c("global", "by_cluster")) {

  direction      <- match.arg(direction)
  p_adjust_scope <- match.arg(p_adjust_scope)

  stopifnot(!is.null(dat$gene.names), !is.null(universe))
  universe <- unique(as.character(universe))

  res_list <- lapply(as.character(clusters), function(cl) {
    bh <- dat$vals[[cl]]$beta.hat
    if (is.null(bh)) return(NULL)
    if (beta_dim > ncol(bh)) stop("beta_dim exceeds ncol(beta.hat) for cluster ", cl)

    rnk <- bh[, beta_dim]
    names(rnk) <- dat$gene.names
    rnk <- rnk[names(rnk) %in% universe]
    if (length(rnk) == 0) return(NULL)

    score <- switch(direction, up = rnk, down = -rnk, abs = abs(rnk))

    if (!is.null(top_n)) {
      genes <- names(sort(score, decreasing = TRUE))[seq_len(min(top_n, length(score)))]
    } else if (!is.null(top_prop)) {
      stopifnot(top_prop > 0 && top_prop < 1)
      k <- max(1, floor(length(score) * top_prop))
      genes <- names(sort(score, decreasing = TRUE))[seq_len(k)]
    } else if (!is.null(dist.threshold)) {
      genes <- names(score[score > dist.threshold])
    } else {
      stop("Provide one of: dist.threshold, top_n, or top_prop.")
    }

    genes <- intersect(unique(genes), universe)
    if (length(genes) < min_genes) return(NULL)

    out <- fgsea::fora(pathways = pathways, genes = genes, universe = universe)
    if (is.null(out) || nrow(out) == 0) return(NULL)

    out$cell_type <- cl
    out$n_genes   <- length(genes)
    out$direction <- direction
    out
  })

  res <- do.call(rbind, res_list)
  if (is.null(res) || nrow(res) == 0) return(res)

  if (p_adjust_scope == "global") {
    res$padj <- p.adjust(res$pval, method = "BH")
  } else {
    res$padj <- ave(res$pval, res$cell_type, FUN = function(p) p.adjust(p, method = "BH"))
  }
  res
}


plot_fora_cluster <- function(fora.res, cluster = "12", padj_thresh = 0.05,
                              top_k = 30) {
  df <- fora.res %>%
    dplyr::filter(cell_type == cluster) %>%
    dplyr::mutate(neglog10 = -log10(padj), sig = padj < padj_thresh) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_k)
  df$pathway <- factor(df$pathway, levels = df$pathway[order(df$neglog10)])

  ggplot(df, aes(x = neglog10, y = pathway)) +
    geom_point(aes(size = overlap, color = sig), alpha = 0.9) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "gray70")) +
    labs(x = expression(-log[10]("BH-FDR")), y = NULL,
         size = "Overlap", color = paste0("FDR < ", padj_thresh),
         title = paste0("ORA (FORA) - cluster ", cluster)) +
    theme_classic(base_size = 12)
}
