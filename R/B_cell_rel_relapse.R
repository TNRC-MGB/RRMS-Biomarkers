
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# B-cell EBNA2 gene-set enrichment - Relapse vs Remission 
# ==============================================================================

rm(list = ls())

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

set.seed(42)

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

# The EBNA2 and control gene sets ship in data/genesets/; see the README for
# the ChIP-seq series they were derived from.
EBNA2_DIR <- "data/genesets"
OBJ_DIR   <- "objects"
OUT_DIR   <- "Intermediate"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# Gene universe
universe_genes <- read.table(file.path(EBNA2_DIR, "genes_universe_bcells_1pct.txt"))[,1]

# scdist results
scdist_bcell_rel <- readRDS(file.path(SCDIST_DIR, "scdist_bcell_relapse_v_remission.RDS"))

# read GMT
pathways_raw <- fgsea::gmtPathways(file.path(EBNA2_DIR, "ebna2_gene_sets.gmt"))


# ---------------------------------------------------------------------------- #
# helpers
clean_symbols <- function(x, force_upper = TRUE) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  if (force_upper) x <- toupper(x)
  unique(x)
}

collapse_duplicate_stats <- function(df) {
  # keep one row per gene using largest absolute statistic
  df %>%
    group_by(gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
}

prepare_pathways <- function(pathways, universe_genes,
                             force_upper = TRUE,
                             min_geneset_size = 10,
                             max_geneset_size = 2000) {
  universe_genes <- clean_symbols(universe_genes, force_upper = force_upper)

  pathways2 <- lapply(pathways, function(gs) {
    gs <- clean_symbols(gs, force_upper = force_upper)
    intersect(gs, universe_genes)
  })

  pathways2 <- pathways2[
    lengths(pathways2) >= min_geneset_size &
      lengths(pathways2) <= max_geneset_size
  ]

  pathways2
}

get_ranked_stats_from_scdist <- function(obj, celltype, universe_genes = NULL,
                                         force_upper = TRUE) {
  stopifnot(celltype %in% names(obj$vals))
  stopifnot("beta.hat" %in% names(obj$vals[[celltype]]))

  bh <- obj$vals[[celltype]]$beta.hat
  stat <- if (is.matrix(bh) || is.data.frame(bh)) bh[, 1] else bh

  if (length(stat) != length(obj$gene.names)) {
    stop(sprintf("Length mismatch in %s: beta.hat=%d, gene.names=%d",
                 celltype, length(stat), length(obj$gene.names)))
  }

  df <- tibble(
    gene = clean_symbols(obj$gene.names, force_upper = force_upper),
    stat = as.numeric(stat)
  ) %>%
    filter(!is.na(stat)) %>%
    collapse_duplicate_stats()

  if (!is.null(universe_genes)) {
    universe_genes <- clean_symbols(universe_genes, force_upper = force_upper)
    df <- df %>% filter(gene %in% universe_genes)
  }

  df <- df %>% arrange(desc(stat))

  ranks <- df$stat
  names(ranks) <- df$gene
  ranks <- sort(ranks, decreasing = TRUE)

  list(
    table = df,
    ranks = ranks
  )
}

run_fgsea_one <- function(ranks, pathways, minSize = 10, maxSize = 2000) {
  # pathways should already be filtered to the same universe as ranks
  fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = minSize,
    maxSize = maxSize,
    eps = 0
  ) %>%
    as_tibble() %>%
    arrange(padj, desc(abs(NES)))
}

run_overlap_tests_one <- function(rank_df,
                                  pathways,
                                  universe_genes,
                                  top_n_up = 300,
                                  top_n_dn = 300,
                                  min_overlap = 3,
                                  min_geneset_size = 10,
                                  force_upper = TRUE) {
  rank_df <- rank_df %>%
    mutate(gene = clean_symbols(gene, force_upper = force_upper),
           stat = as.numeric(stat)) %>%
    filter(!is.na(stat)) %>%
    collapse_duplicate_stats() %>%
    arrange(desc(stat))

  universe_genes <- clean_symbols(universe_genes, force_upper = force_upper)
  rank_df <- rank_df %>% filter(gene %in% universe_genes)

  universe <- unique(universe_genes)
  n_universe <- length(universe)

  top_up <- rank_df %>%
    slice_head(n = min(top_n_up, nrow(rank_df))) %>%
    pull(gene) %>%
    unique()

  top_dn <- rank_df %>%
    slice_tail(n = min(top_n_dn, nrow(rank_df))) %>%
    pull(gene) %>%
    unique()

  test_one <- function(gs_name, gs_genes, tail = c("up", "down")) {
    tail <- match.arg(tail)

    gs <- intersect(clean_symbols(gs_genes, force_upper = force_upper), universe)
    if (length(gs) < min_geneset_size) return(NULL)

    sel <- if (tail == "up") top_up else top_dn
    overlap_genes <- intersect(sel, gs)
    overlap <- length(overlap_genes)

    if (overlap < min_overlap) {
      return(tibble(
        pathway = gs_name,
        tail = tail,
        geneset_size = length(gs),
        selected_size = length(sel),
        overlap = overlap,
        odds_ratio = NA_real_,
        p_value = NA_real_,
        overlap_genes = paste(overlap_genes, collapse = ";")
      ))
    }

    a <- overlap
    b <- length(sel) - a
    c <- length(gs) - a
    d <- n_universe - a - b - c

    if (d < 0) {
      return(tibble(
        pathway = gs_name,
        tail = tail,
        geneset_size = length(gs),
        selected_size = length(sel),
        overlap = overlap,
        odds_ratio = NA_real_,
        p_value = NA_real_,
        overlap_genes = paste(overlap_genes, collapse = ";")
      ))
    }

    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    ft <- fisher.test(mat, alternative = "greater")

    tibble(
      pathway = gs_name,
      tail = tail,
      geneset_size = length(gs),
      selected_size = length(sel),
      overlap = overlap,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      overlap_genes = paste(overlap_genes, collapse = ";")
    )
  }

  bind_rows(
    purrr::map_dfr(names(pathways), ~ test_one(.x, pathways[[.x]], "up")),
    purrr::map_dfr(names(pathways), ~ test_one(.x, pathways[[.x]], "down"))
  ) %>%
    group_by(tail) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
    ungroup() %>%
    arrange(tail, p_adj, desc(overlap))
}

run_fora_one <- function(rank_df, pathways, universe_genes,
                         top_n = 300, tail = c("up", "down"),
                         force_upper = TRUE) {
  tail <- match.arg(tail)

  rank_df <- rank_df %>%
    mutate(gene = clean_symbols(gene, force_upper = force_upper),
           stat = as.numeric(stat)) %>%
    filter(!is.na(stat)) %>%
    collapse_duplicate_stats() %>%
    arrange(desc(stat))

  universe_genes <- clean_symbols(universe_genes, force_upper = force_upper)
  rank_df <- rank_df %>% filter(gene %in% universe_genes)

  fg <- if (tail == "up") {
    rank_df %>% slice_head(n = min(top_n, nrow(rank_df))) %>% pull(gene)
  } else {
    rank_df %>% slice_tail(n = min(top_n, nrow(rank_df))) %>% pull(gene)
  }

  fora(
    genes = fg,
    universe = universe_genes,
    pathways = pathways
  ) %>%
    as_tibble() %>%
    mutate(tail = tail) %>%
    arrange(padj, desc(overlap))
}

make_nes_heatmap <- function(fgsea_all, file_out) {
  if (nrow(fgsea_all) == 0) return(NULL)

  plot_df <- fgsea_all %>%
    mutate(
      sig_label = case_when(
        is.na(padj) ~ "",
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE ~ ""
      )
    )

  p <- ggplot(plot_df, aes(x = celltype, y = pathway, fill = NES)) +
    geom_tile() +
    geom_text(aes(label = sig_label), size = 3) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      title = "EBNA2 gene-set enrichment across B-cell subsets",
      x = "Cell type",
      y = "Gene set",
      fill = "NES"
    )

  ggsave(file_out, p, width = 10,
         height = max(4, 0.35 * dplyr::n_distinct(plot_df$pathway)))
  invisible(p)
}


# ---------------------------------------------------------------------------- #
# RUN

pathways <- prepare_pathways(
  pathways = pathways_raw,
  universe_genes = universe_genes,
  min_geneset_size = 10,
  max_geneset_size = 2000
)

celltypes <- names(scdist_bcell_rel$vals)

fgsea_all <- purrr::map_dfr(celltypes, function(ct) {
  x <- get_ranked_stats_from_scdist(scdist_bcell_rel, ct, universe_genes = universe_genes)
  run_fgsea_one(x$ranks, pathways) %>%
    mutate(celltype = ct)
})

ora_all <- purrr::map_dfr(celltypes, function(ct) {
  x <- get_ranked_stats_from_scdist(scdist_bcell_rel, ct, universe_genes = universe_genes)
  run_overlap_tests_one(
    rank_df = x$table,
    pathways = pathways,
    universe_genes = universe_genes,
    top_n_up = 300,
    top_n_dn = 300,
    min_overlap = 3
  ) %>%
    mutate(celltype = ct)
})

fora_all <- purrr::map_dfr(celltypes, function(ct) {
  x <- get_ranked_stats_from_scdist(scdist_bcell_rel, ct, universe_genes = universe_genes)

  bind_rows(
    run_fora_one(x$table, pathways, universe_genes, top_n = 300, tail = "up"),
    run_fora_one(x$table, pathways, universe_genes, top_n = 300, tail = "down")
  ) %>%
    mutate(celltype = ct)
})

readr::write_csv(fgsea_all, file.path(OUT_DIR, "b_cell_rel_ebna2_fgsea_all.csv"))
readr::write_csv(ora_all,   file.path(OUT_DIR, "b_cell_rel_ebna2_overlap_all.csv"))
readr::write_csv(fora_all,  file.path(OUT_DIR, "b_cell_rel_ebna2_fora_all.csv"))

make_nes_heatmap(fgsea_all, "output/b_cell_rel_ebna2_fgsea_heatmap.pdf")
