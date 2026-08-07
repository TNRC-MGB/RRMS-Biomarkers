
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : B_clean_scdrs.RDS             (B-cell scDRS-scored Seurat object)
#          scdist_bcell_prerelapse_v_remission.RDS
#            (scDist, PreRelapse vs Remission)
#          df_int.RDS                   (integrated B-cell scDist / scDRS table)
#          ebna2_combined_panel_with_bound_broad.gmt   (EBNA2 gene sets)
#          kshv_lana_gse56144_controls.gmt             (KSHV LANA control sets)
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")


set.seed(42)

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

OBJ_DIR   <- "objects"
out_dir   <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# EBNA2 / KSHV gene sets
EBNA2_DIR <- "data/genesets"

# --- Inputs -------------------------------------------------------------------
B_clean <- readRDS(file.path(OBJ_DIR, "B_clean_scdrs.RDS"))      # B-cell scDRS-scored Seurat object
scdist_bcell_pre <- readRDS(file.path(SCDIST_DIR, "scdist_bcell_prerelapse_v_remission.RDS"))  # scDist, PreRelapse vs Remission
df_int  <- readRDS(file.path(OBJ_DIR, "df_int.RDS"))             # integrated B-cell scDist / scDRS table

# --- fgsea helper functions ---------------------------------------------------
clean_symbols <- function(x, force_upper = TRUE) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  if (force_upper) x <- toupper(x)
  unique(x)
}

collapse_duplicate_stats <- function(df) {
  df %>%
    group_by(gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()
}

read_gmt_base <- function(gmt_file, force_upper = TRUE) {
  pathways <- fgsea::gmtPathways(gmt_file)
  lapply(pathways, clean_symbols, force_upper = force_upper)
}

prepare_pathways <- function(pathways,
                             universe_genes = NULL,
                             force_upper = TRUE,
                             min_geneset_size = 2,
                             max_geneset_size = 2000) {
  if (!is.null(universe_genes)) {
    universe_genes <- clean_symbols(universe_genes, force_upper = force_upper)
    pathways <- lapply(pathways, function(gs) intersect(clean_symbols(gs, force_upper), universe_genes))
  } else {
    pathways <- lapply(pathways, clean_symbols, force_upper = force_upper)
  }
  pathways[lengths(pathways) >= min_geneset_size & lengths(pathways) <= max_geneset_size]
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

  list(table = df, ranks = ranks)
}

run_fgsea_one <- function(ranks, pathways, minSize = 2, maxSize = 2000) {
  fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = minSize,
    maxSize = maxSize,
    eps = 0
  ) %>%
    as_tibble() %>%
    arrange(padj, desc(abs(NES)))
}

# --- EBV/KSHV gene sets + fgsea across B-cell subsets -------------------------
pathways_9 <- c(
  "EBNA2_BOUND_GSE246060_SHARED_PROMOTER",
  "EBNA2_BOUND_GSE246060_TYPE1_SPECIFIC_PROMOTER",
  "EBNA2_BOUND_GSE246060_TYPE2_SPECIFIC_PROMOTER",
  "EBNA2_CORE_MULTI_STUDY_UP",
  "EBNA2_CORE_MULTI_STUDY_DN",
  "KSHV_LANA_BOUND_GSE56144_BJAB219_PROMOTER",
  "KSHV_LANA_BOUND_GSE56144_BCBL1_PROMOTER"
)

ebna2_gmt_file <- file.path(EBNA2_DIR, "combined", "ebna2_combined_panel_with_bound_broad.gmt")
kshv_gmt_file  <- file.path(EBNA2_DIR, "combined", "kshv_lana_gse56144_controls.gmt")

pathways_raw <- c(
  read_gmt_base(ebna2_gmt_file),
  read_gmt_base(kshv_gmt_file)
)

# B-cell expressed universe (>=1% detected)
DefaultAssay(B_clean) <- "RNA"
mat <- LayerData(B_clean, assay = "RNA", layer = "counts")
det_frac <- Matrix::rowSums(mat > 0) / ncol(mat)
universe_genes <- clean_symbols(names(det_frac)[det_frac >= 0.01])

pathways_9_list <- pathways_raw[names(pathways_raw) %in% pathways_9] %>%
  prepare_pathways(
    universe_genes = universe_genes,
    min_geneset_size = 2,
    max_geneset_size = 2000
  )

celltypes_use <- intersect(names(scdist_bcell_pre$vals), levels(B_clean$B_annotations))

fgsea_all <- purrr::map_dfr(celltypes_use, function(ct) {
  x <- get_ranked_stats_from_scdist(scdist_bcell_pre, ct, universe_genes = universe_genes)
  run_fgsea_one(
    ranks = x$ranks,
    pathways = pathways_9_list,
    minSize = 2,
    maxSize = 2000
  ) %>%
    mutate(B_annotations = ct)
})

fgsea_plot_in <- fgsea_all
stopifnot(all(c("pathway", "B_annotations", "NES", "padj") %in% colnames(fgsea_plot_in)))


# ---------------------------------------------------------------------------- #
# a - EBV/KSHV NES heatmap across B-cell subtypes
bcell_order <- c(
  "ABC",
  "Atypical/ABC-like Memory B",
  "CD11c++ Activated Memory B",
  "Activated Switched Memory B",
  "Switched Memory B",
  "Plasmablast",
  "Transitional/Immature-like B",
  "TCL1A+ Naive B",
  "Naive B",
  "IFN-stimulated Naive B",
  "Activated Memory B",
  "Early Activated Memory B",
  "Anergic/Naive-leaning B"
)

heat_df <- fgsea_plot_in %>%
  filter(pathway %in% pathways_9) %>%
  mutate(
    pathway = factor(pathway, levels = rev(pathways_9)),
    B_annotations = factor(as.character(B_annotations), levels = bcell_order),
    sig_label = case_when(
      is.na(padj) ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

pathway_labels <- c(
  "EBNA2_BOUND_GSE246060_SHARED_PROMOTER"         = "EBNA2 shared promoter-bound",
  "EBNA2_BOUND_GSE246060_TYPE1_SPECIFIC_PROMOTER" = "EBNA2 type 1-specific promoter-bound",
  "EBNA2_BOUND_GSE246060_TYPE2_SPECIFIC_PROMOTER" = "EBNA2 type 2-specific promoter-bound",
  "EBNA2_CORE_MULTI_STUDY_UP"                     = "EBNA2 multi-study up",
  "EBNA2_CORE_MULTI_STUDY_DN"                     = "EBNA2 multi-study down",
  "KSHV_LANA_BOUND_GSE56144_BJAB219_NEAREST"      = "KSHV LANA BJAB219 nearest",
  "KSHV_LANA_BOUND_GSE56144_BJAB219_PROMOTER"     = "KSHV LANA BJAB219 promoter",
  "KSHV_LANA_BOUND_GSE56144_BCBL1_NEAREST"        = "KSHV LANA BCBL1 nearest",
  "KSHV_LANA_BOUND_GSE56144_BCBL1_PROMOTER"       = "KSHV LANA BCBL1 promoter"
)

p_ebv_kshv_heatmap <- ggplot(
  heat_df,
  aes(x = B_annotations, y = pathway, fill = NES, alpha = padj < 0.05)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig_label), size = 3) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "firebrick",
    midpoint = 0
  ) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  scale_y_discrete(labels = pathway_labels) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(out_dir, "ebv_kshv_9row_bcell_heatmap.pdf"), p_ebv_kshv_heatmap, width = 120, height = 75, units = "mm")
readr::write_csv(heat_df, file.path(out_dir, "ebv_kshv_9row_bcell_heatmap_table.csv"))

p_ebv_kshv_heatmap <- p_ebv_kshv_heatmap + theme(text = element_text(size = 6))
ggsave(file.path(out_dir, "ebv_kshv_heatmap.pdf"), p_ebv_kshv_heatmap, width = 100, height = 50, units = "mm")


# ---------------------------------------------------------------------------- #
# b - Mean EBNA2 supportive NES vs scDist distance (scDRS/MS-GWAS color)
ebv_ebna2_sets <- c(
  "EBNA2_BOUND_GSE246060_SHARED_PROMOTER",
  "EBNA2_BOUND_GSE246060_TYPE1_SPECIFIC_PROMOTER",
  "EBNA2_BOUND_GSE246060_TYPE2_SPECIFIC_PROMOTER",
  "EBNA2_CORE_MULTI_STUDY_UP",
  "EBNA2_CORE_MULTI_STUDY_DN"
)

fgsea_ebv_key <- fgsea_plot_in %>%
  filter(pathway %in% ebv_ebna2_sets) %>%
  mutate(
    direction = case_when(
      grepl("_UP$", pathway) ~ "up",
      grepl("_DN$", pathway) ~ "dn",
      TRUE ~ "bound"
    ),
    supportive_nes = case_when(
      direction == "up" ~ NES,
      direction == "dn" ~ -NES,
      TRUE ~ NES
    )
  )

ebv_subset_summary <- fgsea_ebv_key %>%
  group_by(B_annotations) %>%
  summarize(
    ebv_mean_supportive_nes = mean(supportive_nes, na.rm = TRUE),
    ebv_sd_supportive_nes   = sd(supportive_nes, na.rm = TRUE),
    ebv_n_sets              = sum(is.finite(supportive_nes)),
    ebv_se_supportive_nes   = ebv_sd_supportive_nes / sqrt(ebv_n_sets),
    ebv_low                 = ebv_mean_supportive_nes - 1.96 * ebv_se_supportive_nes,
    ebv_high                = ebv_mean_supportive_nes + 1.96 * ebv_se_supportive_nes,
    ebv_best_fdr            = suppressWarnings(min(padj, na.rm = TRUE)),
    ebv_n_sig               = sum(!is.na(padj) & padj < 0.05),
    .groups = "drop"
  ) %>%
  mutate(
    ebv_best_fdr = ifelse(is.infinite(ebv_best_fdr), NA_real_, ebv_best_fdr),
    ebv_sig = !is.na(ebv_best_fdr) & ebv_best_fdr < 0.05
  )

df_plot_ebv_nes <- df_int %>%
  left_join(ebv_subset_summary, by = "B_annotations")

# match heatmap B-cell ordering
bcell_order <- levels(heat_df$B_annotations)
df_plot_ebv_nes$B_annotations <- factor(as.character(df_plot_ebv_nes$B_annotations), levels = bcell_order)

label_subsets <- c(
  "ABC",
  "Atypical/ABC-like Memory B",
  "CD11c++ Activated Memory B",
  "Activated Switched Memory B",
  "Switched Memory B",
  "Naive B"
)

df_plot_ebv_nes <- df_plot_ebv_nes %>%
  mutate(
    label_plot = ifelse(as.character(B_annotations) %in% label_subsets,
                        as.character(B_annotations), "")
  )

# drop malformed / incomplete rows first
df_plot_ebv_nes2 <- df_plot_ebv_nes %>%
  dplyr::filter(
    !is.na(B_annotations),
    !is.na(scdist),
    !is.na(ebv_mean_supportive_nes),
    !is.na(scdrs_emmean)
  ) %>%
  mutate(
    scdrs_fdr_plot = pmax(scdrs_fdr, 1e-300, na.rm = TRUE),
    neglog10_scdrs_fdr = -log10(scdrs_fdr_plot)
  )

p_b_ebna2nes_vs_scdist <- ggplot(
  df_plot_ebv_nes2,
  aes(
    x = scdist,
    y = ebv_mean_supportive_nes
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray75") +
  geom_errorbarh(
    aes(xmin = scdist_low, xmax = scdist_high),
    height = 0,
    linewidth = 0.25,
    alpha = 0.8
  ) +
  # descriptive spread across EBV/EBNA2 pathways, not a formal CI
  geom_errorbar(
    aes(ymin = ebv_low, ymax = ebv_high),
    width = 0,
    linewidth = 0.25,
    alpha = 0.6
  ) +
  geom_point(size = 5, shape = 21, color = "gray50", aes(fill = neglog10_scdrs_fdr)) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  ggrepel::geom_text_repel(
    aes(label = label_plot),
    size = 2.3,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.15,
    show.legend = FALSE
  ) +
  labs(
    x = "scDist distance (PreRelapse vs Remission)",
    y = "EBV/EBNA2 mean supportive NES",
    fill = expression(-log[10]("scDRS FDR"))
  ) +
  theme_classic(base_size = 9)

ggsave(file.path(out_dir, "ebv_ebna2nes_vs_scdist_scdrscolor.pdf"), p_b_ebna2nes_vs_scdist, width = 150, height = 90, units = "mm")


