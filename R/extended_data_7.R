
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
#
# Inputs : nebula_DE_by_Bannotation/DEG_all_Bannotations_combined.csv
#            (NEBULA DE)
#          scdist_bcell_prerelapse_v_remission.RDS
#            (scDist, PreRelapse vs Remission; preprocessing_bcells.R)
#
# NOTE: per-subset NEBULA DE (raw counts ~ Condition, NBGMM, random intercept =
#       Patient, library size offset, genes expressed in >=10% of cells in
#       either condition) is computed upstream and written to
#       nebula_DE_by_Bannotation/. This script reads those results and renders
#       the volcano grid.
#
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

set.seed(42)

# The NEBULA results and scDist fits are in Zenodo
ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
NEBULA_DIR <- file.path(ZENODO_DIR, "nebula")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")
OBJ_DIR    <- "objects"
out_dir    <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Settings -----------------------------------------------------------------
FDR_THRESH  <- 0.05   # significance cut-off, and the dashed line in the plot
TOP_N_DIST  <- 500    # scDist genes per direction considered for the overlap
N_LABEL     <- 10     # labeled genes per subset
XLIM        <- c(-5, 5)

# --- Inputs -------------------------------------------------------------------
res <- read.csv(file.path(NEBULA_DIR, "DEG_all_Bannotations_combined.csv.gz")) %>%
  mutate(
    B_annotation = as.character(B_annotation),
    gene         = as.character(gene)
  )

scdist_bcell_pre <- readRDS(file.path(SCDIST_DIR, "scdist_bcell_prerelapse_v_remission.RDS"))

message("NEBULA rows: ", nrow(res), " across ",
        length(unique(res$B_annotation)), " B-cell subsets")


# --- Helper -------------------------------------------------------------------
# Genes that are both FDR-significant in NEBULA and among the strongest
# contributors to the scDist distance, in the SAME direction, per B subtype.
find_nebula_scdist_overlap <- function(scdist_obj,
                                       nebula_df,
                                       cluster_col_neb = "B_annotation",
                                       top_n           = TOP_N_DIST,
                                       fdr_thresh      = FDR_THRESH,
                                       concordant_only = TRUE) {

  clusters <- intersect(names(scdist_obj$vals), unique(nebula_df[[cluster_col_neb]]))

  out <- lapply(clusters, function(cl) {

    dist_df <- data.frame(
      gene = scdist_obj$gene.names,
      dist = scdist_obj$vals[[cl]]$beta.hat[, 1],
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        dist_direction = case_when(
          dist > 0 ~ "higher_in_PreRelapse",
          dist < 0 ~ "higher_in_Remission",
          TRUE     ~ "no_difference"
        )
      )

    top_up <- dist_df %>%
      filter(dist > 0) %>%
      arrange(desc(dist)) %>%
      slice_head(n = top_n)

    top_down <- dist_df %>%
      filter(dist < 0) %>%
      arrange(dist) %>%
      slice_head(n = top_n)

    top_dist <- bind_rows(top_up, top_down) %>%
      mutate(B_annotation = cl)

    neb_cl <- nebula_df %>%
      filter(.data[[cluster_col_neb]] == cl, FDR_within_cluster < fdr_thresh) %>%
      mutate(
        neb_direction = case_when(
          logFC > 0 ~ "higher_in_PreRelapse",
          logFC < 0 ~ "higher_in_Remission",
          TRUE      ~ "no_difference"
        )
      )

    ov <- inner_join(top_dist, neb_cl, by = c("gene", "B_annotation"))

    if (concordant_only) ov <- ov %>% filter(dist_direction == neb_direction)

    ov
  })

  bind_rows(out)
}


# ---------------------------------------------------------------------------- #
# NEBULA n scDist overlap
overlap_all <- find_nebula_scdist_overlap(
  scdist_obj      = scdist_bcell_pre,
  nebula_df       = res,
  top_n           = TOP_N_DIST,
  fdr_thresh      = FDR_THRESH,
  concordant_only = TRUE
)

# label positive-logFC genes only (higher in pre-relapse)
overlap_up <- overlap_all[
  overlap_all$dist_direction == "higher_in_PreRelapse" &
  overlap_all$direction      == "higher_in_PreRelapse", ]

label_overlap <- overlap_up %>%
  mutate(overlap_score = abs(dist) * -log10(pmax(FDR_within_cluster, .Machine$double.xmin))) %>%
  group_by(B_annotation) %>%
  arrange(desc(overlap_score), FDR_within_cluster, desc(abs(logFC)), .by_group = TRUE) %>%
  slice_head(n = N_LABEL) %>%
  ungroup() %>%
  select(B_annotation, gene, overlap_score)


# ---------------------------------------------------------------------------- #
# Per-subset volcanos
res_volcano <- res %>%
  mutate(
    neglog10_fdr = -log10(pmax(FDR_within_cluster, .Machine$double.xmin)),
    volcano_group = case_when(
      FDR_within_cluster < FDR_THRESH & logFC > 0 ~ "Higher in PreRelapse",
      FDR_within_cluster < FDR_THRESH & logFC < 0 ~ "Higher in Remission",
      TRUE                                        ~ "Not significant"
    )
  ) %>%
  left_join(label_overlap %>% mutate(label_me = TRUE),
            by = c("B_annotation", "gene")) %>%
  mutate(label_me = ifelse(is.na(label_me), FALSE, label_me))

# shared minimum y so every facet shows the FDR = 0.05 line region, even the
# subsets with no significant genes (CD11c++, Early Activated, Plasmablast)
facet_floor <- res_volcano %>%
  distinct(B_annotation) %>%
  mutate(logFC = 0, neglog10_fdr = 5)

p_volcano <- ggplot(res_volcano, aes(x = logFC, y = neglog10_fdr)) +
  ggrastr::rasterize(
    geom_point(aes(color = volcano_group), alpha = 0.5, size = 0.7),
    dpi = 600
  ) +
  geom_blank(data = facet_floor, aes(x = logFC, y = neglog10_fdr)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_THRESH), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = subset(res_volcano, label_me),
    aes(label = gene),
    size          = 2.8,
    max.overlaps  = Inf,
    box.padding   = 0.25,
    point.padding = 0.15,
    show.legend   = FALSE
  ) +
  xlim(XLIM) +
  facet_wrap(~ B_annotation, scales = "free_y") +
  labs(x = "logFC", y = expression(-log[10]("FDR")), color = NULL) +
  theme_bw()

ggsave(file.path(out_dir, "extended_data_7.pdf"), p_volcano,
       width = 360, height = 180, units = "mm")
