
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Inputs : fgsea_neb_up.RDS   (written by extended_data_9.R)
#
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

set.seed(42)

source("R/enrichment_functions.R")   # calc_fgsea_from_ranks, build_nebula_ranks,
                                     # build_nebula_universe, EBV lifecycle stages

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
NEBULA_DIR <- file.path(ZENODO_DIR, "nebula")
ENRICH_DIR <- "data/genesets"
OBJ_DIR    <- "objects"
out_dir    <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Settings -----------------------------------------------------------------
GSEA_DIRECTION   <- "up"
PATHWAY_MIN_SIZE <- 10
PATHWAY_MAX_SIZE <- 600
P_ADJUST_SCOPE   <- "by_cluster"
PADJ_CAP         <- 5            # -log10(padj_custom) color cap

group_levels <- c("latent", "early", "leaky late", "true late", "other")

# --- Inputs -------------------------------------------------------------------
# can reuse the enrichment computed by extended_data_9.R when available
fgsea_cache <- file.path(OBJ_DIR, "fgsea_neb_up.RDS")

if (file.exists(fgsea_cache)) {
  fgsea_neb_up <- readRDS(fgsea_cache)
  message("Loaded cached fgsea result from ", fgsea_cache)
} else {
  res      <- read.csv(file.path(NEBULA_DIR, "DEG_all_Bannotations_combined.csv.gz"))
  pathways <- readRDS(file.path(ENRICH_DIR, "pebv500.RDS"))

  fgsea_neb_up <- calc_fgsea_from_ranks(
    ranks            = build_nebula_ranks(res),
    pathways         = pathways,
    universe         = build_nebula_universe(res),
    direction        = GSEA_DIRECTION,
    pathway_min_size = PATHWAY_MIN_SIZE,
    pathway_max_size = PATHWAY_MAX_SIZE,
    p_adjust_scope   = P_ADJUST_SCOPE
  )

  saveRDS(fgsea_neb_up, fgsea_cache)
}


# ---------------------------------------------------------------------------- #
# NES by EBV lifecycle stage, faceted by B-cell subset
# Only positively enriched gene sets are shown
# GSEA: positive NES = upregulated in pre-relapse relative to remission.
plot_df <- fgsea_neb_up %>%
  mutate(
    pathway_group        = factor(ebv_lifecycle_stage(pathway), levels = group_levels),
    neglog10_padj_capped = pmin(-log10(pmax(padj_custom, .Machine$double.xmin)), PADJ_CAP)
  ) %>%
  filter(pathway_group %in% group_levels) %>%
  filter(NES > 0)

message("Gene sets plotted: ", nrow(plot_df), " across ",
        length(unique(plot_df$cluster)), " B-cell subsets")

p <- ggplot(plot_df, aes(x = pathway_group, y = NES)) +
  geom_violin(fill = NA, color = "black", linewidth = 0.5,
              trim = FALSE, width = 0.9) +
  geom_point(aes(color = neglog10_padj_capped),
             position = position_jitter(width = 0.22, height = 0),
             size = 1.1, alpha = 0.9) +
  facet_wrap(~ cluster) +
  scale_color_gradientn(
    colors = c("gray35", "gray70", "white", "#F4A582", "#D73027"),
    limits  = c(0, PADJ_CAP),
    name    = expression(-log[10]("padj"))
  ) +
  labs(x = NULL, y = "Normalized enrichment score (NES)", title = "Pre-relapse") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file.path(out_dir, "extended_data_10.pdf"), p,
       width = 240, height = 200, units = "mm")
