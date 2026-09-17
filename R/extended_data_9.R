
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : nebula_DE_by_Bannotation/DEG_all_Bannotations_combined.csv
#            (NEBULA DE)
#          Enrichment Analysis/pebv500.RDS   (EBV factor-derived host gene sets:
#            top 500 genes per factor, Arvey et al. 2012 EBV-human atlas)
#
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

set.seed(42)

source("R/enrichment_functions.R")   # calc_fgsea_from_ranks, build_nebula_ranks,
                                     # build_nebula_universe, EBV lifecycle stages,
                                     # ebv_lifecycle_stage, ebv_stage_colors

# The NEBULA results are in the Zenodo deposit; the EBV gene sets
# are in this repo
ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
NEBULA_DIR <- file.path(ZENODO_DIR, "nebula")
ENRICH_DIR <- "data/genesets"
OBJ_DIR    <- "objects"        # pipeline .RDS objects
out_dir    <- "Intermediate"   # rendered panels
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Settings (match the 3-16-2026 source script) -----------------------------
GSEA_DIRECTION   <- "up"         # positive direction only: up in pre-relapse
PATHWAY_MIN_SIZE <- 10
PATHWAY_MAX_SIZE <- 600
P_ADJUST_SCOPE   <- "by_cluster"
PADJ_CAP         <- 10           # -log10(padj) color cap for panel d

# --- Inputs -------------------------------------------------------------------
res      <- read.csv(file.path(NEBULA_DIR, "DEG_all_Bannotations_combined.csv.gz"))
pathways <- readRDS(file.path(ENRICH_DIR, "pebv500.RDS"))

message("NEBULA rows: ", nrow(res), " across ",
        length(unique(res$B_annotation)), " B-cell subsets")
message("EBV gene sets: ", length(pathways))


# ---------------------------------------------------------------------------- #
# c - EBV gene lifecycle stage reference key
key_df <- tibble::tibble(
  pathway = factor(rev(custom_pathway_order), levels = rev(custom_pathway_order))
) %>%
  mutate(stage = factor(ebv_lifecycle_stage(as.character(pathway)),
                        levels = names(ebv_stage_colors)))

stopifnot(!any(is.na(key_df$stage)))

p_c <- ggplot(key_df, aes(x = 1, y = pathway, fill = stage)) +
  geom_tile(width = 1, height = 0.95) +
  scale_fill_manual(values = ebv_stage_colors,
                    labels = ebv_stage_labels,
                    name   = "EBV lifecycle stage") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = "EBV gene") +
  theme_minimal() +
  theme(
    axis.text.y      = element_text(size = 5, color = "black", hjust = 1),
    axis.text.x      = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "right",
    legend.key.size  = unit(0.3, "cm"),
    legend.title     = element_text(size = 6),
    legend.text      = element_text(size = 5)
  )

ggsave(file.path(out_dir, "extended_data_9c.pdf"), p_c,
       width = 60, height = 200, units = "mm")


# ---------------------------------------------------------------------------- #
# Preranked GSEA of the EBV gene sets (input to panel d and to Fig. 10)
neb_rank_list     <- build_nebula_ranks(res)
neb_universe_list <- build_nebula_universe(res)

fgsea_neb_up <- calc_fgsea_from_ranks(
  ranks            = neb_rank_list,
  pathways         = pathways,
  universe         = neb_universe_list,
  direction        = GSEA_DIRECTION,
  pathway_min_size = PATHWAY_MIN_SIZE,
  pathway_max_size = PATHWAY_MAX_SIZE,
  p_adjust_scope   = P_ADJUST_SCOPE
)


# ---------------------------------------------------------------------------- #
# d - Per-subset GSEA dot plots, ordered by EBV lifecycle stage
# keep only pathways actually tested
path_levels_present <- rev(custom_pathway_order[custom_pathway_order %in% unique(fgsea_neb_up$pathway)])

plot_df <- fgsea_neb_up %>%
  filter(pathway %in% path_levels_present) %>%
  mutate(
    neglog10_padj_capped = pmin(-log10(pmax(padj, .Machine$double.xmin)), PADJ_CAP),
    pathway_group        = ebv_lifecycle_stage(pathway),
    pathway              = factor(pathway, levels = path_levels_present)
  )

make_cluster_pair <- function(cluster_name, plot_df, path_levels_present) {

  df_cl <- plot_df %>%
    filter(cluster == cluster_name) %>%
    mutate(pathway = factor(pathway, levels = path_levels_present))

  # keep every level, including untested pathways, so the strips stay aligned
  strip_df <- tibble::tibble(
    pathway = factor(path_levels_present, levels = path_levels_present)
  ) %>%
    mutate(pathway_group = ebv_lifecycle_stage(as.character(pathway)))

  p_strip <- ggplot(strip_df, aes(x = 1, y = pathway, fill = pathway_group)) +
    geom_tile(width = 1, height = 0.95) +
    scale_fill_manual(values = ebv_stage_colors, guide = "none") +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(axis.text.y = element_blank(),
          plot.margin = margin(5.5, 2, 5.5, 5.5))

  p_main <- ggplot(df_cl, aes(x = NES, y = pathway, color = neglog10_padj_capped)) +
    geom_point() +
    scale_y_discrete(drop = FALSE) +
    # scale_color_distiller(palette = "Reds", direction = 1, limits = c(0, PADJ_CAP),
    #                       name = expression(-log[10]("padj"))) +
    scale_color_viridis_c(limits=c(0,PADJ_CAP)) +
    theme_minimal() +
    labs(title = cluster_name, y = NULL) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
          plot.margin = margin(5.5, 5.5, 5.5, 0))

  p_strip + p_main + plot_layout(widths = c(0.08, 1)) &
    theme(axis.title = element_blank(), axis.text = element_blank())
}

plot_list <- lapply(unique(plot_df$cluster), function(cl) {
  make_cluster_pair(cl, plot_df, path_levels_present)
})

p_d <- wrap_plots(plot_list, ncol = 3, guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "extended_data_9d.pdf"), p_d,
       width = 200, height = 200, units = "mm")

saveRDS(fgsea_neb_up, file.path(OBJ_DIR, "fgsea_neb_up.RDS"))
