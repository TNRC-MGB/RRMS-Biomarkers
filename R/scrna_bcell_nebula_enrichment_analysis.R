
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : nebula_DE_by_Bannotation/DEG_all_Bannotations_combined.csv
#            (NEBULA differential expression, pre-relapse vs remission,
#             one row per gene per B-cell subset)
#          Enrichment Analysis/pebv500.RDS
#            (EBV factor-derived host gene sets: top 500 genes per EBV factor,
#             from the Arvey et al. 2012 EBV-human expression atlas)
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

# fgseaMultilevel estimates p-values by sampling
set.seed(42)

source("R/enrichment_functions.R")     # calc_fora_from_ranks, calc_fgsea_from_ranks,
                                       # build_nebula_ranks, EBV lifecycle stages


ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
NEBULA_DIR <- file.path(ZENODO_DIR, "nebula")
ENRICH_DIR <- "data/genesets"
OBJ_DIR    <- "objects"
out_dir    <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# Inputs
res      <- read.csv(file.path(NEBULA_DIR, "DEG_all_Bannotations_combined.csv.gz"))
pathways <- readRDS(file.path(ENRICH_DIR, "pebv500.RDS"))
message("NEBULA rows: ", nrow(res), " across ",
        length(unique(res$B_annotation)), " B-cell subsets")
message("EBV gene sets: ", length(pathways))


# Ranked gene lists and gene universe, per B-cell subset
neb_rank_list     <- build_nebula_ranks(res)
neb_universe_list <- build_nebula_universe(res)



# Over-representation of EBV gene sets in the top-ranked genes (Figure 5d input)
ora_neb_up <- calc_fora_from_ranks(
  ranks = neb_rank_list,
  pathways = pathways,
  universe = neb_universe_list,
  direction = "up",
  top_n = 200,
  min_genes = 10,
  pathway_min_size = 10,
  pathway_max_size = 500,
  p_adjust_scope = "by_cluster"
)

ora_neb_down <- calc_fora_from_ranks(
  ranks = neb_rank_list,
  pathways = pathways,
  universe = neb_universe_list,
  direction = "down",
  top_n = 200,
  min_genes = 10,
  pathway_min_size = 10,
  pathway_max_size = 500,
  p_adjust_scope = "by_cluster"
)


saveRDS(ora_neb_up,   file.path(OBJ_DIR, "ora_neb_up.RDS"))
saveRDS(ora_neb_down, file.path(OBJ_DIR, "ora_neb_down.RDS"))



# Preranked GSEA of the same gene sets (Extended Data Fig. 9d / Fig. 10 input)
fgsea_neb_std <- calc_fgsea_from_ranks(
  ranks = neb_rank_list,
  pathways = pathways,
  universe = neb_universe_list,
  direction = "std",
  pathway_min_size = 10,
  pathway_max_size = 500,
  p_adjust_scope = "by_cluster"
)

saveRDS(fgsea_neb_std, file.path(OBJ_DIR, "fgsea_neb_std.RDS"))

fgsea_neb_up <- calc_fgsea_from_ranks(
  ranks = neb_rank_list,
  pathways = pathways,
  universe = neb_universe_list,
  direction = "up",
  pathway_min_size = 10,
  pathway_max_size = 600,
  p_adjust_scope = "by_cluster"
)


ggplot(fgsea_neb_up,aes(x=NES,y=pathway,color=-log10(padj))) +
  geom_point() +
  facet_wrap(~cluster)



# Extended Data Fig. 9d / Fig. 10
plot_df <- fgsea_neb_up %>%
  mutate(
    neglog10_padj_capped = pmin(
      -log10(pmax(padj, .Machine$double.xmin)),
      10
    ),
    pathway_group = case_when(
      pathway %in% latent ~ "latent",
      pathway %in% early ~ "early",
      pathway %in% leaky_late ~ "leaky late",
      pathway %in% true_late ~ "true late",
      pathway %in% other ~ "other",
      TRUE ~ "unassigned"
    ),
    pathway = factor(pathway, levels = rev(custom_pathway_order))
  )

ggplot(plot_df, aes(x = NES, y = pathway, color = neglog10_padj_capped)) +
  geom_point() +
  facet_wrap(~ cluster) +
  scale_color_viridis_c(
    name = expression(-log[10]("padj")),
    limits = c(0, 5),
    oob=scales::squish
  ) +
  theme_bw()


group_colors <- c(
  "latent" = "#E41A1C",
  "early" = "#377EB8",
  "leaky late" = "#4DAF4A",
  "true late" = "#984EA3",
  "other" = "#FF7F00"
)

# keep only pathways actually present
path_levels_present <- rev(custom_pathway_order[custom_pathway_order %in% unique(fgsea_neb_up$pathway)])

plot_df <- fgsea_neb_up %>%
  filter(pathway %in% path_levels_present) %>%
  mutate(
    neglog10_padj_capped = pmin(
      -log10(pmax(padj, .Machine$double.xmin)),
      10
    ),
    pathway_group = case_when(
      pathway %in% latent ~ "latent",
      pathway %in% early ~ "early",
      pathway %in% leaky_late ~ "leaky late",
      pathway %in% true_late ~ "true late",
      pathway %in% other ~ "other",
      TRUE ~ "other"
    ),
    pathway = factor(pathway, levels = path_levels_present)
  )

strip_df <- plot_df %>%
  distinct(pathway, pathway_group)

p_strip <- ggplot(strip_df, aes(x = 1, y = pathway, fill = pathway_group)) +
  geom_tile(width = 1, height = 0.95) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.y.left = element_text(hjust = 1),
    plot.margin = margin(5.5, 2, 5.5, 5.5)
  )

p_main <- ggplot(plot_df, aes(x = NES, y = pathway, color = neglog10_padj_capped)) +
  geom_point() +
  facet_wrap(~ cluster) +
  scale_y_discrete(drop = FALSE) +
  scale_color_viridis_c(
    name = expression(-log[10]("padj")),
    limits = c(0, 10)
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

p_strip + p_main + plot_layout(widths = c(0.28, 1))


# from earlier:
#   plot_df
#   path_levels_present
#   group_colors

make_cluster_pair <- function(cluster_name, plot_df, path_levels_present, group_colors) {

  df_cl <- plot_df %>%
    filter(cluster == cluster_name) %>%
    mutate(pathway = factor(pathway, levels = path_levels_present))

  # strip data: one row per pathway, keep all levels for alignment
  strip_df <- tibble(
    pathway = factor(path_levels_present, levels = path_levels_present)
  ) %>%
    mutate(
      pathway_group = case_when(
        as.character(pathway) %in% latent ~ "latent",
        as.character(pathway) %in% early ~ "early",
        as.character(pathway) %in% leaky_late ~ "leaky late",
        as.character(pathway) %in% true_late ~ "true late",
        as.character(pathway) %in% other ~ "other",
        TRUE ~ "other"
      )
    )

  p_strip <- ggplot(strip_df, aes(x = 1, y = pathway, fill = pathway_group)) +
    geom_tile(width = 1, height = 0.95) +
    scale_fill_manual(values = group_colors, guide = "none") +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      axis.text.y = element_blank(),
      plot.margin = margin(5.5, 2, 5.5, 5.5)
    )

  p_main <- ggplot(df_cl, aes(x = NES, y = pathway, color = neglog10_padj_capped)) +
    geom_point() +
    scale_y_discrete(drop = FALSE) +
    # scale_color_viridis_c(
    #  option = 'viridis',
    #  name = expression(-log[10]("padj")),
    #  limits = c(0, 10)
    # ) +
    scale_color_distiller(
      palette="Reds",
      direction = 1,
      limits=c(0,10)
    ) +
    theme_minimal() +
    labs(title = cluster_name, y = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      #axis.text.y = element_blank(),
      #axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )

  p_strip + p_main + plot_layout(widths = c(0.08, 1)) & theme(axis.title = element_blank(), axis.text = element_blank())
}

clusters_use <- unique(plot_df$cluster)

plot_list <- lapply(clusters_use, function(cl) {
  make_cluster_pair(
    cluster_name = cl,
    plot_df = plot_df,
    path_levels_present = path_levels_present,
    group_colors = group_colors
  )
})

ex8_1 <- wrap_plots(plot_list, ncol = 3, guides = "collect") &
  theme(legend.position = "right")


ggsave(file.path(out_dir, "ex8_1.pdf"),
       ex8_1,
       width=200,
       height=200,
       units="mm")


group_levels <- c("latent", "early", "leaky late", "true late", "other")

plot_df <- fgsea_neb_up %>%
  mutate(
    pathway_group = case_when(
      pathway %in% latent ~ "latent",
      pathway %in% early ~ "early",
      pathway %in% leaky_late ~ "leaky late",
      pathway %in% true_late ~ "true late",
      pathway %in% other ~ "other",
      TRUE ~ "other"
    ),
    pathway_group = factor(pathway_group, levels = group_levels),
    neglog10_padj_capped = pmin(
      -log10(pmax(padj_custom, .Machine$double.xmin)),
      5
    )
  ) %>%
  filter(pathway_group %in% group_levels)

# plot_df %>%
#   group_by(pathway_group) %>%
#   summarize(
#     median_NES = median(NES, na.rm = TRUE),
#     .groups = "drop"
#   ) -> plot_df

ggplot(plot_df, aes(x = pathway_group, y = NES)) +
  geom_violin(
    fill = NA,
    color = "black",
    linewidth = 0.5,
    trim = FALSE,
    width = 0.9
  ) +
  geom_point(
    aes(color = neglog10_padj_capped),
    position = position_jitter(width = 0.22, height = 0),
    size = 1.2,
    alpha = 0.9
  ) +
  scale_color_gradientn(
    colors = c("gray35", "gray70", "white", "#F4A582", "#D73027"),
    limits = c(0, 5),
    name = expression(-log[10]("padj"))
  ) +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    title = "Pre-relapse"
  ) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

group_levels <- c("latent", "early", "leaky late", "true late", "other")

plot_df <- fgsea_neb_up %>%
  mutate(
    pathway_group = case_when(
      pathway %in% latent ~ "latent",
      pathway %in% early ~ "early",
      pathway %in% leaky_late ~ "leaky late",
      pathway %in% true_late ~ "true late",
      pathway %in% other ~ "other",
      TRUE ~ "other"
    ),
    pathway_group = factor(pathway_group, levels = group_levels),
    neglog10_padj_capped = pmin(
      -log10(pmax(padj_custom, .Machine$double.xmin)),
      5
    )
  ) %>%
  filter(pathway_group %in% group_levels) %>%
  filter(NES > 0)

ex8_2 <- ggplot(plot_df, aes(x = pathway_group, y = NES)) +
  geom_violin(
    fill = NA,
    color = "black",
    linewidth = 0.5,
    trim = FALSE,
    width = 0.9
  ) +
  geom_point(
    aes(color = neglog10_padj_capped),
    position = position_jitter(width = 0.22, height = 0),
    size = 1.1,
    alpha = 0.9
  ) +
  facet_wrap(~ cluster) +
  scale_color_gradientn(
    colors = c("gray35", "gray70", "white", "#F4A582", "#D73027"),
    limits = c(0, 5),
    name = expression(-log[10]("padj"))
  ) +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    title = "Pre-relapse"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(file.path(out_dir, "ex8_2.pdf"),
       ex8_2,
       width=240,
       height=200,
       units="mm")
