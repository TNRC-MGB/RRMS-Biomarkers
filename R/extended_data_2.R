
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Inputs : pbmc_final.RDS   (preprocessing_pbmc.R)
#
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

set.seed(42)

source("R/scrna_preprocessing_functions.R")

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Clusters with top Azimuth label accounting for at least this fraction of their
# cells are assigned that label directly; the rest are annotated manually.
AZIMUTH_THRESHOLD <- 0.80

# --- Inputs -------------------------------------------------------------------
seu <- readRDS(file.path(OBJ_DIR, "pbmc_final.RDS"))


# ---------------------------------------------------------------------------- #
# a - Top Azimuth-label fraction per Harmony cluster
# NOTE: the per-cluster Azimuth summary is not saved as a separate object;
#       recompute it from pbmc_final with the same call preprocessing_pbmc.R
#       used to assign the labels.
res <- assign_cluster_labels_from_azimuth(
  seu         = seu,
  cluster_col = "clusters.harmony",
  azimuth_col = "predicted.celltype.l2",
  threshold   = AZIMUTH_THRESHOLD,
  new_col     = "clusters.harmony.annot"
)

p_a <- ggplot(res$summary,
              aes(x = reorder(cluster, top_fraction), y = top_fraction,
                  fill = top_azimuth_label)) +
  geom_col() +
  geom_hline(yintercept = AZIMUTH_THRESHOLD, linetype = "dashed") +
  scale_fill_discrete_qualitative(palette = "Dynamic") +
  coord_flip() +
  theme_bw() +
  labs(x = "Harmony cluster", y = "Top Azimuth label fraction",
       fill = "Top Azimuth label")

ggsave(file.path(out_dir, "extended_data_2a.pdf"), p_a,
       width = 80, height = 60, units = "mm")


# ---------------------------------------------------------------------------- #
# b - Manual PBMC annotation vs Azimuth label confusion matrix
df_comp <- seu[[]] %>%
  dplyr::select(predicted.celltype.l2, pbmc_annotations) %>%
  dplyr::filter(!is.na(predicted.celltype.l2), !is.na(pbmc_annotations)) %>%
  dplyr::count(predicted.celltype.l2, pbmc_annotations, name = "n") %>%
  dplyr::group_by(pbmc_annotations) %>%
  dplyr::mutate(frac_within_manual = n / sum(n)) %>%
  dplyr::ungroup()

p_b <- ggplot(df_comp, aes(x = predicted.celltype.l2, y = pbmc_annotations,
                           fill = frac_within_manual)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "steelblue", labels = scales::percent) +
  theme_bw() +
  labs(x = "Azimuth label", y = "Manual PBMC annotation",
       fill = 'Fraction within\nmanual label') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid  = element_blank())

ggsave(file.path(out_dir, "extended_data_2b.pdf"), p_b,
       width = 120, height = 120, units = "mm")


# ---------------------------------------------------------------------------- #
# c - Annotated PBMC UMAP (final manual labels)
col_pbmc <- colorspace::qualitative_hcl(n = length(levels(seu$pbmc_annotations)),
                                        palette = "Set 3")
col_pbmc <- setNames(col_pbmc, levels(seu$pbmc_annotations))

p_c <- DimPlot(seu, label = TRUE, label.size = 4, repel = TRUE,
               reduction = "umap.harmony", group.by = "pbmc_annotations") +
  scale_color_manual(values = col_pbmc) +
  coord_fixed() +
  NoLegend()

ggsave(file.path(out_dir, "extended_data_2c.pdf"), p_c,
       width = 600, height = 600, units = "mm", dpi = 600)
