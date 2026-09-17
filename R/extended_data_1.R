
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : seu_minQC_scDblFinder_filt_firstPass.RDS       (preprocessing_pbmc.R)
#          seu_minQC_scDblFinder_filt_secondPass.RDS      (preprocessing_pbmc.R)
#          seu_minQC_scDblFinder_filt_secondPass_harmony.RDS
#            (preprocessing_pbmc.R)
#
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

set.seed(42)

source("R/scrna_preprocessing_functions.R")

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

font_sizeS <- 12

# Clusters designated for removal at the first-pass QC step. Must stay in sync
# preprocessing_pbmc.R ... panel d reports the resulting % of cells removed.
clusters_to_remove <- c(31, 36, 37, 45, 47, 49)


# ---------------------------------------------------------------------------- #
# First-pass object (panels b, c, d, e)
seu <- readRDS(file.path(OBJ_DIR, "seu_minQC_scDblFinder_filt_firstPass.RDS"))

qc1 <- plot_cluster_qc_panel(
  seu,
  cluster_col       = "clusters.first_pass",
  reduction         = "umap.first_pass",
  label_size        = 4,
  font_sizeS        = 4,
  rasterize_dotplot = TRUE
)


# ---------------------------------------------------------------------------- #
# b - First-pass per-cluster QC dot plot
p_b <- qc1$panel
ggsave(file.path(out_dir, "extended_data_1b.pdf"), p_b,
       width = 500, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# c - First-pass clusters UMAP
Idents(seu) <- "clusters.first_pass"
col1 <- colorspace::qualitative_hcl(n = length(unique(Idents(seu))), palette = "Dynamic")

p_c <- DimPlot(seu, label = TRUE, repel = TRUE, label.size = 4,
               reduction = "umap.first_pass") +
  scale_color_manual(values = col1) +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS) +
  coord_fixed() +
  NoLegend()

ggsave(file.path(out_dir, "extended_data_1c.pdf"), p_c,
       width = 300, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# d - Cells flagged for removal (% removed)
seu$first_pass_remove <- seu$clusters.first_pass %in% clusters_to_remove
percent_removed <- round((sum(seu$first_pass_remove) / ncol(seu)) * 100, 2)
message("Panel d - cells removed at first-pass QC: ", percent_removed, "%")

p_d <- FeaturePlot(seu, features = "first_pass_remove", reduction = "umap.first_pass") +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS, main = 0) +
  coord_fixed() +
  NoLegend() +
  labs(subtitle = paste(percent_removed, "% removed"))

ggsave(file.path(out_dir, "extended_data_1d.pdf"), p_d,
       width = 300, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# e - % removed per library x condition
md <- seu[[]] %>%
  mutate(
    first_pass_remove = as.logical(first_pass_remove),
    Library_ID        = as.character(Library_ID),
    Condition         = as.character(Condition)
  )

by_lib_cond <- md %>%
  group_by(Library_ID, Condition) %>%
  summarize(
    n_total    = n(),
    n_remove   = sum(first_pass_remove, na.rm = TRUE),
    pct_remove = 100 * n_remove / n_total,
    .groups    = "drop"
  ) %>%
  mutate(Library_ID = factor(Library_ID, levels = str_sort(unique(Library_ID), numeric = TRUE)))

p_e <- ggplot(by_lib_cond, aes(x = Condition, y = Library_ID, fill = pct_remove)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(x = NULL, y = NULL, fill = "% removed", title = "% removed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "extended_data_1e.pdf"), p_e,
       width = 100, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# Second-pass object (panels f, g)
seu <- readRDS(file.path(OBJ_DIR, "seu_minQC_scDblFinder_filt_secondPass.RDS"))


# ---------------------------------------------------------------------------- #
# f - Second-pass clusters UMAP
Idents(seu) <- "clusters.second_pass"
col2 <- colorspace::qualitative_hcl(n = length(unique(Idents(seu))), palette = "Pastel")

p_f <- DimPlot(seu, label = TRUE, label.size = 4, repel = TRUE,
               reduction = "umap.second_pass") +
  scale_color_manual(values = col2) +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS) +
  coord_fixed() +
  NoLegend()

ggsave(file.path(out_dir, "extended_data_1f.pdf"), p_f,
       width = 300, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# g - Second-pass UMAP by batch (pre-correction batch assessment)
p_g <- DimPlot(seu, label = TRUE, label.size = 4, repel = TRUE, shuffle = TRUE,
               reduction = "umap.second_pass", group.by = "Batch") +
  scale_color_brewer(palette = "Set2") +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS) +
  coord_fixed() +
  NoLegend()

ggsave(file.path(out_dir, "extended_data_1g.pdf"), p_g,
       width = 300, height = 300, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# h - Harmony-integrated clusters UMAP
# Read the harmonized checkpoint rather than re-running Harmony, so the cluster
# numbering here matches pbmc_final.RDS and every downstream
# figure.
seu <- readRDS(file.path(OBJ_DIR, "seu_minQC_scDblFinder_filt_secondPass_harmony.RDS"))

Idents(seu) <- "clusters.harmony"
col_harmony <- colorspace::qualitative_hcl(n = length(unique(Idents(seu))), palette = "Pastel")

p_h <- DimPlot(seu, label = TRUE, label.size = 4, repel = TRUE,
               reduction = "umap.harmony") +
  scale_color_manual(values = col_harmony) +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS) +
  coord_fixed() +
  NoLegend()

ggsave(file.path(out_dir, "extended_data_1h.pdf"), p_h,
       width = 300, height = 300, units = "mm", dpi = 600)

message("Final atlas: ", ncol(seu), " cells from ",
        length(unique(seu$Sample_ID)), " PBMC vials.")
