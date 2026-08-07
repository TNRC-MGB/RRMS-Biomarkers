
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Inputs : pbmc_final.RDS   (preprocessing_pbmc.R)
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

font_sizeS <- 12

# Minimum target cells per library at microfluidics loading (panel a reference line)
MIN_CELLS_TARGET <- 1500

cond_levels <- c("Healthy", "PreRelapse", "Relapse", "Remission")

# --- Inputs -------------------------------------------------------------------
seu <- readRDS(file.path(OBJ_DIR, "pbmc_final.RDS"))

col_pbmc <- colorspace::qualitative_hcl(n = length(levels(seu$pbmc_annotations)),
                                        palette = "Set 3")
col_pbmc <- setNames(col_pbmc, levels(seu$pbmc_annotations))


# ---------------------------------------------------------------------------- #
# a - Cells recovered per sample, faceted by condition
df_ncells <- seu[[]] %>%
  group_by(Sample_ID, Condition) %>%
  summarize(n_cells = n(), .groups = "drop")

# order Sample_ID within each Condition by descending n_cells
sample_order <- df_ncells %>%
  arrange(Condition, desc(n_cells)) %>%
  pull(Sample_ID)

df_ncells <- df_ncells %>%
  mutate(
    Sample_ID = factor(Sample_ID, levels = unique(sample_order)),
    Condition = factor(Condition, levels = cond_levels)
  )

p_a <- ggplot(df_ncells, aes(x = Sample_ID, y = n_cells, fill = Condition)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = MIN_CELLS_TARGET, linetype = "dashed") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_log10(labels = scales::comma) +
  facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
  theme_bw() +
  labs(x = "Sample ID", y = "Total cells (log10 scale)", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 90))

ggsave(file.path(out_dir, "extended_data_3a.pdf"), p_a,
       width = 200, height = 80, units = "mm")


# ---------------------------------------------------------------------------- #
# b - Cell-type composition per sample
df_comp <- seu[[]] %>%
  dplyr::count(Sample_ID, pbmc_annotations, name = "n") %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::mutate(frac = n / sum(n)) %>%
  dplyr::ungroup()

p_b <- ggplot(df_comp, aes(x = Sample_ID, y = n, fill = pbmc_annotations)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = col_pbmc) +
  theme_bw() +
  labs(x = "Sample ID", y = "Cell fraction", fill = "Cell annotation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file.path(out_dir, "extended_data_3b.pdf"), p_b,
       width = 200, height = 80, units = "mm")


# ---------------------------------------------------------------------------- #
# c - Harmony UMAP by experimental batch
p_c <- DimPlot(seu, label = FALSE, shuffle = TRUE,
               reduction = "umap.harmony", group.by = "Batch") +
  scale_color_brewer(palette = "Set2") +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS) +
  coord_fixed()

ggsave(file.path(out_dir, "extended_data_3c.pdf"), p_c,
       width = 500, height = 500, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# d - Harmony UMAP by donor
patient_cols <- colorspace::qualitative_hcl(n = length(unique(seu$Patient)),
                                            palette = "Set 2")

p_d <- DimPlot(seu, label = FALSE, shuffle = TRUE,
               reduction = "umap.harmony", group.by = "Patient") +
  scale_color_manual(values = patient_cols) +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS) +
  coord_fixed() +
  NoLegend()

ggsave(file.path(out_dir, "extended_data_3d.pdf"), p_d,
       width = 500, height = 500, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# e - Harmony UMAP by disease-modifying therapy
# Both healthy donors and untreated MS donors fall in the "Untreated" category.
p_e <- DimPlot(seu, label = FALSE, shuffle = TRUE,
               reduction = "umap.harmony", group.by = "Patient_Treatment") +
  scale_color_brewer(palette = "Set3", direction = -1) +
  FontSize(x.text = font_sizeS, y.text = font_sizeS,
           x.title = font_sizeS, y.title = font_sizeS) +
  coord_fixed()

ggsave(file.path(out_dir, "extended_data_3e.pdf"), p_e,
       width = 500, height = 500, units = "mm", dpi = 600)


# ---------------------------------------------------------------------------- #
# f - Harmony UMAP by age
p_f <- FeaturePlot(seu, features = "Patient_Age", reduction = "umap.harmony") +
  scale_color_viridis_c(option = "turbo") +
  coord_fixed()

ggsave(file.path(out_dir, "extended_data_3f.pdf"), p_f,
       width = 500, height = 500, units = "mm", dpi = 600)
