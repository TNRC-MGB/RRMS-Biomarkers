
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-22-2026
#
#
# Inputs : pbmc_final.RDS                (scrna preprocessing pipeline)
#          scdist_pbmc_pre/scdist_pbmc_rel/scdist_pbmc_ms_v_healthy.RDS
#            (scrna_pbmc_scdist.R)
#          scdrs/output/pbmc/MS.score.gz (scDRS, MS GWAS - see scdrs/README.md)
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms") 

source("R/packages.R")

set.seed(42)

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

source("R/scrna_preprocessing_functions.R")   # make_impadt_heatmap_df, plot_impadt_heatmap
source("R/scdist_functions.R")                 # DistPlot2, add_scdist_to_seurat

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# scDRS output for the MS GWAS produced by the pipeline in scdrs/
SCDRS_SCORE <- "scdrs/output/pbmc/MS.score.gz"

fs <- 5   # axis text size (pt)
fl <- 5   # axis title size (pt)

# --- Inputs -------------------------------------------------------------------
pbmc <- readRDS(file.path(OBJ_DIR, "pbmc_final.RDS"))

scdist_pbmc_ms_hc   <- readRDS(file.path(SCDIST_DIR, "scdist_pbmc_ms_v_healthy.RDS"))     # MS vs Healthy
scdist_pbmc_rel <- readRDS(file.path(SCDIST_DIR, "scdist_pbmc_relapse_v_remission.RDS"))   # Relapse vs Remission
scdist_pbmc_pre <- readRDS(file.path(SCDIST_DIR, "scdist_pbmc_prerelapse_v_remission.RDS"))   # PreRelapse vs Remission

# scDist distance scored onto each cell (for the panel c UMAPs)
pbmc <- add_scdist_to_seurat(pbmc, scdist_pbmc_ms_hc,   cluster_col = "pbmc_annotations", new_col = "scdist.ms_v_healthy")
pbmc <- add_scdist_to_seurat(pbmc, scdist_pbmc_rel, cluster_col = "pbmc_annotations", new_col = "scdist.relapse")
pbmc <- add_scdist_to_seurat(pbmc, scdist_pbmc_pre, cluster_col = "pbmc_annotations", new_col = "scdist.pre_relapse")

# scDRS scores into metadata (for panels e/f/g)
sc <- fread(SCDRS_SCORE)
setnames(sc, 1, "cell")
stopifnot(identical(sc$cell, colnames(pbmc)))
pbmc <- AddMetaData(pbmc, metadata = as.data.frame(sc), col.name = NULL)


# ---------------------------------------------------------------------------- #
# a - PBMC UMAP atlas
p_a <- DimPlot(
  pbmc, reduction = "umap.harmony", group.by = "pbmc_annotations",
  label = TRUE, repel = TRUE, raster = TRUE, label.size = 2,
  raster.dpi = c(512, 512), pt.size = 1
) +
  scale_color_discrete_qualitative(palette = "Dynamic") +
  theme_bw() + coord_fixed() + NoLegend() + NoGrid() +
  FontSize(x.text = fs, y.text = fs, x.title = fl, y.title = fl, main = 0)

ggsave(file.path(out_dir,"figure_2a.pdf"), p_a, width = 60, height = 60, units = "mm")


# ---------------------------------------------------------------------------- #
# b - Imputed-ADT marker heatmaps (lymphoid | B/myeloid)
adt_markers_lymphoid <- c(
  "CD3-1", "CD4-1", "CD8", "CD25", "CD152", "TIGIT",
  "CD45RA", "CD45RO", "CD27", "CD161", "TCR-V-7.2", "TCR-1",
  "CD56-1", "CD16", "CD57", "CD244", "CD335", "CD122", "CD127"
)
adt_markers_myeloid_b <- c(
  "CD19", "CD20", "CD22", "IgD", "IgM", "CD38-1", "CD138-1",
  "CD123", "CD303", "CD304", "CD11c", "CD1c", "HLA-DR", "CD86", "CD83",
  "CD14", "CD64", "CX3CR1", "CD163"
)

df_lymphoid   <- make_impadt_heatmap_df(pbmc, adt_markers_lymphoid,   group_col = "pbmc_annotations", assay = "impADT", layer = "data")
df_myeloid_b  <- make_impadt_heatmap_df(pbmc, adt_markers_myeloid_b,  group_col = "pbmc_annotations", assay = "impADT", layer = "data")

p_lymphoid <- plot_impadt_heatmap(df_lymphoid, title = "lymphoid markers") +
  coord_fixed() +
  FontSize(x.text = fs, y.text = fs, x.title = fl, y.title = fl, main = fl) +
  NoLegend()

p_myeloid_b <- plot_impadt_heatmap(df_myeloid_b, title = "B/myeloid markers") +
  coord_fixed() +
  FontSize(x.text = fs, y.text = fs, x.title = fl, y.title = fl, main = fl) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.title = element_text(size = fl), legend.text = element_text(size = fs))

p_b <- p_lymphoid + p_myeloid_b
ggsave(file.path(out_dir, "figure_2b.pdf"), p_b, width = 120, height = 60, units = "mm")


# ---------------------------------------------------------------------------- #
# c - scDist-distance UMAPs (MS v Healthy | Relapse v Remission | PreRelapse v Remission)
scdist_umap <- function(feature) {
  FeaturePlot(pbmc, features = feature, reduction = "umap.harmony") +
    coord_fixed() +
    scale_color_viridis_c(limits = c(1, 10), oob = scales::squish) +
    FontSize(x.text = fs, y.text = fs, x.title = fl, y.title = fl, main = fl) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          legend.title = element_text(size = fl), legend.text = element_text(size = fs)) +
    NoLegend()
}

p_c <- scdist_umap("scdist.ms_v_healthy") +
  scdist_umap("scdist.relapse") +
  scdist_umap("scdist.pre_relapse") +
  patchwork::plot_layout(ncol = 3)

ggsave(file.path(out_dir, "figure_2c.pdf"), p_c, width = 174, height = 70, units = "mm")


# ---------------------------------------------------------------------------- #
# d - scDist-distance dotplots (same three contrasts)
YLIM <- c(0, 12)

scdist_dot <- function(res) {
  DistPlot2(res, point_size = 1) +
    ylim(YLIM) +
    FontSize(x.text = fs, y.text = fs, x.title = fl, y.title = fl, main = fl) +
    theme(axis.title.y = element_blank(),
          legend.title = element_text(size = fl), legend.text = element_text(size = fs)) +
    NoLegend()
}

p_d <- scdist_dot(scdist_pbmc_ms_hc) + scdist_dot(scdist_pbmc_rel) + scdist_dot(scdist_pbmc_pre)
ggsave(file.path(out_dir, "figure_2d.pdf"), p_d, width = 174, height = 50, units = "mm")


# ---------------------------------------------------------------------------- #
# e - Patient-aware scDRS (emmeans) per cell type
md <- as.data.table(pbmc[[]])

df_psc <- md[, .(n_cells = .N, mean_norm = mean(norm_score, na.rm = TRUE)),
             by = .(Patient, Sample_ID, pbmc_annotations)]
df_psc <- df_psc[n_cells >= 5]
df_psc[, Patient := factor(Patient)]
df_psc[, pbmc_annotations := factor(pbmc_annotations, levels = levels(pbmc$pbmc_annotations))]

fit_scdrs <- lmer(mean_norm ~ 0 + pbmc_annotations + (1 | Patient), data = df_psc, weights = n_cells)

tab_scdrs <- as.data.table(summary(emmeans(fit_scdrs, ~ pbmc_annotations), infer = c(TRUE, TRUE)))
tab_scdrs[, p_adj := p.adjust(p.value, method = "BH")]
tab_scdrs[, Label := as.character(pbmc_annotations)]
tab_scdrs[, sig := p_adj < 0.05]
tab_scdrs[, ylab := factor(Label, levels = tab_scdrs[order(emmean)]$Label)]

p_e <- ggplot(tab_scdrs, aes(x = emmean, y = ylab)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_segment(aes(x = lower.CL, xend = upper.CL, yend = ylab), linewidth = 0.5) +
  geom_point(aes(color = sig), size = 1) +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "gray55")) +
  labs(x = "Patient-aware mean scDRS norm score ± 95% CI", y = NULL, color = "FDR < 0.05") +
  theme_classic() + theme(text = element_text(size = 5))

ggsave(file.path(out_dir, "figure_2e.pdf"), p_e, width = 65, height = 45, units = "mm")


# ---------------------------------------------------------------------------- #
# f - Mean scDRS score per cluster on UMAP
avg_map <- tapply(pbmc$norm_score, pbmc$pbmc_annotations, mean, na.rm = TRUE)
pbmc$avg_scdrs_by_cluster <- unname(avg_map[as.character(pbmc$pbmc_annotations)])

p_f <- FeaturePlot(pbmc, features = "avg_scdrs_by_cluster", reduction = "umap.harmony", order = FALSE) +
  coord_fixed() +
  scale_color_viridis_c(option = "inferno", limits = c(0, 1), oob = scales::squish) +
  theme_classic() + theme(text = element_text(size = 5))

ggsave(file.path(out_dir, "figure_2f.pdf"), p_f, width = 60, height = 60, units = "mm")


# ---------------------------------------------------------------------------- #
# g - scDist distance vs PreRelapse-Remission scDRS change (by lineage)
res <- scdist_pbmc_pre$results

annotation_to_lineage <- c(
  "CD4 Naive/TCM" = "CD4 T", "CD4 TCM" = "CD4 T", "Activated Treg" = "CD4 T",
  "CD8 Naive" = "CD8 T", "CD8 Naive-like" = "CD8 T", "CD8 Naive/TCM" = "CD8 T", "CD8 Tem (GZMK+)" = "CD8 T",
  "MAIT" = "Innate-like T/NK", "Cytotoxic gdT" = "Innate-like T/NK", "FGFBP2+ Cytotoxic gdT" = "Innate-like T/NK",
  "CD56bright NK" = "Innate-like T/NK", "NK" = "Innate-like T/NK", "Mature NK" = "Innate-like T/NK",
  "Cycling Cytotoxic Lymphocytes" = "Innate-like T/NK",
  "TCL1A+ Naive B" = "B cell", "Naive B" = "B cell", "Intermediate B" = "B cell",
  "Switched Memory B" = "B cell", "Plasmablast" = "B cell",
  "pDC" = "DC", "Activated cDC2" = "DC", "CD14 Mono" = "Monocyte", "CD16 Mono" = "Monocyte"
)
lineage_levels <- c("CD4 T", "CD8 T", "Innate-like T/NK", "B cell", "DC", "Monocyte")

# within-patient PreRelapse - Remission change in mean scDRS, per cell type
df_g <- pbmc[[]] %>%
  dplyr::filter(Condition %in% c("Remission", "PreRelapse")) %>%
  dplyr::group_by(Patient, Condition, pbmc_annotations) %>%
  dplyr::summarize(mean_scdrs = mean(norm_score, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(id_cols = c(Patient, pbmc_annotations),
                     names_from = Condition, values_from = mean_scdrs) %>%
  dplyr::filter(!is.na(Remission), !is.na(PreRelapse)) %>%
  dplyr::mutate(delta_scdrs = PreRelapse - Remission) %>%
  dplyr::group_by(pbmc_annotations) %>%
  dplyr::summarize(delta_scdrs = mean(delta_scdrs), .groups = "drop") %>%
  dplyr::left_join(
    data.frame(pbmc_annotations = rownames(res),
               scdist      = res[["Dist."]],
               scdist_low  = res[["95% CI (low)"]],
               scdist_high = res[["95% CI (upper)"]]),
    by = "pbmc_annotations"
  ) %>%
  dplyr::mutate(lineage = factor(annotation_to_lineage[pbmc_annotations], levels = lineage_levels))

p_g <- ggplot(df_g, aes(x = scdist, y = delta_scdrs, color = lineage, label = pbmc_annotations)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = scdist_low, xmax = scdist_high), height = 0, linewidth = 0.25, alpha = 0.8) +
  geom_point(size = 1) +
  ggrepel::geom_text_repel(size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "scDist distance", y = "PreRelapse - Remission mean scDRS", color = "Lineage") +
  theme_classic() + theme(text = element_text(size = 5))

ggsave(file.path(out_dir, "figure_2g.pdf"), p_g, width = 70, height = 60, units = "mm")
