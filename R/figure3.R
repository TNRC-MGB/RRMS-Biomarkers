
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-22-2026
#
#
# Inputs : B_clean_scdrs.RDS  (B-cell object + B_annotations + scDRS norm_score)
#          scdist_bcell_prerelapse_v_remission.RDS
#            (scDist, PreRelapse vs Remission)
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms") 

source("R/packages.R")

set.seed(42)

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

source("R/scdist_functions.R")         # DistPlot2

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Inputs -------------------------------------------------------------------
B_clean <- readRDS(file.path(OBJ_DIR, "B_clean_scdrs.RDS"))
scdist_bcell_pre <- readRDS(file.path(SCDIST_DIR, "scdist_bcell_prerelapse_v_remission.RDS"))


# ---------------------------------------------------------------------------- #
# a - B-cell UMAP (B_annotations)
B_annotation_cols <- c(
  "Transitional/Immature-like B" = "#80B1D3",
  "TCL1A+ Naive B"               = "#CCEBC5",
  "Naive B"                      = "#B3DE69",
  "Anergic/Naive-leaning B"      = "#FDB462",
  "IFN-stimulated Naive B"       = "#FB8072",
  "Early Activated Memory B"     = "#8DD3C7",
  "Switched Memory B"            = "#FFFF33",
  "Activated Switched Memory B"  = "#D9D9D9",
  "Activated Memory B"           = "#4DAF4A",
  "CD11c++ Activated Memory B"   = "#984EA3",
  "Atypical/ABC-like Memory B"   = "#377EB8",
  "ABC"                          = "#E41A1C",
  "Plasmablast"                  = "#FF7F00",
  "Cycling Plasmablast"          = "#D94801"
)

B_clean$B_annotations <- factor(B_clean$B_annotations, levels = names(B_annotation_cols))

p_b_um <- DimPlot(
  B_clean,
  reduction = "umap.B.harmony",
  group.by  = "B_annotations",
  cols      = B_annotation_cols,
  pt.size   = 2,
  raster    = TRUE,
  label     = TRUE,
  repel     = TRUE
) +
  coord_fixed() +
  theme_minimal() +
  theme(text = element_text(size = 5)) +
  NoLegend()

ggsave(file.path(out_dir, "p_b_um.pdf"), p_b_um, width = 90, height = 60, units = "mm")


# ---------------------------------------------------------------------------- #
# b - RNA marker dot plot
rna_dotplot_markers_main <- c(
  "TCL1A", "IGHD", "FCER2", "MME",
  "CD27", "TNFRSF13B", "FCRL1",
  "CD83", "NR4A1", "AREG",
  "AIM2", "ITGAX", "CD80",
  "PRDM1", "XBP1", "MZB1", "TNFRSF17", "IGHA1",
  "IFI44L", "ISG15"
)

DefaultAssay(B_clean) <- "RNA"
p_b_dp <- DotPlot(
  B_clean,
  features  = rna_dotplot_markers_main,
  group.by  = "B_annotations",
  dot.scale = 2,
  dot.min   = 0.2
) +
  RotatedAxis() +
  coord_fixed() +
  scale_color_viridis_c(option = "turbo") +
  theme(text = element_text(size = 5), axis.title = element_blank()) +
  FontSize(x.text = 5, y.text = 5)

ggsave(file.path(out_dir, "p_b_dp.pdf"), p_b_dp, width = 100, height = 60, units = "mm")


# ---------------------------------------------------------------------------- #
# c - ADT marker UMAPs
markers_for_umap_grid <- c("CD19","CD20","CD21","CD27","CD38-1","CD38-2","HLA-DR","CD11c","IgD","IgM")

DefaultAssay(B_clean) <- "impADT"
p_mk_umap <- FeaturePlot(B_clean, features = markers_for_umap_grid, ncol = 5, reduction = "umap.B.harmony",
                         raster = TRUE, raster.dpi = c(300, 300)) &
  scale_color_viridis_c(option = "inferno") &
  coord_fixed() &
  theme_bw() &
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text  = element_blank(),
        text = element_text(size = 5))

ggsave(file.path(out_dir, "p_mk_umap.pdf"), p_mk_umap, width = 80, height = 80, units = "mm")


# ---------------------------------------------------------------------------- #
# d - scDist distance dotplot (PreRelapse vs Remission)
YLIM <- c(0, 12)

p_d <- DistPlot2(scdist_bcell_pre, point_size = 1) +
  ylim(YLIM) +
  FontSize(x.text = 5, y.text = 5, x.title = 5, y.title = 5, main = 5) +
  theme(axis.title.y = element_blank(),
        legend.title = element_text(size = 5), legend.text = element_text(size = 5)) +
  NoLegend()

ggsave(file.path(out_dir, "figure_3d.pdf"), p_d, width = 58, height = 50, units = "mm")


# ---------------------------------------------------------------------------- #
# e - MS-GWAS scDRS volcano (patient-aware emmeans)
md <- as.data.table(B_clean[[]])

df_b_scdrs <- md[, .(n_cells = .N, mean_norm = mean(norm_score, na.rm = TRUE)),
                 by = .(Patient, B_annotations)]
df_b_scdrs <- df_b_scdrs[n_cells >= 5]
df_b_scdrs[, Patient := factor(Patient)]
df_b_scdrs[, B_annotations := factor(B_annotations, levels = levels(B_clean$B_annotations))]

fit_b_scdrs <- lmer(mean_norm ~ 0 + B_annotations + (1 | Patient), data = df_b_scdrs, weights = n_cells)

tab_b_scdrs <- as.data.table(summary(emmeans(fit_b_scdrs, ~ B_annotations), infer = c(TRUE, TRUE)))
tab_b_scdrs[, p_adj := p.adjust(p.value, method = "BH")]
tab_b_scdrs[, Label := as.character(B_annotations)]
tab_b_scdrs[, sig := p_adj < 0.05]
tab_b_scdrs[, neglog10_fdr := -log10(pmax(p_adj, 1e-300))]

p_b_scdrs_volcano <- ggplot(tab_b_scdrs, aes(x = emmean, y = neglog10_fdr, label = Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") +
  geom_point(aes(color = sig), size = 1.8) +
  ggrepel::geom_text_repel(aes(color = sig), size = 2.2, max.overlaps = Inf) +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "gray55")) +
  labs(x = "Patient-aware mean scDRS norm score (emmean)",
       y = expression(-log[10]("FDR")),
       color = "FDR < 0.05") +
  theme_classic() +
  theme(text = element_text(size = 5))

ggsave(file.path(out_dir, "p_b_scdrs_volcano.pdf"), p_b_scdrs_volcano, width = 76, height = 50, units = "mm")


# ---------------------------------------------------------------------------- #
# f - scDRS change vs scDist distance, prioritized B subsets
md <- as.data.table(B_clean[[]])

df_scdrs_b <- md[Condition %in% c("Remission", "PreRelapse"),
                 .(n_cells = .N, mean_scdrs = mean(norm_score, na.rm = TRUE)),
                 by = .(Patient, Condition, B_annotations)]
df_scdrs_b <- df_scdrs_b[n_cells >= 5]

df_delta_b <- df_scdrs_b %>%
  tidyr::pivot_wider(id_cols = c(Patient, B_annotations),
                     names_from = Condition, values_from = mean_scdrs) %>%
  dplyr::filter(!is.na(Remission), !is.na(PreRelapse)) %>%
  dplyr::mutate(delta_scdrs = PreRelapse - Remission) %>%
  dplyr::group_by(B_annotations) %>%
  dplyr::summarize(delta_scdrs = mean(delta_scdrs), .groups = "drop")

res_b <- scdist_bcell_pre$results
df_scdist_b <- data.frame(
  B_annotations = rownames(res_b),
  scdist      = res_b[["Dist."]],
  scdist_low  = res_b[["95% CI (low)"]],
  scdist_high = res_b[["95% CI (upper)"]],
  stringsAsFactors = FALSE
)

df_int <- df_delta_b %>%
  dplyr::left_join(df_scdist_b, by = "B_annotations") %>%
  dplyr::mutate(B_annotations = factor(B_annotations, levels = levels(B_clean$B_annotations)))

x_cut <- quantile(df_int$scdist, 0.50, na.rm = TRUE)
y_cut <- quantile(df_int$delta_scdrs, 0.50, na.rm = TRUE)
df_int$prioritized <- df_int$scdist >= x_cut & df_int$delta_scdrs >= y_cut

p_b_integrated_delta <- ggplot(df_int, aes(x = scdist, y = delta_scdrs, label = B_annotations)) +
  geom_vline(xintercept = x_cut, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = y_cut, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray75") +
  geom_errorbarh(aes(xmin = scdist_low, xmax = scdist_high), height = 0, linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(color = prioritized), size = 2.3) +
  ggrepel::geom_text_repel(aes(color = prioritized), size = 2.3, max.overlaps = Inf,
                           box.padding = 0.2, point.padding = 0.15) +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "gray50")) +
  labs(x = "scDist distance (PreRelapse vs Remission)",
       y = 'Mean within-patient change in scDRS score\n(PreRelapse vs Remission)',
       color = "Prioritized") +
  theme_classic() +
  theme(text = element_text(size = 5))

ggsave(file.path(out_dir, "figure_3f_bcell_integrated_delta.pdf"), p_b_integrated_delta, width = 70, height = 40, units = "mm")


# ---------------------------------------------------------------------------- #
# Integrated scDist / scDRS table (input to figure6.R panel b)

scdist_stats <- data.frame(
  B_annotations = rownames(res_b),
  scdist_p      = res_b[[if ("p.val" %in% colnames(res_b)) "p.val" else "p.sum"]],
  stringsAsFactors = FALSE
)

df_int_out <- df_int %>%
  dplyr::left_join(
    tab_b_scdrs[, .(B_annotations = Label,
                    scdrs_emmean  = emmean,
                    scdrs_low     = lower.CL,
                    scdrs_high    = upper.CL,
                    scdrs_p       = p.value,
                    scdrs_fdr     = p_adj)],
    by = "B_annotations"
  ) %>%
  dplyr::left_join(scdist_stats, by = "B_annotations") %>%
  dplyr::select(B_annotations, scdrs_emmean, scdrs_low, scdrs_high, scdrs_p,
                scdrs_fdr, scdist, scdist_low, scdist_high, scdist_p,
                delta_scdrs)

saveRDS(df_int_out, file.path(OBJ_DIR, "df_int_regenerated.RDS"))

