
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : R/facs_flowsom_pipeline.R
#            (FCS import -> FlowSOM -> UMAP -> gp350+ -> GLMM;
#             shared with figure4.R)
#
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

set.seed(42)

# Shared with figure4.R
source("R/facs_flowsom_pipeline.R")

# --- Settings -----------------------------------------------------------------
MARKER_LIMS <- c(0, 7)   # shared arcsinh color scale across panel b


# ---------------------------------------------------------------------------- #
# b - Ten-marker UMAP grid
markers_to_plot <- c(
  "FJComp-BV605-A_CD5_asinh",
  "FJComp-PE-Cy5-A_CD11C_asinh",
  "FJComp-BUV661-A_CD18_asinh",
  "FJComp-BB700-A_CD19_asinh",
  "FJComp-BB515-A_CD20_asinh",
  "FJComp-BV711-A_CD21_asinh",
  "FJComp-APC-Cy7-A_CD27_asinh",
  "FJComp-BV786-A_CXCR4_asinh",
  "FJComp-BUV496-A_CD38_asinh",
  "FJComp-BUV563-A_CD172_asinh"
)
stopifnot(all(markers_to_plot %in% names(b.m)))

marker_plots <- lapply(markers_to_plot, function(i) {
  cols_to_select <- c(i, "UMAP_X", "UMAP_Y")
  b.i <- as.data.frame(b.m[, ..cols_to_select])
  ggplot(b.i, aes(x = UMAP_X, y = UMAP_Y, color = get(i))) +
    ggrastr::geom_point_rast(size = 0.1, raster.dpi = 600) +
    scale_color_viridis_c(option = "inferno", limits = MARKER_LIMS,
                          oob = scales::squish) +
    coord_fixed() +
    labs(color = pretty_marker(i)) +
    theme_void()
})

p_b <- wrap_plots(marker_plots, ncol = 5)
ggsave(file.path(out_dir, "extended_data_8b.pdf"), p_b,
       width = 360, height = 360, units = "mm")


# ---------------------------------------------------------------------------- #
# c - Annotated FlowSOM metacluster UMAP
centroids <- b.m[, .(
  UMAP_X = mean(UMAP_X, na.rm = TRUE),
  UMAP_Y = mean(UMAP_Y, na.rm = TRUE)
), by = FlowSOM_metacluster]
centroids$annotation <- ca[as.character(centroids$FlowSOM_metacluster)]

p_c <- ggplot(b.m, aes(x = UMAP_X, y = UMAP_Y,
                       color = as.character(FlowSOM_metacluster))) +
  ggrastr::geom_point_rast(size = 0.1, raster.dpi = 600) +
  scale_color_manual(values = colors) +
  coord_fixed() +
  ggrepel::geom_text_repel(data = centroids, aes(label = annotation),
                           color = "black", size = 4) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "extended_data_8c.pdf"), p_c,
       width = 120, height = 120, units = "mm")


# ---------------------------------------------------------------------------- #
# d - Model-predicted % gp350+ by condition, stratified by run day
# Estimated marginal means from the binomial GLMM fitted in the pipeline
# (Group + exp_day fixed, donor random intercept). gp350 positivity thresholds
# were set per run day at a 1% false-positive rate.
emm_day <- emmeans(fit, ~ Group | exp_day, type = "response")

emm_day_df <- as.data.frame(emm_day) %>%
  mutate(pct = 100 * prob, lo = 100 * asymp.LCL, hi = 100 * asymp.UCL)

p_d <- ggplot(emm_day_df, aes(x = Group, y = pct)) +
  geom_point() +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15) +
  facet_wrap(~ exp_day) +
  labs(x = NULL, y = "Predicted % gp350+")

ggsave(file.path(out_dir, "extended_data_8d.pdf"), p_d,
       width = 60, height = 60, units = "mm")


# ---------------------------------------------------------------------------- #
# e - Odds-ratio forest, pre-relapse binned by proximity to relapse onset
per_sample2 <- per_sample %>%
  mutate(
    Timing_days = as.numeric(Timing),
    Timing_bin = case_when(
      Group == "pre_relapse" & !is.na(Timing_days) & Timing_days >= -90  ~ '0\u201390d',
      Group == "pre_relapse" & !is.na(Timing_days) & Timing_days <  -90 &
                                                     Timing_days >= -150 ~ '90\u2013150d',
      Group == "remission"                                               ~ "remission",
      TRUE                                                               ~ NA_character_
    ),
    Timing_bin = factor(Timing_bin, levels = c("remission", '90\u2013150d', '0\u201390d'))
  )

fit_bin <- glmer(
  cbind(n_gp350, n_B - n_gp350) ~ Timing_bin + exp_day + (1 | Patient_ID),
  data   = per_sample2,
  family = binomial
)

or_tbl_bin <- broom.mixed::tidy(fit_bin, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    OR    = exp(estimate),
    OR_lo = exp(conf.low),
    OR_hi = exp(conf.high),
    term_pretty = dplyr::recode(
      term,
      'Timing_bin0\u201390d'   = 'Pre-relapse 0\u201390d vs remission',
      'Timing_bin90\u2013150d' = 'Pre-relapse 90\u2013150d vs remission',
      "exp_dayexp2"            = "Run day: exp2 vs exp1",
      "exp_dayexp3"            = "Run day: exp3 vs exp1"
    )
  )

p_e <- ggplot(or_tbl_bin, aes(x = term_pretty, y = OR, ymin = OR_lo, ymax = OR_hi)) +
  geom_pointrange() +
  coord_flip() +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(x = NULL, y = "Odds ratio (log scale)")

ggsave(file.path(out_dir, "extended_data_8e.pdf"), p_e,
       width = 60, height = 60, units = "mm")
