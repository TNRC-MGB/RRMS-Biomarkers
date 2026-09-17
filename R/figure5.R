
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : zenodo/qpcr/qpcr_master.csv
#          data/genesets/pebv.RDS         # Arvey 2012 EBV factor atlas
#          data/genesets/Bcell_gene_universe.RDS
#          zenodo/scdist/scdist_bcell_prerelapse_v_remission.RDS
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")
source("R/enrichment_functions.R")   # calc_fora, plot_fora_cluster, ebv_stage_*

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")
ENRICH_DIR <- "data/genesets"
out_dir <- "Intermediate"
dir.create(out_dir, FALSE, TRUE)


# ---------------------------------------------------------------------------- #
# Settings and gene panel

genes <- c("balf1", "bglf4", "bglf5", "bnlf2a", "brlf1", "bzlf1", "eber1",
           "eber2", "ebna1", "ebna2-1", "ebna2-2", "ebna3a", "gp350", "gp42",
           "lmp1")
cols <- c(Remission = "#377EB8", `Pre-relapse` = "#E41A1C")
ct_matrix <- "zenodo/qpcr/qpcr_master.csv"
out_dir <- file.path("intermediate", "qpcr_output")
dir.create(out_dir, FALSE, TRUE)


# ---------------------------------------------------------------------------- #
# Average technical replicates, then B2M-normalize to dCt and expression

master <- read.csv(ct_matrix, stringsAsFactors = FALSE)
stopifnot("b2m" %in% master$Gene, all(genes %in% master$Gene))

cells <- master |>
  group_by(Donor_ID, Condition, Gene) |>
  summarize(Ct = mean(Ct), .groups = "drop")

dct <- cells |>
  group_by(Donor_ID, Condition) |>
  mutate(B2M_Ct = Ct[Gene == "b2m"],
         dCt = Ct - B2M_Ct,
         expression = 40 - dCt) |>
  ungroup() |>
  filter(Gene != "b2m") |>
  mutate(Condition = factor(Condition, c("Remission", "Pre-relapse")))
a <- dct |> filter(Gene %in% genes) |> mutate(Gene = factor(Gene, genes))

a_sum <- a |>
  group_by(Gene, Condition) |>
  summarize(n = n(), median = median(expression),
            q1 = quantile(expression, 0.25), q3 = quantile(expression, 0.75),
            .groups = "drop")


# ---------------------------------------------------------------------------- #
# Univariable logistic regression per target. Odds of being a pre-relapse sample
# per unit increase in 40 - dCt

fit <- a |>
  mutate(case = as.integer(Condition == "Pre-relapse")) |>
  group_by(Gene) |>
  group_modify(\(d, g) {
    m <- glm(case ~ expression, binomial(), data = d)
    co <- summary(m)$coefficients["expression", ]
    ci <- suppressMessages(confint(m, "expression"))
    tibble(
      n_donors = n_distinct(d$Donor_ID), n_samples = nrow(d),
      OR = unname(exp(co["Estimate"])),
      CI_low = unname(exp(ci[1])), CI_high = unname(exp(ci[2])),
      z = unname(co["z value"]),
      P = unname(co["Pr(>|z|)"])
    )
  }) |>
  ungroup() |>
  mutate(
    P_BH = p.adjust(P, "BH"),
    sig = case_when(
      P_BH < 0.001 ~ "***",
      P_BH < 0.01 ~ "**",
      P_BH < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


# ---------------------------------------------------------------------------- #
# PCA of donor dCt profiles

x <- dct |>
  select(Donor_ID, Condition, Gene, dCt) |>
  pivot_wider(names_from = Gene, values_from = dCt)
pca <- prcomp(select(x, -Donor_ID, -Condition), center = TRUE, scale. = FALSE)
scores <- bind_cols(select(x, Donor_ID, Condition),
                    as_tibble(pca$x[, 1:2, drop = FALSE]))
ve <- summary(pca)$importance[2, 1:2] * 100


# ---------------------------------------------------------------------------- #
# Panels

th <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))

pa <- ggplot(a_sum, aes(Gene, median, color = Condition, group = Condition)) +
  geom_pointrange(aes(ymin = q1, ymax = q3), position = position_dodge(0.55),
                  linewidth = 0.4) +
  geom_text(data = fit, aes(Gene, 39, label = sig), inherit.aes = FALSE,
            color = "black") +
  scale_color_manual(values = cols, breaks = c("Remission", "Pre-relapse")) +
  scale_x_discrete(position = "top") +
  coord_cartesian(ylim = c(22, 40)) +
  labs(title = "EBV targets in paired MS pre-relapse and remission samples",
       x = NULL, y = "40 - dCt", color = "Condition") +
  th +
  theme(axis.ticks.x = element_blank(), legend.position = "right")

pb <- ggplot(scores, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.75) +
  scale_color_manual(values = cols, breaks = c("Pre-relapse", "Remission"),
                     name = NULL) +
  labs(x = sprintf("PC1: %.2f%% variance explained", ve[1]),
       y = sprintf("PC2: %.2f%% variance explained", ve[2])) +
  th +
  theme(legend.position = "top")

pc <- ggplot(fit, aes(OR, Gene, color = P)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_color_gradient2(low = "#B40426", mid = "#F7F7F7", high = "#3B4A5A",
                        midpoint = 0.05, limits = c(0, 0.1),
                        oob = scales::squish, name = "p value") +
  labs(x = "Odds ratio", y = NULL) +
  th


# ---------------------------------------------------------------------------- #
# Save

ggsave(file.path(out_dir, "p_median_iqr_goi.pdf"), pa,
       width = 170, height = 60, units = "mm")
ggsave(file.path(out_dir, "pca_sample.pdf"), pb,
       width = 80, height = 70, units = "mm")
ggsave(file.path(out_dir, "p_logr.pdf"), pc,
       width = 80, height = 70, units = "mm")
write_csv(left_join(dct, scores, by = c("Donor_ID", "Condition")),
          file.path(out_dir, "figure5_abc_source_data.csv"))
write_csv(select(fit, -sig), file.path(out_dir, "figure5_abc_statistics.csv"))


# ---------------------------------------------------------------------------- #
# d  ABC host response to EBV genes
# Fold-over-representation of EBV factor-derived host gene sets among the genes
# that distinguish pre-relapse ABC cells.
#

ABC_CLUSTER <- "ABC"

ebv_atlas <- readRDS(file.path(ENRICH_DIR, "pebv.RDS"))   # Arvey 2012 atlas
bcell_pre <- readRDS(file.path(SCDIST_DIR,
                               "scdist_bcell_prerelapse_v_remission.RDS"))
universe  <- readRDS(file.path(ENRICH_DIR, "Bcell_gene_universe.RDS"))

bcell_clusters <- names(bcell_pre$vals)
if (!ABC_CLUSTER %in% bcell_clusters) {
  stop("Cluster '", ABC_CLUSTER, "' not in the scDist object. Clusters found: ",
       paste(bcell_clusters, collapse = ", "),
       "\n  If these are numeric, this is the wrong scDist object; the panel ",
       "needs the annotation-labeled B-cell fit ",
       "(scdist_bcell_prerelapse_v_remission.RDS).")
}

fora.res <- calc_fora(
  bcell_pre,
  pathways       = ebv_atlas,
  clusters       = bcell_clusters,
  universe       = universe,
  dist.threshold = NULL,
  direction      = "abs",
  top_prop       = 0.05,
  min_genes      = 25
)

d_abc <- fora.res %>% dplyr::filter(cell_type == ABC_CLUSTER)
message("Panel d: ", nrow(d_abc), " EBV gene sets tested in ", ABC_CLUSTER,
        "; ", sum(d_abc$padj < 0.05), " significant at FDR < 0.05; ",
        "top 30 plotted")

p_d <- plot_fora_cluster(fora.res, cluster = ABC_CLUSTER) +
  theme(text = element_text(size = 5))

# EBV lifecycle-stage coloring on the y-axis labels
stage_levels <- names(ebv_stage_colors)
axis_paths   <- levels(p_d$data$pathway)
axis_stages  <- ebv_lifecycle_stage(axis_paths)

if (anyNA(axis_stages)) {
  warning("No EBV lifecycle stage for: ",
          paste(axis_paths[is.na(axis_stages)], collapse = ", "),
          ". Those labels render gray; add them to the stage vectors in ",
          "R/enrichment_functions.R.")
}
axis_cols <- unname(ebv_stage_colors[axis_stages])
axis_cols[is.na(axis_cols)] <- "gray30"

# Zero-height dummy layer to generate the stage key
stage_key <- data.frame(
  neglog10 = NA_real_,
  pathway  = factor(axis_paths[1], levels = axis_paths),
  stage    = factor(stage_levels, levels = stage_levels)
)

p_d <- p_d +
  geom_point(data = stage_key,
             aes(x = neglog10, y = pathway, fill = stage),
             shape = 22, size = 2, color = NA, na.rm = TRUE,
             inherit.aes = FALSE) +
  scale_fill_manual(name   = "EBV lifecycle stage",
                    values = ebv_stage_colors,
                    breaks = stage_levels,
                    labels = ebv_stage_labels[stage_levels],
                    guide  = guide_legend(order = 1,
                                          override.aes = list(size = 3))) +
  theme(axis.text.y = element_text(color = axis_cols, face = "bold"))

ggsave(file.path(out_dir, "figure_5d_ABC_EBV_fora.pdf"), p_d,
       width = 80, height = 100, units = "mm")
readr::write_csv(
  d_abc %>% dplyr::mutate(stage = ebv_lifecycle_stage(pathway)),
  file.path(out_dir, "figure_5d_ABC_EBV_fora_table.csv")
)


# ---------------------------------------------------------------------------- #
# e  CD19+ bulk RNA-seq: EBV LMP-1 host signature
# see rrms/bulk