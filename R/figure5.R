
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : data/qpcr/qpcr_well_map.csv    # 384-well plate map, de-identified
#          zenodo/qpcr_ampdata/20251218_RelapseBio_P1.csv .. _P4.csv
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
# qPCR Settings, gene panel and QC

qpcr <- file.path(ZENODO_DIR, "qpcr_ampdata")
genes <- c("balf1", "bglf4", "bglf5", "bnlf2a", "brlf1", "bzlf1", "eber1",
           "eber2", "ebna1", "ebna2-1", "ebna2-2", "ebna3a", "gp350", "gp42",
           "lmp1")
drop <- c("MS050", "MS085", "MS047", "MS030", "MS086")   # failed QC
cols <- c(Remission = "#377EB8", `Pre-relapse` = "#E41A1C")

# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Well map and amplification data
wm <- read.csv("data/qpcr/qpcr_well_map.csv", stringsAsFactors = FALSE) |>
  filter(sample_type != "NTC")
stopifnot(all(drop %in% wm$Donor_ID))

amp <- map_dfr(paste0("P", 1:4), \(p) read.csv(
  file.path(qpcr, paste0("20251218_RelapseBio_", p, ".csv")),
  stringsAsFactors = FALSE))

ct_col <- intersect(c("CT", "Ct", "Cq", "cq", "ct"), names(amp))[1]
stopifnot(!is.na(ct_col))


# ---------------------------------------------------------------------------- #
# Reliable observations, paired donors only
raw <- amp |>
  left_join(select(wm, well, Donor_ID, category),
            by = c("Well.Position" = "well")) |>
  transmute(Donor_ID,
            Condition = recode(category, relapse = "Pre-relapse",
                               remission = "Remission",
                               .default = NA_character_),
            Gene = str_to_lower(Target.Name),
            Ct_text = as.character(.data[[ct_col]]),
            Ct = suppressWarnings(as.numeric(Ct_text))) |>
  filter(!is.na(Donor_ID), !is.na(Condition), !is.na(Gene)) |>
  mutate(reliable = !is.na(Ct) & Ct <= 35 & Ct_text != "Undetermined") |>
  group_by(Donor_ID) |>
  filter(n_distinct(Condition) == 2) |>
  ungroup()


# ---------------------------------------------------------------------------- #
# Donor-level and target-level QC exclusions
sample_qc <- raw |>
  group_by(Donor_ID) |>
  summarize(unreliable_fraction = mean(!reliable), .groups = "drop")
gene_qc <- raw |>
  group_by(Condition, Gene) |>
  summarize(unreliable_fraction = mean(!reliable), .groups = "drop")
bad_donors <- union(drop, filter(sample_qc, unreliable_fraction > .5)$Donor_ID)
bad_genes <- unique(filter(gene_qc, unreliable_fraction > .8)$Gene)


# ---------------------------------------------------------------------------- #
# Technical-well means and condition-mean imputation
dat <- raw |> filter(reliable, !Donor_ID %in% bad_donors, !Gene %in% bad_genes)
pairs <- dat |>
  distinct(Donor_ID, Condition) |>
  count(Donor_ID) |>
  filter(n == 2) |>
  pull(Donor_ID)
obs <- dat |>
  filter(Donor_ID %in% pairs) |>
  group_by(Donor_ID, Condition, Gene) |>
  summarize(Ct = mean(Ct), technical_wells = n(), .groups = "drop")
wide <- obs |>
  select(-technical_wells) |>
  pivot_wider(names_from = Gene, values_from = Ct) |>
  group_by(Condition) |>
  mutate(across(where(is.numeric),
                \(x) replace_na(x, mean(x, na.rm = TRUE)))) |>
  ungroup()
stopifnot("b2m" %in% names(wide), all(genes %in% names(wide)))


# ---------------------------------------------------------------------------- #
# dCt against B2M and expression on the 40 - dCt scale
dct <- wide |>
  pivot_longer(-c(Donor_ID, Condition), names_to = "Gene",
               values_to = "Ct_mean") |>
  left_join(select(obs, Donor_ID, Condition, Gene, technical_wells),
            by = c("Donor_ID", "Condition", "Gene")) |>
  mutate(imputed = is.na(technical_wells),
         technical_wells = replace_na(technical_wells, 0L)) |>
  group_by(Donor_ID, Condition) |>
  mutate(B2M_Ct = Ct_mean[Gene == "b2m"],
         B2M_imputed = imputed[Gene == "b2m"],
         dCt = Ct_mean - B2M_Ct,
         expression = 40 - dCt) |>
  ungroup() |>
  filter(Gene != "b2m") |>
  mutate(Condition = factor(Condition, c("Remission", "Pre-relapse")))
a <- dct |> filter(Gene %in% genes) |> mutate(Gene = factor(Gene, genes))

a_sum <- a |>
  group_by(Gene, Condition) |>
  summarize(n = n(), median = median(expression),
            q1 = quantile(expression, .25), q3 = quantile(expression, .75),
            .groups = "drop")


# ---------------------------------------------------------------------------- #
# Univariable logistic regression per target
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
      P_BH < .001 ~ "***",
      P_BH < .01 ~ "**",
      P_BH < .05 ~ "*",
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
  geom_pointrange(aes(ymin = q1, ymax = q3), position = position_dodge(.55),
                  linewidth = .4) +
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
  geom_point(size = 3, alpha = .75) +
  scale_color_manual(values = cols, breaks = c("Pre-relapse", "Remission"),
                     name = NULL) +
  labs(x = sprintf("PC1: %.2f%% variance explained", ve[1]),
       y = sprintf("PC2: %.2f%% variance explained", ve[2])) +
  th +
  theme(legend.position = "top")

pc <- ggplot(fit, aes(OR, Gene, color = P)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = .15) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_color_gradient2(low = "#B40426", mid = "#F7F7F7", high = "#3B4A5A",
                        midpoint = .05, limits = c(0, .1),
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
write_csv(bind_rows(
  transmute(sample_qc, level = "donor", item = Donor_ID,
            Condition = NA_character_, unreliable_fraction,
            excluded = Donor_ID %in% bad_donors),
  transmute(gene_qc, level = "gene", item = Gene, Condition,
            unreliable_fraction, excluded = Gene %in% bad_genes)),
  file.path(out_dir, "figure5_qpcr_qc.csv"))
n <- n_distinct(a$Donor_ID)
capture.output(sessionInfo(),
               file = file.path(out_dir, "figure5_abc_sessionInfo.txt"))


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