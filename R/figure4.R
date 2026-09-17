
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : fcs/cd19+/samples_*_CD19.fcs         (per-sample CD19+ FCS exports)
#          ../metadata/sample_details.csv
#          (pipeline in R/facs_flowsom_pipeline.R)
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")


set.seed(42)

source("R/facs_flowsom_pipeline.R")


# ---------------------------------------------------------------------------- #
# b - FlowSOM metacluster UMAP
centroids <- b.m[, .(
  UMAP_X = mean(UMAP_X, na.rm = TRUE),
  UMAP_Y = mean(UMAP_Y, na.rm = TRUE)
), by = FlowSOM_metacluster]

p_umap <- ggplot(b.m, aes(x = UMAP_X, y = UMAP_Y, color = factor(FlowSOM_metacluster))) +
  ggrastr::geom_point_rast(size = 0.1, raster.dpi = 600) +
  scale_color_manual(values = colors) +
  coord_fixed() +
  ggrepel::geom_text_repel(data = centroids, aes(label = FlowSOM_metacluster),
                           color = "black", size = 4) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "umap.pdf"), p_umap, width = 101.6, height = 101.6, units = "mm")


# ---------------------------------------------------------------------------- #
# c - Core / identity / overlay marker heatmaps
# overlay markers = all asinh markers except core (train) + lineage/QC channels
all_markers <- grep("_asinh$", names(b.dat), value = TRUE)
all_markers <- setdiff(all_markers, c(
  "FJComp-BV480-A_Live_Dead_asinh",
  "FJComp-BUV395-A_CD3_asinh",
  "FJComp-Alexa Fluor 647-A_GP350_asinh"
))
overlay_markers <- setdiff(all_markers, train_markers)

# core-marker heatmap (long)
cmatl <- as.data.table(as.table(cmat))
setnames(cmatl, c("V1", "V2", "N"), c("cluster", "marker", "z"))
cmatl[, cluster := factor(cluster, levels = row_ord)]
cmatl[, marker := factor(marker, levels = col_ord)]
cmatl[, marker_pretty := factor(pretty_marker(marker),
                                levels = pretty_marker(levels(marker)))]

# overlay-marker heatmap (long)
cmo <- b.m[, lapply(.SD, median), by = FlowSOM_metacluster, .SDcols = overlay_markers]
cmol <- melt(cmo, id.vars = "FlowSOM_metacluster", variable.name = "marker", value.name = "median_asinh")
cmol[, z := scale_mad(median_asinh), by = marker]
cmol[, cluster := factor(FlowSOM_metacluster, levels = row_ord)]
cmol[, marker_pretty := factor(pretty_marker(marker),
                               levels = pretty_marker(levels(marker)))]

p_core <- ggplot(cmatl, aes(x = marker_pretty, y = cluster, fill = z)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D", limits = c(-1.6, 1.6), oob = scales::squish) +
  labs(x = NULL, y = "Cluster", fill = "z", title = "Core markers") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(out_dir, "core.pdf"), p_core, width = 55.88, height = 101.6, units = "mm")

p_overlay <- ggplot(cmol, aes(x = marker_pretty, y = cluster, fill = z)) +
  geom_tile() +
  scale_fill_viridis_c(option = "A", limits = c(-4, 4), oob = scales::squish) +
  labs(x = NULL, y = "Cluster", fill = "z (median)", title = "Overlay markers") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(out_dir, "overlay.pdf"), p_overlay, width = 55.88, height = 101.6, units = "mm")

ca.df <- data.frame(Row = 1:length(ca), Cluster = factor(names(ca), levels = row_ord), Label = ca)
p_ca <- ggplot(ca.df, aes(x = 1, y = Cluster, fill = Cluster)) +
  geom_tile() +
  geom_label(aes(label = Label), fill = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
ggsave(file.path(out_dir, "clusters.pdf"), p_ca, width = 88.9, height = 101.6, units = "mm")

p_combined_B <- (p_core | p_ca | p_overlay) + plot_layout(widths = c(1, 2.5, 1))
ggsave(file.path(out_dir, "combined_panel_B.pdf"), p_combined_B, width = 203.2, height = 101.6, units = "mm")


# ---------------------------------------------------------------------------- #
# d - ABC-proportion (cl13 + cl15) log2-FC waterfall
# Within-patient change in ABC / ABC-like (metaclusters 13 + 15) frequency
# between the paired pre-relapse and remission samples, as log2(fold change),
# ordered ascending across patients.
b.m$cl13cl15 <- ifelse(as.character(b.m$FlowSOM_metacluster) %in% c("13", "15"), TRUE, FALSE)
per_sample_cl13cl15 <- b.m[, .(
  n_B = .N,
  n_cl13cl15 = sum(cl13cl15, na.rm = TRUE),
  pct_cl13cl15 = 100 * mean(cl13cl15, na.rm = TRUE)
), by = .(Patient_ID, FileName, Group, exp_day, Timing)]

# --- paired pre-relapse vs remission, per patient -----------------------------
CAP_LOG2FC <- 3.2    # the published panel is annotated "capped at ±3.2"
PSEUDO     <- 0.01   # percentage points; keeps log2 finite when a patient has
                     # zero ABC events in one of the two samples

# Cast on Group which carries the 'pre_relapse' / 'remission' labels this panel
# needs.
d_wide <- data.table::dcast(
  per_sample_cl13cl15,
  Patient_ID ~ Group,
  value.var = "pct_cl13cl15",
  fun.aggregate = function(x) mean(x, na.rm = TRUE)
)
message("Panel d: Group levels found: ",
        paste(setdiff(names(d_wide), "Patient_ID"), collapse = ", "))

# Identify the two group columns without hard-coding their spelling.
.grp  <- setdiff(names(d_wide), "Patient_ID")
.pre  <- grep("pre",  .grp, ignore.case = TRUE, value = TRUE)[1]
.rem  <- grep("^rem", .grp, ignore.case = TRUE, value = TRUE)[1]
if (is.na(.pre) || is.na(.rem))
  stop("Could not identify pre-relapse / remission Group columns: ",
       paste(.grp, collapse = ", "))

d_wide <- d_wide[is.finite(get(.pre)) & is.finite(get(.rem))]
d_wide[, log2fc := log2((get(.pre) + PSEUDO) / (get(.rem) + PSEUDO))]
d_wide[, log2fc_capped := pmax(pmin(log2fc, CAP_LOG2FC), -CAP_LOG2FC)]
d_wide[, direction := ifelse(log2fc_capped >= 0, "Increase", "Decrease")]
data.table::setorder(d_wide, log2fc_capped)
d_wide[, patient_rank := seq_len(.N)]

message("Panel d: ", nrow(d_wide), " paired patients; ",
        sum(d_wide$direction == "Increase"), " increase, ",
        sum(d_wide$direction == "Decrease"), " decrease")
print(stats::wilcox.test(d_wide[[.pre]], d_wide[[.rem]], paired = TRUE))

p_abc_waterfall <- ggplot(d_wide, aes(x = patient_rank, y = log2fc_capped,
                                      color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = c(-CAP_LOG2FC, CAP_LOG2FC),
             linetype = "dotted", linewidth = 0.4) +
  geom_segment(aes(xend = patient_rank, y = 0, yend = log2fc_capped),
               linetype = "dotted", linewidth = 0.5) +
  geom_point(aes(size = abs(log2fc_capped))) +
  scale_color_manual(values = c(Decrease = "#D7191C", Increase = "#2C7BB6"),
                     name = NULL) +
  scale_size_continuous(range = c(1, 3.5), guide = "none") +
  scale_x_continuous(breaks = NULL) +
  labs(title = "ABC/ABC-like (cl13+cl15)",
       x = "Patient",
       y = sprintf('log2-FC in ABC proportion\n(pre-relapse / remission) (capped at ±%.1f)',
                   CAP_LOG2FC)) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position    = "top",
        legend.key.size    = unit(0.3, "cm"),
        plot.title         = element_text(size = 8))

ggsave(file.path(out_dir, "panelD_ABC_log2fc_waterfall.pdf"), p_abc_waterfall,
       width = 50.8, height = 101.6, units = "mm")
data.table::fwrite(d_wide, file.path(out_dir, "panelD_ABC_log2fc_table.csv"))


# ---------------------------------------------------------------------------- #
# e - EBV gp350+ cells on UMAP
cpg <- min(table(b.m$Group))            # down-sample groups to the smaller n
set.seed(42)                            # which cells are drawn decides panel e
b.ds <- b.m[, .SD[sample(.N, cpg)], by = Group]
b.ds <- b.ds[order(`FJComp-Alexa Fluor 647-A_GP350_asinh`)]

p_ebv <- ggplot(b.ds, aes(x = UMAP_X, y = UMAP_Y, color = gp350_pos, size = gp350_pos)) +
  ggrastr::geom_point_rast(size = 0.1, raster.dpi = 600) +
  scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 1)) +
  scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "gray80")) +
  coord_fixed() +
  theme_few() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "EBV_umap.pdf"), p_ebv, width = 101.6, height = 101.6, units = "mm")


# ---------------------------------------------------------------------------- #
# f - gp350+ fraction per cluster
gp350_by_cluster <- b.m[, .(
  n_cells   = .N,
  n_gp350   = sum(gp350_pos, na.rm = TRUE),
  pct_gp350 = 100 * mean(gp350_pos, na.rm = TRUE)
), by = FlowSOM_metacluster][order(FlowSOM_metacluster)]
gp350_by_cluster$FlowSOM_metacluster <- factor(gp350_by_cluster$FlowSOM_metacluster, levels = row_ord)

p_gpc <- ggplot(gp350_by_cluster, aes(x = "gp350", y = FlowSOM_metacluster, fill = pct_gp350)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "inferno") +
  labs(x = NULL, y = "FlowSOM metacluster", fill = "% gp350+") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
ggsave(file.path(out_dir, "gp_by_cluster.pdf"), p_gpc, width = 50.8, height = 101.6, units = "mm")

p_panelF <- (p_gpc | p_ca) + plot_layout(widths = c(1, 4))
ggsave(file.path(out_dir, "panelF.pdf"), p_panelF, width = 127, height = 101.6, units = "mm")


# ---------------------------------------------------------------------------- #
# g - gp350+ by disease state
# Model-predicted % gp350+ (marginal over patients) + paired contrast p-value
emm <- emmeans(fit, ~ Group, type = "response")
p_df <- as.data.frame(pairs(emm))
p_txt <- sprintf("two-sided p = %.2g", p_df$p.value[1])
emm_df <- as.data.frame(emm) %>%
  mutate(pct = 100 * prob, pct_lo = 100 * asymp.LCL, pct_hi = 100 * asymp.UCL)

p_model <- ggplot(emm_df, aes(x = Group, y = pct)) +
  geom_errorbar(aes(ymin = pct_lo, ymax = pct_hi), width = 0.15) +
  geom_point(size = 4, aes(color = Group)) +
  annotate("text", x = 1.5, y = max(emm_df$pct_hi) * 1.08, label = p_txt, size = 3.4) +
  scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, 0.45)) +
  labs(y = "Model-predicted % gp350+", x = NULL) +
  theme_classic2() +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "model-pred perc gp350.pdf"), p_model, width = 45.72, height = 101.6, units = "mm")

# Observed per-patient paired %gp350+ (pre_relapse vs remission)
file_col    <- "FileName"
patient_col <- "Patient_ID"
state_col   <- "Group"
freq_by_sample <- b.m[, .(
  n_bcells   = .N,
  n_gp350pos = sum(gp350_pos, na.rm = TRUE),
  pct_gp350  = 100 * mean(gp350_pos, na.rm = TRUE)
), by = c(file_col, patient_col, state_col)]
freq_by_sample <- freq_by_sample[order(get(patient_col), get(state_col))]
paired <- dcast(
  freq_by_sample,
  formula = as.formula(paste(patient_col, "~", state_col)),
  value.var = "pct_gp350"
)
plot_dt <- melt(
  paired,
  id.vars = patient_col,
  measure.vars = c("remission", "pre_relapse"),
  variable.name = "Group",
  value.name = "pct_gp350"
)

p_gpp <- ggplot(plot_dt, aes(x = Group, y = pct_gp350, group = .data[[patient_col]])) +
  geom_line(alpha = 0.35) +
  geom_point(size = 2, aes(color = Group)) +
  scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, 1.1)) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "% gp350+ (of CD19+ B cells)",
       title = "gp350+ B cells per patient (paired pre_relapse vs remission)") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "gp_by_patient.pdf"), p_gpp, width = 50.8, height = 101.6, units = "mm")

p_panel_FG <- (p_gpp | p_model) +
  plot_layout(widths = c(2, 2), guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "panel_FG.pdf"), p_panel_FG, width = 101.6, height = 101.6, units = "mm")
