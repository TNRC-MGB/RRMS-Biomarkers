
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
#
# Inputs : bulk/results/dge.rds
#          Supplementary Information/S7 B Cell Bulk RNA-seq Libraries.xlsx
#          bulk/reference/genesets/pebv500.RDS
# Outputs: bulk/results/figure5e_bigpoints.pdf
#          bulk/results/figure5e_scores.tsv
#          bulk/results/figure5e_qc_metrics.tsv
#
# ==============================================================================


setwd("C:/Users/devin/Desktop/rrms")   # <- repository root

source("R/packages.R")

OUT_DIR     <- "bulk/results"
S7_XLSX     <- "Supplementary Information/S7 B Cell Bulk RNA-seq Libraries.xlsx"
P_MODULES   <- "bulk/reference/genesets/pebv500.RDS"
MODULE_NAME <- "LMP-1"

N_TOP  <- 1000L   # most variable genes used for the PCA
N_PC   <- 10L
N_NULL <- 500L    # expression-matched random gene sets
PERI   <- 21      # peri-relapse half-width (d), three weeks either side
SEED   <- 42L
STOP_ON_SEX_MISMATCH <- TRUE
set.seed(SEED)

fs <- 6   # axis text size (pt)
fl <- 6   # axis title size (pt)
OK   <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")   # Okabe-Ito
SET2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854")   #

PT_SIZE <- 3.2    
FIG_W   <- 75     
FIG_H   <- 88     



# --- Inputs -------------------------------------------------------------------
dge <- readRDS(file.path(OUT_DIR, "dge.rds"))
s7  <- as.data.table(read_excel(S7_XLSX))
setnames(s7, "Time to relapse (d)", "time_to_relapse_d")
s7[, time_to_relapse_d := as.numeric(time_to_relapse_d)]

message("quantified samples: ", ncol(dge), "   S7 rows: ", nrow(s7))
no_s7 <- setdiff(colnames(dge), s7$Sample_ID)
no_q  <- setdiff(s7$Sample_ID, colnames(dge))

meta <- s7[match(colnames(dge), Sample_ID)]
for (v in c("Donor_ID", "Condition", "time_to_relapse_d", "Sex")) {
  old <- as.character(dge$samples[[v]]); new <- as.character(meta[[v]])
  n_diff <- sum(!is.na(old) & !is.na(new) & old != new)
  if (n_diff) message("  NOTE: ", v, " differs from the imported metadata at ",
                  n_diff, " sample(s); S7 is used.")
  dge$samples[[v]] <- meta[[v]]
}
dge$samples$group <- factor(dge$samples$Condition)

sm   <- as.data.table(dge$samples)
lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
message("samples carried forward: ", ncol(dge), "   genes: ", nrow(dge),
    "   donors: ", uniqueN(sm$Donor_ID))


# per-sample metrics
qc <- data.table(
  Sample_ID  = colnames(dge),
  Donor_ID   = as.character(sm$Donor_ID),
  Condition  = as.character(sm$Condition),
  Sex_table  = as.character(sm$Sex),
  time_to_relapse_d = as.numeric(sm$time_to_relapse_d),
  lib_size   = colSums(dge$counts),
  norm_factor = dge$samples$norm.factors,
  n_detected = colSums(dge$counts > 0),
  frac_top50 = apply(dge$counts, 2, function(v) sum(sort(v, decreasing = TRUE)[1:50]) / sum(v)))


message("---- per-sample metrics ------------------------------------------------")
message(paste(capture.output(print(
  qc[, .(metric = c("lib_size", "n_detected", "norm_factor", "frac_top50"),
         min = c(min(lib_size), min(n_detected), round(min(norm_factor), 3),
                 round(min(frac_top50), 3)),
         median = c(median(lib_size), median(n_detected),
                    round(median(norm_factor), 3), round(median(frac_top50), 3)),
         max = c(max(lib_size), max(n_detected), round(max(norm_factor), 3),
                 round(max(frac_top50), 3)))])), collapse = '\n'))


# PCA
row_var <- function(M) { m <- rowMeans(M); rowSums((M - m)^2) / (ncol(M) - 1) }
vg  <- head(order(row_var(lcpm), decreasing = TRUE), N_TOP)
pc  <- prcomp(t(lcpm[vg, ]), center = TRUE, scale. = FALSE)
ve  <- 100 * pc$sdev^2 / sum(pc$sdev^2)
npc <- min(N_PC, ncol(pc$x))

message("---- PCA on the ", N_TOP, " most variable genes ------------------------")
message("variance explained: ", paste0("PC", seq_len(min(5, npc)), " ",
                                   round(ve[seq_len(min(5, npc))], 1), "%", collapse = "   "))
pcd <- data.table(qc, pc$x[, seq_len(npc), drop = FALSE])


message("---- PCA after removing the donor effect -------------------------------")
rep_don <- qc[, .N, by = Donor_ID][N > 1]$Donor_ID
kd  <- qc$Donor_ID %in% rep_don
pc2 <- NULL
if (sum(kd) >= 6 && uniqueN(qc$Condition[kd]) > 1) {
  adj <- limma::removeBatchEffect(lcpm[vg, kd], batch = factor(qc$Donor_ID[kd]),
                                  design = model.matrix(~ factor(qc$Condition[kd])))
  pc2 <- prcomp(t(adj), center = TRUE, scale. = FALSE)
  ve2 <- 100 * pc2$sdev^2 / sum(pc2$sdev^2)
} else message("too few repeat donors for a donor-adjusted PCA.")


# qc panels
th_qc <- theme_bw(base_size = fs) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")
lab    <- function(k, v) paste0("PC", k, " (", round(v[k], 1), "%)")
p_base <- ggplot(pcd, aes(PC1, PC2)) + labs(x = lab(1, ve), y = lab(2, ve)) + th_qc

q_a <- p_base + geom_point(aes(colour = Condition), size = 1.4) +
  scale_colour_manual(values = OK)

q_b <- p_base + geom_point(aes(colour = Sex_table), size = 1.4) +
  scale_colour_manual(values = OK)

q_c <- p_base +
  geom_line(aes(group = Donor_ID), colour = "grey65", linewidth = .3) +
  geom_point(aes(colour = Condition), size = 1.2) +
  scale_colour_manual(values = OK)

q_d <- ggplot(pcd, aes(PC1, log10(lib_size))) +
  geom_point(aes(colour = Condition), size = 1.2) +
  scale_colour_manual(values = OK) +
  labs(x = lab(1, ve), y = "log10 library size") + th_qc

# the gene set and the score
mods <- readRDS(P_MODULES)
if (!MODULE_NAME %in% names(mods))
  stop(MODULE_NAME, " not in ", P_MODULES)
sig   <- unique(mods[[MODULE_NAME]])
genes <- rownames(dge)[dge$genes$symbol %in% sig]

message("---- gene set ----------------------------------------------------------")
message('"', MODULE_NAME, '": ', length(sig), " symbols   matched to expressed genes: ",
    length(genes), " (", round(100 * length(genes) / length(sig)), "%)")

if (!requireNamespace("singscore", quietly = TRUE))
  stop("singscore is required (Bioconductor).")

score_singscore <- function(R, gs)
  setNames(singscore::simpleScore(R, upSet = gs)$TotalScore, colnames(R))

rk      <- singscore::rankGenes(lcpm)
primary <- score_singscore(rk, genes)

lraw  <- log2(t(t(dge$counts) / colSums(dge$counts)) * 1e6 + 1)
delta <- max(abs(primary - score_singscore(singscore::rankGenes(lraw), genes)[names(primary)]))


# expression-matched null
bin  <- cut(rowMeans(lcpm), quantile(rowMeans(lcpm), seq(0, 1, .05)),
            include.lowest = TRUE, labels = FALSE)
tabg <- table(bin[match(genes, rownames(lcpm))])
draw <- function() unlist(lapply(names(tabg), function(b)
  sample(rownames(lcpm)[bin == as.integer(b)], tabg[[b]], replace = FALSE)))

icc_of <- function(v) {
  a <- anova(lm(v ~ factor(qc$Donor_ID)))
  max(0, (a$`Mean Sq`[1] - a$`Mean Sq`[2]) /
         (a$`Mean Sq`[1] + (nrow(qc) / uniqueN(qc$Donor_ID) - 1) * a$`Mean Sq`[2]))
}
sd_obs  <- sd(primary); icc_obs <- icc_of(primary)
nullm   <- vapply(seq_len(N_NULL), function(i) score_singscore(rk, draw()), numeric(ncol(dge)))
sd_null <- apply(nullm, 2, sd)
icc_null <- apply(nullm, 2, icc_of)

message("---- expression-matched null (", N_NULL, " sets) -----------------------")
message("across-sample SD : observed ", round(sd_obs, 4),
    "   null median ", round(median(sd_null), 4),
    "   p = ", signif((1 + sum(sd_null >= sd_obs)) / (N_NULL + 1), 3))
message("donor ICC        : observed ", round(icc_obs, 4),
    "   null median ", round(median(icc_null), 4),
    "   p = ", signif((1 + sum(icc_null >= icc_obs, na.rm = TRUE)) / (N_NULL + 1), 3))


# the score table
d <- data.table(Sample_ID = names(primary), score = unname(primary),
                null_z = unname((primary - rowMeans(nullm)) / apply(nullm, 1, sd)))
d <- merge(d, qc[, .(Sample_ID, Donor_ID, Condition, time_to_relapse_d, Sex_table)],
           by = "Sample_ID")

WIN <- c("Remission (pre)", "pre-relapse", "peri-relapse", "post-relapse",
         "Remission (post)")
names(SET2) <- WIN
d[, win := factor(fifelse(
  Condition == "Remission", fifelse(time_to_relapse_d < 0, WIN[1], WIN[5]),
  fifelse(time_to_relapse_d < -PERI, WIN[2],
  fifelse(abs(time_to_relapse_d) <= PERI, WIN[3], WIN[4]))), levels = WIN)]


message("---- colour groups (peri band +/-", PERI, " d) --------------------------")
message(paste(capture.output(print(
  d[, .(n = .N, t_min = min(time_to_relapse_d), t_max = max(time_to_relapse_d)),
    by = win][order(t_min)])), collapse = '\n'))
fwrite(d, file.path(OUT_DIR, "figure5e_scores.tsv"), sep = '\t')



# LMP1 signature score against time to relapse
XLIM <- c(-400, 415)
YBREAKS <- pretty(d$score)
YLIM    <- range(YBREAKS)
n_out   <- sum(d$score < YLIM[1] | d$score > YLIM[2])

message("score range ", sprintf("%.4f .. %.4f", min(d$score), max(d$score)),
    "   axis ", sprintf("%.3f .. %.3f", YLIM[1], YLIM[2]),
    "   points outside: ", n_out)
if (n_out) stop(n_out, " point(s) fall outside the axis and would be dropped.")

th_e <- theme_classic(base_size = fs) +
  theme(plot.title   = element_text(size = fl, hjust = .5),
        axis.title.y = element_blank(),
        axis.line    = element_line(linewidth = .4, colour = "black"),
        axis.ticks   = element_line(linewidth = .4, colour = "black"),
        axis.ticks.length = unit(2, "pt"),
        axis.text    = element_text(colour = "black", size = fs),
        legend.position = "none",
        plot.margin  = margin(4, 22, 2, 4))

p_e <- ggplot(d, aes(time_to_relapse_d, score)) +
  geom_point(aes(colour = win), size = PT_SIZE, stroke = 0) +
  scale_colour_manual(values = SET2, drop = FALSE) +
  scale_x_continuous(limits = XLIM, breaks = c(-200, 0, 200, 400),
                     expand = expansion(0)) +
  scale_y_continuous(limits = YLIM, breaks = YBREAKS, expand = expansion(0)) +
  annotate("text", x = XLIM[2] + 18, y = YLIM[1], hjust = 0, vjust = .5,
           label = 'Time\n(days)', size = fs / .pt, lineheight = .95) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL,
       title = "CD19+ bulk RNA-seq: EBV LMP-1 host signature") + th_e

bar <- data.table(x1 = c(-163, -56, 54), x2 = c(-56, 54, 161),
                  lab = c("pre-", "peri-", "post-"), col = SET2[c(2, 3, 4)])

p_strip <- ggplot() +
  annotate("point", x = 0, y = 0.86, shape = 17, size = 4.2, colour = "black") +
  annotate("text", x = XLIM[1] + 8, y = 0.62, label = "Remission",
           hjust = 0, size = fs / .pt) +
  annotate("text", x = XLIM[2] - 8, y = 0.62, label = "Remission",
           hjust = 1, size = fs / .pt) +
  geom_segment(data = bar, aes(x = x1, xend = x2, y = 0.40, yend = 0.40,
                               colour = I(col)), linewidth = 1.6) +
  annotate("text", x = 0, y = 0.20, label = "Relapse window", size = fs / .pt) +
  geom_text(data = bar, aes(x = (x1 + x2) / 2, y = 0.04, label = lab),
            size = fs / .pt) +
  scale_x_continuous(limits = XLIM, expand = expansion(0)) +
  scale_y_continuous(limits = c(-0.05, 1), expand = expansion(0)) +
  theme_void() +
  theme(plot.margin = margin(0, 22, 2, 4), legend.position = "none")

p_fig5e <- wrap_plots(p_e, p_strip, ncol = 1, heights = c(1, 0.30))
ggsave(file.path(OUT_DIR, "figure5e_bigpoints.pdf"), p_fig5e,
       width = FIG_W, height = FIG_H, units = "mm")


message("wrote figure5e_bigpoints.pdf, figure5e_scores.tsv,")
message("      figure5e_qc_metrics.tsv to ", OUT_DIR)


# ---------------------------------------------------------------------------- #
# LMP-1 score, pre-relapse (0-90 d before onset) vs remission
# Donor-random-intercept LMM

stopifnot(all(c("PreRelapse", "Remission") %in% d$Condition))

dd <- d[Condition == "Remission" |
          (Condition == "PreRelapse" & time_to_relapse_d >= -90 &
             time_to_relapse_d < 0)]
dd[, group := factor(fifelse(Condition == "Remission", "Remission", "PreRelapse"),
                     levels = c("Remission", "PreRelapse"))]

fit_lmp1 <- lmerTest::lmer(score ~ group + (1 | Donor_ID), data = dd)

co <- summary(fit_lmp1)$coefficients["groupPreRelapse", ]
ci <- confint(fit_lmp1, "groupPreRelapse")

message("n = ", nrow(dd), " samples (",
        sum(dd$group == "PreRelapse"), " pre-relapse, ",
        sum(dd$group == "Remission"), " remission) from ",
        dd[, uniqueN(Donor_ID)], " donors")
message("PreRelapse - Remission: ", signif(co[["Estimate"]], 3),
        " (95% CI ", signif(ci[1], 3), " to ", signif(ci[2], 3),
        "), t = ", signif(co[["t value"]], 3),
        ", df = ", signif(co[["df"]], 4),
        ", p = ", signif(co[["Pr(>|t|)"]], 3))

fwrite(data.table(contrast = "PreRelapse_0to90d_vs_Remission",
                  n_samples = nrow(dd), n_donors = dd[, uniqueN(Donor_ID)],
                  estimate = co[["Estimate"]], ci_low = ci[1], ci_high = ci[2],
                  t = co[["t value"]], df = co[["df"]], p = co[["Pr(>|t|)"]]),
       file.path(OUT_DIR, "figure5e_lmp1_lmm.tsv"), sep = "\t")
