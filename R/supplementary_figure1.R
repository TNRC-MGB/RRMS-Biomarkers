
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-22-2026
#
# Supplementary Figure 1 - MS PBMC cohort overview
#   a  Time-to-relapse distribution
#   b  Assay coverage UpSet (SCRNA / FLOW / QPCR / BULK)
#   c  Samples per patient
#
# Inputs : PBMC Samples.xlsx
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

filter    <- dplyr::filter
mutate    <- dplyr::mutate
count     <- dplyr::count
arrange   <- dplyr::arrange
summarize <- dplyr::summarize
rename    <- dplyr::rename

data_dir <- "Supplementary Information"
out_dir  <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

condition_colors <- c(
  "PreRelapse" = "#E69F00",
  "Relapse"    = "#D55E00",
  "Remission"  = "#0072B2",
  "Healthy"    = "#009E73"
)
cond_levels <- c("PreRelapse", "Relapse", "Remission", "Healthy")

assay_colors <- c(
  "SCRNA" = "#7570B3",
  "BULK"  = "#D95F02",
  "FLOW"  = "#1B9E77",
  "QPCR"  = "#E7298A"
)

theme_publication <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "gray95", color = NA),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, color = "gray40"),
    axis.title         = element_text(size = 10)
  )

# --- Inputs -------------------------------------------------------------------
pbmc <- read_excel(file.path(data_dir, "S3 PBMC Samples.xlsx"))

pbmc <- pbmc %>%
  mutate(
    Condition = factor(Condition, levels = cond_levels),
    has_SCRNA = str_detect(Assays, "SCRNA"),
    has_BULK  = str_detect(Assays, "BULK"),
    has_FLOW  = str_detect(Assays, "FLOW"),
    has_QPCR  = str_detect(Assays, "QPCR")
  )

patients <- pbmc %>% distinct(Donor_ID, .keep_all = TRUE)


# a - Time to Relapse Distribution
p4 <- pbmc %>%
  filter(Condition != "Healthy", !is.na(`Time to relapse (d)`)) %>%
  ggplot(aes(x = `Time to relapse (d)`, fill = Condition)) +
  geom_histogram(bins = 40, alpha = 0.7, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#D55E00", linewidth = 0.5) +
  annotate("text", x = 0, y = Inf, label = "Relapse event", vjust = 2,
           size = 3, color = "#D55E00", fontface = "bold") +
  scale_fill_manual(values = condition_colors) +
  labs(title = "Time to Relapse Distribution",
       subtitle = "Negative = before relapse, positive = after relapse",
       x = "Time to Relapse (days)", y = "Count") +
  theme_publication

ggsave(file.path(out_dir, "Fig4_TTR_Distribution.pdf"), p4,
       width = 9, height = 5, dpi = 300)


# b - Assay coverage UpSet (SCRNA / FLOW / QPCR / BULK)
upset_data <- pbmc %>%
  mutate(
    SCRNA = as.integer(has_SCRNA),
    BULK  = as.integer(has_BULK),
    FLOW  = as.integer(has_FLOW),
    QPCR  = as.integer(has_QPCR)
  ) %>%
  select(SCRNA, BULK, FLOW, QPCR) %>%
  as.data.frame()

pdf(file.path(out_dir, "Fig2A_Assay_UpSet.pdf"), width = 8, height = 5)
print(upset(upset_data,
            sets = c("SCRNA", "BULK", "FLOW", "QPCR"),
            sets.bar.color = unname(assay_colors[c("SCRNA", "BULK", "FLOW", "QPCR")]),
            order.by = "freq",
            main.bar.color = "#4575B4",
            text.scale = c(1.3, 1.2, 1.1, 1.1, 1.3, 1.1),
            mb.ratio = c(0.6, 0.4),
            point.size = 3.5,
            line.size = 1.2))
dev.off()


# c - Samples per Patient
samples_per_patient <- pbmc %>%
  count(Donor_ID, name = "N_Samples") %>%
  left_join(patients %>% select(Donor_ID, Condition), by = "Donor_ID") %>%
  mutate(Group = ifelse(Condition == "Healthy", "Healthy", "MS"))

p3 <- samples_per_patient %>%
  ggplot(aes(x = factor(N_Samples), fill = Group)) +
  geom_bar(position = "stack", color = "white", linewidth = 0.3) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("MS" = "#4575B4", "Healthy" = "#009E73")) +
  labs(title = "Samples per Patient",
       subtitle = sprintf("%d patients with a total of %d samples",
                          nrow(patients), nrow(pbmc)),
       x = "Number of Samples", y = "Number of Patients") +
  theme_publication

ggsave(file.path(out_dir, "Fig3_Samples_per_Patient.pdf"), p3,
       width = 7, height = 5, dpi = 300)
