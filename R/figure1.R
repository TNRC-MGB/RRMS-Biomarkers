
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-22-2026
#
# Figure 1
#
# Input: PBMC Samples.xlsx
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Colorblind-safe palette (Okabe-Ito derived)
cond_colors <- c(
  "PreRelapse" = "#E69F00",
  "Relapse"    = "#D55E00",
  "Remission"  = "#0072B2",
  "Healthy"    = "#009E73"
)
cond_levels <- c("PreRelapse", "Relapse", "Remission", "Healthy")

assay_colors <- c(
  "scRNA"  = "#7570B3",
  "Bulk"   = "#D95F02",
  "Flow"   = "#1B9E77",
  "qPCR"   = "#E7298A"
)

sex_colors <- c("F" = "#CC79A7", "M" = "#56B4E9")

theme_fig <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill = "gray95", color = NA),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle      = element_text(size = 8, color = "gray40"),
    axis.title         = element_text(size = 9),
    axis.text          = element_text(size = 8)
  )


# ---------------------------------------------------------------------------- #
# Load data
pbmc <- read_excel(file.path("Supplementary Information/S3 PBMC Samples.xlsx"))

pbmc <- pbmc %>%
  mutate(
    Condition = factor(Condition, levels = cond_levels),
    has_SCRNA = str_detect(Assays, "SCRNA"),
    has_BULK  = str_detect(Assays, "BULK"),
    has_FLOW  = str_detect(Assays, "FLOW"),
    has_QPCR  = str_detect(Assays, "QPCR")
  )

patients    <- pbmc %>% distinct(Donor_ID, .keep_all = TRUE)
ms_patients <- patients %>% filter(Condition != "Healthy")

n_ms    <- ms_patients %>% nrow()
n_ms_s  <- pbmc %>% filter(Condition != "Healthy") %>% nrow()
n_hc    <- patients %>% filter(Condition == "Healthy") %>% nrow()
n_hc_s  <- pbmc %>% filter(Condition == "Healthy") %>% nrow()
n_scrna <- pbmc %>% filter(has_SCRNA) %>% nrow() + 10 # 10 patient-timepoints had a second vial
n_bulk  <- pbmc %>% filter(has_BULK) %>% nrow()
n_flow  <- pbmc %>% filter(has_FLOW) %>% nrow()
n_qpcr  <- pbmc %>% filter(has_QPCR) %>% nrow()


# ---------------------------------------------------------------------------- #
# PANEL A: Study design overview
study_text <- tibble(
  x = 0, y = rev(seq_len(13)),
  label = c(
    "bold('BWH CLIMB Study')",
    sprintf("paste('MS patients: ', bold('%d'), ' (', bold('%d'), ' PBMC samples)')", n_ms, n_ms_s),
    sprintf("paste('Healthy donors: ', bold('%d'), ' (', bold('%d'), ' PBMC samples)')", n_hc, n_hc_s),
    "",
    "bold('Assays')",
    sprintf("paste('scRNA-seq: ', bold('%d'), ' samples (+ 1 WGS per subject)')", n_scrna),
    sprintf("paste('CD19+ Bulk RNA-seq: ', bold('%d'), ' samples')", n_bulk),
    sprintf("paste('Flow cytometry: ', bold('%d'), ' samples')", n_flow),
    sprintf("paste('qPCR: ', bold('%d'), ' samples')", n_qpcr),
    "",
    "bold('Conditions')",
    "'PreRelapse \u00b7 Relapse \u00b7 Remission \u00b7 Healthy'",
    ""
  )
)

p_a <- ggplot(study_text, aes(x = x, y = y)) +
  geom_text(aes(label = label), parse = TRUE, hjust = 0, size = 3) +
  xlim(-0.1, 5) + ylim(0, 14) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 10)) +
  labs(title = "A")


# ---------------------------------------------------------------------------- #
# PANEL B: Sample distribution by Condition (bar chart)
cond_counts <- pbmc %>%
  dplyr::count(Condition) %>%
  mutate(label = paste0("n=", n))

p_b <- ggplot(cond_counts, aes(x = Condition, y = n, fill = Condition)) +
  geom_col(color = "white", linewidth = 0.4, width = 0.7) +
  geom_text(aes(label = label), vjust = -0.4, size = 3) +
  scale_fill_manual(values = cond_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "B  Sample Distribution", x = NULL, y = "Samples") +
  theme_fig +
  theme(legend.position = "none", panel.grid.major.x = element_blank())


# ---------------------------------------------------------------------------- #
# PANEL C: Age violin + sex bar (side by side)
# C left: Age by condition
p_c_age <- pbmc %>%
  ggplot(aes(x = Condition, y = `Age (yrs)`, fill = Condition)) +
  geom_violin(alpha = 0.35, color = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.8, alpha = 0.9,
               fill = "white", color = "gray30") +
  scale_fill_manual(values = cond_colors) +
  labs(title = "C  Age by Condition", x = NULL, y = "Age (years)") +
  theme_fig +
  theme(legend.position = "none")

# C middle: Sex ratio - patient level, MS vs Healthy
sex_summary <- patients %>%
  mutate(Group = ifelse(Condition == "Healthy", "Healthy", "MS")) %>%
  dplyr::count(Group, Sex) %>%
  group_by(Group) %>%
  mutate(pct = n / sum(n) * 100,
         label = paste0(n, " (", round(pct), "%)"))

p_c_sex <- ggplot(sex_summary, aes(x = Group, y = pct, fill = Sex)) +
  geom_col(position = "stack", color = "white", linewidth = 0.4, width = 0.6) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 2.8, lineheight = 0.9) +
  scale_fill_manual(values = sex_colors) +
  labs(title = "    Sex (patient-level)", x = NULL, y = "Percentage (%)") +
  theme_fig +
  theme(legend.position = "right",
        legend.key.size = unit(0.35, "cm"))

# C right: Treatment pie chart - MS patients only
tx_summary <- ms_patients %>%
  dplyr::count(Treatment) %>%
  arrange(desc(n)) %>%
  mutate(
    pct = n / sum(n) * 100,
    label = ifelse(pct >= 5, paste0(Treatment, "\n", round(pct), "%"), "")
  )

# Muted qualitative palette for treatments
tx_pal <- c("#4575B4", "#91BFDB", "#FEE090", "#FC8D59", "#D73027",
            "#ABD9E9", "#E0F3F8", "#FDAE61", "#F46D43", "#A50026",
            "#74ADD1", "#313695", "#FFFFBF", "#66BD63")

p_c_tx <- ggplot(tx_summary, aes(x = "", y = n, fill = Treatment)) +
  geom_col(width = 1, color = "white", linewidth = 0.4) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 2, lineheight = 0.85) +
  scale_fill_manual(values = tx_pal[seq_len(nrow(tx_summary))]) +
  labs(title = "    Treatment (MS)") +
  theme_void(base_size = 10) +
  theme(
    plot.title   = element_text(face = "bold", size = 11, hjust = 0),
    legend.position  = "right",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))

# Combine side by side: age | sex | treatment pie
p_c <- p_c_age + p_c_sex + p_c_tx + plot_layout(widths = c(1.3, 0.8, 1.2))


# ---------------------------------------------------------------------------- #
# PANEL D: Time-to-relapse swim lane (MS patients only)
swim <- pbmc %>%
  filter(Condition != "Healthy") %>%
  arrange(Donor_ID, `Time to relapse (d)`)

# Assign y-position by patient (ordered by median TTR)
patient_order <- swim %>%
  group_by(Donor_ID) %>%
  summarize(med_ttr = median(`Time to relapse (d)`, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_ttr) %>%
  mutate(y_pos = row_number())

swim <- swim %>%
  left_join(patient_order, by = "Donor_ID")

# Expand to one row per sample-assay for shape mapping
swim_long <- swim %>%
  select(Sample_ID, Donor_ID, `Time to relapse (d)`, Condition, y_pos,
         has_SCRNA, has_BULK, has_FLOW, has_QPCR) %>%
  pivot_longer(cols = starts_with("has_"),
               names_to = "Assay", values_to = "present",
               names_prefix = "has_") %>%
  filter(present) %>%
  mutate(Assay = recode(Assay, "SCRNA" = "scRNA", "BULK" = "Bulk",
                        "FLOW" = "Flow", "QPCR" = "qPCR"),
         Assay = factor(Assay, levels = c("scRNA", "Bulk", "Flow", "qPCR")))

p_d <- ggplot() +
  # Shaded ±90-day relapse window
  annotate("rect", xmin = -90, xmax = 90, ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.06) +
  # Patient timelines (connect samples from same patient)
  geom_line(data = swim %>% distinct(Sample_ID, .keep_all = TRUE),
            aes(x = `Time to relapse (d)`, y = y_pos, group = Donor_ID),
            color = "gray80", linewidth = 0.3) +
  # Threshold reference lines
  geom_vline(xintercept = c(-90, 90), linetype = "dotted",
             color = "gray50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#D55E00", linewidth = 0.5) +
  annotate("text", x = -90, y = max(patient_order$y_pos) + 3, label = "\u221290d",
           size = 2, color = "gray50") +
  annotate("text", x = 0, y = max(patient_order$y_pos) + 3, label = "Relapse",
           size = 2.2, color = "#D55E00", fontface = "bold") +
  annotate("text", x = 90, y = max(patient_order$y_pos) + 3, label = "+90d",
           size = 2, color = "gray50") +
  # Sample points: color = Condition, shape = Assay
  geom_point(data = swim_long,
             aes(x = `Time to relapse (d)`, y = y_pos,
                 color = Condition, shape = Assay),
             size = 1.4, alpha = 0.8,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = cond_colors) +
  scale_shape_manual(values = c("scRNA" = 16, "Bulk" = 17, "Flow" = 15, "qPCR" = 18)) +
  scale_x_continuous(limits = c(-500, 500), breaks = seq(-400, 400, 200)) +
  labs(title = "D  Temporal Sampling Design",
       subtitle = "MS samples ordered by median time to relapse; dashed line = relapse event",
       x = "Time to Nearest Relapse (days)", y = "Patient (ordered)") +
  theme_fig +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position  = "right",
    legend.key.size  = unit(0.35, "cm"),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7)
  )


# ---------------------------------------------------------------------------- #
# PANEL E: Condition-Assay heatmap
cond_assay <- pbmc %>%
  pivot_longer(cols = c(has_SCRNA, has_BULK, has_FLOW, has_QPCR),
               names_to = "Assay", values_to = "Present",
               names_prefix = "has_") %>%
  filter(Present) %>%
  dplyr::count(Condition, Assay) %>%
  mutate(
    Assay = recode(Assay, "SCRNA" = "scRNA", "BULK" = "Bulk",
                   "FLOW" = "Flow", "QPCR" = "qPCR"),
    Assay = factor(Assay, levels = c("scRNA", "Bulk", "Flow", "qPCR"))
  )

p_e <- ggplot(cond_assay, aes(x = Assay, y = Condition, fill = n)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = n), size = 3.5, fontface = "bold") +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", name = "N",
                      guide = guide_colorbar(barwidth = 0.5, barheight = 3)) +
  labs(title = "E  Assay Coverage", x = NULL, y = NULL) +
  theme_fig +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 6)) +
  coord_equal()


# ---------------------------------------------------------------------------- #
# Figure 1
# Layout: top row = A | B | C;  bottom row = D (wide) | E
top_row    <- p_a + p_b + p_c + plot_layout(widths = c(0.8, 0.8, 2.2))
bottom_row <- p_d + p_e + plot_layout(widths = c(2.5, 1))

n_total   <- nrow(pbmc)
n_donors  <- n_distinct(pbmc$Donor_ID)

fig1 <- top_row / bottom_row +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    title    = "Figure 1. Cohort Overview",
    subtitle = sprintf("BWH CLIMB Study: %d participants, %d PBMC samples, 4 assay platforms",
                        n_donors, n_total),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  )

# Save
ggsave(file.path(out_dir, "Fig1_Cohort_Overview.pdf"), fig1,
       width = 14, height = 10, dpi = 600, device = cairo_pdf)

ggsave(file.path(out_dir, "Fig1_Cohort_Overview.tiff"), fig1,
       width = 14, height = 10, dpi = 600, compression = "lzw")
