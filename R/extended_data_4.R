
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : healthy_split_null.RDS     (scrna_pbmc_scdist.R, "Controls" section)
#          remission_split_null.RDS   (scrna_pbmc_scdist.R, "Controls" section)
#          scdist_pbmc_prerelapse_v_remission.RDS          (scrna_pbmc_scdist.R)
#
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

ZENODO_DIR <- Sys.getenv("ZENODO_DIR", "zenodo")
SCDIST_DIR <- file.path(ZENODO_DIR, "scdist")

source("R/scdist_functions.R")   # DistPlot2

OBJ_DIR <- "objects"
out_dir <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_SHOW <- 10          # null iterations displayed in panels a and b
YLIM   <- c(0, 12.5)  # shared x-range across every null subplot

# --- Inputs -------------------------------------------------------------------
healthy_split_null   <- readRDS(file.path(OBJ_DIR, "healthy_split_null.RDS"))
remission_split_null <- readRDS(file.path(OBJ_DIR, "remission_split_null.RDS"))
scdist_pbmc_pre              <- readRDS(file.path(SCDIST_DIR, "scdist_pbmc_prerelapse_v_remission.RDS"))


# --- Helper -------------------------------------------------------------------
# Flatten a split-null object into one row per (iteration x cell type).
extract_null_dist_table <- function(null_obj) {
  purrr::map_dfr(null_obj$results, function(x) {
    if (inherits(x$scdist, "error")) return(NULL)
    res <- x$scdist$results
    tibble::tibble(
      run_id    = x$run_id,
      null_type = x$null_type,
      cluster   = rownames(res),
      Dist      = res[["Dist."]],
      CI_low    = res[["95% CI (low)"]],
      CI_high   = res[["95% CI (upper)"]],
      p_val     = if ("p.val" %in% colnames(res)) res[["p.val"]] else NA_real_
    )
  })
}

# One DistPlot2 subplot per null iteration, on a common x-range.
null_iteration_grid <- function(null_obj, n_show = N_SHOW) {
  plots <- lapply(seq_len(n_show), function(i) {
    DistPlot2(null_obj$results[[i]]$scdist) +
      ylim(YLIM) +
      NoLegend() +
      theme(axis.title = element_blank())
  })
  wrap_plots(plots, ncol = 5)
}


# ---------------------------------------------------------------------------- #
# a - Healthy vs healthy split null (first 10 of 100 iterations)
p_a <- null_iteration_grid(healthy_split_null)
ggsave(file.path(out_dir, "extended_data_4a.pdf"), p_a,
       width = 500, height = 200, units = "mm")


# ---------------------------------------------------------------------------- #
# b - MS remission vs remission split null (first 10 of 100 iterations)
p_b <- null_iteration_grid(remission_split_null)
ggsave(file.path(out_dir, "extended_data_4b.pdf"), p_b,
       width = 500, height = 200, units = "mm")


# ---------------------------------------------------------------------------- #
# c - Observed maximum pre-relapse scDist vs the healthy split null
healthy_null_df <- extract_null_dist_table(healthy_split_null)

max_null_healthy <- healthy_null_df %>%
  group_by(run_id) %>%
  summarize(max_Dist = max(Dist, na.rm = TRUE), .groups = "drop")

obs_pre_max <- max(scdist_pbmc_pre$results[["Dist."]], na.rm = TRUE)

message("Observed max pre-relapse scDist: ", round(obs_pre_max, 3),
        "  |  null max range: ",
        round(min(max_null_healthy$max_Dist), 3), "-",
        round(max(max_null_healthy$max_Dist), 3),
        "  (", nrow(max_null_healthy), " iterations)")

p_c <- ggplot(max_null_healthy, aes(x = max_Dist)) +
  geom_histogram(bins = 30, fill = "gray80", color = "white") +
  geom_vline(xintercept = obs_pre_max, color = "firebrick", linewidth = 1) +
  theme_classic() +
  labs(x = "Maximum scDist distance across PBMC cell types",
       y = "Null iterations",
       title = "Observed maximum pre-relapse scDist vs healthy split null")

ggsave(file.path(out_dir, "extended_data_4c.pdf"), p_c,
       width = 150, height = 50, units = "mm")
