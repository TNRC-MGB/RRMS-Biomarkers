
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-28-2026
#
# CD19+ B-cell flow cytometry: Spectre/FlowSOM pipeline
#
#
# note: source('R/packages.R') first. Spectre::package.load() is called after 
# R/packages.R, so Spectre's versions of any shared names take effect last.
#
# ==============================================================================

FLOW_DIR <- Sys.getenv("FLOW_DIR", "zenodo/flow_fcs")
out_dir  <- "Intermediate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

Spectre::package.check()
Spectre::package.load()

InputDirectory <- FLOW_DIR
MetaDirectory  <- "data/flow"
metadata_file  <- "flow_sample_details.csv"
s5_file        <- "Supplementary Information/S5 Flow Cytometry Samples.xlsx"
asinh_cofactor <- 150L


# ---------------------------------------------------------------------------- #
# functions
# robust (median/MAD) scaling used before FlowSOM
fit_scale_mad <- function(x) {
  med <- median(x, na.rm = TRUE)
  sc  <- mad(x, center = med, constant = 1, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- 1
  list(center = med, scale = sc)
}
apply_scale <- function(x, fit) (x - fit$center) / fit$scale

scale_mad <- function(x) {
  x <- as.numeric(x)
  m <- stats::median(x, na.rm = TRUE)
  s <- stats::mad(x, constant = 1, na.rm = TRUE)  # constant=1 => raw MAD
  if (!is.finite(s) || s == 0) return(x - m)      # fallback
  (x - m) / s
}

pretty_marker <- function(x) {
  if (!is.character(x)) x <- as.character(x)
  sub(".*_([^_]+)_asinh$", "\\1", x)
}

# Spectre::read.files copied but edited to pass extra args to flowCore::read.FCS
# (e.g. truncate_max_range=FALSE)
read.files2 <- function (file.loc = getwd(), file.type = c(".csv", ".fcs"),
                         files = NULL, nrows = NULL, do.embed.file.names = TRUE, header = TRUE, ...)
{
  file.type <- tryCatch(match.arg(file.type), error = function(e) {
    stop("Invalid value for 'file.type': must be either '.csv' or '.fcs'.",
         call. = FALSE)
  })
  abs.file.loc <- fs::path_abs(file.loc)
  if (!fs::dir_exists(abs.file.loc)) {
    stop(paste0("Directory not found: '", file.loc, "'\nAre you sure that location exists?"))
  }
  if (length(list.files(path = abs.file.loc, pattern = file.type)) ==
      0) {
    stop(paste0("We did not find any ", file.type, " files in ",
                file.loc, ".\nAre you sure this is the right place?"))
  }
  if (is.null(files)) {
    file.names <- list.files(path = abs.file.loc, pattern = file.type)
  }
  else {
    file.names <- files
  }
  if (file.type == ".csv") {
    message("Reading CSV files...")
    data.list <- lapply(seq_len(length(file.names)), function(i) {
      file.name <- file.names[i]
      message(paste("Reading", file.name))
      if (is.null(nrows)) {
        tempdata <- data.table::fread(file.path(abs.file.loc,
                                                file.name), check.names = FALSE, header = header)
      }
      else {
        message(paste0("Reading ", nrows, " rows (cells) per file"))
        tempdata <- data.table::fread(file.path(abs.file.loc,
                                                file.name), check.names = FALSE, header = header,
                                      nrows = nrows)
      }
      if (do.embed.file.names) {
        tempdata$FileName <- gsub(".csv", "", file.name)
        tempdata$FileNo <- i
      }
      return(tempdata)
    })
    names(data.list) <- gsub(".csv", "", file.names)
  }
  else if (file.type == ".fcs") {
    data.list <- lapply(seq_len(length(file.names)), function(i) {
      file.name <- file.names[i]
      message(paste("Reading", file.name))
      if (is.null(nrows)) {
        x <- flowCore::read.FCS(file.path(abs.file.loc,
                                          file.name), ...)
      }
      else {
        message(paste0("Reading ", nrows, " rows (cells) per file"))
        x <- flowCore::read.FCS(file.path(abs.file.loc,
                                          file.name), which.lines = nrows, ...)
      }
      channel_name <- x@parameters@data$name
      antibody_name <- x@parameters@data$desc
      tempdata_colnames <- ifelse(is.na(antibody_name),
                                  channel_name, paste0(channel_name, "_", antibody_name))
      tempdata <- data.table::as.data.table(Biobase::exprs(x))
      names(tempdata) <- tempdata_colnames
      if (do.embed.file.names) {
        tempdata$FileName <- gsub(".fcs", "", file.name)
        tempdata$FileNo <- i
      }
      return(tempdata)
    })
    names(data.list) <- gsub(".fcs", "", file.names)
  }
  return(data.list)
}


# ---------------------------------------------------------------------------- #
# Pipeline: FCS import -> FlowSOM metaclusters -> UMAP -> gp350+ -> mixed model

# Import FCS files and merge
#
# READ ORDER IS IMPORTANT
meta_order <- data.table::fread(file.path(MetaDirectory, metadata_file))
data.table::setorder(meta_order, Read_Order)
sample_files <- paste0(meta_order$FileName, ".fcs")
stopifnot(all(file.exists(file.path(InputDirectory, sample_files))))
data.list <- read.files2(
  file.loc = InputDirectory,
  file.type = ".fcs",
  files = sample_files,
  do.embed.file.name = TRUE,
  transformation = FALSE,
  truncate_max_range = FALSE
)
cell.dat <- Spectre::do.merge.files(dat = data.list)
cell.dat$FileName <- factor(cell.dat$FileName, levels = meta_order$FileName)

# Metadata
meta.dat <- data.table::fread(file.path(MetaDirectory, metadata_file))

required_cols <- c("FileName", "Sample_ID", "Group", "Patient_ID", "Batch", "Timing")
stopifnot(all(required_cols %in% names(meta.dat)))

# Trim character keys before anything is factored
chr_cols <- intersect(c("FileName", "Group", "Patient_ID", "Batch"), names(meta.dat))
meta.dat[, (chr_cols) := lapply(.SD, function(x) trimws(as.character(x))), .SDcols = chr_cols]

stopifnot(setequal(meta.dat$FileName, as.character(cell.dat$FileName)))   # 1:1 join, both ways
stopifnot(!anyDuplicated(meta.dat$FileName))
stopifnot(setequal(unique(meta.dat$Group), c("pre_relapse", "remission")))

# Supplementary Table 5 sanity check
s5 <- readxl::read_excel(s5_file)
s5_group <- ifelse(s5$Condition == "PreRelapse", "pre_relapse", "remission")
names(s5_group) <- s5$Sample_ID
s5_ttr <- s5$`Time to relapse (d)`[s5$Condition == "PreRelapse"]
names(s5_ttr) <- s5$Donor_ID[s5$Condition == "PreRelapse"]

stopifnot(setequal(meta.dat$Sample_ID, s5$Sample_ID))
stopifnot(all(meta.dat$Group == s5_group[meta.dat$Sample_ID]))
stopifnot(all(as.numeric(meta.dat$Timing) == s5_ttr[meta.dat$Patient_ID]))

meta.dat$FileName <- factor(meta.dat$FileName, levels = levels(cell.dat$FileName))
cell.dat <- do.add.cols(cell.dat, base.col = "FileName", add.dat = meta.dat, add.by = "FileName", rmv.ext = TRUE)
cell.dat$FileName <- droplevels(cell.dat$FileName)

# Drop unused detectors
drop_cols <- c("FJComp-BV570-A", "FJComp-BV750-A", "FJComp-BUV805-A")
stopifnot(all(drop_cols %in% names(cell.dat)))
cell.dat <- cell.dat[, !..drop_cols]

# Arcsinh transform of compensated fluorescence channels
fluoro <- grep("^FJComp-", names(cell.dat), value = TRUE)
stopifnot(all(fluoro %in% names(cell.dat)))
b.dat <- Spectre::do.asinh(dat = cell.dat, use.cols = fluoro, cofactor = asinh_cofactor)

# Downsample per file so high-count samples do not dominate FlowSOM
set.seed(1)
target_n <- 1000
b.sub <- b.dat[, .SD[sample(.N, min(.N, target_n))], by = FileName]

# Core markers + robust z-score (median/MAD) for FlowSOM training
train_markers <- c(
  "FJComp-BB515-A_CD20_asinh",
  "FJComp-APC-Cy7-A_CD27_asinh",
  "FJComp-BV711-A_CD21_asinh",
  "FJComp-BUV496-A_CD38_asinh",
  "FJComp-PE-Cy5-A_CD11C_asinh",
  "FJComp-BV605-A_CD5_asinh",
  "FJComp-BV786-A_CXCR4_asinh",
  "FJComp-BUV563-A_CD172_asinh",
  "FJComp-BUV661-A_CD18_asinh"
)
stopifnot(all(train_markers %in% names(b.sub)))
scaled_markers <- paste0(train_markers, "_scaled")
scale_fits <- lapply(train_markers, \(m) fit_scale_mad(b.sub[[m]]))
names(scale_fits) <- train_markers
b.sub[, (scaled_markers) := Map(\(m, fit) apply_scale(get(m), fit), train_markers, scale_fits)]
b.dat[, (scaled_markers) := Map(\(m, fit) apply_scale(get(m), fit), train_markers, scale_fits)]
stopifnot(all(scaled_markers %in% names(b.sub)))
stopifnot(all(scaled_markers %in% names(b.dat)))

# Train FlowSOM (10x10 SOM -> 16 metaclusters) on the down-sampled cells
meta_k <- 16
X_sub <- as.matrix(b.sub[, ..scaled_markers])
fsom <- FlowSOM::FlowSOM(
  X_sub,
  colsToUse = seq_len(ncol(X_sub)),
  xdim = 10, ydim = 10,
  scale = FALSE,  # already scaled
  silent = TRUE,
  seed = 1
)
set.seed(1)
fsom$metaclustering <- FlowSOM::metaClustering_consensus(fsom$map$codes, k = meta_k)

# Project the trained SOM onto the full dataset
X_full <- as.matrix(b.dat[, ..scaled_markers])
node_full <- FlowSOM:::MapDataToCodes(fsom$map$codes, X_full)
b.dat[, FlowSOM_cluster := as.integer(node_full[, 1])]
b.dat[, FlowSOM_dist := as.numeric(node_full[, 2])]   # node distance (QC)
b.dat[, FlowSOM_metacluster := as.integer(fsom$metaclustering[FlowSOM_cluster])]

# UMAP on the full dataset
set.seed(1)
b.dat <- Spectre::run.umap(dat = b.dat, use.cols = scaled_markers)

# Metacluster palette
colors <- c(
  "6"="#A65628", "1"="#FB8072", "3"="#FFFFB3", "4"="#8DD3C7", "7"="#CCEBC5",
  "10"="#FFFF33", "14"="#984EA3", "16"="#D9D9D9", "13"="#E41A1C", "9"="#377EB8",
  "8"="#FF7F00", "15"="#BEBADA", "2"="#B3DE69", "11"="#4DAF4A", "12"="#80B1D3",
  "5"="#FDB462"
)

# Drop metacluster 6 (high-gp350 debris/contaminant that would bias downstream).
# This is the "one cluster removed as likely non B cell debris" of the
# Extended Data Fig. 8c legend; 15 annotated metaclusters remain.
b.m <- b.dat
setDT(b.m)
b.m <- b.m[FlowSOM_metacluster != 6]
colors <- colors[names(colors) != "6"]

# Cluster ordering (from core-marker median z-scores) + identities
cm <- b.m[, lapply(.SD, median), by = FlowSOM_metacluster, .SDcols = train_markers]
cmat <- as.matrix(cm[, ..train_markers])
row.names(cmat) <- as.character(cm$FlowSOM_metacluster)
cmat <- scale(cmat)
row_ord <- row.names(cmat)[hclust(dist(cmat))$order]
col_ord <- colnames(cmat)[hclust(dist(t(cmat)))$order]

ca <- c(
  "1"="1. Classical memory",
  "3"="3. Resting memory",
  "4"="4. Resting/intermediate memory",
  "7"="7. Activated memory",
  "10"="10. CD11c+ activated memory",
  "14"="14. CD11c++ activated memory",
  "16"="16. CD21lo memory",
  "13"="13. Atypical (ABC)",
  "9"="9. Pre-ABC",
  "8"="8. Plasmablast",
  "15"="15. ABC-like",
  "2"="2. Naive",
  "11"="11. CD11c+ memory",
  "12"="12. Transitional/immature-like",
  "5"="5. Anergic/naive-leaning"
)
stopifnot(setdiff(names(ca), row_ord) == 0)

# Run day
gp350_col <- "FJComp-Alexa Fluor 647-A_GP350_asinh"
stopifnot(gp350_col %in% names(b.m))
b.m[, exp_day := Batch]
stopifnot(identical(
  b.m$exp_day,
  local({
    ro <- meta_order$Read_Order[match(b.m$FileName, meta_order$FileName)]
    stopifnot(!anyNA(ro))
    fifelse(ro <= 14, "exp1", fifelse(ro <= 30, "exp2", "exp3"))
  })
))

# FMO-derived gp350+ threshold per run day (active FPR = 0.01)
thr_tbl <- data.table(
  exp_day = c("exp1", "exp2", "exp3"),
  thr_gp350 = c(4.601678, 4.884224, 4.209576)
)
b.m <- merge(b.m, thr_tbl, by = "exp_day", all.x = TRUE)
b.m[, gp350_pos := get(gp350_col) > thr_gp350]

# Per-sample gp350+ counts (feed the mixed model, Fig. 4g and ED Fig. 8d/e)
per_sample <- b.m[, .(
  n_B = .N,
  n_gp350 = sum(gp350_pos, na.rm = TRUE),
  pct_gp350 = 100 * mean(gp350_pos, na.rm = TRUE)
), by = .(Patient_ID, FileName, Group, exp_day, Timing)]

# Batch-adjusted binomial mixed model (Group + run day, patient random effect)
per_sample$Group <- factor(per_sample$Group)
per_sample$Group <- relevel(per_sample$Group, ref = "remission")
per_sample$exp_day <- factor(per_sample$exp_day)
fit <- glmer(
  cbind(n_gp350, n_B - n_gp350) ~ Group + exp_day + (1 | Patient_ID),
  data = per_sample,
  family = binomial()
)

message("R/facs_flowsom_pipeline.R: ", nrow(b.m), " CD19+ events across ",
        length(unique(b.m$FileName)), " samples; ",
        length(unique(b.m$FlowSOM_metacluster)), " metaclusters retained.")
