
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-22-2026
#
# scRNA import / demultiplexing helper functions
# ==============================================================================


get_pool_id <- function(x,
                        prefix = c("BW","W"),
                        normalize = TRUE) {
  stopifnot(is.character(x))

  # normalize path separators + trim trailing slashes
  x2 <- gsub("\\\\", "/", x)
  x2 <- sub("/+$", "", x2)

  pat <- "(?<![A-Za-z0-9])((?:BW\\d+)|(?:W\\d{1,2}))(?![A-Za-z0-9])"
  m <- regmatches(x2, regexpr(pat, x2, perl = TRUE))

  m <- ifelse(is.na(m) | m == "", NA_character_, m)

  if (normalize) {
    # BW01 -> BW1
    bw <- grepl("^BW\\d+$", m) & !is.na(m)
    m[bw] <- sub("^BW0+([0-9]+)$", "BW\\1", m[bw])

    # W1 -> W01, W01 stays W01
    w <- grepl("^W\\d{1,2}$", m) & !is.na(m)
    m[w] <- paste0("W", sprintf("%02d", as.integer(sub("^W", "", m[w]))))
  }

  m
}



infer_unpooled_donor <- function(pool_id, md) {
  sub <- md[md$Pool_ID == pool_id, , drop = FALSE]
  pats <- unique(sub$Patient_ID)
  if (length(pats) == 1) return(pats[1])
  NA_character_
}


# ---------------------------------------------------------------------------- #
# Super important
# Two library-metadata schemes exist for this study; both describe the same 52
# 10x libraries. The first was used for the pilot and early stages when this was
# scRNA only. The later is when we incorporated all remaining validation
# datasets
#
#   internal    the working sheet used during the analysis
#               (metadata/scrna_metadata.xlsx)
#   published   Supplementary Table S4
#               ('Supplementary Information/S4 scRNA Libraries.xlsx')
#
# Note the published table renames several columns:
#
#   published                internal              note
#   -----------------------  -------------------   ---------------------------
#   Sample_ID  ('S077')      Public_Sample_ID      de-identified publication ID
#   Prior_Sample_ID          Sample_ID             analysis ID, MS01_PreRelapse
#   Donor_ID   ('MS001')     Patient_ID ('MS01')   zero padding removed
#   Time to relapse (d)      Time_to_Relapse
#   Sex                      Patient_Sex
#   Age (yrs)                Patient_Age
#   Treatment                Patient_Treatment
#   Original_ID, Diagnosis   (carried through unchanged)
#
# ---------------------------------------------------------------------------- #

normalize_library_meta <- function(library_meta, verbose = TRUE) {
  stopifnot(is.data.frame(library_meta))
  df <- as.data.frame(library_meta, stringsAsFactors = FALSE, check.names = FALSE)

  published <- all(c("Donor_ID", "Prior_Sample_ID") %in% colnames(df))
  if (!published) {
    if (verbose) message("library metadata: internal schema detected; used as-is.")
    return(df)
  }
  if (verbose) message("library metadata: Supplementary Table S4 schema detected; ",
                       "mapping to internal column names.")

  # Publication sample ID kept for traceability; analysis ID restored.
  df$Public_Sample_ID <- as.character(df$Sample_ID)
  df$Sample_ID        <- as.character(df$Prior_Sample_ID)

  # Donor_ID 'MS001' -> Patient_ID 'MS01'
  stem <- sub("^([A-Za-z]+)0*([0-9]+)$", "\\1", as.character(df$Donor_ID))
  num  <- suppressWarnings(as.integer(sub("^([A-Za-z]+)0*([0-9]+)$", "\\2",
                                          as.character(df$Donor_ID))))
  if (anyNA(num) || any(!nzchar(stem))) {
    stop("Could not parse Donor_ID into letter prefix + number: ",
         paste(unique(df$Donor_ID[is.na(num) | !nzchar(stem)]), collapse = ", "))
  }
  df$Patient_ID <- sprintf("%s%02d", stem, num)

  # Cross-check against the patient prefix of the analysis Sample_ID.
  from_sample <- sub("_.*$", "", df$Sample_ID)
  bad <- which(df$Patient_ID != from_sample)
  if (length(bad) > 0) {
    warning("Donor_ID and Prior_Sample_ID disagree on the patient for ",
            length(bad), " row(s): ",
            paste(sprintf("%s (Donor_ID->%s, Sample_ID->%s)",
                          df$Library_ID[bad], df$Patient_ID[bad],
                          from_sample[bad]), collapse = "; "))
  }

  renames <- c("Time to relapse (d)" = "Time_to_Relapse",
               "Sex"                 = "Patient_Sex",
               "Age (yrs)"           = "Patient_Age",
               "Treatment"           = "Patient_Treatment")
  for (from in names(renames)) {
    to <- unname(renames[[from]])
    if (from %in% colnames(df)) {
      if (to %in% colnames(df)) {
        warning("both '", from, "' and '", to, "' present; keeping '", to, "'.")
      } else {
        df[[to]] <- df[[from]]
      }
      df[[from]] <- NULL
    }
  }

  # Drop the now-redundant staging column; keep Donor_ID / Original_ID /
  # Diagnosis as additional descriptive metadata.
  df$Prior_Sample_ID <- NULL

  if (verbose) {
    message("  mapped ", nrow(df), " libraries; ",
            length(unique(df$Patient_ID)), " donors; ",
            length(unique(df$Sample_ID)), " samples.")
  }
  df
}




# Validate library-level metadata for pooled/unpooled 10x design.
#
# Accepts either the internal sheet or Supplementary Table S4 - the published
# schema is mapped onto the internal one by normalize_library_meta() first.
#
# Required (post-normalization) columns:
# Library_ID, Sample_ID, Pool_ID, Patient_ID, Condition, Time_to_Relapse, Pooled, Batch,
# Replicate, Patient_Sex, Patient_Age, Patient_Treatment
#
# Returns the normalized data frame invisibly, so call it as:
#     md <- validate_library_meta(md)

validate_library_meta <- function(library_meta,
                                  pools = NULL,
                                  allow_multiple_libraries_per_sample = TRUE,
                                  allow_multiple_patients_per_pool = TRUE,
                                  require_condition_levels = c("PreRelapse","Remission","Relapse","Healthy"),
                                  normalize = TRUE,
                                  verbose = TRUE) {
  stopifnot(is.data.frame(library_meta))

  if (normalize) library_meta <- normalize_library_meta(library_meta, verbose = verbose)

  req <- c("Library_ID","Sample_ID","Pool_ID","Patient_ID","Condition","Time_to_Relapse",
           "Pooled","Batch","Replicate","Patient_Sex","Patient_Age","Patient_Treatment")
  miss <- setdiff(req, colnames(library_meta))
  if (length(miss) > 0) stop("library_meta missing required columns: ", paste(miss, collapse=", "))

  df <- library_meta

  # Optional pool subset
  if (!is.null(pools)) {
    pools <- unique(pools)
    df <- df[df$Pool_ID %in% pools, , drop = FALSE]
    if (nrow(df) == 0) stop("No rows in library_meta after filtering to pools: ", paste(pools, collapse=", "))
  }

  # Basic uniqueness checks
  if (anyDuplicated(df$Library_ID) > 0) {
    dups <- unique(df$Library_ID[duplicated(df$Library_ID)])
    stop("Library_ID must be unique. Duplicates: ", paste(dups, collapse=", "))
  }

  if (any(is.na(df$Pool_ID)) || any(df$Pool_ID == "")) stop("Pool_ID contains NA/empty.")
  if (any(is.na(df$Patient_ID)) || any(df$Patient_ID == "")) stop("Patient_ID contains NA/empty.")
  if (any(is.na(df$Sample_ID)) || any(df$Sample_ID == "")) stop("Sample_ID contains NA/empty.")
  if (any(is.na(df$Condition)) || any(df$Condition == "")) stop("Condition contains NA/empty.")

  # Condition sanity (do not hard fail on unexpected, but warn)
  conds <- unique(as.character(df$Condition))
  unexpected <- setdiff(conds, require_condition_levels)
  if (length(unexpected) > 0) {
    warning("Unexpected Condition level(s) found: ", paste(unexpected, collapse=", "),
            ". Expected subset of: ", paste(require_condition_levels, collapse=", "))
  }

  # Sample_ID should map to exactly one Patient_ID (note, biospecimen is patient+timepoint instance)
  tab_sample_patient <- aggregate(Patient_ID ~ Sample_ID, df, function(x) length(unique(x)))
  bad1 <- tab_sample_patient$Sample_ID[tab_sample_patient$Patient_ID != 1]
  if (length(bad1) > 0) {
    ex <- df[df$Sample_ID %in% bad1, c("Sample_ID","Patient_ID","Pool_ID","Library_ID","Condition","Replicate"), drop=FALSE]
    stop("Sample_ID must map to exactly 1 Patient_ID. Violations for Sample_ID: ",
         paste(unique(bad1), collapse=", "), "\nExamples:\n", paste(capture.output(print(ex)), collapse="\n"))
  }

  # Sample_ID should map to exactly one Condition
  tab_sample_cond <- aggregate(Condition ~ Sample_ID, df, function(x) length(unique(x)))
  bad2 <- tab_sample_cond$Sample_ID[tab_sample_cond$Condition != 1]
  if (length(bad2) > 0) {
    ex <- df[df$Sample_ID %in% bad2, c("Sample_ID","Patient_ID","Pool_ID","Library_ID","Condition","Replicate"), drop=FALSE]
    stop("Sample_ID must map to exactly 1 Condition. Violations for Sample_ID: ",
         paste(unique(bad2), collapse=", "), "\nExamples:\n", paste(capture.output(print(ex)), collapse="\n"))
  }

  # Within a Pool_ID, a Patient_ID should map to exactly one Sample_ID
  # makes (Pool_ID, Patient_ID) join unambiguous.)
  tab_pool_patient_sample <- aggregate(Sample_ID ~ Pool_ID + Patient_ID, df, function(x) length(unique(x)))
  bad3 <- tab_pool_patient_sample[tab_pool_patient_sample$Sample_ID != 1, c("Pool_ID","Patient_ID")]
  if (nrow(bad3) > 0) {
    ex <- merge(df, bad3, by=c("Pool_ID","Patient_ID"))
    ex <- ex[, c("Pool_ID","Patient_ID","Sample_ID","Library_ID","Condition","Replicate"), drop=FALSE]
    stop("Within a pool, each Patient_ID must map to exactly 1 Sample_ID. Violations:\n",
         paste(capture.output(print(unique(ex))), collapse="\n"),
         "\nIf you truly have multiple samples for the same patient in the same pool, you must include a more specific join key (e.g., Library_ID).")
  }

  # Replicate behavior: allow multiple libraries per Sample_ID only if allow_multiple_libraries_per_sample
  tab_sample_libs <- aggregate(Library_ID ~ Sample_ID, df, function(x) length(unique(x)))
  bad4 <- tab_sample_libs$Sample_ID[tab_sample_libs$Library_ID > 1]
  if (length(bad4) > 0 && !allow_multiple_libraries_per_sample) {
    ex <- df[df$Sample_ID %in% bad4, c("Sample_ID","Library_ID","Pool_ID","Patient_ID","Condition","Replicate"), drop=FALSE]
    stop("Multiple libraries per Sample_ID detected but allow_multiple_libraries_per_sample=FALSE. Examples:\n",
         paste(capture.output(print(ex)), collapse="\n"))
  }

  # Pooled flag sanity: if Pooled==FALSE, pool should contain exactly 1 patient (unless allow_multiple_patients_per_pool)
  # This catches metadata errors where a solo library is incorrectly labeled.
  pool_n_pat <- aggregate(Patient_ID ~ Pool_ID, df, function(x) length(unique(x)))
  names(pool_n_pat)[2] <- "n_patients"
  pool_pooled <- aggregate(Pooled ~ Pool_ID, df, function(x) length(unique(as.character(x))))
  # if Pool_ID has mixed Pooled values, that's a metadata problem
  mixed_pooled <- pool_pooled$Pool_ID[pool_pooled$Pooled != 1]
  if (length(mixed_pooled) > 0) {
    ex <- df[df$Pool_ID %in% mixed_pooled, c("Pool_ID","Library_ID","Sample_ID","Patient_ID","Pooled"), drop=FALSE]
    stop("Pool_ID has inconsistent Pooled values across rows for pool(s): ",
         paste(mixed_pooled, collapse=", "), "\nExamples:\n", paste(capture.output(print(ex)), collapse="\n"))
  }

  pooled_status <- aggregate(Pooled ~ Pool_ID, df, function(x) unique(as.character(x)))
  pooled_status <- merge(pooled_status, pool_n_pat, by="Pool_ID", all.x=TRUE)
  solo_bad <- pooled_status$Pool_ID[pooled_status$Pooled == "FALSE" & pooled_status$n_patients != 1]
  if (length(solo_bad) > 0 && !allow_multiple_patients_per_pool) {
    ex <- df[df$Pool_ID %in% solo_bad, c("Pool_ID","Library_ID","Sample_ID","Patient_ID","Pooled"), drop=FALSE]
    stop("Pools labeled Pooled==FALSE must have exactly 1 patient (unless allow_multiple_patients_per_pool=TRUE). Bad pools: ",
         paste(solo_bad, collapse=", "), "\nExamples:\n", paste(capture.output(print(ex)), collapse="\n"))
  }

  # Time_to_Relapse 
  tab_sample_ttr <- aggregate(Time_to_Relapse ~ Sample_ID, df, function(x) length(unique(x[!is.na(x)])))
  bad6 <- tab_sample_ttr$Sample_ID[tab_sample_ttr$Time_to_Relapse > 1]
  if (length(bad6) > 0) {
    ex <- df[df$Sample_ID %in% bad6, c("Sample_ID","Patient_ID","Library_ID","Time_to_Relapse"), drop=FALSE]
    stop("Time_to_Relapse must be consistent within Sample_ID. Violations for Sample_ID: ",
         paste(unique(bad6), collapse=", "), "\nExamples:\n", paste(capture.output(print(ex)), collapse="\n"))
  }

  if (verbose) {
    message("library_meta validation passed.",
            "\n  n_rows: ", nrow(df),
            "\n  n_pools: ", length(unique(df$Pool_ID)),
            "\n  n_patients: ", length(unique(df$Patient_ID)),
            "\n  n_samples (patient-timepoint): ", length(unique(df$Sample_ID)),
            "\n  Condition counts:\n",
            paste(capture.output(print(sort(table(df$Condition), decreasing=TRUE))), collapse="\n"),
            "\n  Treatment counts:\n",
            paste(capture.output(print(sort(table(df$Patient_Treatment), decreasing=TRUE))), collapse="\n"))
  }

  invisible(df)
}




import_counts_vireo <- function(
    data_dir10x,
    pool_id,
    library_meta = NULL,
    vireo_dir = NULL,
    expected_map = NULL,
    donor_id_unpooled = NULL,
    library_id = NULL,
    barcode_suffix = "-1",
    allow_partial_intersection = FALSE,
    allow_missing_vireo = TRUE,
    verbose = TRUE,
    vireo_filename = "donor_ids.tsv"
) {
  stopifnot(dir.exists(data_dir10x))
  if (is.null(library_id)) library_id <- pool_id
  if (verbose) message(pool_id)

  # ---- normalize data_dir10x if given ".../outs" ----
  data_dir10x_norm <- normalizePath(data_dir10x, winslash = "/", mustWork = FALSE)
  if (basename(data_dir10x_norm) == "outs") {
    cand <- c(
      file.path(data_dir10x, "filtered_feature_bc_matrix"),
      file.path(data_dir10x, "filtered_feature_bc_matrix.h5"),
      file.path(data_dir10x, "raw_feature_bc_matrix")
    )
    if (dir.exists(cand[1])) {
      data_dir10x <- cand[1]
    } else if (file.exists(cand[2])) {
      data_dir10x <- cand[2]
    } else if (dir.exists(cand[3])) {
      data_dir10x <- cand[3]
    } else {
      stop("data_dir10x points to an 'outs' directory but no matrix found under it. Tried: ",
           paste(cand, collapse = ", "))
    }
  }

  # ---- read counts ----
  counts <- Seurat::Read10X(data.dir = data_dir10x)
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }
  cb_counts <- colnames(counts)
  if (is.null(cb_counts) || length(cb_counts) == 0) {
    stop("Read10X produced an empty matrix for: ", data_dir10x)
  }

  # ---- expected donors for this pool ----
  expected <- NULL
  if (!is.null(expected_map)) {
    if (!(pool_id %in% names(expected_map))) stop("pool_id not found in expected_map: ", pool_id)
    expected <- unique(as.character(expected_map[[pool_id]]))
    expected <- expected[!is.na(expected) & nzchar(expected)]
    if (length(expected) == 0) expected <- NULL
  }

  # ---- read vireo (if present) ----
  vireo_fl  <- if (!is.null(vireo_dir)) file.path(vireo_dir, vireo_filename) else NA_character_
  has_vireo <- !is.null(vireo_dir) && file.exists(vireo_fl)

  if (!has_vireo && !allow_missing_vireo) {
    stop("Vireo missing for pool_id=", pool_id, ". Set allow_missing_vireo=TRUE to import counts only.")
  }

  if (has_vireo) {
    demux <- utils::read.table(
      vireo_fl, header = TRUE, stringsAsFactors = FALSE,
      sep = "\t", quote = "", comment.char = ""
    )
    req <- c("cell","donor_id","prob_max","prob_doublet","n_vars",
             "best_singlet","best_doublet","doublet_logLikRatio")
    miss <- setdiff(req, colnames(demux))
    if (length(miss) > 0) stop("Missing required Vireo columns: ", paste(miss, collapse=", "))

    # barcode suffix harmonization
    cb_vireo <- demux$cell
    if (any(grepl(paste0(barcode_suffix, "$"), cb_counts)) &&
        !any(grepl(paste0(barcode_suffix, "$"), cb_vireo))) {
      demux$cell <- paste0(demux$cell, barcode_suffix)
    } else if (!any(grepl(paste0(barcode_suffix, "$"), cb_counts)) &&
               any(grepl(paste0(barcode_suffix, "$"), cb_vireo))) {
      demux$cell <- sub(paste0(barcode_suffix, "$"), "", demux$cell)
    }

    in_both <- intersect(colnames(counts), demux$cell)
    if (length(in_both) == 0) stop("No overlapping barcodes between counts and vireo after suffix normalization.")

    if (!allow_partial_intersection) {
      if (length(setdiff(colnames(counts), demux$cell)) > 0 ||
          length(setdiff(demux$cell, colnames(counts))) > 0) {
        stop("Barcode sets not identical. Investigate vireo vs cellranger filtering.")
      }
    }

    counts <- counts[, colnames(counts) %in% in_both, drop = FALSE]
    demux  <- demux[match(colnames(counts), demux$cell), , drop = FALSE]
    stopifnot(identical(demux$cell, colnames(counts)))

    demux$donor_id_raw      <- demux$donor_id
    demux$Pool_ID           <- pool_id                # ALWAYS set for every cell
    demux$library_id_import <- library_id
    demux$has_vireo         <- TRUE
    demux$doublet_method    <- "vireo"

    # flag unexpected donors (preserve 'doublet'/'unassigned')
    if (!is.null(expected)) {
      demux$donor_id <- ifelse(
        demux$donor_id %in% c(expected, "doublet", "unassigned"),
        demux$donor_id,
        paste0("misassigned_", demux$donor_id)
      )
    }

    md_cell <- demux
    rownames(md_cell) <- md_cell$cell
    md_cell$cell <- NULL

  } else {
    if (verbose) message("Vireo missing for pool_id=", pool_id, "; importing counts only.")
    n <- length(cb_counts)
    md_cell <- data.frame(
      donor_id            = rep(NA_character_, n),
      donor_id_raw        = rep(NA_character_, n),
      prob_max            = rep(NA_real_, n),
      prob_doublet        = rep(NA_real_, n),
      n_vars              = rep(NA_integer_, n),
      best_singlet        = rep(NA_character_, n),
      best_doublet        = rep(NA_character_, n),
      doublet_logLikRatio = rep(NA_real_, n),
      Pool_ID             = rep(pool_id, n),          # ALWAYS set
      library_id_import   = rep(library_id, n),
      has_vireo           = rep(FALSE, n),
      doublet_method      = rep("scrna", n),
      stringsAsFactors    = FALSE
    )
    rownames(md_cell) <- cb_counts
  }

  # ---- join library_meta (ONLY for mappable donors; NEVER overwrite Pool_ID) ----
  if (!is.null(library_meta)) {
    req_meta <- c("Library_ID","Sample_ID","Pool_ID","Patient_ID","Condition","Time_to_Relapse",
                  "Pooled","Batch","Replicate","Patient_Sex","Patient_Age","Patient_Treatment")
    miss_meta <- setdiff(req_meta, colnames(library_meta))
    if (length(miss_meta) > 0) stop("library_meta missing required columns: ", paste(miss_meta, collapse=", "))

    lib_sub <- library_meta[library_meta$Pool_ID == pool_id, , drop = FALSE]
    if (nrow(lib_sub) == 0) stop("library_meta has no rows for Pool_ID=", pool_id)

    # if no vireo: assign the single patient for this pool
    if (!md_cell$has_vireo[1]) {
      if (is.null(donor_id_unpooled)) {
        if (nrow(lib_sub) != 1) {
          stop("No Vireo and library_meta has ", nrow(lib_sub), " rows for Pool_ID=", pool_id,
               ". Supply donor_id_unpooled.")
        }
        donor_id_unpooled <- lib_sub$Patient_ID[1]
      }
      md_cell$donor_id     <- donor_id_unpooled
      md_cell$donor_id_raw <- donor_id_unpooled
      md_cell$best_singlet <- donor_id_unpooled
    }

    # donors that should map to library_meta
    mappable <- !(md_cell$donor_id %in% c("doublet","unassigned")) &
      !grepl("^misassigned_", md_cell$donor_id) &
      !is.na(md_cell$donor_id) & nzchar(md_cell$donor_id)

    key_md  <- paste(pool_id, md_cell$donor_id, sep="|")
    key_lib <- paste(lib_sub$Pool_ID, lib_sub$Patient_ID, sep="|")
    m <- match(key_md, key_lib)

    # only fill metadata for mappable cells, and do not overwrite existing columns (esp Pool_ID)
    cols_to_add <- setdiff(colnames(lib_sub), intersect(colnames(lib_sub), colnames(md_cell)))
    cols_to_add <- setdiff(cols_to_add, "Pool_ID")  # defensive

    for (cc in cols_to_add) {
      md_cell[[cc]] <- NA
      md_cell[[cc]][mappable] <- lib_sub[[cc]][m[mappable]]
    }

    # also allow filling existing-but-empty cols safely (common ones)
    safe_fill <- setdiff(colnames(lib_sub), "Pool_ID")
    safe_fill <- intersect(safe_fill, colnames(md_cell))
    safe_fill <- setdiff(safe_fill, c("Pool_ID")) # keep Pool_ID untouched

    for (cc in safe_fill) {
      # only fill where currently NA AND mappable
      is_na <- is.na(md_cell[[cc]])
      md_cell[[cc]][is_na & mappable] <- lib_sub[[cc]][m[is_na & mappable]]
    }

    # validate mapping for mappable cells
    unmapped <- mappable & is.na(md_cell$Sample_ID)
    if (any(unmapped)) {
      bad <- unique(md_cell$donor_id[unmapped])
      stop("Unmapped mappable donors in pool ", pool_id, ": ", paste(bad, collapse=", "),
           ". Check library_meta Pool_ID/Patient_ID entries.")
    }
  }

  # ---- scDist-friendly columns ----
  if (!("Patient_ID" %in% colnames(md_cell)) && "donor_id" %in% colnames(md_cell)) {
    md_cell$Patient_ID <- md_cell$donor_id
  }
  md_cell$Patient <- md_cell$Patient_ID

  if (verbose) {
    message("donor_id table (", pool_id, "):")
    print(sort(table(md_cell$donor_id), decreasing = TRUE))
    if ("Condition" %in% colnames(md_cell)) {
      message("Condition table (", pool_id, "):")
      print(sort(table(md_cell$Condition), decreasing = TRUE))
    }
    message("Pool_ID NA count: ", sum(is.na(md_cell$Pool_ID)))
  }

  seu <- Seurat::CreateSeuratObject(counts = counts, meta.data = md_cell, min.cells = 0, min.features = 0)
  seu <- Seurat::RenameCells(seu, add.cell.id = pool_id)
  seu
}
