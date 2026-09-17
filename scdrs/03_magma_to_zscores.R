# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
# Updated: 7-25-2026
#
# scDRS step 03 - MAGMA Entrez IDs -> gene symbols, in scDRS z-score format
#
# MAGMA reports gene-level statistics against Entrez IDs; scDRS matches gene
# sets to an AnnData object by symbol. Where several Entrez IDs map to the same
# symbol, the one with the largest |Z| is kept.
#
# Inputs : scdrs/work/MS_GWAS.genes.out         (from 02_magma.sh)
# Outputs: scdrs/work/MS_magma_symbol_zscore.tsv
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

for (p in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("'", p, "' is required for the Entrez -> symbol mapping. ",
         "BiocManager::install('", p, "')")
}

WORK_DIR <- "scdrs/work"
in_file  <- file.path(WORK_DIR, "MS_GWAS.genes.out")
out_file <- file.path(WORK_DIR, "MS_magma_symbol_zscore.tsv")

stopifnot(file.exists(in_file))

# MAGMA gene-level output
go <- data.table::fread(in_file)
stopifnot(all(c("GENE", "ZSTAT") %in% names(go)))

dt <- go[, .(ENTREZID = as.character(GENE), MS = ZSTAT)]

# Entrez -> SYMBOL
map <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys     = unique(dt$ENTREZID),
  keytype  = "ENTREZID",
  columns  = "SYMBOL"
)
map <- data.table::as.data.table(map)[!is.na(SYMBOL)]

# Merge, and where a symbol has several Entrez IDs keep the strongest |Z|
dt2 <- merge(dt, map, by = "ENTREZID")
data.table::setorder(dt2, -abs(MS))
dt2 <- dt2[!duplicated(SYMBOL)]

out <- dt2[, .(GENE = SYMBOL, MS = MS)]
data.table::fwrite(out, out_file, sep = '\t')

message("wrote ", out_file, "  (", nrow(out), " genes)")
message("  strongest: ",
        paste(utils::head(out[order(-MS)]$GENE, 10), collapse = ", "))
