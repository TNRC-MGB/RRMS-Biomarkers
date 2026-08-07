
# ==============================================================================
# Author: Devin A. King, PhD
#   Translational Neuroimmunology Research Center (TNRC) | Chitnis Lab
#   Mass General Brigham
#   Harvard Medical School
#
#
# Inputs : zenodo/bulk/counts/<Sample_ID>/quant.sf # salmon quantification
#          Supplementary Information/S3 PBMC Samples.xlsx
#          Supplementary Information/S7 B Cell Bulk RNA-seq Libraries.xlsx
#          bulk/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf
# Outputs: bulk/results/dge.rds            filtered, TMM-normalized DGEList
#          bulk/results/txi.rds            tximport objects
#          bulk/results/gene_counts.tsv
#
# ==============================================================================

setwd("C:/Users/devin/Desktop/rrms")

source("R/packages.R")

QUANT    <- "zenodo/bulk/counts"
OUT_DIR  <- "bulk/results"
GTF      <- "bulk/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
S3_XLSX  <- "Supplementary Information/S3 PBMC Samples.xlsx"
S7_XLSX  <- "Supplementary Information/S7 B Cell Bulk RNA-seq Libraries.xlsx"

EXCLUDE_IG_TR <- TRUE

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (p in c("tximport", "edgeR", "limma"))
  if (!requireNamespace(p, quietly = TRUE))
    stop(p, " is required (Bioconductor).")

# ==============================================================================



files <- Sys.glob(file.path(QUANT, "*", "quant.sf"))
if (!length(files)) stop("no quant.sf under ", QUANT, "  run 02_quantify.sh first.")
sids  <- basename(dirname(files))
files <- setNames(files[order(sids)], sort(sids))
message("quantifications found under ", QUANT, ": ", length(files))

s3 <- as.data.table(read_excel(S3_XLSX))
s7 <- as.data.table(read_excel(S7_XLSX))
need <- c("Sample_ID", "Donor_ID", "Condition", "Time to relapse (d)", "Sex")
miss <- setdiff(need, names(s3))
if (length(miss)) stop("S3 is missing: ", paste(miss, collapse = ", "))
setnames(s3, "Time to relapse (d)", "time_to_relapse_d")
s3[, time_to_relapse_d := as.numeric(time_to_relapse_d)]
message("S3 vials: ", nrow(s3), "    S7 libraries: ", nrow(s7))

meta <- s3[match(names(files), Sample_ID)]

# tx2gene
if (!file.exists(GTF)) stop("GTF not found: ", GTF)

message("reading tx2gene from ", GTF)
gl    <- fread(cmd = paste('grep -P "\\ttranscript\\t"', shQuote(GTF)),
               sep = '\t', header = FALSE, showProgress = FALSE)
attr9 <- gl$V9
tx2gene <- data.table(
  TXNAME = sub('.*transcript_id "([^"]+)".*', '\\1', attr9),
  GENEID = sub('.*gene_id "([^"]+)".*',       '\\1', attr9),
  SYMBOL = ifelse(grepl('gene_name "', attr9),
                  sub('.*gene_name "([^"]+)".*', '\\1', attr9), NA_character_))
tx2gene <- unique(tx2gene[TXNAME != attr9])
message("  transcripts: ", nrow(tx2gene), "   genes: ", uniqueN(tx2gene$GENEID))
if (!nrow(tx2gene)) stop("tx2gene came out empty  check the GTF attribute format.")

# Salmon's Name column must match TXNAME or every transcript is dropped
hit <- mean(fread(files[[1]], nrows = 5)$Name %in% tx2gene$TXNAME)
message("  transcript IDs in quant.sf matching the GTF: ", round(100 * hit), "%")
if (hit < 0.5)
  stop("quant.sf transcript names do not match the GTF (", round(100 * hit),
       "%). Usually the --gencode flag in 01_build_index.sh, or a version ",
       "suffix such as ENST00000456328.2 on one side only.")


# import
message("tximport (countsFromAbundance = lengthScaledTPM) ...")
txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene[, .(TXNAME, GENEID)],
                          countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
txi_raw <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene[, .(TXNAME, GENEID)],
                              countsFromAbundance = "no", ignoreTxVersion = TRUE)
saveRDS(list(lengthScaledTPM = txi, raw = txi_raw), file.path(OUT_DIR, "txi.rds"))

cts <- txi$counts
message("gene x sample matrix: ", nrow(cts), " x ", ncol(cts))
message("library size (summed counts): min ", format(round(min(colSums(cts))), big.mark = ","),
    "  median ", format(round(median(colSums(cts))), big.mark = ","),
    "  max ", format(round(max(colSums(cts))), big.mark = ","))


# immunoglobulin / T-cell receptor share
gsym  <- tx2gene[, .(SYMBOL = SYMBOL[1]), by = GENEID]
gsym  <- gsym[match(rownames(cts), GENEID)]
is_ig <- !is.na(gsym$SYMBOL) & grepl("^(IGH|IGK|IGL|TRA|TRB|TRD|TRG)[VDJC]", gsym$SYMBOL)
frac  <- colSums(cts[is_ig, , drop = FALSE]) / colSums(cts)

message("---- immunoglobulin / TCR segment genes --------------------------------")
message("genes matched: ", sum(is_ig))
message("share of counts: min ", round(100 * min(frac), 1), "%  median ",
    round(100 * median(frac), 1), "%  max ", round(100 * max(frac), 1), "%")
if (diff(range(frac)) > 0.15)
  message("  NOTE: the share varies by more than 15 points across samples. Left in,",
      ' that alone would move every other gene\'s CPM.')
if (EXCLUDE_IG_TR) {
  cts <- cts[!is_ig, , drop = FALSE]
  message("EXCLUDED before normalisation (EXCLUDE_IG_TR = TRUE). ", nrow(cts), " genes remain.")
} else {
  message("RETAINED (EXCLUDE_IG_TR = FALSE).")
}


# 5  filter and normalize
# filterByExpr is given the Condition grouping, so a gene expressed in only one
# condition is kept rather than filtered as low-abundance overall.
dge <- edgeR::DGEList(counts  = cts,
                      samples = data.frame(meta, row.names = meta$Sample_ID,
                                           check.names = FALSE))
dge$samples$group <- factor(meta$Condition)
keep_gene <- edgeR::filterByExpr(dge, group = dge$samples$group)

message("filterByExpr: ", sum(keep_gene), " of ", nrow(dge), " genes retained (",
    round(100 * mean(keep_gene)), "%)")
dge <- edgeR::calcNormFactors(dge[keep_gene, , keep.lib.sizes = FALSE])
message("TMM factors: ", round(min(dge$samples$norm.factors), 3), " .. ",
    round(max(dge$samples$norm.factors), 3))
if (max(dge$samples$norm.factors) / min(dge$samples$norm.factors) > 2)
  message("  NOTE: >2x spread in TMM factors. Worth checking whether it tracks the",
      " IG share or a submission batch.")

dge$genes <- data.frame(gene_id = rownames(dge),
                        symbol  = gsym$SYMBOL[match(rownames(dge), gsym$GENEID)],
                        row.names = rownames(dge))
saveRDS(dge, file.path(OUT_DIR, "dge.rds"))
fwrite(data.table(gene_id = rownames(dge), symbol = dge$genes$symbol,
                  as.data.frame(dge$counts)),
       file.path(OUT_DIR, "gene_counts.tsv"), sep = '\t')


message("wrote dge.rds, txi.rds, gene_counts.tsv to ", OUT_DIR)

message("R ", as.character(getRversion()), "   ",
    paste(vapply(c("tximport", "edgeR", "limma"),
                 function(p) paste0(p, " ", as.character(packageVersion(p))), ""),
          collapse = "   "))
