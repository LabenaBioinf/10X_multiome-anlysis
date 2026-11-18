# Ambient RNA correction with SoupX on the raw Cell Ranger ARC droplet matrices, then pseudobulk.
#
# For each of the 20 libraries: the ambient profile is estimated from the empty droplets of the raw matrix, the
# retained nuclei of the analysis object (same barcodes, same annotation) are the cells, the annotated populations are
# the clusters, the contamination fraction is estimated with autoEstCont (if that fails, from immunoglobulin genes in
# non-B populations), and counts are corrected with adjustCounts (integers). Corrected counts are summed per
# population and library exactly as in pseudobulk_edgeR.R.
# Check: uncorrected counts taken from the raw matrices must reproduce the cached pseudobulk counts of the object.
#
# usage: Rscript soupx_pseudobulk.R <harmonized_seurat.rds> <Cell Ranger results dir> <cached pseudobulk_counts.rds> <out dir> <R lib>
args <- commandArgs(trailingOnly = TRUE)
OBJ <- args[1]; CR <- args[2]; PB0 <- args[3]; OUT <- args[4]
if (length(args) > 4) .libPaths(c(args[5], .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(SoupX) })
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
mdf <- file.path(OUT, "nuclei_metadata.rds")
if (!file.exists(mdf)) {
  obj <- readRDS(OBJ); md <- obj@meta.data[, c("orig.ident", "Batch", "Condition", "predicted.id")]
  md$cell <- colnames(obj); rm(obj); gc(); saveRDS(md, mdf)
}
md <- readRDS(mdf)
md$lib_code <- sub("_.*$", "", md$cell)                       # add.cell.ids = Cell Ranger sample (c..)
md$barcode  <- sub("_[0-9]+$", "", sub("^[^_]+_", "", md$cell))
cat("nuclei:", nrow(md), "| library codes:", length(unique(md$lib_code)), "\n")
print(table(md$lib_code, md$orig.ident)[, 1:4])
stopifnot(all(tapply(md$orig.ident, md$lib_code, function(x) length(unique(x))) == 1))
A <- readRDS(PB0); pb0 <- A$counts; info <- A$info
ig_pat <- "^IG[HKL][VDJC]|^IGH[AGMDE][0-9]?$|^IGKC$|^IGLC[0-9]$|^JCHAIN$"
BPOP <- c("B naive", "B intermediate", "B memory", "Plasmablast")
pb_raw <- list(); pb_cor <- list(); rho <- list()
for (lc in sort(unique(md$lib_code))) {
  m <- md[md$lib_code == lc, ]
  f <- file.path(CR, lc, "outs", "raw_feature_bc_matrix.h5")
  tod <- Read10X_h5(f); if (is.list(tod)) tod <- tod[["Gene Expression"]]
  stopifnot(all(m$barcode %in% colnames(tod)))
  toc <- tod[, m$barcode]
  sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
  cl <- setNames(as.character(m$predicted.id), m$barcode)
  sc <- setClusters(sc, cl)
  est <- tryCatch({ s2 <- autoEstCont(sc, doPlot = FALSE, verbose = FALSE); list(sc = s2, method = "autoEstCont") },
                  error = function(e) { cat(lc, "autoEstCont failed:", conditionMessage(e), "\n"); NULL })
  # immunoglobulin-based estimate (non-B populations do not express these genes), always computed for comparison
  ig <- grep(ig_pat, rownames(toc), value = TRUE)
  nonB <- !(cl %in% BPOP)
  use <- matrix(nonB, ncol = 1, dimnames = list(names(cl), "IG"))   # nuclei x gene sets, as SoupX expects
  sig <- tryCatch(calculateContaminationFraction(sc, list(IG = ig), use, forceAccept = TRUE, verbose = FALSE), error = function(e) NULL)
  rho_ig <- if (is.null(sig)) NA else unique(sig$metaData$rho)[1]
  if (is.null(est)) { stopifnot(!is.null(sig)); est <- list(sc = sig, method = "immunoglobulin genes") }
  r <- unique(est$sc$metaData$rho)[1]
  soup_ig <- sum(sc$soupProfile[ig, "est"])
  set.seed(20261118)                                                  # adjustCounts(roundToInt = TRUE) rounds stochastically: seeded for reproducibility
  out <- adjustCounts(est$sc, roundToInt = TRUE, verbose = 0)
  agg <- function(x) { g <- factor(cl[colnames(x)]); M <- sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1, dims = c(length(g), nlevels(g)), dimnames = list(NULL, levels(g))); as.matrix(x %*% M) }
  a <- agg(toc); b <- agg(out); lib <- as.character(m$orig.ident[1])
  colnames(a) <- colnames(b) <- paste(colnames(a), lib, sep = "|")
  pb_raw[[lc]] <- a; pb_cor[[lc]] <- b
  rho[[lc]] <- data.frame(library = lib, condition = as.character(m$Condition[1]), nuclei = nrow(m), method = est$method, rho = r,
                          rho_immunoglobulin = rho_ig, soup_fraction_immunoglobulin = soup_ig,
                          counts_removed_fraction = 1 - sum(out) / sum(toc))
  cat(lc, lib, "nuclei", nrow(m), est$method, "rho", round(r, 3), "rho_IG", round(rho_ig, 3), "\n")
  rm(tod, toc, sc, out); gc()
}
genes <- Reduce(intersect, lapply(pb_raw, rownames))
RAW <- do.call(cbind, lapply(pb_raw, function(x) x[genes, , drop = FALSE])); COR <- do.call(cbind, lapply(pb_cor, function(x) x[genes, , drop = FALSE]))
# consistency with the object: same columns, same counts
common_g <- intersect(rownames(pb0), genes); common_c <- intersect(colnames(pb0), colnames(RAW))
cat("columns object/raw/common:", ncol(pb0), ncol(RAW), length(common_c), "| genes object/raw/common:", nrow(pb0), length(genes), length(common_g), "\n")
d <- abs(pb0[common_g, common_c] - RAW[common_g, common_c])
cat("uncorrected counts from raw matrices vs object: max abs difference", max(d), "| total counts", sum(pb0[common_g, common_c]), sum(RAW[common_g, common_c]), "\n")
info2 <- info[match(colnames(COR), info$column), ]; stopifnot(!any(is.na(info2$column)))
saveRDS(list(counts = COR, info = info2), file.path(OUT, "pseudobulk_counts_soupx.rds"))
saveRDS(list(counts = RAW, info = info2), file.path(OUT, "pseudobulk_counts_raw_check.rds"))
R <- do.call(rbind, rho)
R$pair <- NULL
write.table(R, file.path(OUT, "soupx_contamination_per_library_internal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (v in c("rho", "rho_immunoglobulin", "soup_fraction_immunoglobulin")) { cat(v, "median by condition:\n"); print(tapply(R[[v]], R$condition, median, na.rm = TRUE)) }
print(R[, c("library", "condition", "method", "rho", "rho_immunoglobulin")])
