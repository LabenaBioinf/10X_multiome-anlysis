# Object-level analyses on the harmonised multiome object (flare vs remission).
#
# B1  TRBV11-2: is the gene in the RNA feature set; fraction of T-cell nuclei with
#     TRBV11-2 UMIs per sample, flare vs remission (paired, exact signed-rank);
#     detection of all TRBV genes for context (3' GEX captures V segments poorly).
# B2  γδ T identity: RNA detection of TRDC, TRGC1, TRGC2, TRDV2, TRGV9 per cell type,
#     and ATAC accessibility of peaks in the TRG locus (chr7) and at TRDC (chr14),
#     γδ T vs other T-cell populations (per sample, CPM of reads in locus peaks).
# B3  ATAC counts summed per cell type x sample (pseudobulk) are saved here; motif
#     deviations and the paired flare vs remission test are in b3_chromvar.R.
# Also: per-sample nuclei and detection fractions for key genes (cell-level context).

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix); library(GenomicRanges)
  library(SummarizedExperiment)
})
args <- commandArgs(trailingOnly = TRUE)
OBJ <- args[1]; OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261016)

obj <- readRDS(OBJ)
md <- obj@meta.data
md$cond <- ifelse(md$Condition == "MISC", "Flare", "Remission")
md$sample <- paste(md$Batch, md$cond, sep = "_")
DefaultAssay(obj) <- "RNA"
if (length(grep("^counts", Layers(obj[["RNA"]]))) > 1) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
rna <- GetAssayData(obj, assay = "RNA", layer = "counts")
atac <- GetAssayData(obj, assay = "ATAC", layer = "counts")
peaks <- granges(obj[["ATAC"]])
rm(obj); gc()
ct <- as.character(md$predicted.id)
t_types <- c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating", "Treg",
             "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "MAIT", "gdT", "dnT")
paired_exact <- function(x, y) {
  d <- x - y; d <- d[d != 0 & !is.na(d)]
  if (length(d) < 2) return(NA_real_)
  r <- rank(abs(d)); obs <- sum(r[d > 0]); mu <- sum(r) / 2
  s <- as.matrix(expand.grid(rep(list(c(0, 1)), length(d))))
  mean(abs(s %*% r - mu) >= abs(obs - mu) - 1e-9)
}

## ---- B1: TRBV11-2 -----------------------------------------------------------
trbv <- grep("^TRBV", rownames(rna), value = TRUE)
cat("TRBV genes in RNA features:", length(trbv), "| TRBV11-2 present:", "TRBV11-2" %in% rownames(rna), "\n")
is_t <- ct %in% t_types
det_trbv <- Matrix::colSums(rna[trbv, , drop = FALSE] > 0) > 0
cat("T-cell nuclei with any TRBV UMI:", round(100 * mean(det_trbv[is_t]), 2), "%\n")
b1 <- NULL
if ("TRBV11-2" %in% rownames(rna)) {
  pos <- rna["TRBV11-2", ] > 0
  b1 <- do.call(rbind, lapply(split(seq_len(nrow(md))[is_t], md$sample[is_t]), function(i) {
    data.frame(sample = md$sample[i[1]], patient = md$Batch[i[1]], cond = md$cond[i[1]],
               t_nuclei = length(i), trbv_any = sum(det_trbv[i]), trbv11_2_pos = sum(pos[i]),
               frac_of_T = mean(pos[i]), frac_of_TRBV_detected = sum(pos[i]) / max(1, sum(det_trbv[i])))
  }))
  w <- reshape(b1[, c("patient", "cond", "frac_of_T", "frac_of_TRBV_detected")], idvar = "patient",
               timevar = "cond", direction = "wide")
  cat("TRBV11-2+ fraction of T nuclei: median flare", signif(median(w$frac_of_T.Flare), 3),
      "remission", signif(median(w$frac_of_T.Remission), 3),
      "| exact paired p =", signif(paired_exact(w$frac_of_T.Flare, w$frac_of_T.Remission), 3), "\n")
  cat("TRBV11-2 share of TRBV-detected T nuclei: median flare", signif(median(w$frac_of_TRBV_detected.Flare), 3),
      "remission", signif(median(w$frac_of_TRBV_detected.Remission), 3),
      "| exact paired p =", signif(paired_exact(w$frac_of_TRBV_detected.Flare, w$frac_of_TRBV_detected.Remission), 3), "\n")
  write.table(b1, file.path(OUT, "B1_TRBV11-2_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # per T-cell subtype detection (flare vs remission pooled counts, descriptive)
  sub <- aggregate(cbind(n = 1, pos = as.numeric(pos[is_t]), any = as.numeric(det_trbv[is_t])) ~ ct[is_t] + md$cond[is_t], FUN = sum)
  names(sub)[1:2] <- c("cell_type", "cond")
  write.table(sub, file.path(OUT, "B1_TRBV11-2_by_celltype.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

## ---- B2: γδ T identity --------------------------------------------------------
gd_genes <- intersect(c("TRDC", "TRGC1", "TRGC2", "TRDV2", "TRGV9", "TRAC", "TRBC1", "TRBC2", "CD8A", "CD4"), rownames(rna))
det <- sapply(gd_genes, function(g) tapply(rna[g, ] > 0, ct, mean))
write.table(data.frame(cell_type = rownames(det), round(det, 4), check.names = FALSE),
            file.path(OUT, "B2_TR_gene_detection_by_celltype.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
loci <- GRanges(c("chr7", "chr14"), IRanges(c(38240024, 22462932), c(38368055, 22465787)),
                locus = c("TRG_locus", "TRDC"))
total_atac <- Matrix::colSums(atac)
b2 <- NULL
for (k in seq_along(loci)) {
  hit <- unique(queryHits(findOverlaps(peaks, loci[k])))
  cat(loci$locus[k], ": peaks overlapping =", length(hit), "\n")
  if (!length(hit)) next
  reads <- Matrix::colSums(atac[hit, , drop = FALSE])
  keep <- ct %in% c("gdT", "CD8 TEM", "CD8 Naive", "CD4 Naive", "CD4 TCM", "MAIT", "NK")
  agg <- aggregate(cbind(reads = reads[keep], total = total_atac[keep]) ~ ct[keep] + md$sample[keep], FUN = sum)
  names(agg)[1:2] <- c("cell_type", "sample")
  agg$cpm <- 1e6 * agg$reads / agg$total; agg$locus <- loci$locus[k]
  b2 <- rbind(b2, agg)
}
write.table(b2, file.path(OUT, "B2_TR_locus_accessibility_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(b2)) {
  for (L in unique(b2$locus)) {
    x <- b2[b2$locus == L, ]
    cat("\n", L, "median CPM by cell type:\n"); print(round(tapply(x$cpm, x$cell_type, median), 1))
    g <- x[x$cell_type == "gdT", c("sample", "cpm")]
    for (other in c("CD8 TEM", "MAIT", "CD4 Naive")) {
      o <- x[x$cell_type == other, c("sample", "cpm")]; m <- merge(g, o, by = "sample")
      cat("  gdT vs", other, ": higher in", sum(m$cpm.x > m$cpm.y), "of", nrow(m),
          "samples; exact paired p =", signif(paired_exact(m$cpm.x, m$cpm.y), 3), "\n")
    }
  }
}

## ---- B3: ATAC pseudobulk (input for b3_chromvar.R) --------------------------------
key <- factor(paste(ct, md$sample, sep = "|"))
M <- sparseMatrix(i = seq_along(key), j = as.integer(key), x = 1, dims = c(length(key), nlevels(key)),
                  dimnames = list(NULL, levels(key)))
apb <- atac %*% M
ncell <- as.integer(table(key)[colnames(apb)])
colinfo <- do.call(rbind, strsplit(colnames(apb), "|", fixed = TRUE))
cd <- DataFrame(cell_type = colinfo[, 1], sample = colinfo[, 2], n_cells = ncell)
cd$patient <- sub("_(Flare|Remission)$", "", cd$sample); cd$cond <- sub("^.*_", "", cd$sample)
saveRDS(list(counts = apb, info = as.data.frame(cd)), file.path(OUT, "atac_pseudobulk_counts.rds"))
# Motif deviation analysis on these pseudobulk counts: b3_chromvar.R
