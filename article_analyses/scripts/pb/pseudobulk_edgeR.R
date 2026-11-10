# Paired pseudobulk differential expression, MIS-C flare vs remission, per cell type.
#
# Unit of analysis = one patient's sample (not one cell). For each cell type:
#   1. sum raw RNA counts over all cells of that type in each library
#   2. keep libraries with >= MIN_CELLS cells; keep patients with both samples kept
#   3. edgeR quasi-likelihood GLM, design ~ patient + condition (paired)
#   4. BH over all expressed genes (genome-wide), not a curated panel
#   5. Hallmark gene sets with camera (competitive, accounts for inter-gene correlation)
# Then compare with the manuscript's two statistics (Table S6a: harmony cell-level
# FindMarkers adj.p and the per-patient paired Wilcoxon padj; Table S5: GSEA).
#
# Object: harmonized_seurat.rds (89,154 nuclei). Metadata: Batch = patient, Condition MISC = flare /
# Control = remission, orig.ident = library, predicted.id = cell type.

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(edgeR); library(limma); library(msigdbr)
})
args      <- commandArgs(trailingOnly = TRUE)
OBJ       <- args[1]
OUT       <- args[2]
TABLES    <- args[3]
MIN_CELLS <- 10   # per library and cell type
MIN_PAIRS <- 4
dir.create(file.path(OUT, "de"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "camera"), recursive = TRUE, showWarnings = FALSE)

agg_file <- file.path(OUT, "pseudobulk_counts.rds")
if (!file.exists(agg_file)) {
  cat("reading", OBJ, "\n"); t0 <- Sys.time()
  obj <- readRDS(OBJ)
  cat("read in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min,", ncol(obj), "cells\n")
  DefaultAssay(obj) <- "RNA"
  if (length(grep("^counts", Layers(obj[["RNA"]]))) > 1) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  cts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  md  <- obj@meta.data
  rm(obj); gc()
  key <- paste(md$predicted.id, md$orig.ident, sep = "|")
  f   <- factor(key)
  M   <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                      dims = c(length(f), nlevels(f)), dimnames = list(NULL, levels(f)))
  pb  <- as.matrix(cts %*% M)
  cols <- do.call(rbind, strsplit(colnames(pb), "|", fixed = TRUE))
  lib_info <- unique(data.frame(library = as.character(md$orig.ident),
                                patient = as.character(md$Batch),
                                condition = as.character(md$Condition)))
  info <- data.frame(column = colnames(pb), cell_type = cols[, 1], library = cols[, 2],
                     n_cells = as.integer(table(f)[colnames(pb)]))
  info <- merge(info, lib_info, by = "library")
  info <- info[match(colnames(pb), info$column), ]
  saveRDS(list(counts = pb, info = info), agg_file)
  rm(cts, M); gc()
}
A    <- readRDS(agg_file)
pb   <- A$counts
info <- A$info
info$condition <- factor(ifelse(info$condition == "MISC", "Flare", "Remission"),
                         levels = c("Remission", "Flare"))
write.table(info[, c("cell_type", "patient", "condition", "n_cells")],
            file.path(OUT, "cells_per_celltype_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Hallmark sets (msigdbr is bundled; no download)
hm <- msigdbr(species = "Homo sapiens", category = "H")
hm_sets <- split(hm$gene_symbol, hm$gs_name)

summ <- list(); allde <- list(); allcam <- list()
for (ct in sort(unique(info$cell_type))) {
  s <- info[info$cell_type == ct & info$n_cells >= MIN_CELLS, ]
  both <- names(which(table(s$patient) == 2))
  s <- s[s$patient %in% both, ]
  n_pairs <- length(both)
  if (n_pairs < MIN_PAIRS) {
    summ[[ct]] <- data.frame(cell_type = ct, n_pairs = n_pairs, cells_flare = NA, cells_remission = NA,
                             genes_tested = NA, n_fdr05 = NA, n_up = NA, n_down = NA,
                             note = paste0("skipped: < ", MIN_PAIRS, " pairs with >= ", MIN_CELLS, " cells"))
    next
  }
  s <- s[order(s$patient, s$condition), ]
  y <- DGEList(pb[, s$column], samples = s)
  patient   <- factor(s$patient)
  condition <- s$condition
  design <- model.matrix(~ patient + condition)
  keep <- filterByExpr(y, group = condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "conditionFlare")
  tt  <- topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt); tt$cell_type <- ct; tt$n_pairs <- n_pairs
  tt <- tt[, c("cell_type", "gene", "logFC", "logCPM", "F", "PValue", "FDR", "n_pairs")]
  write.table(tt, file.path(OUT, "de", paste0(gsub("[ /]", "_", ct), "_edgeR_QL.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  allde[[ct]] <- tt

  idx <- ids2indices(hm_sets, rownames(y))
  cam <- camera(y, idx, design, contrast = ncol(design))
  cam$set <- rownames(cam); cam$cell_type <- ct
  write.table(cam, file.path(OUT, "camera", paste0(gsub("[ /]", "_", ct), "_camera_hallmark.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  allcam[[ct]] <- cam

  summ[[ct]] <- data.frame(cell_type = ct, n_pairs = n_pairs,
                           cells_flare = sum(s$n_cells[s$condition == "Flare"]),
                           cells_remission = sum(s$n_cells[s$condition == "Remission"]),
                           genes_tested = nrow(tt), n_fdr05 = sum(tt$FDR < 0.05),
                           n_up = sum(tt$FDR < 0.05 & tt$logFC > 0),
                           n_down = sum(tt$FDR < 0.05 & tt$logFC < 0), note = "")
  cat(ct, ": pairs", n_pairs, "| genes", nrow(tt), "| FDR<0.05", sum(tt$FDR < 0.05), "\n")
}
summ <- do.call(rbind, summ)
write.table(summ, file.path(OUT, "summary_per_celltype.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
de  <- do.call(rbind, allde)
cam <- do.call(rbind, allcam)
write.table(cam, file.path(OUT, "camera_hallmark_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- comparison with the manuscript (Table S6a) ----
s6 <- read.delim(file.path(TABLES, "AF2__12_Table_S6a.tsv"), skip = 1, check.names = FALSE,
                 stringsAsFactors = FALSE)
s6$cell_type_pb <- gsub("_", " ", s6$cell_type)
s6$cell_type_pb[s6$cell_type == "gdT"] <- "gdT"
cmp <- merge(s6[, c("gene", "cell_type_pb", "claim", "n_patients", "harmony_avg_log2FC",
                    "harmony_p_val_adj", "median_log2FC", "paired_p", "paired_padj")],
             de[, c("cell_type", "gene", "logFC", "PValue", "FDR", "n_pairs")],
             by.x = c("gene", "cell_type_pb"), by.y = c("gene", "cell_type"), all.x = TRUE)
names(cmp)[names(cmp) == "cell_type_pb"] <- "cell_type"
names(cmp)[names(cmp) == "logFC"] <- "pb_logFC"
names(cmp)[names(cmp) == "PValue"] <- "pb_p"
names(cmp)[names(cmp) == "FDR"] <- "pb_FDR_genomewide"
cmp$harmony_sig <- !is.na(cmp$harmony_p_val_adj) & cmp$harmony_p_val_adj < 0.05
cmp$paired_sig  <- !is.na(cmp$paired_padj) & cmp$paired_padj < 0.05
cmp$pb_sig      <- !is.na(cmp$pb_FDR_genomewide) & cmp$pb_FDR_genomewide < 0.05
cmp$same_direction_pb_vs_harmony <- sign(cmp$pb_logFC) == sign(cmp$harmony_avg_log2FC)
write.table(cmp, file.path(OUT, "comparison_S6_vs_pseudobulk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n== Manuscript key-gene rows (Table S6a):", nrow(cmp), "| not testable in pseudobulk:",
    sum(is.na(cmp$pb_FDR_genomewide)), "\n")
print(with(cmp, table(harmony = harmony_sig, pseudobulk = pb_sig)))
print(with(cmp, table(paired_wilcoxon = paired_sig, pseudobulk = pb_sig)))
cat("direction agreement pb vs harmony:",
    round(mean(cmp$same_direction_pb_vs_harmony, na.rm = TRUE), 3), "\n")

focus <- c("IFNG", "IFNG-AS1", "STAT4", "STAT1", "BATF", "PRDM1", "TBX21", "IL12RB1", "IL12RB2",
           "PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT", "GZMB", "CXCL10")
f <- cmp[cmp$gene %in% focus & cmp$cell_type %in% c("CD8 TEM", "NK", "gdT"), ]
f <- f[order(f$gene, f$cell_type), c("gene", "cell_type", "harmony_avg_log2FC", "harmony_p_val_adj",
                                     "n_patients", "paired_padj", "pb_logFC", "pb_FDR_genomewide", "n_pairs")]
options(width = 220)
cat("\n== Focus genes in CD8 TEM / NK / gdT\n"); print(f, row.names = FALSE, digits = 3)

cat("\n== Hallmark IFN-gamma and TGF-beta (camera, FDR within cell type)\n")
print(cam[cam$set %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TGF_BETA_SIGNALING"),
          c("cell_type", "set", "NGenes", "Direction", "PValue", "FDR")], row.names = FALSE, digits = 3)
cat("\n"); print(summ, row.names = FALSE)
