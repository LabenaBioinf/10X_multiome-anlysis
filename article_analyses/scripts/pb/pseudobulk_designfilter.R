# Sensitivity analysis: paired pseudobulk edgeR QL as in pseudobulk_edgeR.R, but expression filtering with the
# full design (filterByExpr(y, design = design)) instead of the condition groups. Writes key-gene results for the
# populations shown in Figure 4B and camera Hallmark results for all testable populations.
# usage: Rscript pseudobulk_designfilter.R <pseudobulk_counts.rds> <out dir>
suppressPackageStartupMessages({ library(edgeR); library(limma); library(msigdbr) })
args <- commandArgs(trailingOnly = TRUE)
A <- readRDS(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MIN_CELLS <- 10; MIN_PAIRS <- 4
pb <- A$counts; info <- A$info
info$condition <- factor(ifelse(info$condition == "MISC", "Flare", "Remission"), levels = c("Remission", "Flare"))
hm <- msigdbr(species = "Homo sapiens", category = "H")
hm_sets <- split(hm$gene_symbol, hm$gs_name)
genes <- c("STAT1", "PRDM1", "BATF", "IRF4", "IL6ST", "STAT3", "STAT4", "TBX21", "IFNG", "IFNG-AS1", "GZMB", "PRF1",
           "NKG7", "CTLA4", "LAG3", "HAVCR2", "PDCD1", "TOX", "NR4A2")
key <- list(); cams <- list()
for (ct in sort(unique(info$cell_type))) {
  s <- info[info$cell_type == ct & info$n_cells >= MIN_CELLS, ]
  both <- names(which(table(s$patient) == 2)); s <- s[s$patient %in% both, ]
  if (length(both) < MIN_PAIRS) next
  s <- s[order(s$patient, s$condition), ]
  y <- DGEList(pb[, s$column], samples = s)
  design <- model.matrix(~ factor(s$patient) + s$condition)
  keep <- filterByExpr(y, design = design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, coef = ncol(design)), n = Inf)$table
  g <- intersect(genes, rownames(tt))
  key[[ct]] <- data.frame(cell_type = ct, gene = g, logFC = tt[g, "logFC"], PValue = tt[g, "PValue"],
                          FDR_genomewide = tt[g, "FDR"], n_pairs = length(both), genes_tested = nrow(tt))
  cam <- camera(y, ids2indices(hm_sets, rownames(y)), design, contrast = ncol(design))
  cam$set <- rownames(cam); cam$cell_type <- ct; cams[[ct]] <- cam
}
write.table(do.call(rbind, key), file.path(OUT, "key_genes_designfilter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, cams), file.path(OUT, "camera_hallmark_designfilter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("populations:", length(key), "\n")
