# Sensitivity analysis over the full targeted panel (Table S6): paired pseudobulk edgeR QL with expression filtering
# by the full design (filterByExpr(y, design = design)), as in pseudobulk_designfilter.R, but reporting every panel
# gene in every tested population. Genes removed by the design-aware filter are listed with passed_design_filter = FALSE.
# usage: Rscript pseudobulk_designfilter_panel.R <pseudobulk_counts.rds> <key_gene_panel_pseudobulk.tsv> <out dir>
suppressPackageStartupMessages({ library(edgeR); library(limma) })
args <- commandArgs(trailingOnly = TRUE)
A <- readRDS(args[1]); panel <- read.delim(args[2], check.names = FALSE); OUT <- args[3]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MIN_CELLS <- 10; MIN_PAIRS <- 4
pb <- A$counts; info <- A$info
info$condition <- factor(ifelse(info$condition == "MISC", "Flare", "Remission"), levels = c("Remission", "Flare"))
genes <- sort(unique(panel$gene))
key <- list()
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
  g <- genes[genes %in% rownames(pb)]
  ok <- g %in% rownames(tt)
  key[[ct]] <- data.frame(cell_type = ct, gene = g, passed_design_filter = ok,
                          logFC = ifelse(ok, tt[g, "logFC"], NA), PValue = ifelse(ok, tt[g, "PValue"], NA),
                          FDR_genomewide = ifelse(ok, tt[g, "FDR"], NA), n_pairs = length(both), genes_tested = nrow(tt))
}
res <- do.call(rbind, key)
write.table(res, file.path(OUT, "key_genes_designfilter_panel.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("populations:", length(key), " panel genes:", length(genes), " rows:", nrow(res),
    " passed:", sum(res$passed_design_filter), "\n")
