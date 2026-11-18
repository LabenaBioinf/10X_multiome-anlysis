# Article (final version): only log2FC, p value and genome-wide FDR are used. The panel FDR columns are still written
# but no longer used (genes below the primary filter are reported from pseudobulk_designfilter_panel.R instead).
# Corrected key-gene and targeted gene-set analysis (article version, v2.1 of this repository):
#  - Primary statistics (logFC, PValue, genome-wide FDR) come from the standard-filtered fit (same model as
#    pseudobulk_edgeR.R, filterByExpr(group = condition)); forced low-expression genes have no genome-wide FDR.
#  - The panel model (standard + low-expression panel genes with at least five reads) is fitted separately and
#    used only for the post hoc panel FDR (BH over all 120 panel genes) and for the estimates of forced genes.
#  - Targeted gene sets are de-duplicated and tested with camera on the standard-filtered model.
# Output columns of key_gene_panel_pseudobulk.tsv: logFC/PValue = primary estimate if the gene passes the
# filter, otherwise the panel-model estimate (estimate_source says which; not_estimable = fewer than five reads);
# *_primary and *_panel_model give both; programme = category of the key-gene table, or additional genes.
# Follow-up analyses on the paired pseudobulk counts (flare vs remission, per cell type).
#
#  1. Key-gene panel (all genes of Table S6 plus exhaustion / IL-27 / IL-12 genes):
#     edgeR QL with the panel genes forced into the model even when they fail
#     filterByExpr if they have at least five reads in total (others are reported as not estimable).
#     Genome-wide FDR is taken from the standard (filtered) fit. Panel FDR: BH with n = all 120 panel genes
#     per population (fixed family; not-estimable genes count as untested hypotheses).
#  2. Gene-set tests (camera) for targeted signatures chosen after the initial analysis (post hoc):
#       IL-12 / STAT4 : PID_IL12_STAT4_PATHWAY, PID_IL12_2PATHWAY,
#                       IL-12-induced genes in activated CD8 T cells (GSE15930 24/48/72 h,
#                       GSE13173), IL-12-treated PBMC (GSE12839),
#                       STAT4-dependent genes in CD8 T cells (GSE40666 WT vs STAT4 KO)
#       IL-27         : PID_IL27_PATHWAY, IL-27-induced genes (GSE44732 DC, GSE44955 macrophage)
#       gp130 / STAT3 : HALLMARK_IL6_JAK_STAT3_SIGNALING
#       reference     : HALLMARK_INTERFERON_GAMMA_RESPONSE, HALLMARK_TGF_BETA_SIGNALING
#     MSigDB C7 naming: "A_VS_B_DN" = genes higher in B. Signature genes are therefore
#     taken from the side that is higher after cytokine exposure.
#
# Input : results/pseudobulk/pseudobulk_counts.rds (from pseudobulk_edgeR.R)
# Output: results/pseudobulk_followup/ (key_gene_panel_pseudobulk.tsv, signature_tests_camera.tsv, gene_sets_used.tsv)

suppressPackageStartupMessages({ library(edgeR); library(limma); library(msigdbr) })
args   <- commandArgs(trailingOnly = TRUE)
PB     <- args[1]
OUT    <- args[2]
TABLES <- args[3]
MIN_CELLS <- 10
MIN_PAIRS <- 4
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

A <- readRDS(PB); pb <- A$counts; info <- A$info
info$condition <- factor(ifelse(info$condition == "MISC", "Flare", "Remission"), levels = c("Remission", "Flare"))

s6 <- read.delim(file.path(TABLES, "AF2__12_Table_S6a.tsv"), skip = 1, check.names = FALSE)
extra <- c("TOX", "TOX2", "PDCD1", "TCF7", "EOMES", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1",
           "CXCL13", "NR4A1", "NR4A2", "NR4A3", "IL27", "EBI3", "IL27RA", "IL6ST", "IL6R", "IL12A",
           "IL12B", "IL12RB1", "IL12RB2", "STAT4", "TBX21", "IFNG", "IFNG-AS1", "TRBV11-2", "TRDC",
           "TRGC1", "TRGC2", "IL2RA", "IL2RB", "IL2RG", "IL15RA", "IL10", "SOCS1", "SOCS3")
panel <- sort(unique(c(s6$gene, extra)))
panel <- panel[panel %in% rownames(pb)]
# programme of each panel gene: category of the key-gene table, or one of the ten genes added to the panel in the second step
prog <- tapply(as.character(s6$claim), s6$gene, function(x) x[1])
prog <- setNames(ifelse(panel %in% names(prog), prog[panel], "additional genes"), panel)
cat("panel genes present:", length(panel), "\n")

msig <- function(cat, sub = NULL) {
  m <- if (is.null(sub)) msigdbr("Homo sapiens", category = cat) else msigdbr("Homo sapiens", category = cat, subcategory = sub)
  split(m$gene_symbol, m$gs_name)
}
h  <- msig("H"); c2 <- msig("C2"); c7 <- msig("C7")
sets <- list(
  IL12_STAT4_PID           = c2[["PID_IL12_STAT4_PATHWAY"]],
  IL12_PID                 = c2[["PID_IL12_2PATHWAY"]],
  IL12_induced_CD8_24h     = c7[["GSE15930_STIM_VS_STIM_AND_IL12_24H_CD8_T_CELL_DN"]],
  IL12_induced_CD8_48h     = c7[["GSE15930_STIM_VS_STIM_AND_IL12_48H_CD8_T_CELL_DN"]],
  IL12_induced_CD8_72h     = c7[["GSE15930_STIM_VS_STIM_AND_IL12_72H_CD8_T_CELL_DN"]],
  IL12_induced_actCD8      = c7[["GSE13173_UNTREATED_VS_IL12_TREATED_ACT_CD8_TCELL_DN"]],
  IL12_induced_PBMC        = c7[["GSE12839_CTRL_VS_IL12_TREATED_PBMC_DN"]],
  STAT4_dependent_CD8      = c7[["GSE40666_WT_VS_STAT4_KO_CD8_TCELL_UP"]],
  IL27_PID                 = c2[["PID_IL27_PATHWAY"]],
  IL27_induced_DC          = c7[["GSE44732_UNSTIM_VS_IL27_STIM_IMATURE_DC_DN"]],
  IL27_induced_macrophage  = c7[["GSE44955_MCSF_VS_MCSF_AND_IL27_STIM_MACROPHAGE_DN"]],
  IL6_JAK_STAT3_hallmark   = h[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]],
  IFNG_response_hallmark   = h[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
  TGFB_signalling_hallmark = h[["HALLMARK_TGF_BETA_SIGNALING"]]
)
sets <- lapply(sets, function(g) unique(g[!is.na(g) & nzchar(g)]))
stopifnot(all(lengths(sets) > 0))
write.table(data.frame(set = names(sets), n_unique_genes = lengths(sets),
                       genes = sapply(sets, paste, collapse = ",")),
            file.path(OUT, "gene_sets_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

panel_res <- list(); cam_res <- list()
for (ct in sort(unique(info$cell_type))) {
  s <- info[info$cell_type == ct & info$n_cells >= MIN_CELLS, ]
  both <- names(which(table(s$patient) == 2)); s <- s[s$patient %in% both, ]
  if (length(both) < MIN_PAIRS) next
  s <- s[order(s$patient, s$condition), ]
  patient <- factor(s$patient); condition <- s$condition
  design <- model.matrix(~ patient + condition)
  y0 <- DGEList(pb[, s$column])
  keep <- filterByExpr(y0, group = condition)
  forced <- !keep & rownames(y0) %in% panel & rowSums(y0$counts) >= 5
  # primary model: standard filter only
  ys <- calcNormFactors(y0[keep, , keep.lib.sizes = FALSE])
  ys <- estimateDisp(ys, design, robust = TRUE)
  tp <- topTags(glmQLFTest(glmQLFit(ys, design, robust = TRUE), coef = "conditionFlare"), n = Inf, sort.by = "none")$table
  # panel model: standard + forced low-expression panel genes (post hoc panel only)
  y <- calcNormFactors(y0[keep | forced, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  tt <- topTags(glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = "conditionFlare"), n = Inf, sort.by = "none")$table
  pr <- tt[rownames(tt) %in% panel, ]
  pr$FDR_panel <- p.adjust(pr$PValue, "BH", n = length(panel))   # fixed family: all panel genes
  g <- rownames(pr); inprim <- g %in% rownames(tp)
  lfc_p <- ifelse(inprim, tp[g, "logFC"], NA); p_p <- ifelse(inprim, tp[g, "PValue"], NA); fdr_gw <- ifelse(inprim, tp[g, "FDR"], NA)
  cpm_f <- rowMeans(cpm(y)[g, condition == "Flare", drop = FALSE])
  cpm_r <- rowMeans(cpm(y)[g, condition == "Remission", drop = FALSE])
  panel_res[[ct]] <- data.frame(cell_type = ct, gene = g, n_pairs = length(both),
                                logFC = ifelse(inprim, lfc_p, pr$logFC), mean_CPM_flare = cpm_f, mean_CPM_remission = cpm_r,
                                PValue = ifelse(inprim, p_p, pr$PValue), FDR_genomewide = fdr_gw, FDR_panel = pr$FDR_panel,
                                low_expression_forced = g %in% rownames(y0)[forced],
                                estimate_source = ifelse(inprim, "primary", "panel_model"),
                                logFC_primary = lfc_p, PValue_primary = p_p,
                                logFC_panel_model = pr$logFC, PValue_panel_model = pr$PValue)
  miss <- setdiff(panel, g)   # panel genes not estimable in this population (fewer than five reads in total)
  if (length(miss)) {
    na <- panel_res[[ct]][rep(1, length(miss)), ]; na[] <- NA
    na$cell_type <- ct; na$gene <- miss; na$n_pairs <- length(both); na$low_expression_forced <- FALSE
    na$estimate_source <- "not_estimable"
    panel_res[[ct]] <- rbind(panel_res[[ct]], na)
  }
  panel_res[[ct]]$programme <- prog[panel_res[[ct]]$gene]
  idx <- ids2indices(sets, rownames(ys))
  cam <- camera(ys, idx, design, contrast = ncol(design), sort = FALSE)
  cam$set <- rownames(cam); cam$cell_type <- ct; cam$n_pairs <- length(both)
  cam_res[[ct]] <- cam
}
panel_res <- do.call(rbind, panel_res)
cam_res <- do.call(rbind, cam_res)
cam_res$FDR_all_tests <- p.adjust(cam_res$PValue, "BH")
write.table(panel_res, file.path(OUT, "key_gene_panel_pseudobulk.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cam_res[, c("cell_type", "n_pairs", "set", "NGenes", "Direction", "PValue", "FDR", "FDR_all_tests")],
            file.path(OUT, "signature_tests_camera.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

options(width = 220)
cat("\n== Signatures (camera; FDR within cell type over 14 sets) in effector and naive populations\n")
show <- c("CD8 TEM", "NK", "gdT", "CD4 TCM", "CD4 Naive", "CD8 Naive", "CD14 Mono", "CD16 Mono", "NK_CD56bright")
print(cam_res[cam_res$cell_type %in% show, c("cell_type", "set", "NGenes", "Direction", "PValue", "FDR")],
      row.names = FALSE, digits = 3)
cat("\n== Exhaustion / checkpoint panel (B8)\n")
ex <- c("TOX", "PDCD1", "TCF7", "EOMES", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1", "NR4A2", "NR4A3")
print(panel_res[panel_res$gene %in% ex & panel_res$cell_type %in% c("CD8 TEM", "NK", "gdT", "CD8 TCM"),
                c("cell_type", "gene", "logFC", "mean_CPM_flare", "mean_CPM_remission", "PValue", "FDR_panel", "low_expression_forced")],
      row.names = FALSE, digits = 3)
cat("\n== IL-27 / gp130 / IL-12 genes\n")
g2 <- c("IL27", "EBI3", "IL27RA", "IL6ST", "IL12RB1", "IL12RB2", "STAT4", "STAT3", "TBX21", "BATF", "IRF4")
print(panel_res[panel_res$gene %in% g2 & panel_res$cell_type %in% c("CD8 TEM", "NK", "gdT", "CD14 Mono"),
                c("cell_type", "gene", "logFC", "PValue", "FDR_panel", "low_expression_forced")],
      row.names = FALSE, digits = 3)
