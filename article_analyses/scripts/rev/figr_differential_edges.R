# Differential TF–DORC association, flare vs remission, in CD8 TEM and γδ T nuclei.
#
# Replaces the "significant in flare and not significant in remission" logic used for
# Table S8 with a direct test of the difference.
#
#   1. Nuclei of one cell type from both conditions are analysed together, so every
#      TF–DORC pair is defined once and tested in both conditions.
#   2. ATAC LSI (TF-IDF + SVD) on the cell type; kNN smoothing (k <= 30) is done within
#      each library, so neighbours never mix patients or conditions.
#   3. Peak–gene correlations (FigR::runGenePeakcorr, n_bg = 100, pvalZ <= 0.05) for the
#      DORC genes of the published panels (Table S8 gene sets of that cell type).
#   4. Motif enrichment of each TF in each DORC's peaks (FigR::runFigRGRN on all nuclei,
#      all expressed TFs with a cisBP motif). Edges are carried forward only if the TF
#      motif is enriched (Enrichment.P < 0.05).
#   5. For each edge: Spearman rho between smoothed TF RNA and smoothed DORC score in
#      flare and in remission nuclei; statistic = atanh(rho_flare) - atanh(rho_rem).
#      Null distribution: all 2^P relabellings that swap flare/remission within patients
#      (exact; P = number of patients). Two-sided p, BH over edges per cell type.
#   6. Patient-level consistency: rho per library (>= 30 nuclei), sign of the per-patient
#      difference.
#   7. Comparison with the published flare-exclusive edges (Table S8).

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix); library(SummarizedExperiment)
  library(GenomicRanges); library(FigR); library(FNN); library(irlba)
})
args <- commandArgs(trailingOnly = TRUE)
OBJ <- args[1]; OUT <- args[2]; S8 <- args[3]; NCORES <- as.integer(args[4])
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261016)

s8 <- read.delim(S8, skip = 1, stringsAsFactors = FALSE)
obj <- readRDS(OBJ)
md0 <- obj@meta.data

tf_idf_lsi <- function(m, dims = 2:30) {
  m <- m[Matrix::rowSums(m) > 0, ]
  tf <- t(t(m) / Matrix::colSums(m))
  idf <- log(1 + ncol(m) / Matrix::rowSums(m > 0))
  x <- tf * idf
  x@x <- log1p(x@x * 1e4)
  sv <- irlba(x, nv = max(dims))
  emb <- sv$v %*% diag(sv$d)
  emb[, dims]
}
smooth_within <- function(emb, lib, k = 30) {
  n <- nrow(emb); ii <- integer(0); jj <- integer(0); ww <- numeric(0)
  for (L in unique(lib)) {
    idx <- which(lib == L); kk <- min(k, length(idx) - 1)
    nn <- if (kk >= 1) get.knn(emb[idx, , drop = FALSE], k = kk)$nn.index else matrix(integer(0), length(idx), 0)
    for (a in seq_along(idx)) {
      nb <- c(idx[a], idx[nn[a, ]])
      ii <- c(ii, rep(idx[a], length(nb))); jj <- c(jj, nb); ww <- c(ww, rep(1 / length(nb), length(nb)))
    }
  }
  sparseMatrix(i = ii, j = jj, x = ww, dims = c(n, n))
}
rank_rows <- function(m) t(apply(m, 1, rank))
spearman_xy <- function(rx, ry) {
  # rx, ry: row-ranked matrices over the same cells; returns cor(rx_i, ry_j)
  zx <- rx - rowMeans(rx); zx <- zx / sqrt(rowSums(zx^2))
  zy <- ry - rowMeans(ry); zy <- zy / sqrt(rowSums(zy^2))
  zx %*% t(zy)
}
fz <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))

summary_all <- list()
for (CT in c("CD8 TEM", "gdT")) {
  s8ct <- s8[s8$cell_type == sub(" ", "_", CT), ]
  panel_genes <- unique(s8ct$gene)
  cells <- rownames(md0)[md0$predicted.id == CT]
  sub <- subset(obj, cells = cells)
  md <- sub@meta.data
  lib <- as.character(md$orig.ident)
  patient <- as.character(md$Batch)
  cond <- ifelse(md$Condition == "MISC", "Flare", "Remission")
  DefaultAssay(sub) <- "RNA"
  if (length(grep("^counts", Layers(sub[["RNA"]]))) > 1) sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
  sub <- NormalizeData(sub, verbose = FALSE)
  rna <- GetAssayData(sub, assay = "RNA", layer = "data")
  atac <- GetAssayData(sub, assay = "ATAC", layer = "counts")
  pk <- granges(sub[["ATAC"]])
  std_chr <- paste0("chr", c(1:22, "X", "Y"))
  keep_pk <- Matrix::rowSums(atac > 0) >= 10 & as.character(seqnames(pk)) %in% std_chr
  atac <- atac[keep_pk, ]; pk <- pk[keep_pk]
  seqlevels(pk) <- std_chr
  cat("\n=====", CT, ": nuclei", ncol(atac), "| peaks", nrow(atac), "| panel genes", length(panel_genes), "\n")

  emb <- tf_idf_lsi(atac)
  W <- smooth_within(emb, lib, k = 30)

  ATAC.se <- SummarizedExperiment(assays = list(counts = atac), rowRanges = pk)
  genes_in <- intersect(panel_genes, rownames(rna))
  cis <- runGenePeakcorr(ATAC.se = ATAC.se, RNAmat = rna, genome = "hg38", geneList = genes_in,
                         nCores = NCORES, n_bg = 100, p.cut = NULL)
  cis.f <- cis[cis$pvalZ <= 0.05, ]
  dorc_genes <- unique(cis.f$Gene)
  cat("panel genes with >=1 significant peak (DORCs):", length(dorc_genes), "of", length(genes_in), "\n")
  dorc <- getDORCScores(ATAC.se, dorcTab = cis.f, geneList = dorc_genes, nCores = NCORES)
  dorc[is.na(dorc)] <- 0
  dorc.s <- as.matrix(dorc %*% t(W))
  colnames(dorc.s) <- colnames(atac)

  pfm <- readRDS(file.path(find.package("FigR"), "data", "cisBP_human_pfms_2021.rds"))
  if (all(grepl("_", names(pfm), fixed = TRUE))) names(pfm) <- FigR::extractTFNames(names(pfm))
  tfs <- intersect(names(pfm), rownames(rna))
  tfs <- tfs[Matrix::rowSums(rna[tfs, ] > 0) >= 0.01 * ncol(rna)]
  rna.s <- as.matrix(rna[union(tfs, dorc_genes), ] %*% t(W))
  colnames(rna.s) <- colnames(atac)

  grn <- runFigRGRN(ATAC.se = ATAC.se, dorcTab = cis.f, genome = "hg38", dorcMat = dorc.s,
                    rnaMat = rna.s, dorcGenes = dorc_genes, nCores = NCORES)
  grn <- grn[grn$Motif %in% tfs, ]
  cat("TF motifs tested:", length(unique(grn$Motif)), "| DORCs:", length(unique(grn$DORC)),
      "| edges:", nrow(grn), "| motif-enriched edges:", sum(grn$Enrichment.P < 0.05), "\n")
  write.table(grn, file.path(OUT, paste0(gsub(" ", "_", CT), "_figR_pooled.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  E <- grn[grn$Enrichment.P < 0.05, c("Motif", "DORC", "Enrichment.P", "Corr", "Score")]
  tfE <- unique(E$Motif); dE <- unique(E$DORC)
  X <- rna.s[tfE, , drop = FALSE]; Y <- dorc.s[dE, , drop = FALSE]
  pts <- sort(unique(patient))
  stat_for <- function(flare_mask) {
    rf <- spearman_xy(rank_rows(X[, flare_mask, drop = FALSE]), rank_rows(Y[, flare_mask, drop = FALSE]))
    rr <- spearman_xy(rank_rows(X[, !flare_mask, drop = FALSE]), rank_rows(Y[, !flare_mask, drop = FALSE]))
    list(rf = rf, rr = rr, d = fz(rf) - fz(rr))
  }
  is_flare <- cond == "Flare"
  obs <- stat_for(is_flare)
  key <- cbind(match(E$Motif, tfE), match(E$DORC, dE))
  d_obs <- obs$d[key]
  # exact within-patient label swaps (2^P); flips[p] = TRUE swaps that patient's labels
  flips <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(pts))))
  exceed <- numeric(nrow(E))
  for (b in seq_len(nrow(flips))) {
    sw <- pts[flips[b, ]]
    mask <- ifelse(patient %in% sw, !is_flare, is_flare)
    db <- stat_for(mask)$d[key]
    exceed <- exceed + (abs(db) >= abs(d_obs) - 1e-12)
  }
  p_perm <- exceed / nrow(flips)

  # patient-level consistency: rho per library with >= 30 nuclei
  libs_ok <- names(which(table(lib) >= 30))
  per_lib <- lapply(libs_ok, function(L) {
    m <- lib == L
    fz(spearman_xy(rank_rows(X[, m, drop = FALSE]), rank_rows(Y[, m, drop = FALSE])))[key]
  })
  names(per_lib) <- libs_ok
  lib_info <- unique(data.frame(lib = lib, patient = patient, cond = cond))
  lib_info <- lib_info[lib_info$lib %in% libs_ok, ]
  pp <- names(which(table(lib_info$patient) == 2))
  dpat <- sapply(pp, function(p) {
    lf <- lib_info$lib[lib_info$patient == p & lib_info$cond == "Flare"]
    lr <- lib_info$lib[lib_info$patient == p & lib_info$cond == "Remission"]
    per_lib[[lf]] - per_lib[[lr]]
  })
  if (is.null(dim(dpat))) dpat <- matrix(dpat, ncol = length(pp))
  same_sign <- rowSums(sign(dpat) == sign(d_obs), na.rm = TRUE)

  res <- data.frame(cell_type = CT, Motif = E$Motif, DORC = E$DORC, Enrichment.P = E$Enrichment.P,
                    rho_flare = obs$rf[key], rho_remission = obs$rr[key], delta_z = d_obs,
                    p_perm = p_perm, FDR = p.adjust(p_perm, "BH"),
                    patients_evaluable = length(pp), patients_same_direction = same_sign)
  pub <- s8ct[, c("motif", "gene", "corr", "gene_set")]
  pub$published <- TRUE
  res <- merge(res, pub, by.x = c("Motif", "DORC"), by.y = c("motif", "gene"), all.x = TRUE)
  res$published[is.na(res$published)] <- FALSE
  write.table(res, file.path(OUT, paste0(gsub(" ", "_", CT), "_differential_edges.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  # status of every published S8 edge
  st <- merge(pub, grn[, c("Motif", "DORC", "Enrichment.P", "Corr")], by.x = c("motif", "gene"),
              by.y = c("Motif", "DORC"), all.x = TRUE)
  st <- merge(st, res[, c("Motif", "DORC", "delta_z", "p_perm", "FDR")], by.x = c("motif", "gene"),
              by.y = c("Motif", "DORC"), all.x = TRUE)
  st$status <- ifelse(is.na(st$Enrichment.P), "not testable (TF or DORC absent in re-analysis)",
               ifelse(st$Enrichment.P >= 0.05, "no motif enrichment",
               ifelse(st$FDR < 0.05 & sign(st$delta_z) == sign(st$corr), "flare-specific difference (FDR<0.05)",
               ifelse(st$p_perm < 0.05, "nominal difference only (p<0.05, FDR>=0.05)", "no flare-remission difference"))))
  write.table(st, file.path(OUT, paste0(gsub(" ", "_", CT), "_published_S8_edge_status.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nPublished S8 edges by status:\n"); print(table(st$status))
  cat("Differential edges (FDR < 0.05) among motif-enriched edges:", sum(res$FDR < 0.05), "of", nrow(res),
      "| minimum attainable p =", signif(2 / nrow(flips), 3), "\n")
  top <- res[order(res$p_perm, -abs(res$delta_z)), ][1:min(25, nrow(res)), ]
  print(top[, c("Motif", "DORC", "rho_flare", "rho_remission", "delta_z", "p_perm", "FDR",
                "patients_same_direction", "patients_evaluable", "published")], row.names = FALSE, digits = 3)
  focus <- res[res$Motif %in% c("TBX21", "STAT1", "BATF", "PRDM1", "STAT4", "EOMES", "FOSB", "NFIL3") &
               res$DORC %in% c("IFNG", "IFNG-AS1", "GZMB", "STAT1", "STAT4", "CTLA4", "HAVCR2", "IL6ST", "TGFB1", "SMAD7", "BHLHE40"), ]
  cat("\nFocus edges named in the manuscript:\n")
  print(focus[, c("Motif", "DORC", "Enrichment.P", "rho_flare", "rho_remission", "p_perm", "FDR")], row.names = FALSE, digits = 3)
  summary_all[[CT]] <- table(st$status)
  rm(sub, rna, atac, ATAC.se, dorc, dorc.s, rna.s); gc()
}
