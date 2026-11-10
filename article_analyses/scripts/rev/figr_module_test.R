# Module-level (TF programme) flare vs remission test on the joint FigR re-analysis, CD8 TEM and γδ T.
#
# Why: the edge-wise exact patient permutation (figr_differential_edges.R) has 2^10 relabellings,
# so its smallest two-sided p is 0.00195; with 2,594 / 1,534 edges, BH needs >= 101 / 60 edges at
# that floor before any edge can reach FDR < 0.05. A programme claim ("TBX21 drives the effector
# programme more strongly in flare") is one hypothesis and is tested as one.
#
# Steps 1-4 are copied unchanged from figr_differential_edges.R (same seed), so the motif-supported
# edge set E is the same; the job checks it against results/figr_differential.
#
# Modules, declared before the first run (23 Oct 2026; results/figr_module/MODULES.md):
#   CD8 TEM  TBX21_all        all motif-supported TBX21 pairs of the joint run
#            TBX21_published  the subset that is in Table S8
#            published_all    all motif-supported pairs that are in Table S8
#   γδ T     NFKB1, KLF3, FOSL1, JUNB, FOSL2, BACH2 (all motif-supported pairs of each TF)
#            published_all    all motif-supported pairs that are in Table S8
# Statistics per module:
#   (a) primary: mean over pairs of atanh(rho_flare) - atanh(rho_rem) on pooled nuclei; exact
#       two-sided p over all 2^P within-patient relabellings; BH over the 10 modules.
#   (b) per patient: mean over pairs of the per-library difference (libraries >= 30 nuclei);
#       patients with flare > remission; exact Wilcoxon signed-rank over patients.
#   (c) descriptive: median rho flare / remission, pairs with rho_flare > rho_rem.


suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix); library(SummarizedExperiment)
  library(GenomicRanges); library(FigR); library(FNN); library(irlba)
})
args <- commandArgs(trailingOnly = TRUE)
OBJ <- args[1]; OUT <- args[2]; S8 <- args[3]; NCORES <- as.integer(args[4]); PREV <- args[5]
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

res_all <- list(); pat_all <- list(); edges_all <- list()
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
  key <- cbind(match(E$Motif, tfE), match(E$DORC, dE))

  # check: same motif-supported edge set as the edge-wise run
  prev <- read.delim(file.path(PREV, paste0(gsub(" ", "_", CT), "_differential_edges.tsv")))
  same <- setequal(paste(prev$Motif, prev$DORC), paste(E$Motif, E$DORC))
  cat("motif-supported pairs:", nrow(E), "| identical to results/figr_differential:", same, "\n")

  pub <- paste(s8ct$motif, s8ct$gene)
  is_pub <- paste(E$Motif, E$DORC) %in% pub
  mods <- if (CT == "CD8 TEM") {
    list(TBX21_all = E$Motif == "TBX21", TBX21_published = E$Motif == "TBX21" & is_pub,
         published_all = is_pub)
  } else {
    c(setNames(lapply(c("NFKB1", "KLF3", "FOSL1", "JUNB", "FOSL2", "BACH2"), function(t) E$Motif == t),
               c("NFKB1", "KLF3", "FOSL1", "JUNB", "FOSL2", "BACH2")),
      list(published_all = is_pub))
  }
  mods <- mods[sapply(mods, sum) > 0]

  obs <- stat_for(is_flare)
  d_obs <- obs$d[key]
  m_obs <- sapply(mods, function(m) mean(d_obs[m]))
  flips <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), length(pts))))
  exceed <- numeric(length(mods))
  for (b in seq_len(nrow(flips))) {
    sw <- pts[flips[b, ]]
    mask <- ifelse(patient %in% sw, !is_flare, is_flare)
    db <- stat_for(mask)$d[key]
    mb <- sapply(mods, function(m) mean(db[m]))
    exceed <- exceed + (abs(mb) >= abs(m_obs) - 1e-12)
  }
  p_mod <- exceed / nrow(flips)

  # per-patient differences, libraries with >= 30 nuclei
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
  colnames(dpat) <- pp

  for (nm in names(mods)) {
    m <- mods[[nm]]
    pm <- colMeans(dpat[m, , drop = FALSE], na.rm = TRUE)
    wx <- suppressWarnings(wilcox.test(pm, exact = TRUE))
    res_all[[length(res_all) + 1]] <- data.frame(
      cell_type = CT, module = nm, n_pairs = sum(m),
      median_rho_flare = median(obs$rf[key][m]), median_rho_remission = median(obs$rr[key][m]),
      pairs_flare_higher = sum(d_obs[m] > 0), mean_delta_z = m_obs[[nm]],
      p_perm_module = p_mod[match(nm, names(mods))],
      patients_evaluable = length(pm), patients_flare_higher = sum(pm > 0),
      median_patient_delta_z = median(pm), p_wilcoxon_patients = wx$p.value)
    pat_all[[length(pat_all) + 1]] <- data.frame(cell_type = CT, module = nm, patient = names(pm),
                                                 mean_delta_z = pm)
  }
  edges_all[[CT]] <- data.frame(cell_type = CT, Motif = E$Motif, DORC = E$DORC, published = is_pub,
                                rho_flare = obs$rf[key], rho_remission = obs$rr[key], delta_z = d_obs,
                                patients_flare_higher = rowSums(dpat > 0, na.rm = TRUE),
                                patients_evaluable = rowSums(!is.na(dpat)))
  rm(sub, rna, atac, ATAC.se, dorc, dorc.s, rna.s); gc()
}

res <- do.call(rbind, res_all)
res$FDR_module <- p.adjust(res$p_perm_module, "BH")
write.table(res, file.path(OUT, "module_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, pat_all), file.path(OUT, "module_per_patient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, edges_all), file.path(OUT, "module_edges.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 220)
print(res, row.names = FALSE, digits = 3)
