# chromVAR sequencing-depth sensitivity analysis (Supplementary Table S8g).
#
# Same input, peaks, peak filter, motifs, background peaks, populations, pairing, model (limma ~ patient + condition) and
# test families (BH over motifs within population, BH over families within population) as chromvar_naive_extension.R.
# Only the depth handling changes:
#   A  z        deviationScores(dev)         = published analysis (must reproduce results/chromvar_naive_extension)
#   B  dev      deviations(dev)              bias-corrected accessibility deviations (not scaled by their SD)
#   C  down_z / down_dev, seeds 1-5         within each population and patient, the deeper of the flare and remission
#                                            pseudobulks is binomially thinned per peak to the fragment total of the
#                                            shallower one, then chromVAR is rerun with the background peaks and the
#                                            expectation of the original data, so only depth differs
#   D  depth report                          fragments and nuclei per pseudobulk, paired log2 depth ratio, Spearman of
#                                            STAT1/STAT4 z with log fragments, and of the paired change with the depth ratio
# Equal numbers of nuclei per phase cannot be drawn here: the input is aggregated per population and sample.
#
# usage: Rscript chromvar_depth_sensitivity.R <results/object_level> <out dir>
suppressPackageStartupMessages({
  library(Matrix); library(GenomicRanges); library(SummarizedExperiment); library(chromVAR)
  library(motifmatchr); library(BSgenome.Hsapiens.UCSC.hg38); library(BiocParallel); library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
IN <- args[1]; OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
register(MulticoreParam(8))
set.seed(20261016)                                   # as the published run, so the background peaks are identical
A <- readRDS(file.path(IN, "atac_pseudobulk_counts.rds"))
apb <- A$counts; cd <- DataFrame(A$info)
pp <- do.call(rbind, strsplit(rownames(apb), "-", fixed = TRUE))
peaks <- GRanges(pp[, 1], IRanges(as.integer(pp[, 2]), as.integer(pp[, 3])))
std <- as.character(seqnames(peaks)) %in% paste0("chr", c(1:22, "X", "Y"))
apb <- apb[std, ]; peaks <- peaks[std]
seqlevels(peaks) <- paste0("chr", c(1:22, "X", "Y"))
use <- cd$n_cells >= 20
se <- SummarizedExperiment(assays = list(counts = apb[, use]), rowRanges = peaks, colData = cd[use, ])
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
se <- filterPeaks(se, min_fragments_per_peak = 10, non_overlapping = TRUE)
pfm <- readRDS(file.path(find.package("FigR"), "data", "cisBP_human_pfms_2021.rds"))
if (all(grepl("_", names(pfm), fixed = TRUE))) names(pfm) <- FigR::extractTFNames(names(pfm))
fam <- list(AP1 = c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "BATF", "BATF3"),
            NFKB = c("NFKB1", "NFKB2", "RELA", "REL", "RELB"),
            MAF = c("MAF", "MAFB", "MAFF", "MAFG", "MAFK"),
            CEBP = c("CEBPA", "CEBPB", "CEBPD", "CEBPG"), KLF3 = "KLF3", TBET = c("TBX21", "EOMES"))
tf_want <- unique(c(unlist(fam), "BACH2", "BCL11B", "SMAD3", "ZNF711", "IKZF3", "STAT1", "STAT4",
                    "BHLHE40", "RORC", "ZBTB16"))
nm <- names(pfm)
sel <- unique(unlist(lapply(tf_want, function(t) which(nm == t | sub("_.*$", "", nm) == t))))
tf_of <- sub("_.*$", "", nm[sel])
mm <- matchMotifs(pfm[sel], se, genome = BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(se, niterations = 200)
expct <- computeExpectations(se)
POPS <- c("gdT", "CD8 TEM", "CD4 TCM", "NK", "MAIT", "CD8 TCM", "CD4 Naive", "CD8 Naive", "B naive", "Treg")

paired <- function(M, j, info) {
  patient <- factor(info$patient); cnd <- factor(info$cond, levels = c("Remission", "Flare"))
  fit <- eBayes(lmFit(M[, j, drop = FALSE], model.matrix(~ patient + cnd)))
  tt <- topTable(fit, coef = "cndFlare", number = Inf, sort.by = "none")
  d <- sapply(levels(patient), function(p) {
    f <- j[info$patient == p & info$cond == "Flare"]; r <- j[info$patient == p & info$cond == "Remission"]
    M[, f] - M[, r]
  })
  d <- matrix(d, nrow = nrow(M))
  data.frame(feature = rownames(tt), delta = tt$logFC, P = tt$P.Value, FDR = tt$adj.P.Val,
             patients_flare_higher = rowSums(d > 0), stringsAsFactors = FALSE)
}
test_all <- function(M, cdat, mode) {
  rownames(M) <- make.unique(nm[sel])
  famM <- t(sapply(fam, function(f) colMeans(M[tf_of %in% f, , drop = FALSE])))
  res <- NULL; fres <- NULL
  for (cty in POPS) {
    j <- which(cdat$cell_type == cty)
    if (length(j) < 8) next
    info <- as.data.frame(cdat[j, ])
    pts <- names(which(table(info$patient) == 2)); if (length(pts) < 4) next
    j <- j[info$patient %in% pts]; info <- info[info$patient %in% pts, ]
    r <- paired(M, j, info); r$TF <- tf_of[match(r$feature, rownames(M))]
    res <- rbind(res, data.frame(mode = mode, cell_type = cty, n_pairs = length(pts), r))
    fres <- rbind(fres, data.frame(mode = mode, cell_type = cty, n_pairs = length(pts), paired(famM, j, info)))
  }
  list(motif = res, family = fres)
}

dev0 <- computeDeviations(se, annotations = mm, background_peaks = bg, expectation = expct)
cdat <- colData(dev0)
R <- list(test_all(deviationScores(dev0), cdat, "z"), test_all(deviations(dev0), cdat, "dev"))

# ---- D: depth report
frag <- colSums(assay(se, "counts"))
info0 <- as.data.frame(colData(se)); info0$fragments <- frag
z0 <- deviationScores(dev0); rownames(z0) <- make.unique(nm[sel])
for (tf in c("STAT1", "STAT4")) info0[[paste0("z_", tf)]] <- z0[tf, ]
write.table(info0[, c("cell_type", "patient", "cond", "n_cells", "fragments", "z_STAT1", "z_STAT4")],
            file.path(OUT, "pseudobulk_depth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dep <- NULL
for (cty in POPS) {
  x <- info0[info0$cell_type == cty, ]
  pts <- names(which(table(x$patient) == 2)); if (length(pts) < 4) next
  x <- x[x$patient %in% pts, ]
  f <- x[x$cond == "Flare", ][order(x$patient[x$cond == "Flare"]), ]
  r <- x[x$cond == "Remission", ][order(x$patient[x$cond == "Remission"]), ]
  lr <- log2(f$fragments / r$fragments)
  row <- data.frame(cell_type = cty, n_pairs = length(pts),
                    median_fragments_flare = median(f$fragments), median_fragments_remission = median(r$fragments),
                    median_nuclei_flare = median(f$n_cells), median_nuclei_remission = median(r$n_cells),
                    pairs_flare_deeper = sum(lr > 0), median_log2_depth_ratio = median(lr))
  for (tf in c("STAT1", "STAT4")) {
    v <- paste0("z_", tf)
    row[[paste0("rho_", tf, "_z_vs_log_fragments")]] <- suppressWarnings(cor(x[[v]], log(x$fragments), method = "spearman"))
    row[[paste0("rho_", tf, "_dz_vs_log2_depth_ratio")]] <- suppressWarnings(cor(f[[v]] - r[[v]], lr, method = "spearman"))
  }
  dep <- rbind(dep, row)
}
write.table(dep, file.path(OUT, "depth_by_population.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- C: equal-depth downsampling within population x patient
cnt <- assay(se, "counts")
for (s in 1:5) {
  set.seed(20261118 + s)
  C2 <- cnt
  for (cty in POPS) {
    j <- which(info0$cell_type == cty)
    for (p in unique(info0$patient[j])) {
      k <- j[info0$patient[j] == p]
      if (length(k) != 2) next
      tot <- frag[k]; hi <- k[which.max(tot)]; target <- min(tot)
      if (tot[hi == k] == target) next
      x <- C2[, hi]
      nz <- which(x > 0)
      x[nz] <- rbinom(length(nz), size = as.integer(x[nz]), prob = target / sum(x))
      C2[, hi] <- x
    }
  }
  se2 <- se; assay(se2, "counts") <- C2
  dv <- computeDeviations(se2, annotations = mm, background_peaks = bg, expectation = expct)
  R[[length(R) + 1]] <- test_all(deviationScores(dv), colData(dv), paste0("down_z_seed", s))
  R[[length(R) + 1]] <- test_all(deviations(dv), colData(dv), paste0("down_dev_seed", s))
  cat("seed", s, "done | fragments after thinning, max pair ratio:",
      round(max(sapply(split(colSums(C2), paste(info0$cell_type, info0$patient)), function(v) if (length(v) == 2) max(v) / min(v) else 1)), 4), "\n")
}
M <- do.call(rbind, lapply(R, `[[`, "motif")); Fm <- do.call(rbind, lapply(R, `[[`, "family"))
write.table(M, file.path(OUT, "motif_results_all_modes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Fm, file.path(OUT, "family_results_all_modes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# reproduction check against the published table
pub <- read.delim(file.path(dirname(OUT), "chromvar_naive_extension", "motif_deviation_flare_vs_remission.tsv"))
a <- M[M$mode == "z", ]; key <- function(d) paste(d$cell_type, d$feature)
m <- match(key(pub), key(a))
cat("published rows matched:", sum(!is.na(m)), "of", nrow(pub), "| max |delta_z difference|:",
    signif(max(abs(pub$delta_z - a$delta[m]), na.rm = TRUE), 3), "| max |FDR difference|:", signif(max(abs(pub$FDR - a$FDR[m]), na.rm = TRUE), 3), "\n")

# summary for STAT1 / STAT4
st <- M[M$TF %in% c("STAT1", "STAT4") & !duplicated(paste(M$mode, M$cell_type, M$TF)), ]
st$kind <- sub("_seed[0-9]+$", "", st$mode)
S <- do.call(rbind, lapply(split(st, list(st$cell_type, st$TF, st$kind), drop = TRUE), function(d) data.frame(
  cell_type = d$cell_type[1], TF = d$TF[1], mode = d$kind[1], n_runs = nrow(d), n_pairs = d$n_pairs[1],
  delta_median = median(d$delta), P_median = median(d$P), FDR_median = median(d$FDR), FDR_max = max(d$FDR),
  runs_FDR_lt_0.05 = sum(d$FDR < 0.05), patients_flare_higher_min = min(d$patients_flare_higher))))
S <- S[order(S$TF, S$cell_type, match(S$mode, c("z", "dev", "down_z", "down_dev"))), ]
write.table(S, file.path(OUT, "STAT_summary_by_mode.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 220)
print(dep, row.names = FALSE, digits = 3)
print(S, row.names = FALSE, digits = 3)
