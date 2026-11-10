# chromVAR motif accessibility, flare vs remission, paired by patient, on the ATAC pseudobulk
# (cell type x sample; results/object_level/atac_pseudobulk_counts.rds), extending b3_chromvar.R to
# the TFs of the FigR analysis. Unlike FigR edges this is a per-patient measurement, so the paired
# design (limma ~ patient + condition on deviation z) applies directly.
#
# Declared before the first run (23 Oct 2026; results/gdt_chromvar/MOTIFS.md):
#   primary cell types: gdT, CD8 TEM; other T/NK types reported as context.
#   motifs: γδ T FigR hubs with motif support (NFKB1, KLF3, JUNB, FOSL2, BACH2, FOSL1), TFs named in
#   the submitted γδ T text (FOSB, MAF, CEBPD, BCL11B, FOS, SMAD3, TBX21), CD8 TEM hubs (TBX21,
#   ZNF711, IKZF3, STAT1, STAT4, BHLHE40, EOMES), NF-κB (NFKB2, RELA, REL, RELB), AP-1 (JUN, JUND,
#   BATF, BATF3), MAF/CEBP relatives (MAFB, MAFF, MAFG, MAFK, CEBPA, CEBPB, CEBPG), RORC, ZBTB16.
#   family scores (mean deviation z of the family's motifs; one test each, BH over families):
#   AP1 = FOS, FOSB, FOSL1, FOSL2, JUN, JUNB, JUND, BATF, BATF3; NFKB = NFKB1, NFKB2, RELA, REL, RELB;
#   MAF = MAF, MAFB, MAFF, MAFG, MAFK; CEBP = CEBPA, CEBPB, CEBPD, CEBPG; KLF3; TBET = TBX21, EOMES.
#   BH within cell type over single motifs; AP-1 motifs are near-identical in sequence, so single
#   AP-1 members are not interpreted separately.
suppressPackageStartupMessages({
  library(Matrix); library(GenomicRanges); library(SummarizedExperiment); library(chromVAR)
  library(motifmatchr); library(BSgenome.Hsapiens.UCSC.hg38); library(BiocParallel); library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
IN <- args[1]; OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
register(MulticoreParam(8))
set.seed(20261016)
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
cat("motifs selected:", length(sel), "| TFs without a cisBP motif:",
    paste(setdiff(tf_want, tf_of), collapse = ", "), "\n")
mm <- matchMotifs(pfm[sel], se, genome = BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(se, niterations = 200)
dev <- computeDeviations(se, annotations = mm, background_peaks = bg)
z <- deviationScores(dev)
rownames(z) <- make.unique(nm[sel])
famz <- t(sapply(fam, function(f) colMeans(z[tf_of %in% f, , drop = FALSE])))
if (length(fam) == 1) famz <- matrix(famz, nrow = 1, dimnames = list(names(fam), colnames(z)))

paired <- function(M, j, info) {
  patient <- factor(info$patient); cnd <- factor(info$cond, levels = c("Remission", "Flare"))
  fit <- eBayes(lmFit(M[, j, drop = FALSE], model.matrix(~ patient + cnd)))
  tt <- topTable(fit, coef = "cndFlare", number = Inf, sort.by = "none")
  # per-patient direction
  d <- sapply(levels(patient), function(p) {
    f <- j[info$patient == p & info$cond == "Flare"]; r <- j[info$patient == p & info$cond == "Remission"]
    M[, f] - M[, r]
  })
  d <- matrix(d, nrow = nrow(M))
  data.frame(feature = rownames(tt), delta_z = tt$logFC, P = tt$P.Value, FDR = tt$adj.P.Val,
             patients_flare_higher = rowSums(d > 0), stringsAsFactors = FALSE)
}
res <- NULL; fres <- NULL
for (cty in c("gdT", "CD8 TEM", "CD4 TCM", "NK", "MAIT", "CD8 TCM")) {
  j <- which(colData(dev)$cell_type == cty)
  if (length(j) < 8) next
  info <- as.data.frame(colData(dev)[j, ])
  pts <- names(which(table(info$patient) == 2)); if (length(pts) < 4) next
  j <- j[info$patient %in% pts]; info <- info[info$patient %in% pts, ]
  r <- paired(z, j, info); r$TF <- tf_of[match(r$feature, rownames(z))]
  res <- rbind(res, data.frame(cell_type = cty, n_pairs = length(pts), r))
  f <- paired(famz, j, info)
  fres <- rbind(fres, data.frame(cell_type = cty, n_pairs = length(pts), f))
}
write.table(res, file.path(OUT, "motif_deviation_flare_vs_remission.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fres, file.path(OUT, "family_deviation_flare_vs_remission.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
cat("\n== Family scores (primary)\n"); print(fres, row.names = FALSE, digits = 3)
cat("\n== Single motifs, gdT and CD8 TEM\n")
print(res[res$cell_type %in% c("gdT", "CD8 TEM"), ], row.names = FALSE, digits = 3)
