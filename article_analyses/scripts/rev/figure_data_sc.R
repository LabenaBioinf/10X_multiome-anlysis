# Per-sample values for the revised figures (single-nucleus data):
#   umap.tsv           UMAP coordinates, cell type and phase per nucleus (Figure 3A-B)
#   pb_ifng.tsv        pseudobulk log2 CPM of IFNG and IFNG-AS1 per cell type and sample (Figure 5B)
#   chromvar_stat.tsv  chromVAR deviation z of STAT1/STAT4 motifs per cell type and sample (Figure 6A);
#                      same peaks, motifs, seed and background as gdt_chromvar.R, so values match its tests
# usage: Rscript figure_data_sc.R <harmonized_seurat.rds> <results dir> <out dir>
suppressPackageStartupMessages({
  library(Seurat); library(edgeR); library(chromVAR); library(motifmatchr); library(SummarizedExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38); library(BiocParallel)
})
args <- commandArgs(trailingOnly = TRUE)
OBJ <- args[1]; RES <- args[2]; OUT <- args[3]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Pseudobulk IFNG / IFNG-AS1
P <- readRDS(file.path(RES, "pseudobulk", "pseudobulk_counts.rds"))
lcpm <- cpm(DGEList(P$counts), log = TRUE, prior.count = 1)
g <- intersect(c("IFNG", "IFNG-AS1"), rownames(lcpm))
pb <- do.call(rbind, lapply(g, function(x) data.frame(P$info[, c("cell_type", "patient", "condition", "n_cells")],
                                                       gene = x, log2cpm = lcpm[x, P$info$column])))
write.table(pb, file.path(OUT, "pb_ifng.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# chromVAR STAT motifs, as in gdt_chromvar.R
register(MulticoreParam(8))
set.seed(20261016)
A <- readRDS(file.path(RES, "object_level", "atac_pseudobulk_counts.rds"))
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
nm <- names(pfm)
sel <- which(sub("_.*$", "", nm) %in% c("STAT1", "STAT4"))
mm <- matchMotifs(pfm[sel], se, genome = BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(se, niterations = 200)
z <- deviationScores(computeDeviations(se, annotations = mm, background_peaks = bg))
info <- as.data.frame(colData(se))
cv <- do.call(rbind, lapply(seq_len(nrow(z)), function(i)
  data.frame(info[, intersect(c("cell_type", "patient", "cond", "condition", "n_cells"), names(info))],
             motif = sub("_.*$", "", nm[sel][i]), z = z[i, ])))
write.table(cv, file.path(OUT, "chromvar_stat.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# UMAP
obj <- readRDS(OBJ)
red <- grep("umap", Reductions(obj), ignore.case = TRUE, value = TRUE)
cat("reductions:", paste(Reductions(obj), collapse = ", "), "| using", red[1], "\n")
emb <- Embeddings(obj, red[1])
md <- obj@meta.data
um <- data.frame(umap1 = emb[, 1], umap2 = emb[, 2], cell_type = md$predicted.id, condition = md$Condition,
                 patient = md$Batch)
write.table(um, file.path(OUT, "umap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("done\n")
