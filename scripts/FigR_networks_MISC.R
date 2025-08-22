# Bla.


# needed packages:
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(cisTopic)
library(umap)
library(FNN)
library(Matrix)
library(BuenColors)
library(ComplexHeatmap)
library(FigR)
library(ggrepel)
library(doParallel)
library(ggrastr)
library(scales)
library(grid)
library(networkD3)
library(r2d3)


################################################################################


# set up input data folders
data_folder = "/home/data"

# set up genes of special interest
# genes_of_interest <- read.table("./genes_of_interest_all_pat_6.txt")[[1]]

# set the number of threads to use for scDblFinder, cisTopic and FigR functions
nCores = 50

# create output directory
# output_dir <- "./FigR_results_4_cisCorr_genelist"
# output_dir <- "./FigR_results_4_cisCorr_cutoff_6"
output_dir <- "./FigR_results_4_all"


dir.create(output_dir, recursive = TRUE)

# setup log file
#log_file <- file(paste0(output_dir, "/FigR.log"), open="wt")
#sink(log_file, type="message")


################################################################################


# import needed data
pbmc <- readRDS(paste0(data_folder, "/harmonized_seurat.rds"))
pbmc <- SeuratObject::UpdateSeuratObject(pbmc)

# split Seurat object according to sample condition
pbmc_split <- Seurat::SplitObject(pbmc, split.by = "Condition")
rm(pbmc)

# infer gene regulatory networks
sample = "MISC"

# create a FigR output directory
dir.create(paste0(output_dir, "/", sample, "_FigR"), recursive = TRUE)
  
# extract normalized RNA counts as a sparse matrix
RNAmat <- as(pbmc_split[[sample]][["SCT"]]$counts, 'sparseMatrix')
  
# extract ATAC counts in sparse matrix format needed for cisTopic
ATACmat <- as(pbmc_split[[sample]][["ATAC"]]$counts, 'sparseMatrix')
rownames(ATACmat) <- sub("-", ":", rownames(ATACmat))
  
# create cisTopic object
#cisTopicObject <- cisTopic::createcisTopicObject(ATACmat, project.name = sample)
  
# apply cisTopic to define a number of “topics” that are best able to summarize the variability in ATAC counts
#cisTopicObject <- cisTopic::runCGSModels(cisTopicObject, topic = 1:50,
#                                         seed = 987, nCores = nCores, burnin = 90,
#                                         iterations = 100, addModels = FALSE)
  
# select the best model, i.e. the optimal number of topics
#par(mfrow = c(2,1))
#cisTopicObject <- cisTopic::selectModel(cisTopicObject, type = 'maximum')

# retrieve and save topic-cell and region-topic assignments
#cisAssign <- t(cisTopic::modelMatSelection(cisTopicObject, target = 'cell', method = 'Probability'))
#saveRDS(cisAssign, paste0(output_dir, "/", sample, "_FigR/", sample, "_cisTopic_assignments.rds"))
cisAssign <- readRDS(paste0(output_dir, "/", sample, "_FigR/", sample, "_cisTopic_assignments.rds"))
  
# calculate a kNN-graph using the topics matrix from cisTopic
set.seed(123)
cellkNN <- FNN::get.knn(cisAssign, k = 30)$nn.index

# prepare ATAC data in SummarizedExperiment object format needed for runGenePeakcorr
ATAC.se <- SummarizedExperiment::makeSummarizedExperimentFromDataFrame(
  as.data.frame(ATACmat) %>%
    dplyr::mutate(seqnames =  sapply(strsplit(rownames(.), ":|-"), `[`, 1)) %>%
    dplyr::mutate(start = sapply(strsplit(rownames(.), ":|-"), `[`, 2)) %>%
    dplyr::mutate(end = sapply(strsplit(rownames(.), ":|-"), `[`, 3)) %>%
    dplyr::filter(grepl('chr', rownames(.)))
  )
assay(ATAC.se) <- as(assay(ATAC.se), 'sparseMatrix')
counts(ATAC.se) <- assay(ATAC.se)

# compute correlation between RNA expression and peak accessibility
# cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se, RNAmat = RNAmat,
#                                 genome = "hg38", nCores = nCores, n_bg = 100)
# saveRDS(cisCorr, paste0(output_dir, "/", sample, "_FigR/", sample, "_FigR_cisCorr.rds"))
# cisCorr <- readRDS(paste0(output_dir, "/", sample, "_FigR/", sample, "_FigR_cisCorr.rds"))
cisCorr <- readRDS(paste0(output_dir, "/", sample, "_FigR/", sample, "_only_FigR_cisCorr.rds"))  
  
##############################################################################
  
  
# filter peak-gene correlations by p-value     
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
cat('All associations:', nrow(cisCorr))
cat('Pval <= 0.05 associations:', nrow(cisCorr.filt))
  
# get the number of significant associations for each gene
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
write.table(numDorcs, paste0(output_dir, "/", sample, "_FigR/", sample, "_numDorcs.txt"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
# determine DORC (domains of regulatory chromatin) genes

# filter DOCR genes to include only those of interest

# dorcGenes.filt <- cisCorr.filt %>% FigR::dorcJPlot(cutoff = 6, returnGeneList = TRUE)

# dorcGenes <-  unique(cisCorr.filt$Gene)
# dorcGenes.filt <- dorcGenes[dorcGenes %in% genes_of_interest]
# print(dorcGenes.filt)
# write.table(dorcGenes.filt, paste0(output_dir, "/", sample, "_FigR/", sample, "_filteredDorc.txt"), 
#            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

png(paste0(output_dir, "/", sample, "_FigR/", sample, "_dorcJPlot.png"), width = 840, height = 480)
cisCorr.filt %>% FigR::dorcJPlot(cutoff = 2)
dev.off()  

# get DORC scores
# dorcMat <- FigR::getDORCScores(ATAC.se, dorcTab = cisCorr.filt, geneList = dorcGenes.filt, nCores = nCores)

dorcMat <- FigR::getDORCScores(ATAC.se, dorcTab = cisCorr.filt, nCores = nCores)
rownames(cellkNN) <- colnames(dorcMat)
write.table(dorcMat, paste0(output_dir, "/", sample, "_FigR/", sample, "_dorcMat.txt"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  

# smooth DORC scores using cell KNNs
dorcMat.s <- FigR::smoothScoresNN(NNmat = cellkNN, mat = dorcMat, nCores = nCores)
  
# smooth RNA matrix using cell KNNs
RNAmat.s <- FigR::smoothScoresNN(NNmat = cellkNN, mat = RNAmat, nCores = nCores)
  
# visualize on pre-computed UMAP
# UMAP <- as.data.frame(pbmc_split[[sample]][["umap"]]@cell.embeddings)
# dir.create(paste0(output_dir, "/", sample, "_FigR/UMAP_plots"), recursive = TRUE)
# for (gene in dorcGenes.filt) {
#   dorcg <- plotMarker2D(UMAP, dorcMat.s, markers = gene, maxCutoff = "q0.99", colorPalette = "brewer_heat") + 
#     ggtitle(paste0(gene, " DORC"))
#   rnag <- plotMarker2D(UMAP, RNAmat.s, markers = gene, maxCutoff = "q0.99", colorPalette = "brewer_purple") + 
#     ggtitle(paste0(gene, " RNA"))
#   png(paste0(output_dir, "/", sample, "_FigR/UMAP_plots/", sample, "_UMAP_", gene, ".png"), width = 840, height = 480)
#   print(dorcg + rnag)
#   dev.off()
# }
  
# run TF motif-to-gene associations using reference DORC peak-gene mappings and TF RNA expression levels
figR.d <- FigR::runFigRGRN(ATAC.se = ATAC.se,
                           dorcTab = cisCorr.filt,
                           genome = "hg38",
                           dorcMat = dorcMat.s,
                           rnaMat = RNAmat.s,
                           nCores = nCores)

write.table(figR.d, paste0(output_dir, "/", sample, "_FigR/", sample, "_figR_d.txt"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
                             
# plot TF-DORC regulation scores scatter plot
figR.d %>%
  ggplot2::ggplot(aes(Corr.log10P, Enrichment.log10P, color = Score)) +
  ggrastr::geom_point_rast(size = 0.05, shape = 16) +
  ggplot2::theme_classic() +
  ggplot2::scale_color_gradientn(colours = BuenColors::jdb_palette("solar_extra"), limits = c(-3, 3), 
                                 oob = scales::squish, breaks = scales::breaks_pretty(n = 3))
ggplot2::ggsave(paste0(output_dir, "/", sample, "_FigR/", sample, "_TF-DORC_regulation_scatter_plot.jpg"), width = 12, height = 8, dpi = 300, units = "in")
  
# plot a rank based visualization of driver TFs
png(paste0(output_dir, "/", sample, "_FigR/", sample, "_rankDrivers_meanScore.png"), width = 840, height = 480)
FigR::rankDrivers(figR.d, rankBy = "meanScore", interactive = FALSE)
dev.off()

# plot the relevant annotated motifs by their associated number of activated or repressed target genes
png(paste0(output_dir, "/", sample, "_FigR/", sample, "_rankDrivers_nTargets.png"), width = 840, height = 480)
FigR::rankDrivers(figR.d, score.cut = 1, rankBy = "nTargets", interactive = FALSE)
dev.off()
  
# plot heatmap of TFs (columns) x DORCs (rows), colored by their associated FigR regulation score
png(paste0(output_dir, "/", sample, "_FigR/", sample, "_heatmap.png"), width = 840, height = 840)
FigR::plotfigRHeatmap(figR.d = figR.d,
                      score.cut = 1,
                      TFs = unique(figR.d$Motif),
                      column_names_gp = gpar(fontsize = 10), # from ComplexHeatmap
                      show_row_dend = FALSE) # from ComplexHeatmap
dev.off()
  
# plot network of associations
network <- FigR::plotfigRNetwork(figR.d,
                                 score.cut = 1,
                                 TFs = unique(figR.d$Motif),
                                 weight.edges = TRUE)
r2d3::save_d3_html(network, paste0(output_dir, "/", sample, "_FigR/", sample, "_network.html"))
saveRDS(network, paste0(output_dir, "/", sample, "_FigR/", sample, "_network.rds"))


