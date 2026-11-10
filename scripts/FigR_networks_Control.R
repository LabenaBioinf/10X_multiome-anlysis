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
# output_dir <- "./FigR_results_4_all"
# output_dir <- "./FigR_results_4_cisCorr_genelist"
# output_dir <- "./FigR_results_4_cisCorr_cutoff_6"
output_dir <- "./FigR_results_gdT"

dir.create(output_dir, recursive = TRUE)

# setup log file
#log_file <- file(paste0(output_dir, "/FigR.log"), open="wt")
#sink(log_file, type="message")

# Function to clean sparse matrices efficiently by removing rows and columns with only zeros
clean_sparse_matrix_efficient <- function(mat) {
  # Remove rows with all zeros
  non_zero_rows <- rowSums(mat != 0) > 0
  mat <- mat[non_zero_rows, , drop = FALSE]
  
  # Remove columns with all zeros
  non_zero_cols <- colSums(mat != 0) > 0
  mat <- mat[, non_zero_cols, drop = FALSE]
  
  return(mat)
}

run_analysis <- function(ATAC.se, dorcMat.s, RNAmat.s, cisCorr.filt, output_dir, sample, nCores) {

  # Function to clean sparse matrices efficiently by removing rows and columns with only zeros
  clean_sparse_matrix_efficient <- function(mat) {
    # Remove rows with all zeros
    non_zero_rows <- rowSums(mat != 0) > 0
    mat <- mat[non_zero_rows, , drop = FALSE]
    
    # Remove columns with all zeros
    non_zero_cols <- colSums(mat != 0) > 0
    mat <- mat[, non_zero_cols, drop = FALSE]
    
    return(mat)
  }

  # Process and align matrices function
  process_and_align_matrices <- function(ATAC.se, dorcMat.s, RNAmat.s) {
    # Clean and process the ATAC matrix
    if (is(ATAC.se, "SummarizedExperiment")) {
      ATAC_clean <- assay(ATAC.se)
      if (is(ATAC_clean, "dgCMatrix")) {
        ATAC_clean <- clean_sparse_matrix_efficient(ATAC_clean)
        
        # Automatically align row and column names
        common_rows <- intersect(rownames(ATAC_clean), rownames(ATAC.se))
        common_cols <- intersect(colnames(ATAC_clean), colnames(ATAC.se))
        if (length(common_rows) == 0 || length(common_cols) == 0) {
          stop("No matching rows or columns between ATAC_clean and ATAC.se.")
        }
        
        # Subset the ATAC.se object to the common rows and columns
        ATAC.se <- ATAC.se[common_rows, common_cols]
        
        # Update all assays in ATAC.se to match the cleaned dimensions
        for (i in seq_along(assays(ATAC.se))) {
          assay_data <- assays(ATAC.se)[[i]]
          
          # Clean and subset each assay
          if (is(assay_data, "dgCMatrix")) {
            assay_data <- assay_data[common_rows, common_cols, drop = FALSE]
            assay_data <- clean_sparse_matrix_efficient(assay_data)
          } else {
            assay_data <- assay_data[common_rows, common_cols, drop = FALSE]
          }
          
          # Assign the cleaned assay back to ATAC.se
          assays(ATAC.se)[[i]] <- assay_data
        }
        
        print(paste("Cleaned ATAC.se dimensions:", dim(assays(ATAC.se)[[1]])))
      } else {
        message("ATAC_clean is not a sparse matrix.")
      }
    } else {
      message("ATAC.se is not a SummarizedExperiment object.")
    }

    # Clean and process the dorcMat.s matrix
    if (!is.null(dorcMat.s)) {
      if (is(dorcMat.s, "dgCMatrix")) {
        dorcMat.s <- clean_sparse_matrix_efficient(dorcMat.s)
        print(paste("Cleaned dorcMat.s dimensions:", dim(dorcMat.s)))
      } else {
        message("dorcMat.s is not a sparse matrix.")
      }
    } else {
      message("dorcMat.s is NULL or invalid.")
    }

    # Clean and process the RNAmat.s matrix
    if (!is.null(RNAmat.s)) {
      if (is(RNAmat.s, "dgCMatrix")) {
        RNAmat.s <- clean_sparse_matrix_efficient(RNAmat.s)
        print(paste("Cleaned RNAmat.s dimensions:", dim(RNAmat.s)))
      } else {
        message("RNAmat.s is not a sparse matrix.")
      }
    } else {
      message("RNAmat.s is NULL or invalid.")
    }

    # Align rows between dorcMat.s and RNAmat.s
    common_rows <- intersect(rownames(dorcMat.s), rownames(RNAmat.s))
    if (length(common_rows) > 0) {
      dorcMat.s <- dorcMat.s[common_rows, , drop = FALSE]
      RNAmat.s <- RNAmat.s[common_rows, , drop = FALSE]
      message("Row names between dorcMat.s and RNAmat.s are now aligned.")
    } else {
      stop("No common row names between dorcMat.s and RNAmat.s.")
    }

    # Align columns between ATAC_clean and dorcMat.s
    common_cols <- intersect(colnames(assays(ATAC.se)[[1]]), colnames(dorcMat.s))
    if (length(common_cols) == 0) {
      stop("No common column names between ATAC_clean and dorcMat.s.")
    } else {
      message("Column names between ATAC_clean and dorcMat.s are aligned.")
    }

    return(list(ATAC_clean = assays(ATAC.se)[[1]], dorcMat.s = dorcMat.s, RNAmat.s = RNAmat.s, ATAC.se = ATAC.se))
  }

  # Main try-catch for the analysis
  tryCatch({
    # Process matrices
    matrices <- process_and_align_matrices(ATAC.se, dorcMat.s, RNAmat.s)
    ATAC_clean <- matrices$ATAC_clean
    dorcMat.s <- matrices$dorcMat.s
    RNAmat.s <- matrices$RNAmat.s
    ATAC.se <- matrices$ATAC.se  # Updated ATAC.se object

    # Subset cisCorr.filt to only include peaks in ATAC_clean
    cisCorr.filt <- cisCorr.filt[cisCorr.filt$Peak %in% rownames(ATAC_clean), ]
    
    # Check for NA values
    lapply(list(ATAC_clean, dorcMat.s, RNAmat.s), function(mat) {
      if (anyNA(mat)) {
        message("NAs detected in matrix: ", deparse(substitute(mat)))
      }
    })
    
    # Proceed with FigRGRN calculation if no NAs
    figR.d <- tryCatch({
      result <- FigR::runFigRGRN(ATAC.se = ATAC.se,
                                 dorcTab = cisCorr.filt,
                                 genome = "hg38",
                                 dorcMat = dorcMat.s,
                                 rnaMat = RNAmat.s,
                                 nCores = nCores)
      if (is.null(result)) {
        message("FigRGRN returned NULL.")
      }
      result
    }, error = function(e) {
      message("Error in FigRGRN for sample ", sample, ": ", e$message)
      NULL
    })

    # Write result if not NULL and not empty
    if (!is.null(figR.d) && nrow(figR.d) > 0) {
      tryCatch({
        write.table(figR.d, file = paste0(output_dir, "/", sample, "_FigR/", sample, "_figR_d.txt"), 
                    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
        message("Successfully wrote FigRGRN output to file for sample: ", sample)
      }, error = function(e) {
        message("Error writing table for sample ", sample, ": ", e$message)
      })
    } else {
      message("No valid FigRGRN output for sample: ", sample)
    }
  }, error = function(e) {
    message("Error encountered during processing of sample ", sample, ": ", e$message)
  })
}


################################################################################


# import needed data
pbmc <- readRDS(paste0(data_folder, "/harmonized_seurat.rds"))
pbmc <- SeuratObject::UpdateSeuratObject(pbmc)

# split Seurat object according to sample condition
# pbmc_split <- Seurat::SplitObject(pbmc, split.by = "Condition")
pbmc_split <- Seurat::SplitObject(pbmc, split.by = "celltype.condition")
rm(pbmc)

# infer gene regulatory networks
# sample = "Control"
celltype_conditions = list(
                          # "CD8 TCM_Control",
                          # "CD4 Naive_Control",
                          # "B intermediate_Control",
                          # "B memory_Control",
                          # "B naive_Control",
                          # "CD14 Mono_Control",
                          # "CD4 TCM_Control",
                          # "CD8 Proliferating_Control",
                          # "CD8 TEM_Control",
                          # "NK_Control",
                          #"NK Proliferating_Control",
                          "gdT_Control"
                          )

for (sample in celltype_conditions) {

  print(sample)
 # sample = "B naive_Control"

    # create a FigR output directory
    # dir.create(paste0(output_dir, "/", sample, "_FigR"), recursive = TRUE)
  output_path <- paste0(output_dir, "/", sample, "_FigR")
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  # extract normalized RNA counts as a sparse matrix
  RNAmat <- as(pbmc_split[[sample]][["SCT"]]$counts, 'sparseMatrix')
    
  # extract ATAC counts in sparse matrix format needed for cisTopic
  ATACmat <- as(pbmc_split[[sample]][["ATAC"]]$counts, 'sparseMatrix')
  rownames(ATACmat) <- sub("-", ":", rownames(ATACmat))
    
  # create cisTopic object
  cisTopicObject <- cisTopic::createcisTopicObject(ATACmat, project.name = sample)
    
  # # apply cisTopic to define a number of “topics” that are best able to summarize the variability in ATAC counts
  cisTopicObject <- cisTopic::runCGSModels(cisTopicObject, topic = 1:50,
                                            seed = 987, nCores = nCores, burnin = 90,
                                            iterations = 100, addModels = FALSE)
    
  # # select the best model, i.e. the optimal number of topics
  par(mfrow = c(2,1))
  cisTopicObject <- cisTopic::selectModel(cisTopicObject, type = 'maximum')

  # # retrieve and save topic-cell and region-topic assignments
  cisAssign <- t(cisTopic::modelMatSelection(cisTopicObject, target = 'cell', method = 'Probability'))
  saveRDS(cisAssign, paste0(output_dir, "/", sample, "_FigR/", sample, "_cisTopic_assignments.rds"))
  
  # cisAssign <- readRDS(paste0(output_dir, "/", sample, "_FigR/", sample, "_cisTopic_assignments.rds"))
    
  # # calculate a kNN-graph using the topics matrix from cisTopic
  set.seed(123)
  cellkNN <- FNN::get.knn(cisAssign, k = 40)$nn.index

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

  # # compute correlation between RNA expression and peak accessibility
  cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC.se, RNAmat = RNAmat,
                                    genome = "hg38", nCores = nCores, n_bg = 100)
  saveRDS(cisCorr, paste0(output_dir, "/", sample, "_FigR/", sample, "_FigR_cisCorr.rds"))

  # cisCorr <- readRDS(paste0(output_dir, "/", sample, "_FigR/", sample, "_FigR_cisCorr.rds"))    
    
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
  #             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


  png(paste0(output_dir, "/", sample, "_FigR/", sample, "_dorcJPlot.png"), width = 840, height = 480)
  cisCorr.filt %>% FigR::dorcJPlot(cutoff = 2)
  dev.off()
    
  # get DORC scores
  dorcMat <- FigR::getDORCScores(ATAC.se, dorcTab = cisCorr.filt, nCores = nCores)
  rownames(cellkNN) <- colnames(dorcMat)
  write.table(dorcMat, paste0(output_dir, "/", sample, "_FigR/", sample, "_dorcMat.txt"), 
          sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


  # dorcMat <- FigR::getDORCScores(ATAC.se, dorcTab = cisCorr.filt, geneList = dorcGenes.filt, nCores = nCores)
  # rownames(cellkNN) <- colnames(dorcMat)

  if (anyNA(dorcMat)) {
    message("Warning: dorcMat contains missing values (NAs). Handling them before proceeding.")
    dorcMat[is.na(dorcMat)] <- 0  # Or choose another method like mean imputation
  }

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
  print("figR.d")
  print(dim(ATAC.se))
  print(dim(dorcMat.s))
  print(dim(RNAmat.s))
  print(head(cisCorr.filt))

  run_analysis(ATAC.se, dorcMat.s, RNAmat.s, cisCorr.filt, output_dir, sample, nCores)


  # # plot TF-DORC regulation scores scatter plot
  # figR.d %>%
  #   ggplot2::ggplot(aes(Corr.log10P, Enrichment.log10P, color = Score)) +
  #   ggrastr::geom_point_rast(size = 0.05, shape = 16) +
  #   ggplot2::theme_classic() +
  #   ggplot2::scale_color_gradientn(colours = BuenColors::jdb_palette("solar_extra"), limits = c(-3, 3), 
  #                                  oob = scales::squish, breaks = scales::breaks_pretty(n = 3))
  # ggplot2::ggsave(paste0(output_dir, "/", sample, "_FigR/", sample, "_TF-DORC_regulation_scatter_plot.jpg"), width = 12, height = 8, dpi = 300, units = "in")
    
  # # plot a rank based visualization of driver TFs
  # png(paste0(output_dir, "/", sample, "_FigR/", sample, "_rankDrivers_meanScore.png"), width = 840, height = 480)
  # FigR::rankDrivers(figR.d, rankBy = "meanScore", interactive = FALSE)
  # dev.off()

  # # plot the relevant annotated motifs by their associated number of activated or repressed target genes
  # png(paste0(output_dir, "/", sample, "_FigR/", sample, "_rankDrivers_nTargets.png"), width = 840, height = 480)
  # FigR::rankDrivers(figR.d, score.cut = 1, rankBy = "nTargets", interactive = FALSE)
  # dev.off()
    
  # # plot heatmap of TFs (columns) x DORCs (rows), colored by their associated FigR regulation score
  # png(paste0(output_dir, "/", sample, "_FigR/", sample, "_heatmap.png"), width = 840, height = 840)
  # FigR::plotfigRHeatmap(figR.d = figR.d,
  #                       score.cut = 1,
  #                       TFs = unique(figR.d$Motif),
  #                       column_names_gp = gpar(fontsize = 10), # from ComplexHeatmap
  #                       show_row_dend = FALSE) # from ComplexHeatmap
  # dev.off()
    
  # # plot network of associations
  # network <- FigR::plotfigRNetwork(figR.d,
  #                                  score.cut = 1,
  #                                  TFs = unique(figR.d$Motif),
  #                                  weight.edges = TRUE)
  # r2d3::save_d3_html(network, paste0(output_dir, "/", sample, "_FigR/", sample, "_network.html"))
  # saveRDS(network, paste0(output_dir, "/", sample, "_FigR/", sample, "_network.rds"))

}
