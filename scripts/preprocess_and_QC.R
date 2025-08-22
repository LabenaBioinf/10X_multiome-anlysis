# needed packages
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(IRanges)
library(biovizBase)
library(scDblFinder)
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(ggpubr)


################################################################################


# Function to plot QC metrics

plot_qc <- function(data, patient_id, output_dir) {
  
  # supress scientific notation
  options(scipen=999)
  
  # create metadata dataframe
  metadata <- data@meta.data
  
  # visualize the number of cells per sample
  cells_per_sample <- metadata %>% 
    ggplot2::ggplot(aes(x = Condition, fill = Condition)) + 
    ggplot2::geom_bar() +
    ggplot2::geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.25, size = 3) +
    ggplot2::ggtitle("No. cells") +
    ggplot2::ylab("No. cells") +
    ggplot2::xlab("Sample") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
                   axis.text.x = element_text(size = 10), legend.position = "top", legend.title = element_blank(),
                   plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  # visualize the number of doublets per sample
  dublets_per_sample <- metadata %>% 
    ggplot2::ggplot(aes(Condition, fill = scDblFinder.class)) + 
    ggplot2::geom_bar(position = "stack", stat = "count") +
    ggplot2::geom_text(stat = 'count', aes(label = after_stat(count)), position = position_stack(vjust = 0.5), size = 3) +
    ggplot2::ggtitle("No. dublets") +
    ggplot2::ylab("No. cells") +
    ggplot2::xlab("Sample") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
                   axis.text.x = element_text(size = 10), legend.position = "top", legend.title = element_blank(),
                   plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  # merge the number of cells/doublets per sample plots
  per_sample_plots <- ggarrange(cells_per_sample, dublets_per_sample,
                                nrow = 1, ncol = 2)
  
  # visualize the  the distribution of no. UMIs/transcripts detected per cell
  per_cell_umi <- metadata %>% 
    ggplot2::ggplot(aes(color = Condition, x = nCount_RNA, fill = Condition)) + 
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept = 500) +
    ggplot2::ggtitle("Per cell UMI distribution")  +
    ggplot2::ylab("Cell density") +
    ggplot2::xlab("No. UMIs") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
                   axis.text = element_text(size = 10), legend.title = element_blank(), legend.position = "top",
                   legend.text = element_text(size = 10), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  # visualize the distribution of no. genes detected per cell
  per_cell_genes <- metadata %>% 
    ggplot2::ggplot(aes(color = Condition, x = nFeature_RNA, fill = Condition)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::scale_x_log10() + 
    ggplot2::geom_vline(xintercept = 250) +
    ggplot2::ggtitle("Per cell gene distribution")  +
    ggplot2::ylab("Cell density") +
    ggplot2::xlab("No. genes") +
    ggplot2::theme_classic()  +
    ggplot2::theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
                   axis.text = element_text(size = 10), legend.title = element_blank(), legend.position = "top",
                   legend.text = element_text(size = 10), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  # visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  novelty_score <- metadata %>%
    ggplot2::ggplot(aes(x = log10GenesPerUMI, color = Condition, fill = Condition)) +
    ggplot2::geom_density(alpha = 0.2) +
    ggplot2::geom_vline(xintercept = 0.8) +
    ggplot2::ggtitle("Per UMI gene distribution")  +
    ggplot2::ylab("UMI density") +
    ggplot2::xlab("No. genes") +
    ggplot2::theme_classic()  +
    ggplot2::theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
                   axis.text = element_text(size = 10), legend.title = element_blank(), legend.position = "top",
                   legend.text = element_text(size = 10), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  # visualize the distribution of mitochondrial gene expression detected per cell
  mito_per_cell <- metadata %>% 
    ggplot2::ggplot(aes(color = Condition, x = mitoRatio, fill = Condition)) + 
    ggplot2::geom_density(alpha = 0.2) + 
    ggplot2::scale_x_log10() + 
    ggplot2::geom_vline(xintercept = 0.2) +
    ggplot2::ggtitle("Per cell mitochondrial gene expression")  +
    ggplot2::ylab("Cell density") +
    ggplot2::xlab("Mitochondrial ratio") +
    ggplot2::theme_classic()  +
    ggplot2::theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
                   axis.text = element_text(size = 10), legend.title = element_blank(), legend.position = "top",
                   legend.text = element_text(size = 10), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
  
  # visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
  gene_UMI_mito_corr <- metadata %>% 
    ggplot2::ggplot(aes(x = nCount_RNA, y = nFeature_RNA)) + 
    ggplot2::geom_point(aes(color = mitoRatio), size = 0.5) + 
    ggplot2::scale_colour_gradient(low = "gray90", high = "black") +
    ggplot2::stat_smooth(method = lm) +
    ggplot2::scale_x_log10() + 
    ggplot2::scale_y_log10() + 
    ggplot2::theme_classic() +
    ggplot2::geom_vline(xintercept = 500) +
    ggplot2::geom_hline(yintercept = 250) +
    ggplot2::ylab("No. genes") +
    ggplot2::xlab("No. UMI") +
    ggplot2::labs(color = paste0("Mitochondrial", "\n", "ratio")) +
    ggplot2::theme_classic()  +
    ggplot2::theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12),
                   axis.text = element_text(size = 10), legend.position = "top", legend.title = element_text(size = 10), 
                   legend.spacing.x = unit(0.5, "cm"), strip.text = element_text(size = 12)) +
    ggplot2::facet_wrap(~Condition)
  
  ################################################################################

  # combine plots and export
  ggpubr::ggarrange(per_sample_plots, per_cell_umi, per_cell_genes,
                    gene_UMI_mito_corr, novelty_score, mito_per_cell,
                    labels = c("A)", "B)", "C)", "D)", "E)", "F)"), ncol = 3, nrow = 2)
  ggplot2::ggsave(paste0(output_dir, "/", deparse(substitute(data)), "_QC/", patient_id, "_", deparse(substitute(data)), "_QC_plots.jpg"), 
                  width = 14, height = 8, dpi = 300, units = "in", create.dir = TRUE)
 
  
  # visualize QC metrics as a violin plots
  Seurat::VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "mitoPerc"), ncol = 3, pt.size = 0, split.by = "Condition")
  ggplot2::ggsave(paste0(output_dir, "/", deparse(substitute(data)), "_QC/", patient_id, "_", deparse(substitute(data)), "_QC_RNA_violin.jpg"), 
                  width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  Seurat::VlnPlot(object = data, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 3, pt.size = 0, split.by = "Condition")
  ggplot2::ggsave(paste0(output_dir, "/", deparse(substitute(data)), "_QC/", patient_id, "_", deparse(substitute(data)), "_QC_ATAC_violin.jpg"),
                  width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
  # visualize QC metrics as a feature scatter plot
  feature_plot_1 <- Seurat::FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "mitoPerc")
  feature_plot_2 <- Seurat::FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  Seurat::CombinePlots(plots = list(feature_plot_1, feature_plot_2))
  ggplot2::ggsave(paste0(output_dir, "/", deparse(substitute(data)), "_QC/", patient_id, "_", deparse(substitute(data)), "_QC_feature_scatter.jpg"),
                  width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
}


################################################################################


# Function to preprocess each patient separately

preprocess <- function(patient_id, misc_sample, control_sample, output_dir, nCores) {

  # Import needed data
  
  # load ATAC and RNA data for the MISC sample
  misc_peaks_bed <- read.table(file = paste0(data_dir, misc_sample, "/outs/atac_peaks.bed"), col.names = c("chr", "start", "end"))
  misc_metadata <- read.csv(file = paste0(data_dir, misc_sample, "/outs/per_barcode_metrics.csv"), header = TRUE, row.names = 1) %>% dplyr::filter(is_cell > 0)
  misc_frag_obj <- Signac::CreateFragmentObject(path = paste0(data_dir, misc_sample, "/outs/atac_fragments.tsv.gz"), cells = rownames(misc_metadata))
  misc_RNA <- Seurat::Read10X_h5(paste0(data_dir, misc_sample, "/outs/filtered_feature_bc_matrix.h5"))

  # load ATAC and RNA data for the control sample
  control_peaks_bed <- read.table(file = paste0(data_dir, control_sample, "/outs/atac_peaks.bed"), col.names = c("chr", "start", "end"))
  control_metadata <- read.csv(file = paste0(data_dir, control_sample, "/outs/per_barcode_metrics.csv"), header = TRUE, row.names = 1) %>% dplyr::filter(is_cell > 0)
  control_frag_obj <- Signac::CreateFragmentObject(path = paste0(data_dir, control_sample, "/outs/atac_fragments.tsv.gz"), cells = rownames(control_metadata))
  control_RNA <- Seurat::Read10X_h5(paste0(data_dir, control_sample, "/outs/filtered_feature_bc_matrix.h5"))

  # get gene annotations for hg38
  gene_annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  # edit chromosome names and genome version in annotation data from GRCh38 to hg38
  Signac::seqlevels(gene_annotation) <- paste0('chr', Signac::seqlevels(gene_annotation))
  Signac::genome(gene_annotation) <- "hg38"
  
  ################################################################################

  # Set up Seurat object
  
  # create a common set of peaks across the samples so they can be merged later
  combined_peaks <- GenomicRanges::reduce(x = c(GenomicRanges::makeGRangesFromDataFrame(misc_peaks_bed), GenomicRanges::makeGRangesFromDataFrame(control_peaks_bed)))
  # filter out bad peaks based on length
  combined_peaks <- combined_peaks[IRanges::width(combined_peaks) < 10000 & IRanges::width(combined_peaks) > 20]

  # quantify peaks in each sample
  misc_peak_counts <- Signac::FeatureMatrix(fragments = misc_frag_obj, features = combined_peaks, cells = rownames(misc_metadata))
  control_peak_counts <- Signac::FeatureMatrix(fragments = control_frag_obj, features = combined_peaks, cells = rownames(control_metadata))
  
  # create chromatin assays
  misc_atac_assay <- Signac::CreateChromatinAssay(misc_peak_counts, fragments = misc_frag_obj, annotation = gene_annotation)
  control_atac_assay <- Signac::CreateChromatinAssay(control_peak_counts, fragments = control_frag_obj, annotation = gene_annotation)

  # initialize Seurat objects with chromatin assays
  misc_pbmc_multi <- Seurat::CreateSeuratObject(misc_atac_assay, assay = "ATAC", meta.data = misc_metadata, project = paste0(patient_id, "_MISC"))
  control_pbmc_multi <- Seurat::CreateSeuratObject(control_atac_assay, assay = "ATAC", meta.data = control_metadata, project = paste0(patient_id, "_CONTROL"))
  
  # add RNA counts to initialized Seurat objects
  misc_pbmc_multi[['RNA']] <- Seurat::CreateAssayObject(counts = misc_RNA$`Gene Expression`)
  control_pbmc_multi[['RNA']] <- Seurat::CreateAssayObject(counts = control_RNA$`Gene Expression`)
  
  # add dataset-identifying metadata
  misc_pbmc_multi <- Seurat::AddMetaData(misc_pbmc_multi, metadata = "MISC", col.name = "Condition")
  misc_pbmc_multi <- Seurat::AddMetaData(misc_pbmc_multi, metadata = patient_id, col.name = "Batch")
  
  control_pbmc_multi <- Seurat::AddMetaData(control_pbmc_multi, metadata = "Control", col.name = "Condition")
  control_pbmc_multi <- Seurat::AddMetaData(control_pbmc_multi, metadata = patient_id, col.name = "Batch")
  
  # merge Seurat objects
  pbmc_original <- merge(x = misc_pbmc_multi, y = control_pbmc_multi, add.cell.ids = c(misc_sample, control_sample), project = patient_id)
  
  ################################################################################

  # Perform QC and filter merged dataset
  
  # calculate number of genes detected per UMI (novelty score) and add that info to metadata
  pbmc_original$log10GenesPerUMI <- log10(pbmc_original$nFeature_RNA) / log10(pbmc_original$nCount_RNA)
  
  # calculate percentage and ratio of transcripts mapping to mitochondrial genes
  Seurat::DefaultAssay(pbmc_original) <- "RNA"
  pbmc_original$mitoPerc <- PercentageFeatureSet(object = pbmc_original, pattern = "^MT-")
  pbmc_original$mitoRatio <- pbmc_original@meta.data$mitoPerc / 100
  
  # identify doublets and export the resulting scores back to the Seurat object
  pbmc_original_sce <- scDblFinder::scDblFinder(Seurat::GetAssayData(pbmc_original, slot = "counts"), 
                                                samples = pbmc_original$Condition, 
                                                BPPARAM = BiocParallel::MulticoreParam(nCores, RNGseed=1234))
  pbmc_original$scDblFinder.class <- pbmc_original_sce$scDblFinder.class
  pbmc_original$scDblFinder.score <- pbmc_original_sce$scDblFinder.score
  
  # calculate TTS (transcription start site enrichment score) and the strength of the nucleosome signal per cell
  Seurat::DefaultAssay(pbmc_original) <- "ATAC"
  pbmc_original <- Signac::TSSEnrichment(pbmc_original)
  pbmc_original <- Signac::NucleosomeSignal(pbmc_original)
  
  # plot QC plots on the original dataset
  plot_qc(pbmc_original, patient_id, output_dir)
  
  # filter out low quality cells
  pbmc_filtered <- subset(
    x = pbmc_original,
    subset = nCount_ATAC > 1000 & nCount_ATAC < 100000 &
      TSS.enrichment > 1 & nucleosome_signal < 4 &
      nCount_RNA > 500 & nFeature_RNA > 250 &
      mitoPerc < 20 & scDblFinder.class == "singlet"
  )

  # plot QC plots on the filtered dataset
  plot_qc(pbmc_filtered, patient_id, output_dir)

  pbmc_split <- SplitObject(pbmc_filtered, split.by = "Condition")
  
  return (pbmc_split)
}
