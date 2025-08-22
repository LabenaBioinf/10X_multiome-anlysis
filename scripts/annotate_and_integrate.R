# needed packages
library(SeuratDisk)
library(Seurat)
library(ggplot2)


################################################################################


# Function to perform cell type annotation

annotate <- function(data, cell_reference, reduction) {

  # load and update PBMC reference
  reference <- SeuratDisk::LoadH5Seurat(cell_reference, assays = list("SCT" = "counts"), reductions = 'spca')
  reference <- Seurat::UpdateSeuratObject(reference)

  # find transfer anchors
  Seurat::DefaultAssay(data) <- "SCT"
  transfer_anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = data,
    normalization.method = "SCT",
    reference.reduction = "spca",
    recompute.residuals = FALSE,
    dims = 1:50
  )

  # transfer cell type data across datasets
  predictions <- Seurat::TransferData(
    anchorset = transfer_anchors, 
    refdata = reference$celltype.l2,
    weight.reduction = data[[reduction]],
    dims = 1:50
  )

  # add cell type predictions to the Seurat object
  annotated_seurat <- Seurat::AddMetaData(
    object = data,
    metadata = predictions
  )

  # set a reasonable order for cell types to be displayed when plotting
  Seurat::Idents(annotated_seurat) <- "predicted.id"
  levels(annotated_seurat) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating", 
                              "CD8 Naive", "dnT", "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", 
                              "NK", "NK_CD56bright", "NK Proliferating", "gdT", "Treg", "B naive", 
                              "B intermediate", "B memory", "Plasmablast", "CD14 Mono", "CD16 Mono", 
                              "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

  return(annotated_seurat)

}


################################################################################


# Function to integrate and plot integrated data


integrate <- function(data, genes_of_interest, output) {

  # data integration
  Seurat::DefaultAssay(data) <- "SCT"
  data <- Seurat::IntegrateLayers(object = data, method = CCAIntegration, normalization.method = "SCT", verbose = T)
  data <- Seurat::FindNeighbors(data, reduction = "integrated.dr", dims = 1:30)
  data <- Seurat::FindClusters(data, resolution = 0.6)
  data <- Seurat::RunUMAP(data, dims = 1:30, reduction = "integrated.dr")

  # plot UMAP
  Seurat::DimPlot(data, reduction = "umap", group.by = c("Condition", "predicted.id"), label = TRUE)
  ggplot2::ggsave(paste0(output, "/annotation/DimPlot_integrated_SCT_ANN.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # plot UMAP split
  Seurat::DimPlot(data, reduction = "umap", group.by = c("predicted.id"), split.by = "Condition")
  ggplot2::ggsave(paste0(output, "/annotation/DimPlot_integrated_SCT_ANN_split.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # plot DotPlot
  Seurat::Idents(data) <- "predicted.id"
  genes_of_interest_list <- split(genes_of_interest, cut(seq_along(genes_of_interest), 3, labels = FALSE))
  for (gene in names(genes_of_interest_list)) {
    Seurat::DotPlot(data, features = genes_of_interest_list[[gene]], cols = c("red", "blue"), dot.scale = 8, split.by = "Condition") +
      Seurat::RotatedAxis()
    ggplot2::ggsave(paste0(output, "/annotation/DE_DotPlot_", gene,".jpg"), width = 14, height = 14, dpi = 300, units = "in", create.dir = TRUE)
  }

  return(data)
}
