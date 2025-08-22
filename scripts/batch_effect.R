# needed packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)


################################################################################


# Function to plot UMAP for ATAC and RNA (Batch, Condition) before integration

harmonize <- function(data, output_dir) {

  # Plot UMAP for ATAC and RNA (Batch, Condition) before integration
  
  # preprocess ATAC data
  data_temp <- Seurat::FindNeighbors(data, assay="ATAC", reduction = "lsi", dims = 1:30)
  data_temp <- Seurat::FindClusters(data_temp, resolution = 0.6, graph.name = "ATAC_nn")
  data_temp <- Seurat::RunUMAP(data_temp, dims = 1:30, reduction = "lsi")

  # plot ATAC data
  Seurat::DimPlot(data_temp, reduction = "umap", group.by = c("Condition", "Batch"), label = TRUE)
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_ATAC.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # plot ATAC data split
  Seurat::DimPlot(data_temp, reduction = "umap", group.by = c("Batch"), split.by = "Condition")
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_ATAC_split.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)


  # preprocess RNA data
  data_temp <- Seurat::FindNeighbors(data_temp, assay="SCT", reduction = "pca", dims = 1:30)
  data_temp <- Seurat::FindClusters(data_temp, resolution = 0.6, graph.name = "SCT_nn")
  data_temp <- Seurat::RunUMAP(data_temp, dims = 1:30, reduction = "pca")

  # plot RNA data
  dimplot <- Seurat::DimPlot(data_temp, reduction = "umap", group.by = c("Condition", "Batch"), label = TRUE)
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_SCT.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # plot RNA data split
  dimplot_split <- Seurat::DimPlot(data_temp, reduction = "umap", group.by = c("Batch"), split.by = "Condition")
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_SCT_split.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  ################################################################################

  # remove batch effects in ATAC data
  harmonized_seurat <- harmony::RunHarmony(data, 
                                           group.by.vars = "Batch", 
                                           reduction = "lsi", assay.use = "ATAC", reduction.save = "harmony_lsi", project.dim = F)

  # process harmonized ATAC data
  harmonized_seurat <- Seurat::FindNeighbors(harmonized_seurat, assay="ATAC", reduction = "harmony_lsi", dims = 1:30)
  harmonized_seurat <- Seurat::FindClusters(harmonized_seurat, resolution = 0.6, graph.name = "ATAC_nn")
  harmonized_seurat <- Seurat::RunUMAP(harmonized_seurat, reduction="harmony_lsi", dims=1:30)
  
  # plot harmonized ATAC data
  dimplot <- Seurat::DimPlot(harmonized_seurat, reduction = "umap", group.by = c("Condition", "Batch"), label = TRUE)
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_harmony_ATAC.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
  # plot harmonized ATAC data split
  dimplot_split <- Seurat::DimPlot(harmonized_seurat, reduction = "umap", group.by = c("Batch"), split.by = "Condition")
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_harmony_ATAC_split.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
  # remove batch effects in RNA data
  harmonized_seurat <- harmony::RunHarmony(data, 
                                           group.by.vars = "Batch", 
                                           reduction = "pca", assay.use = "SCT", reduction.save = "harmony_pca")

  # process harmonized RNA data
  harmonized_seurat <- Seurat::FindNeighbors(harmonized_seurat, assay="SCT", reduction = "harmony_pca", dims = 1:30)
  harmonized_seurat <- Seurat::FindClusters(harmonized_seurat, resolution = 0.6, graph.name = "SCT_nn")
  harmonized_seurat <- Seurat::RunUMAP(harmonized_seurat, reduction="harmony_pca", dims=1:30)

  # plot harmonized RNA data
  Seurat::DimPlot(harmonized_seurat, reduction = "umap", group.by = c("Condition", "Batch"), label = TRUE)
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_harmony_SCT.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # plot harmonized RNA data split
  Seurat::DimPlot(harmonized_seurat, reduction = "umap", group.by = c("Batch"), split.by = "Condition")
  ggplot2::ggsave(paste0(output_dir, "/batch_effect/DimPlot_harmony_SCT_split.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  return(harmonized_seurat)
}
