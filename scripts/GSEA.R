# needed packages
library(Seurat)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tibble)
library(ggplot2)
library(irGSEA)
library(RcppML)


################################################################################


# Function to run fgsea

run_fgsea <- function(data, output_dir) {

  # get cell types
  cell_types <- names(data)

  for (cell in cell_types) {
    
    # get Molecular Signatures Database (MSigDB) gene sets
    m_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2")

    # filter the .gmt file and convert it to a format accepted by fgsea
    fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

    # rank DE genes
    ranks <- data[[cell]]$avg_log2FC
    names(ranks) <- rownames(data[[cell]])
  
    # run fgsea
    #fgseaRes <- fgsea::fgsea(fgsea_sets, stats = ranks,
    #                         minSize = 10,
    #                         maxSize = 500,
    #                         nperm = 1000000)
  
    # run fgsea
    fgseaRes <- fgsea::fgsea(fgsea_sets, stats = ranks,
                             minSize = 10,
                             maxSize = 500)
  
    # format results
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))

    # plot top 20 pathways
    ggplot2::ggplot(fgseaResTidy %>% dplyr::filter(padj < 0.05) %>% head(n = 20), aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill = NES < 7.5)) +
      coord_flip() +
      labs(x = "Pathway", 
           y = "Normalized Enrichment Score",
           title = "REACTOME pathways NES from GSEA") + 
      theme_minimal()
    ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/fgsea/fgsea_", sub(" ", "_", cell), "_ggplot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
    # plot GSEA enrichment plot
    fgsea::plotEnrichment(fgsea_sets[["REACTOME_INTERLEUKIN_3_INTERLEUKIN_5_AND_GM_CSF_SIGNALING"]], ranks) + 
      labs(title = "REACTOME_INTERLEUKIN_3_INTERLEUKIN_5_AND_GM_CSF_SIGNALING")
    ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/fgsea/fgsea_", sub(" ", "_", cell), "_enrichmentplot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
    # convert leadingEdge list column to character vector
    fgseaResTidy$leadingEdge <- sapply(fgseaResTidy$leadingEdge, paste, collapse = ",")
  
    # write results to file
    write.table(fgseaResTidy, file = paste0(output_dir, "/all_patients/GSEA/fgsea/fgseaResTidy_", sub(" ", "_", cell), ".txt"), sep = "\t", col.names = TRUE, quote = FALSE)
  
  }
}


################################################################################


# Function to run irGSEA

run_irGSEA <- function(data, output_dir, nCores) {

  Idents(data) <- data$predicted.id
  
  # calculate enrichment scores
  data <- irGSEA::irGSEA.score(object = data, assay = "SCT", 
                               slot = "data", seeds = 123, ncores = nCores,
                               min.cells = 3, min.feature = 0,
                               custom = F, geneset = NULL, msigdb = T, 
                               species = "Homo sapiens", category = "C2",  
                               subcategory = "REACTOME", geneid = "symbol",
                               method = c("UCell", "singscore"),
                               aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                               kcdf = 'Gaussian')
  
  # integrate differential gene set
  result.dge <- irGSEA::irGSEA.integrate(object = data, 
                                         group.by = "predicted.id",
                                         metadata = NULL, col.name = "celltype.condition",
                                         method = c("UCell", "singscore"))
  

  # show co-upregulated or co-downregulated gene sets per cluster in RRA
  irGSEA::irGSEA.heatmap(object = result.dge, method = "RRA", top = 50, show.geneset = NULL)
  ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.heatmap.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # show co-upregulated or co-downregulated gene sets per cluster in RRA.
  irGSEA::irGSEA.bubble(object = result.dge, method = "RRA", top = 50)
  ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.bubble.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # show the intersections of significant gene sets among clusters in RRA
  irGSEA::irGSEA.upset(object = result.dge, method = "RRA")
  ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.upset.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  # show the intersections of significant gene sets among clusters in all methods
  irGSEA::irGSEA.barplot(object = result.dge, method = c("UCell", "singscore"))
  ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/", "irGSEA.barplot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)

  ################################################################################

  # Show the expression and distribution of special gene sets in special gene set enrichment analysis method
  
  # show the expression and distribution of “REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING” in Ucell on UMAP plot
  #irGSEA::irGSEA.density.scatterplot(object = data, method = "singscore", show.geneset = "REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING", reduction = "umap")
  #ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.R-I3-I5-GM-CSF.density.scatterplot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
  # show the expression and distribution of “REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING” in Ucell among clusters
  #irGSEA::irGSEA.halfvlnplot(object = data, method = "UCell", show.geneset = "REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING")
  #ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.R-I3-I5-GM-CSF.halfvlnplot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
  # show the expression and distribution of “REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING” between  UCell and singscore among clusters
  #irGSEA::irGSEA.vlnplot(object = data, method = c("UCell", "singscore"), show.geneset = "REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING")
  #ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.R-I3-I5-GM-CSF.vlnplot.plot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
  # show the expression and distribution of “REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING” in Ucell among clusters
  #irGSEA::irGSEA.ridgeplot(object = data, method = "UCell", show.geneset = "REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING")
  #ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.R-I3-I5-GM-CSF.ridgeplot.plot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
  # show the expression and distribution of “REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING” in Ucell among clusters
  #irGSEA::irGSEA.densityheatmap(object = data, method = "UCell", show.geneset = "REACTOME-INTERLEUKIN-3-INTERLEUKIN-5-AND-GM-CSF-SIGNALING")
  #ggplot2::ggsave(paste0(output_dir, "/all_patients/GSEA/irGSEA/irGSEA.R-I3-I5-GM-CSF.densityheatmap.plot.jpg"), width = 14, height = 7, dpi = 300, units = "in", create.dir = TRUE)
  
}

