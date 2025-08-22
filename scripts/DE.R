# needed packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)


################################################################################


# Function to prepare data for differential expression anaylsis

prepare_DE <- function(data) {

	# prepare the dataset for DE
	data <- Seurat::PrepSCTFindMarkers(data)
	data$celltype.condition <- paste(data$predicted.id, data$Condition, sep = "_")
	Seurat::Idents(data) <- "celltype.condition"

	return(data)
}


################################################################################


# Function to perform differential expression anaylsis

run_DE <- function(data, genes_of_interest, output) {


	# get the number of cells for each cell type
	Seurat::Idents(data) <- "predicted.id"
	cell_types <- levels(Seurat::Idents(data))
	cell_type_count <- as.data.frame.matrix(table(data$predicted.id, data$Condition))
	write.table(cell_type_count, paste0(output, "/cell_type_count.txt"), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)


	# run DE for all cell types
	options(scipen = 0)
	dir.create(paste0(output, "/DE/tables"), recursive = TRUE)
	Seurat::Idents(data) <- "celltype.condition"
	DE_results <- list()
	for (cell in cell_types) {
  		try({
    		DE_results[[cell]] <- Seurat::FindMarkers(data, ident.1 = paste0(cell, "_MISC"), ident.2 = paste0(cell, "_Control"), verbose = FALSE)
    		write.table(DE_results[[cell]], paste0(output, "/DE/tables/DE_", sub(" ", "_", cell), "_results.txt"), sep = "\t", col.names = TRUE, quote = FALSE)
  		})
	}


	# plot FeaturePlot
	genes_of_interest_list <- split(genes_of_interest, cut(seq_along(genes_of_interest), 12, labels = FALSE))
	for (gene in names(genes_of_interest_list)) {
  		Seurat::FeaturePlot(data, features = genes_of_interest_list[[gene]], 
    		        split.by = "Condition", max.cutoff = 3, cols = c("grey", "red"), reduction = "umap")
  		ggplot2::ggsave(paste0(output, "/DE/plots/DE_Feature_plot_", gene, ".jpg"), width = 14, height = 14, dpi = 300, units = "in", create.dir = TRUE)
	}

	# plot Violin plots
	genes_of_interest_list <- split(genes_of_interest, cut(seq_along(genes_of_interest), 8, labels = FALSE))
	for (gene in names(genes_of_interest_list)) {
		vlnplots <- Seurat::VlnPlot(data, genes_of_interest_list[[gene]],
           		                    split.by = "Condition", group.by = "predicted.id", pt.size = 0, combine = FALSE, add.noise = FALSE)
  		vlnplot <- patchwork::wrap_plots(plots = vlnplots, ncol = 1)
  		ggplot2::ggsave(paste0(output, "/DE/plots/DE_VlnPlot_", gene, ".jpg"), width = 14, height = 20, dpi = 300, units = "in", create.dir = TRUE)
	}

	# plot volcano
	for (cell in names(DE_results)) {
  
  		# edit DE results for volcano plots
  		DE_results[[cell]] %>%
    		dplyr::add_rownames(var = "gene") %>%
    		# add up, down or normally regulated labels
    		dplyr::mutate(regulation = case_when(avg_log2FC > 1 & p_val_adj < 0.05 ~ "Up",
            		                             avg_log2FC < -1 & p_val_adj < 0.05 ~ "Down", 
                    		                     .default = "Normal")) %>%
    	# add labels for volcano plots
    	dplyr::mutate(labels = case_when(avg_log2FC < quantile(avg_log2FC, 0.02) ~ gene,
        	                             avg_log2FC > quantile(avg_log2FC, 0.98) ~ gene,
            	                         -log10(p_val_adj) > quantile(-log10(p_val_adj), 0.98) ~ gene)) %>%
    	# plot volcano plots
    	ggplot2::ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    	ggplot2::geom_point(aes(color = regulation)) +
    	ggrepel::geom_text_repel(aes(label = labels)) +
    	ggplot2::scale_color_manual(values = c("blue", "gray", "red")) +
    	ggplot2::ggtitle(cell) +
    	ggplot2::xlab(expression("Log"[2]*" fold change")) +
    	ggplot2::ylab(expression("-Log"[10]*" P")) +
    	ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    	ggplot2::geom_vline(xintercept = 1, linetype="dashed") +
    	ggplot2::geom_vline(xintercept = -1, linetype="dashed") +
    	ggplot2::theme_bw() +
    	ggplot2::theme(legend.title = element_blank(), legend.text = element_text(size = 14),
        	           title = element_text(size = 14, face = "bold"),
            	       axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
                	   axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
  		ggplot2::ggsave(paste0(output, "/DE/volcanos/", sub(" ", "_", cell), "_volcano.png"), width = 7, height = 6,  dpi = 300, units = "in", create.dir = TRUE)
	}

	return(DE_results)
}


################################################################################


# Function to extract DE genes common across 60% samples

get_DE_merged <- function(data, output) {

	# get cell types
	cell_types <- list()
	for (patient in names(data)) {
  		cell_types[[patient]] <- names(data[[patient]])
	}
	cell_types <- unique(unlist(cell_types))


	# get a list of significantly DE genes for each cell type
	DE_results_up_filtered <- list()
	DE_results_down_filtered <- list()
	for (patient in names(data))  {
  		for (cell in cell_types) {
    		if (cell %in% names(data[[patient]])) {
      			DE_results_up_filtered[[cell]][[patient]] <- data[[patient]][[cell]] %>%
        			dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
        			rownames()
      			DE_results_down_filtered[[cell]][[patient]] <- data[[patient]][[cell]] %>%
        			dplyr::filter(avg_log2FC < -1 & p_val_adj < 0.05) %>%
        			rownames()
    		}
  		}
	}

	# for each cell type find DE genes common across 60% samples
	dir.create(output, recursive = TRUE)
	sample_threshold <- length(names(data)) * 0.6 
	DE_up_genes <- list()
	DE_down_genes <- list()
	for(cell in cell_types) {
  		DE_up_genes_table <- table(unlist(DE_results_up_filtered[[cell]]))
  		DE_up_genes[[cell]] <- names(DE_up_genes_table[DE_up_genes_table >= sample_threshold])
  		DE_down_genes_table <- table(unlist(DE_results_down_filtered[[cell]]))
  		DE_down_genes[[cell]] <- names(DE_down_genes_table[DE_down_genes_table >= sample_threshold])
	}

	# output to file
	write.table(sapply(DE_up_genes, "length<-", max(lengths(DE_up_genes))),
    	        paste0(output, "/DE_up_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
	write.table(sapply(DE_down_genes, "length<-", max(lengths(DE_down_genes))),
    	        paste0(output, "/DE_down_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}


################################################################################


# Function to extract DE genes common across all samples

get_DE_harmonized <- function(data, output) {

	# get cell types
	cell_types <- names(data)

	# get a list of significantly DE genes for each cell type
	DE_results_up_filtered <- list()
	DE_results_down_filtered <- list()
	for (cell in cell_types) {
      	DE_results_up_filtered[[cell]] <- data[[cell]] %>%
        	dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>%
        	rownames()
      	DE_results_down_filtered[[cell]] <- data[[cell]] %>%
        	dplyr::filter(avg_log2FC < -1 & p_val_adj < 0.05) %>%
        	rownames()
	}

	# output to file
	write.table(sapply(DE_results_up_filtered, "length<-", max(lengths(DE_results_up_filtered))),
    	        paste0(output, "/DE/DE_up_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
	write.table(sapply(DE_results_down_filtered, "length<-", max(lengths(DE_results_down_filtered))),
    	        paste0(output, "/DE/DE_down_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}
