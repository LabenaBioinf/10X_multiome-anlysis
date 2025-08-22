# needed packages
library(Seurat)
library(Signac)


################################################################################


# Function to normalize ATAC and RNA data and run dimension reduction

normalize <- function(data) {
	
	# normalize ATAC assay
	Seurat::DefaultAssay(data) <- "ATAC"
	data <- Signac::RunTFIDF(data)

	# run dimension reduction
	data <- Signac::FindTopFeatures(data, min.cutoff = 5)
	data <- Signac::RunSVD(data)


	# normalize RNA assay
	Seurat::DefaultAssay(data) <- "RNA"

	# split the RNA measurements into two layers one for MISC sample cells, one for control sample cells
	data[["RNA"]] <- split(data[["RNA"]], f = data$Condition)
	data <- Seurat::SCTransform(data)

	# run dimension reduction
	data <- Seurat::RunPCA(data)

	return(data)

}
