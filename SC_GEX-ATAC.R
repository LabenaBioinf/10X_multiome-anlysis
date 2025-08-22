# needed functions
source("./scripts/preprocess_and_QC.R")
source("./scripts/normalize.R")
source("./scripts/batch_effect.R")
source("./scripts/annotate_and_integrate.R")
source("./scripts/DE.R")
source("./scripts/GSEA.R")
source("./scripts/FigR.R")


################################################################################


# set variables
patients <- c("P2", "P3", "P4", "P5", "P7", "P8", "P10", "P11", "P12", "P13")
misc_samples <- c("c42", "c40", "c38", "c29", "c25", "c16", "c7", "c5", "c6", "c4")
control_samples <- c("c61", "c59", "c64", "c50", "c48", "c49", "c51", "c35", "c45", "c37")

# set input and output directories
data_dir <- "./data/" #'[ For HPC change the path to the one you bind (Example: `/data/SC_RNA_seq_results/` ]
output_dir <- "./results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# set cell type reference
# download an annotated PBMC reference dataset from Hao et al. (2020) from here: 
# https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
cell_reference <- "./resources/pbmc_multimodal.h5seurat" #'[ For HPC change the path to the one you bind (Example: `/home/resources/pbmc_multimodal.h5seurat` ]

# set up genes of special interest
genes_of_interest <- read.table("genes_of_interest.txt")[[1]]

# set the number of thread to use
nCores = 8


################################################################################


# preprocess each patient separately and then merge then into a single Seurat object

# create a list to store Seurat objects for each patient
seurat_objects <- list()
data.type <- c("MISC", "Control")

# iterate over each patient and merge them into one Seurat object
for (i in 1:length(patients)) {
  data <- preprocess(patient_id = patients[i], 
                     misc_sample = misc_samples[i], 
                     control_sample = control_samples[i], 
                     output_dir = output_dir, 
                     nCores = nCores)
  k <- 1
  for (x in data) {
    seurat_objects[[paste0(patients[i], "_", data.type[k])]] <- x
    k <- k + 1
  }
}

# save the list with individual Seurat objects
saveRDS(seurat_objects, file = paste0(output_dir, "/preprocessed_data.rds"))
#seurat_objects <- readRDS(paste0(output_dir, "/preprocessed_data.rds"))

# clean up
rm(data, data.type, x, k, i)
gc()


# merge Seurat for each patient separately
pbmc_merged <- sapply(patients, function(patient) {
  merge(x = seurat_objects[startsWith(names(seurat_objects), patient)][[1]], 
        y = seurat_objects[startsWith(names(seurat_objects), patient)][[2]])
  }, simplify = FALSE, USE.NAMES = TRUE)

# merge all Seurat objects into one
pbmc_combined <- merge(x = seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])

# clean up
rm(seurat_objects)
gc()


################################################################################


# normalize merged datasets
pbmc_merged <- lapply(pbmc_merged, normalize)
pbmc_combined <- normalize(pbmc_combined)


################################################################################


# perform batch correction
harmonized_seurat <- harmonize(data = pbmc_combined, output_dir = output_dir)

# clean up
rm(pbmc_combined)
gc()


################################################################################


# annotate cell types
pbmc_merged <- lapply(pbmc_merged, 
                      function(dataset) annotate(data = dataset, 
                                                 cell_reference = cell_reference, 
                                                 reduction = 'pca'))
harmonized_seurat <- annotate(data = harmonized_seurat, cell_reference = cell_reference, reduction = 'harmony_pca')
gc()

# integrate across conditions
pbmc_merged <- sapply(names(pbmc_merged), 
                      function(patient) integrate(data = pbmc_merged[[patient]], 
                                                  genes_of_interest = genes_of_interest, 
                                                  output = paste0(output_dir, "/per_patient/", patient)),
                      simplify = FALSE, USE.NAMES = TRUE)
harmonized_seurat <- integrate(data = harmonized_seurat, genes_of_interest = genes_of_interest, output = paste0(output_dir, "/all_patients"))

# save final Seurat objects
saveRDS(pbmc_merged, file = paste0(output_dir, "/pbmc_merged_final.rds"))
saveRDS(harmonized_seurat, file = paste0(output_dir, "/harmonized_seurat_final.rds"))
#pbmc_merged <- readRDS(paste0(output_dir, "/pbmc_merged_final.rds"))
#harmonized_seurat <- readRDS(paste0(output_dir, "/harmonized_seurat_final.rds"))


################################################################################


# prepare data for differential expression analysis
pbmc_merged <- lapply(pbmc_merged, prepare_DE)
harmonized_seurat <- prepare_DE(harmonized_seurat)

# perform differential expression analysis
#options(future.globals.maxSize=8500*1024^2)
DE_merged <- sapply(names(pbmc_merged), 
                    function(patient) run_DE(data = pbmc_merged[[patient]], 
                                             genes_of_interest = genes_of_interest, 
                                             output = paste0(output_dir, "/per_patient/", patient)),
                    simplify = FALSE, USE.NAMES = TRUE)
DE_harmonized <- run_DE(data = harmonized_seurat, genes_of_interest = genes_of_interest, output = paste0(output_dir, "/all_patients"))

# extracts DE genes common across 60% samples
get_DE_merged(data = DE_merged, output = paste0(output_dir, "/per_patient/"))
get_DE_harmonized(data = DE_harmonized, output = paste0(output_dir, "/all_patients"))


################################################################################


# run gene set enrichment analysis
run_fgsea(data = DE_harmonized, output_dir = output_dir)
run_irGSEA(data = harmonized_seurat, output_dir = output_dir, nCores = nCores)

# clean up
rm(DE_merged, DE_harmonized)
gc()


################################################################################


# run regulatory network and DORC analysis
sapply(names(pbmc_merged), 
       function(patient) run_FigR(data = pbmc_merged[[patient]], 
                                  output = paste0(output_dir, "/per_patient/", patient),
                                  nCores = nCores,
                                  patient_id = patient,
                                  genes_of_interest = genes_of_interest),
       simplify = FALSE, USE.NAMES = TRUE)
gc()
run_FigR(data = harmonized_seurat,
         output = paste0(output_dir, "/per_patient/", patient),
         nCores = nCores,
         patient_id = "all_patients",
         genes_of_interest = genes_of_interest)
