# Methods (sequencing depth): median raw read pairs per nucleus per library.
# usage: Rscript seq_depth.R <seurat .rds with Cell Ranger ARC per_barcode_metrics in meta.data> <out dir>
# Per-library raw read pairs per nucleus (GEX and ATAC) from Cell Ranger ARC per_barcode_metrics carried in the object metadata
args <- commandArgs(TRUE); obj <- readRDS(args[1]); out <- args[2]
md <- obj@meta.data
writeLines(colnames(md), file.path(out, "metadata_columns.txt"))
lib <- intersect(c("orig.ident","sample","library","dataset"), colnames(md))[1]
cols <- intersect(c("gex_raw_reads","atac_raw_reads","gex_mapped_reads","atac_fragments"), colnames(md))
if (length(cols) == 0) quit(status = 0)
s <- do.call(rbind, lapply(split(md, md[[lib]]), function(d) data.frame(library = d[[lib]][1], nuclei = nrow(d),
      t(sapply(cols, function(cn) median(d[[cn]], na.rm = TRUE))), check.names = FALSE)))
write.table(s, file.path(out, "per_library_median_reads_per_nucleus.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary(s))
