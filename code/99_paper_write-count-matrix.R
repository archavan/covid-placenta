# write count matrix and metadata (to be uploaded to GEO or similar)

setwd("/home/arc78/scratch60/covid-placenta")

# package 
library(Seurat)
library(tibble)

# read seurat object
seur <- readRDS("data/seurat-object_annotated.rds")

# write count matrix
counts <- as.data.frame(seur@assays$RNA@counts)
counts <- rownames_to_column(.data = counts, var = "gene_name")

write.table(counts, 
            "data/count_matrix/counts.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# write metadata
metadata <- seur@meta.data
metadata <- rownames_to_column(.data = metadata, var = "cell_barcode")

write.table(metadata, 
            "data/count_matrix/metadata.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

