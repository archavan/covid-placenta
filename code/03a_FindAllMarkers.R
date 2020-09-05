library(Seurat)

# read seurat-processed data
plac <- readRDS("/home/arc78/scratch60/covid-placenta/data/placenta.rds")

# set idents to seurat clusters
Idents(plac) <- factor(plac@meta.data$seurat_clusters)

# set default assay to SCT
DefaultAssay(object = plac) <- "SCT"

# find markers for all clusters
markers <- FindAllMarkers(object = plac, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  

# write.csv
write.csv(markers, "/home/arc78/scratch60/covid-placenta/markers_sct.csv", row.names = FALSE)
