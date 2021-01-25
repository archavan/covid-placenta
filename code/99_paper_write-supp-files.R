# write associated supplementary files

library(tidyverse)

### Cluster marker genes ------------------------------------------------------
# marker genes are output from FindAllMarkers function, for all clusters, top genes upto logFC of 0.25. 
markers <- read.csv("results/02_annotation/files/markers_sct.csv")

# rename clusters from 0, 1 to clust_00, clust_01 etc.
markers$cluster <- paste0("clust_", markers$cluster)
markers$cluster <- gsub("(clust_)(\\d)$", "\\10\\2", markers$cluster)

# rearrange columns
markers <- markers[, c("cluster", "gene", "pct.1", "pct.2", "avg_logFC", "p_val", "p_val_adj")]

# filter by adjusted p value
markers <- markers %>% dplyr::filter(p_val_adj < 0.05)

# write
write.csv(markers, "results/99_paper-supp-files/cluster_marker-genes.csv", row.names = FALSE)
### end =======================================================================




