# Generate input file for metascape. The input is a single file where each column is a list of DE genes in a celltype. 

### packages ==================================================================
library(tidyverse)

### DE genes data =============================================================
celltypes <- c("APC", "Bcell", "dec.DSC", "dec.Endo", "dec.FB", "Gran", "Mono_1", "Mono_2", "NK_1", "NK_2", "NK_3", "dec.SMC", "Tcell_1", "Tcell_2", "Tcell_3", "vil.Ery", "vil.EVT", "vil.FB", "vil.Hofb", "vil.SCT", "vil.VCT")

de.genes <- list()
for(i in celltypes){
  de.genes[[i]] <- read.csv(
    paste0("results/04_de-genes-by-celltype/logfc_0.40/de/files/de-genes_", i, ".csv")
  )
  de.genes[[i]] <- dplyr::rename(de.genes[[i]], "gene" = "X")
  de.genes[[i]]$celltype <- i
}

### create metascape input file ===============================================
meta <- as.data.frame(
  matrix(
    nrow = lapply(de.genes, nrow) %>% unlist() %>% max(),
    ncol = length(celltypes)
  )
)

names(meta) <- celltypes

for(i in celltypes) {
  meta[[i]] <- c(de.genes[[i]][["gene"]],
                 rep("", nrow(meta) - nrow(de.genes[[i]])))
}

write.csv(meta, "results/07_metascape/metascape_input-file.csv", row.names = FALSE)

# subset of celltypes without B cells, granulocytes etc. 
meta %>% 
  select(-c("Bcell", "Gran", "Tcell_3", "vil.Ery")) %>% 
  write.csv(., 
            "results/07_metascape/metascape_input-file_fewer-celltypes.csv",
            row.names = FALSE)

# subset of celltypes: only immune cells 
meta %>% 
  select("APC", "Mono_1", "Mono_2", "NK_1", "NK_2", "NK_3", "Tcell_1", "Tcell_2", "vil.Hofb") %>% 
  write.csv(.,
            "results/07_metascape/metascape_input-file_immune-celltypes.csv",
            row.names = FALSE)

# subset of celltypes: only innate immune celltypes
meta %>% 
  select("APC", "Mono_1", "Mono_2", "NK_1", "NK_2", "NK_3", "vil.Hofb") %>% 
  write.csv(.,
            "results/07_metascape/metascape_input-file_innate-immune-celltypes.csv",
            row.names = FALSE)

### end  ======================================================================
