# Generate input file for metascape. The input is a single file where each column is a list of DE genes in a celltype. 

### packages ==================================================================
library(tidyverse)

### DE genes data =============================================================
celltypes <- c("dec.APC", "dec.Bcells", "dec.DSC", "dec.Endo", "dec.FB", "dec.Gran", "dec.Mono_1", "dec.Mono_2", "dec.NK_1", "dec.NK_2", "dec.NK_3", "dec.SMC", "dec.Tcell_1", "dec.Tcell_2", "dec.Tcell_3", "vil.Ery", "vil.EVT", "vil.FB", "vil.Hofb", "vil.SCT", "vil.VCT")

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
  select(-c("dec.Bcells", "dec.Gran", "dec.Tcell_3", "vil.Ery")) %>% 
  write.csv(., 
            "results/07_metascape/metascape_input-file_fewer-celltypes.csv",
            row.names = FALSE)

# subset of celltypes: only immune cells 
meta %>% 
  select("dec.APC", "dec.Mono_1", "dec.Mono_2", "dec.NK_1", "dec.NK_2", "dec.NK_3", "dec.Tcell_1", "dec.Tcell_2", "vil.Hofb") %>% 
  write.csv(.,
            "results/07_metascape/metascape_input-file_immune-celltypes.csv",
            row.names = FALSE)

# subset of celltypes: only innate immune celltypes
meta %>% 
  select("dec.APC", "dec.Mono_1", "dec.Mono_2", "dec.NK_1", "dec.NK_2", "dec.NK_3", "vil.Hofb") %>% 
  write.csv(.,
            "results/07_metascape/metascape_input-file_innate-immune-celltypes.csv",
            row.names = FALSE)

### end  ======================================================================
