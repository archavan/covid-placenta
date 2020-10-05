# Repeat GO enrichment for DE genes from bulk RNA-seq with clusterProfiler to be consistent in terms of pacakges used between sc-RNA-seq and bulk. 
# 2020-10-04

### packages ==================================================================
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

### data ======================================================================
up <- readxl::read_excel(
  "info/from-alice/for-go/2020.10.04_bulkseqhits-forGO.xlsx", 
  sheet = "upregulated (COVID vs. ctrl)", 
  col_names = TRUE, trim_ws = TRUE
)

down <- readxl::read_excel(
  "info/from-alice/for-go/2020.10.04_bulkseqhits-forGO.xlsx",
  sheet = "downregulated (COVID vs. ctrl)",
  col_names = TRUE, trim_ws = TRUE
)

bg <- readxl::read_excel(
  "info/from-alice/for-go/2020.10.04_bulkseqhits-forGO.xlsx",
  sheet = "background list",
  col_names = TRUE, trim_ws = TRUE
)

### GO enrichment =============================================================
## function --------------------------------------------------------------------
enrichGO2 <- function(target) {
  
  # coerce to character
  target <- as.character(target)
  
  # convert ids
  target.ids <- bitr(target,
                     fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = "org.Hs.eg.db")
  
  # enrichment
  ego <- list()

  for(i in c("BP", "MF", "CC")) {
    ego[[i]] <-  enrichGO(gene          = target.ids$ENTREZID,
                          #universe      = bg.ids$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          ont           = i,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01, # default 0.01
                          qvalueCutoff  = 0.05, # default 0.05
                          readable = TRUE)
  }

  return(ego)
  
}

## enrichment -----------------------------------------------------------------
go.up <-   enrichGO2(target = up$ens_gene)
go.down <- enrichGO2(target = down$ens_gene)

## write csv ------------------------------------------------------------------
for(i in c("up", "down")) {
  for(j in c("BP", "MF", "CC")) {
    write.csv(
      get(paste0("go.", i))[[j]]@result,
      paste0("results/05_go-enrichment_bulk-rna-seq/bulk_GO-enrichment_", 
             i, "_", j, ".csv"),
      row.names = FALSE
    )
  }
}

### dotplots ------------------------------------------------------------------
for(i in c("up", "down")) {
  for(j in c("BP", "MF", "CC")) {
    
    p <- dotplot(get(paste0("go.", i))[[j]]) + 
      ggtitle(paste0(i, '_', j)) +
      theme(aspect.ratio = 1)
    
    ggsave(
      p, 
      filename = paste0("results/05_go-enrichment_bulk-rna-seq/bulk_GO-dotplot_", 
                        i, '_', j, '.pdf'), 
      width = 8, height = 6, units = 'in'
    )
  }
}

### end  ======================================================================
