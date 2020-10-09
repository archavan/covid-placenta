# Supplementary figure about cell type annotation of sc-RNA-seq clusters.

### packages ==================================================================
library(tidyverse)
library(Seurat)

library(patchwork) # plotting
library(cowplot)

### data ======================================================================
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")

# average by cluster
Idents(seur) <- "seurat_clusters"
plac <- AverageExpression(seur) %>% .[["SCT"]]

### reference data ============================================================
# data from pavlicev et al, suryavanshi et al, and vento-tormo et al, averaged by cell type, and compiled into one dataframe.
ref <- read.csv("results/01_reference-atlas/vento-surya-pavli_joined.csv", stringsAsFactors = FALSE)

### get data in shape =========================================================
plac$gene_name <- rownames(plac)
plac <- dplyr::inner_join(plac, ref, by = "gene_name")
plac <- dplyr::relocate(plac, gene_name)

### annotation by correlation =================================================
refdata <- c("vento", "surya", "pavli")

# build correlation matrix from expression data
cor.matrix <- cor(plac[, names(plac)[names(plac) != "gene_name"]], 
                  method = 'spearman')

# reorder correlation matrix based on clustering
dd <- as.dist((1 - cor.matrix))
hc <- hclust(dd, method = 'complete')
cor.matrix <- cor.matrix[hc$order, hc$order]

# melt correlation matrix
cormat <- reshape2::melt(cor.matrix, na.rm = T)
cormat$Var1_source <- sapply(strsplit(as.character(cormat$Var1), split = "_"), "[[", 1)
cormat$Var2_source <- sapply(strsplit(as.character(cormat$Var2), split = "_"), "[[", 1)

# subset
# cormat <- cormat[which(cormat$Var1_source %in% refdata & 
#                          cormat$Var2_source == "clust"), ]

cormat$top3 <- NA

for(i in unique(cormat$Var2)){
  for(j in refdata){
    # identify top3 match indices
    ind <- which(cormat$Var2 == i & cormat$Var1_source == j)
    val <- cormat$value[ind]
    top3ind <- ind[order(val, decreasing = TRUE)[1:3]]
    
    # assign match ranking to correlation data for plotting
    cormat$top3[top3ind[1]] <- "1"
    cormat$top3[top3ind[2]] <- "2"
    cormat$top3[top3ind[3]] <- "3"
  }
}

### heatmap plotting function =================================================
clustAnnoPlot <- function(dat, query, reference){
  
  dat <- dat[which(dat$Var2_source == query & dat$Var1_source == reference), ]
  
  p <- ggplot(data = dat, 
              aes(Var1, Var2, fill = value)) +
    geom_tile(colour = "white") +
    scale_fill_gradient(low = 'white', high = 'red',
                        name = "Spearman\nCorrelation",
                        limits = c(0.40, 0.95),
                        breaks = c(0.40, 0.60, 0.80),
                        labels = c("0.4", "0.6", "0.8")) +
    geom_point(aes(Var1, Var2, alpha = top3),
               size = 1, shape = 19, stroke  = 0) +
    scale_alpha_manual(values = c(1, 0.5, 0.25), 
                       breaks = c(1, 2, 3),
                       name = "Top 3 match", na.value = 0) +
    coord_fixed(ratio = 1) +
    xlab("reference") +
    ylab("query") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, size = 5.5, 
                                 hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(size = 5.5, color = "black"),
      axis.ticks = element_line(size = 0.25),
      axis.ticks.length = unit(0.15, units = c('lines')),
      axis.title = element_text(size = 7),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.key.size = unit(0.8, units = c('lines')), 
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.tag = element_text(size = 8, face = 2, color = "black"),
      plot.margin = margin(0, 0, 0, 0) # helps alignment by panel area
    )
  
  return(p)
}

### individual plots ==========================================================
surya.vento <- clustAnnoPlot(dat = cormat, query = "surya", reference = "vento") + labs(tag = "A") # between vento and surya

by.pavli <- clustAnnoPlot(dat = cormat, query = "clust", reference = "pavli") + labs(tag = "B")
by.vento <- clustAnnoPlot(dat = cormat, query = "clust", reference = "vento") + labs(tag = "C")
by.surya <- clustAnnoPlot(dat = cormat, query = "clust", reference = "surya") + labs(tag = "D")


### arrange ===================================================================
layout <- c(area(t = 1,  l = 1,  b = 24, r = 32),
            area(t = 1,  l = 33, b = 27, r = 32+8),
            area(t = 3,  l = 41, b = 22, r = 50),
            area(t = 31, l = 1,  b = 70, r = 32),
            area(t = 31, l = 33, b = 70, r = 32+22))

composite <- surya.vento  + by.pavli + guide_area() + by.vento + by.surya + 
  plot_layout(guides = "collect", design = layout)

cowplot::ggsave2(
  composite, device = "pdf", width = 6.75, height = 8.25, units = "in",
  filename = "results/99_paper-figures/supp-fig_annotation/supp-fig_annotation.pdf"
)



### end =======================================================================



