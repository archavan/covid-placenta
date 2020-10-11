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

### heatmap ===================================================================
ann <- unique(seur@meta.data[, c("seurat_clusters", "annotation", "annotation_merged")])

plot.dat <- cormat %>% 
  filter(Var1_source %in% c("pavli", "surya", "vento"), 
         Var2_source %in% c("clust"))

plot.dat$lab <- ann$annotation_merged[match(x = plot.dat$Var2, 
                                            table = ann$seurat_clusters)]
plot.dat$lab[!(plot.dat$Var1 %in% c("vento_HB"))] <- NA # this is a lazy hack to make geom_text add label only once. If this was filtered to Var1_source == vento, it would add 35 (number of vento clusters) labels on top of each other.

facet.names <- c(
  pavli = "Pavlicev et al",
  surya = "Suryawanshi et al",
  vento = "Vento-Tormo et al"
)

cor.plot <- ggplot(plot.dat, aes(Var1, Var2, fill = value)) +
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
  facet_grid(. ~ Var1_source, scales = "free_x", space = "free_x",
             labeller = as_labeller(facet.names)) +
  coord_cartesian(clip = "off") +
  geom_text(aes(label = lab), x = 33, size = 5.25/.pt, 
            hjust = 0, color = "black") +
  labs(x = "Reference", y = "Query", tag = "A") +
  theme_bw() +
  theme(
    aspect.ratio = 35,
    plot.tag = element_text(size = 8, colour = "black", face = 2),
    strip.text = element_text(size = 7, color = "black"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, size = 5.25, 
                               hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 5.25, color = "black"),
    axis.ticks = element_line(size = 0.25),
    axis.ticks.length = unit(0.15, units = c('lines')),
    axis.title = element_text(size = 7),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.box.spacing = unit(2.5, "lines"),
    legend.key.size = unit(0.65, units = c('lines')), 
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6)
  )

### UMAP figure, optional =====================================================
theme.umap <-   theme_classic() +
  theme(axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 2),
        legend.position = "none",
        aspect.ratio = 1)

# all annotations
Idents(seur) <- "annotation"
umap.full <- DimPlot(seur, label = TRUE, label.size = 6/.pt, shuffle = TRUE, repel = TRUE) + 
  labs(tag = "B") +
  theme.umap
umap.full$layers[[2]]$geom_params$segment.size <- 0.25

# merged annotations
Idents(seur) <- "annotation_merged"
umap.merged <- DimPlot(seur, label = TRUE, label.size = 6/.pt, shuffle = TRUE, repel = TRUE) + 
  labs(tag = "C") +
  theme.umap
umap.merged$layers[[2]]$geom_params$segment.size <- 0.25

# test umap plot with custom colours used in the stacked bar plot so that the colours match between plots
stack.colors <- c("#6cb9a0",
                  "#a90090",
                  "#2abb3e",
                  "#4140c3",
                  "#a6d468",
                  "#ff8dee",
                  "#00b56e",
                  "#b70048",
                  "#39ddbe",
                  "#e8591a",
                  "#02aded",
                  "#ba5600",
                  "#68d6e3",
                  "#ff6db4",
                  "#00ba9d",
                  "#ff8f94",
                  "#4e5719",
                  "#ffaf52",
                  "#b27b5e",
                  "#e2c36e",
                  "#836400") # colors from https://medialab.github.io/iwanthue/
Idents(seur) <- "annotation_merged"
umap.merged1 <- DimPlot(seur, label = TRUE, label.size = 6/.pt, shuffle = TRUE, repel = TRUE, cols = stack.colors) + 
  labs(tag = "C") +
  theme.umap
umap.merged1$layers[[2]]$geom_params$segment.size <- 0.25
ggsave(
  filename = "results/99_paper-figures/supp-fig_annotation/umap_colour-test.png",
  umap.merged1, width = 3.5, height = 3.5, units = "in", type = "cairo"
)

### arrange ===================================================================
layout <- c(area(t = 1, l = 1, b = 45, r = 70),
            area(t = 46, l = 1, b = 80, r = 35),
            area(t = 46, l = 36, b = 80, r = 70))

composite <- cor.plot + umap.full + umap.merged + 
  plot_layout(design = layout, guides = "keep")
  
cowplot::ggsave2(
  filename = "results/99_paper-figures/supp-fig_annotation/supp-fig_annotation.pdf",
  composite, 
  device = "pdf", width = 7, height = 7, units = "in"
)

cowplot::ggsave2(
  filename = "results/99_paper-figures/supp-fig_annotation/supp-fig_annotation.png",
  composite, 
  width = 7, height = 7, units = "in", type = "cairo", dpi = 600
)

### write associated supplementary files ======================================
## cluster annotations --------------------------------------------------------
ann <- ann[order(ann$seurat_clusters), ]
write.csv(ann, "results/99_paper-supp-files/cluster_annotations.csv", row.names = FALSE)

## Cluster marker genes -------------------------------------------------------
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

## Number of cells of each celltype by control and covid status ---------------
cellnum <- table(seur$annotation_merged, seur$covid) %>% as.data.frame()
cellnum <- pivot_wider(data = cellnum, id_cols = "Var1", names_from = "Var2", values_from = "Freq")
cellnum <- dplyr::rename(cellnum, "celltype" = "Var1")

write.csv(cellnum, "results/99_paper-supp-files/n-cells_by-celltype-and-covid-status.csv", row.names = FALSE)

### end =======================================================================




