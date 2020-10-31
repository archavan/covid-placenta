# Supplementary figure about cell type annotation of sc-RNA-seq clusters.

### packages ==================================================================
library(tidyverse)
library(Seurat)

library(patchwork) # plotting
library(cowplot)
library(ggrepel)

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

# calculate top 3 matches
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

plot.dat$lab <- ann$annotation[match(x = plot.dat$Var2, 
                                     table = ann$seurat_clusters)]
plot.dat$lab[!(plot.dat$Var1 %in% c("vento_HB"))] <- NA # this is a lazy hack to make geom_text add label only once. If this was filtered to Var1_source == vento, it would add 35 (number of vento clusters) labels on top of each other.

facet.names <- c(
  pavli = "Pavlicev et al",
  surya = "Suryawanshi et al",
  vento = "Vento-Tormo et al"
)

cor.plot <- ggplot(plot.dat, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = 'white', high = '#d94801',
                      name = "Spearman\nCorrelation",
                      limits = c(0.50, 0.95),
                      oob = scales::squish,
                      breaks = c(0.50, 0.7, 0.9),
                      labels = c("<0.5", "0.7", "0.9")) +
  geom_point(aes(Var1, Var2, alpha = top3),
             size = 1, shape = 19, stroke  = 0) +
  scale_alpha_manual(values = c(1, 0.5, 0.25), 
                     breaks = c(1, 2, 3),
                     name = "Top 3 match", na.value = 0) +
  guides(alpha = guide_legend(override.aes = list(size = 2)),
         fill = guide_colorbar(ticks.colour = "black")) +
  facet_grid(. ~ Var1_source, 
             scales = "free_x", 
             space = "free_x",
             labeller = as_labeller(facet.names)) +
  coord_cartesian(clip = "off") +
  geom_text(aes(label = lab), x = 33, size = 5.25/.pt, 
            hjust = 0, color = "black") +
  labs(x = "Reference", y = "Query", tag = "A") +
  theme_bw() +
  theme(
    aspect.ratio = 35,
    plot.tag = element_text(size = 8, colour = "black", face = 2),
    strip.text = element_text(size = 6, color = "black"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, size = 5, 
                               hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.ticks = element_line(size = 0.25),
    axis.ticks.length = unit(0.15, units = c('lines')),
    axis.title = element_text(size = 6),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.25, color = "black", fill = NA),
    legend.box.spacing = unit(2.5, "lines"),
    legend.key.size = unit(0.5, units = c('lines')), 
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )

### UMAP figures ==============================================================
## colors for plotting --------------------------------------------------------
stack.colors.20 <- c("#6cb9a0",
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

stack.colors.35 <- c("#b8bd38",
                     "#a454cc",
                     "#67be40",
                     "#cc3fa5",
                     "#4fc16f",
                     "#dd74d2",
                     "#4c8e2c",
                     "#5b69d8",
                     "#e0a736",
                     "#824b9e",
                     "#91b143",
                     "#b786dd",
                     "#888c2d",
                     "#d9417f",
                     "#5fc7a8",
                     "#d4403f",
                     "#48b1da",
                     "#e1722b",
                     "#4e65a6",
                     "#ac7c21",
                     "#8396de",
                     "#ad522d",
                     "#339580",
                     "#a05084",
                     "#4b9054",
                     "#dc88b5",
                     "#2a6a45",
                     "#d97578",
                     "#566926",
                     "#a14456",
                     "#a5b670",
                     "#e6926c",
                     "#7e7336",
                     "#cea864",
                     "#956534")

## function -------------------------------------------------------------------
plot_umap <- function(object, idents, colors, tag) {
  
  # set idents
  Idents(object) <- idents
  
  # extract plotting data
  gg <- DimPlot(object, label = TRUE, shuffle = TRUE)
  udat <- gg$data # point coordinates
  ulab <- gg$layers[[2]]$data # label coordinates
  
  # new colors
  ulab <- ulab[order(as.character(ulab$ident)), ]
  ulab$newcolors <- colors
  
  # factors
  udat$ident <- factor(udat$ident, levels = ulab$ident)
  
  # plot
  p <- ggplot() +
    geom_point(data = udat, 
               aes(x = UMAP_1, y = UMAP_2, color = ident),
               size = 0.1) +
    scale_color_manual(values = alpha(ulab$newcolors, 1), 
                       breaks = ulab$ident) +
    labs(x = "UMAP 1", y = "UMAP 2", tag = tag) +
    geom_label_repel(data = ulab, 
                     aes(x = UMAP_1, y = UMAP_2, label = ident),
                     size = 5/.pt,
                     segment.size = 0.2, segment.color = "black",
                     force = 1,
                     fontface = 1, 
                     label.padding = unit(0.1, "lines"), 
                     label.r = 0,
                     label.size = NA, # https://stackoverflow.com/questions/43417514/getting-rid-of-border-in-pdf-output-for-geom-label-for-ggplot2-in-r
                     seed = 0, 
                     fill = alpha(c("white"), 0.7)) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      plot.tag = element_text(face = 2, size = 8, color = "black"),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_text(size = 5),
      axis.line = element_line(size = 0.25),
      axis.ticks = element_blank()
    )

}

## plot: all annotations ------------------------------------------------------
umap.full <- plot_umap(object = seur, idents = "annotation", colors = stack.colors.35, tag = "B")

## plot: merged annotations ---------------------------------------------------
umap.merged <- plot_umap(object = seur, idents = "annotation_merged", colors = stack.colors.20, tag = "C")

### arrange ===================================================================
layout <- c(area(t = 1, l = 1, b = 45, r = 70),
            area(t = 46, l = 1, b = 80, r = 35),
            area(t = 46, l = 36, b = 80, r = 70))

composite <- cor.plot + umap.full + umap.merged + 
  plot_layout(design = layout, guides = "keep")
  
cowplot::ggsave2(
  filename = "results/99_paper-figures/supp-fig_annotation/supp-fig_annotation_composite.pdf",
  composite, 
  device = "pdf", width = 7, height = 7, units = "in"
)

cowplot::ggsave2(
  filename = "results/99_paper-figures/supp-fig_annotation/supp-fig_annotation_composite.png",
  composite, 
  width = 7, height = 7, units = "in", type = "cairo", dpi = 600
)

### end  ======================================================================
