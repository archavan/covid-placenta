# Make supplementary figures related to sample QC

### packages ==================================================================
library(tidyverse)

library(Seurat) # single cell

library(ggrepel) # plotting
library(patchwork)
library(ggthemes)

### data ======================================================================
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")

### individual plots ==========================================================
## uUMI per sample ------------------------------------------------------------
n.umi <- ggplot(data = seur@meta.data, aes(x = orig.ident, y = nCount_RNA)) +
  geom_violin(aes(fill = covid)) +
  scale_fill_tableau(palette = "Classic 10 Medium", name = "Sample") +
  labs(tag = "A", x = "Sample", y = "nUMI") 

## nGene per sample -----------------------------------------------------------
n.gene <- ggplot(data = seur@meta.data, aes(x = orig.ident, y = nFeature_RNA)) +
  geom_violin(aes(fill = covid)) +
  scale_fill_tableau(palette = "Classic 10 Medium", name = "Sample") +
  labs(tag = "B", x = "Sample", y = "nGene") 

## percent mt reads -----------------------------------------------------------
pct.mt <- ggplot(data = seur@meta.data, aes(x = orig.ident, y = percent.mt)) +
  geom_violin(aes(fill = covid)) +
  scale_fill_tableau(palette = "Classic 10 Medium", name = "Sample") +
  labs(tag = "C", x = "Sample", y = "% MT UMI") 

## umap plottig function ------------------------------------------------------
## plotting function
plot_umap <- function(object, colourby, ...){
  
  Idents(object) <- colourby
  ggobject <- DimPlot(object, label = FALSE, shuffle = TRUE, ...)
  udat <- ggobject$data # point coordinates
  
  if(colourby == "covid"){
    udat$ident <- factor(udat$ident, levels = c("cntrl", "covid"))
  }

  # plot
  p <- ggplot() +
    geom_point(data = udat, aes(x = UMAP_1, y = UMAP_2, color = ident),
               size = 0.1) +
    labs(x = "", y = "") +
    scale_color_tableau(palette = "Classic 10 Medium", 
                        name = "Sample",
                        guide = guide_legend(override.aes = list(size = 1))) +
    annotate(geom = "text", x = Inf, y = -Inf, hjust = 1, vjust = -0.5,
             label = "UMAP", size = 6/.pt) +
    theme(aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.text = element_blank())

  return(p)
  
}

## umap coloured by sample ----------------------------------------------------
umap.sample <- plot_umap(object = seur, colourby = "orig.ident") + labs(tag = "E")

## colored by covid -----------------------------------------------------------
umap.covid <- plot_umap(object = seur, colourby = "covid") + labs(tag = "G")

## stacked bar: cells/celltype/sample -----------------------------------------
cell.sample <- table(seur@meta.data$annotation_merged, seur@meta.data$orig.ident) %>% 
  as.data.frame()

cell.sample$frac <- NA
for(i in 1:nrow(cell.sample)) {
  cell.sample$frac[i] <- cell.sample$Freq[i] / 
    sum(cell.sample$Freq[cell.sample$Var1 == cell.sample$Var1[i]])
}

bar.sample <- ggplot(data = cell.sample, aes(x = Var1, y = frac, fill = Var2)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_tableau("Classic 10 Medium", name = "Sample") +
  labs(x = "Cell type", y = "Fraction of cells", tag = "D")

## stacked bar: cells/celltype/covid -----------------------------------------
cell.covid <- table(seur@meta.data$annotation_merged, seur@meta.data$covid) %>% 
  as.data.frame()

cell.covid$frac <- NA
for(i in 1:nrow(cell.covid)) {
  cell.covid$frac[i] <- cell.covid$Freq[i] / 
    sum(cell.covid$Freq[cell.covid$Var1 == cell.covid$Var1[i]])
}

bar.covid <- ggplot(data = cell.covid, aes(x = Var1, y = frac, fill = Var2)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_tableau("Classic 10 Medium", name = "Covid") +
  labs(x = "Cell type", y = "Fraction of cells", tag = "F")
  
### arrange ===================================================================
row1 <- n.umi + n.gene + pct.mt + plot_layout(guides = "collect")
row2 <- bar.sample + umap.sample + plot_layout(widths = c(2, 1))
row3 <- bar.covid + umap.covid + plot_layout(widths = c(2, 1))

composite <- row1 / row2 / row3 &
  theme_classic() +
  theme(
    plot.tag = element_text(size = 9, face = 2),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    axis.title = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 7, color = "black",
                               angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.75, "lines")
  )

cowplot::ggsave2(
  composite, 
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC.pdf",
  width = 6.5, height = 6.5, units = "in"
)

cowplot::ggsave2(
  composite, 
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC.png",
  width = 6.5, height = 6.5, units = "in", type = "cairo", dpi = 600
)

