# Make supplementary figures related to sample QC

### packages ==================================================================
library(tidyverse)

library(Seurat) # single cell

library(ggrepel) # plotting
library(patchwork)
library(ggthemes)
library(cowplot)

### data ======================================================================
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")

### individual plots ==========================================================
## themes
theme2 <-   theme_classic() +
  theme(
    plot.tag = element_text(size = 8, face = 2),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    axis.title = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 6, color = "black",
                               angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 6, color = "black"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.75, "lines")
  )

theme_umap <- theme_classic() +
  theme(aspect.ratio = 1,
        plot.tag = element_text(size = 8, face = 2),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(size = 0.25, fill = NA))

## uUMI per sample ------------------------------------------------------------
n.umi <- ggplot(data = seur@meta.data, 
                aes(x = orig.ident_renamed, y = nCount_RNA)) +
  geom_violin(aes(fill = covid)) +
  scale_fill_tableau(palette = "Classic 10 Medium", name = "Sample") +
  labs(tag = "A", x = "Sample", y = "nUMI") +
  theme2

## nGene per sample -----------------------------------------------------------
n.gene <- ggplot(data = seur@meta.data, 
                 aes(x = orig.ident_renamed, y = nFeature_RNA)) +
  geom_violin(aes(fill = covid)) +
  scale_fill_tableau(palette = "Classic 10 Medium", name = "Sample") +
  labs(tag = "B", x = "Sample", y = "nGene") +
  theme2

## stacked bar: cells/celltype/sample -----------------------------------------
cell.sample <- table(seur@meta.data$orig.ident_renamed,
                     seur@meta.data$annotation_merged) %>% 
  as.data.frame()

cell.sample$frac <- NA
for(i in 1:nrow(cell.sample)) {
  cell.sample$frac[i] <- cell.sample$Freq[i] / 
    sum(cell.sample$Freq[cell.sample$Var1 == cell.sample$Var1[i]])
}

bar.sample <- ggplot(data = cell.sample, 
                     aes(x = Var1, y = frac, fill = Var2)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(name = "Cell type", 
                    values = c("#6cb9a0",
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
                               "#836400")) + # colors from https://medialab.github.io/iwanthue/
  labs(x = "Sample", y = "Fraction of cells", tag = "C") +
  theme2

## umap colored by sample -----------------------------------------------------
Idents(seur) <- "orig.ident_renamed"
ggobject <- DimPlot(seur, label = FALSE, shuffle = TRUE) 

umap.sample <- ggplot() +
  geom_point(data = ggobject$data, aes(x = UMAP_1, y = UMAP_2, color = ident),
             size = 0.1) +
  labs(x = "UMAP 1", y = "UMAP 2", tag = "D") +
  scale_color_discrete(name = "Sample",
                       guide = guide_legend(override.aes = list(size = 1))) +
  theme_umap

## UMAP split by covid and control --------------------------------------------
Idents(seur) <- "covid"
ggobject <- DimPlot(seur, shuffle = TRUE, label = FALSE)

umap.cntrl <- ggplot(
  data = ggobject$data[ggobject$data$ident == "cntrl", ],
  aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, color = "grey") +
  annotate(geom = "text", label = "Control", 
           x = Inf, y = -Inf, hjust = 1.25, vjust = -0.5, size = 8/.pt) +
  labs(tag = "E") +
  theme_umap

umap.covid <- ggplot(
  data = ggobject$data[ggobject$data$ident == "covid", ],
  aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, color = "grey") +
  annotate(geom = "text", label = "Covid", 
           x = Inf, y = -Inf, hjust = 1.25, vjust = -0.5, size = 8/.pt) +
  theme_umap

### arrange ===================================================================
row1 <- n.umi + n.gene + plot_layout(guides = "collect")
row2 <- bar.sample 
row3 <- umap.sample + umap.cntrl + umap.covid + plot_layout(widths = c(1, 1, 1))

composite <- row1 / row2 / row3 +
  plot_layout(heights = c(2, 3, 1.75))

cowplot::ggsave2(
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC.pdf",
  composite, 
  width = 6.5, height = 6, units = "in"
)

cowplot::ggsave2(
  composite, 
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC.png",
  width = 6.5, height = 6, units = "in", type = "cairo", dpi = 600
)

