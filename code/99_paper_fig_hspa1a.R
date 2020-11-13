# HSPA1A violin plot

# packages
library(tidyverse)
library(Seurat)

# data
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")

# plot, without splitting by covid status
seur$annotation_merged <- factor(seur$annotation_merged, 
                                 levels = sort(unique(seur$annotation_merged)))
Idents(seur) <- seur$annotation_merged

p <- VlnPlot(object = seur, 
             features = c("HSPA1A"), 
             assay = "SCT", 
             pt.size = 0, 
             cols = rep("grey70", length(unique(seur$annotation_merged)))) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                               size = 5.25, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 5.25, color = "black"),
    axis.text.y = element_text(size = 5.25, color = "black"),
    legend.position = "none"
  )

p$layers[[1]]$aes_params$size = 0.25 # violin stroke

ggsave(
  filename = "results/99_paper-figures/fig_hspa1a/vln_HSPA1A_v1.pdf",
  p,
  width = 3, height = 2, units = "in"
)
