# HSPA1A violin plot

# packages
library(tidyverse)
library(Seurat)

# data
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")

### violin plot ===============================================================
# colours to match umap
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

# plot, without splitting by covid status
seur$annotation_merged <- factor(seur$annotation_merged, 
                                 levels = sort(unique(seur$annotation_merged)))
Idents(seur) <- seur$annotation_merged

p <- VlnPlot(object = seur, 
             features = c("HSPA1A"), 
             assay = "SCT", 
             pt.size = 0, 
             cols = stack.colors) +
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
  filename = "results/99_paper-figures/fig_hspa1a/violinplot_HSPA1A_v1.pdf",
  p,
  width = 3, height = 2, units = "in"
)

### Feature plot ==============================================================

DefaultAssay(seur) <- "SCT"
fplot <- FeaturePlot(seur, 
                     feature = c("HSPA1A"), 
                     pt.size = 0.01, 
                     order = TRUE,
                     min.cutoff = "q1", max.cutoff = "q99") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    legend.key.size = unit(0.25, "line"),
    legend.title = element_text(size = 5.25),
    legend.text = element_text(size = 5),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 5),
    axis.title.y.right = element_blank(),
    axis.text = element_blank(),
    aspect.ratio = 1
  )

ggsave(
  filename = "results/99_paper-figures/fig_hspa1a/featureplot_HSPA1A_v1.pdf",
  fplot,
  width = 2, height = 2, units = "in"
)

### end =======================================================================
