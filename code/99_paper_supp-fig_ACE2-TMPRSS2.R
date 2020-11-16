# ACE2 and TMPRSS2 plots

# packages
library(tidyverse)
library(Seurat)
library(patchwork)

# data
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")

### violin plots ===============================================================
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

VlnPlot2 <- function(object = seur, features) {
  
  Idents(object) <- object$annotation_merged
  object$annotation_merged <- factor(object$annotation_merged, 
                                     levels = sort(unique(object$annotation_merged)))
  
  p <- VlnPlot(object = object, 
               features = features, 
               assay = "SCT", 
               pt.size = 0, 
               cols = stack.colors) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 6, color = "black"),
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
  
  p$layers[[1]]$aes_params$size = 0.20 # violin stroke
  
  return(p)
  
}

vln.ace2 <- VlnPlot2(features = c("ACE2"))
vln.tmprss2 <- VlnPlot2(features = c("TMPRSS2"))

vln <- vln.ace2 + vln.tmprss2 + plot_layout(ncol = 2)

ggsave(
  filename = "results/99_paper-figures/supp-fig_ACE2-TMPRSS2/supp-fig_ACE2-TMPRSS2_violin.jpeg",
  vln,
  width = 5, height = 2, units = "in",
  type = "cairo", dpi = 600
)

ggsave(
  filename = "results/99_paper-figures/supp-fig_ACE2-TMPRSS2/supp-fig_ACE2-TMPRSS2_violin.pdf",
  vln,
  width = 5, height = 2, units = "in"
)

### Feature plot ==============================================================

FeaturePlot2 <- function(object = seur, features,
                         min.cutoff = "q1", max.cutoff = "q99") {
  DefaultAssay(object) <- "SCT"
  p <- FeaturePlot(object, 
                   feature = features, 
                   pt.size = 0.01, 
                   order = TRUE,
                   min.cutoff = min.cutoff, 
                   max.cutoff = max.cutoff) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 7, color = "black"),
      legend.key.size = unit(0.4, "line"),
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
  
  return(p)
}

f.ace2 <- FeaturePlot2(features = "ACE2")
f.tmprss2 <- FeaturePlot2(features = "TMPRSS2")

fplot <- f.ace2 + f.tmprss2 + plot_layout(ncol = 2)

ggsave(
  filename = "results/99_paper-figures/supp-fig_ACE2-TMPRSS2/supp-fig_ACE2-TMPRSS2_feature-plot.jpeg",
  fplot,
  width = 4.5, height = 2.5, units = "in",
  type = "cairo", dpi = 600
)

ggsave(
  filename = "results/99_paper-figures/supp-fig_ACE2-TMPRSS2/supp-fig_ACE2-TMPRSS2_feature-plot.pdf",
  fplot,
  width = 4.5, height = 2.5, units = "in"
)

### end =======================================================================
