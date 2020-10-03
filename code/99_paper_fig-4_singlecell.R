# 2020-10-01
# Arun Chavan

### packages ==================================================================
library(tidyverse)

# single cell
library(Seurat)

# plotting
library(patchwork)
library(ggthemes)
library(ggrepel)
library(ggforce) # for sina plots
library(cowplot)
library(ggtext)
library(ungeviz) # for geom_hpline

### data ======================================================================
## Seurat
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")
celltypes <- seur@meta.data$annotation_merged %>% unique()

## DE genes
de.genes <- list()
for(i in celltypes){
  de.genes[[i]] <- read.csv(
    paste0("results/04_de-genes-by-celltype/logfc_0.40/de/files/de-genes_", i, ".csv")
  )
  de.genes[[i]] <- dplyr::rename(de.genes[[i]], "gene" = "X")
  de.genes[[i]]$celltype <- i
}

## Interferome output
ifome.n <- read.csv("results/04_de-genes-by-celltype/logfc_0.40/interferome/files/number-of-DE-genes-in-interferome.csv")

## GO output
de.go <- list()
for(i in celltypes){
  de.go[[i]] <- read.csv(
    paste0("results/04_de-genes-by-celltype/logfc_0.40/go/files/enrichedGO_", i, ".csv")
  )
  de.go[[i]]$celltype <- i
}

### UMAP ======================================================================

### DE dotplot ================================================================

### Interferome plot ==========================================================
# color by pval
ifome.n$color <- ifelse(test = ifome.n$pval < 0.05, yes = "Black", no = "Grey40")

# tidy
ifome.n.long <- pivot_longer(
  data = ifome.n,
  cols = 3:6,
  names_to = c(".value", "interferome"), 
  names_pattern = "(.+)_(.+)"
  )
ifome.n.long$interferome <- gsub("interferome.", "", ifome.n.long$interferome)

ifome.n.long$celltype <- factor(
  ifome.n.long$celltype, 
  levels = ifome.n$celltype[order(ifome.n$pct_interferome.yes)]
)

# plot
plot.dat <- ifome.n.long
p <- ggplot(data = plot.dat, 
            aes(y = celltype, x = pct)) +
  geom_bar(aes(fill = interferome),
           position = "stack", 
           stat = "identity", 
           width = 0.8) +
  scale_fill_manual(values = c("grey90", "#4e79a7"), 
                    breaks = c("no", "yes"),
                    name = "Interferome", 
                    labels = c("No", "Yes")) +
  scale_y_discrete(expand = c(0.05, 0)) +
  scale_x_continuous(expand = c(0.025, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "DE genes present in Interferome (%)") +
  geom_text(data = plot.dat[plot.dat$interferome == "yes", ], # number of genes
            aes(y = celltype, x = pct, 
                label = paste0(n, "/", n_de)),
            size = 1.8, vjust = 0.5, hjust = 0.0, nudge_x = 0.5) +
  geom_text(data = plot.dat[, c("celltype", "pval")] %>% unique(), # pvalues
            aes(y = celltype, x = 100, 
                label = formatC(pval, format = "e", digits = 0)),
            hjust = 0, size = 1.75, nudge_x = 1,
            color = unique(plot.dat[, c("celltype", "color")])$color) +
  annotate(geom = "richtext", x  = 101, y = 22, size = 2,
           label = "enrichment *p* value", hjust = 0,
           fill = NA, label.color = NA, # remove background and outline
           label.padding = grid::unit(0, "pt")) + # remove padding
  theme_classic() +
  theme(
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 6, 
                               color = "black"),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks.y = element_line(size = 0.25),
    axis.ticks.x = element_line(size = 0.25),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 3))

ggsave(p, filename = "results/99_paper-figures/interferome.pdf", 
       width = 3.5, height = 3, units = "in")

### GO plot ===================================================================

### Putting it all together ===================================================
