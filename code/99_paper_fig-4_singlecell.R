# started: 2020-10-01
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

### Make individual plots =====================================================
fig4 <- list()  # list for saving panels

## Panel A: UMAP --------------------------------------------------------------
Idents(seur) <- seur@meta.data$annotation_merged

# plotting fn
plot_umap <- function(object, title = NA, label, 
                      split = FALSE, split.by, split.plot, ...){
  
  if(split == FALSE){
    ggobject <- DimPlot(object, label = TRUE, shuffle = TRUE, ...)
    
    udat <- ggobject$data # point coordinates
    ulab <- ggobject$layers[[2]]$data # label coordinates
  }
  
  if(split == TRUE){
    ggobject <- DimPlot(object, split.by = as.character(split.by),
                        label = TRUE, shuffle = TRUE, ...)
    
    udat <- ggobject$data[eval(
      expr(`$`(ggobject$data, !!split.by))
    ) == split.plot, ] # geom point coordinates
    
    ulab <- ggobject$layers[[2]]$data[eval(
      expr(`$`(ggobject$layers[[2]]$data, !!split.by))
    ) == split.plot, ] # label coordinates
  }  # see section 19.4 https://adv-r.hadley.nz/quasiquotation.html
  
  # plot
  p <- ggplot() +
    geom_point(data = udat, aes(x = UMAP_1, y = UMAP_2, color = ident),
               size = 0.1) +
    annotate(geom = "text", x = Inf, y = -Inf, label = title,
             hjust = 1, vjust = -1, size = 7/.pt, fontface = 2) +
    theme_classic() +
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
  if(label == TRUE)
    p <- p + 
    geom_text_repel(data = ulab, aes(x = UMAP_1, y = UMAP_2, label = ident),
                    size = 6/.pt)
  
  return(p)
  
}

## individual plots
umap.all <- plot_umap(object = seur, label = TRUE, title = "All cells")
umap.cntrl <- plot_umap(object = seur, label = FALSE, 
                        split = TRUE, split.by = "covid", split.plot = "cntrl",
                        title = "Control")
umap.covid <- plot_umap(object = seur, label = FALSE, 
                        split = TRUE, split.by = "covid", split.plot = "covid",
                        title = "Covid")

## arrange 
umap.v <- plot_grid(plotlist = list(umap.all, umap.cntrl, umap.covid), ncol = 1, axis = 'lr') # vertical
umap.h <- plot_grid(plotlist = list(umap.all, umap.cntrl, umap.covid), nrow = 1, axis = 'tb')

cowplot::ggsave2(
  umap.v, 
  filename = "results/99_paper-figures/fig4_single-cell/panelA_umap_composite_vertical.png", 
  width = 2, height = 6, units = "in", type = "cairo"
)

cowplot::ggsave2(
  umap.h, 
  filename = "results/99_paper-figures/fig4_single-cell/panelA_umap_composite_horizontal.png", 
  width = 6, height = 2, units = "in", type = "cairo"
)

## Panel B: DE dotplot --------------------------------------------------------
## create celltype_covid ident to use for seurat dotplot
seur$celltype_covid <- paste0(seur$annotation_merged, "_", seur$covid)

## theme for DE genes faceted dotplot
theme_dotplot <- theme(
  panel.background = element_blank(),
  text = element_text(color = "Black"),
  panel.grid = element_blank(),
  panel.border = element_rect(fill = NA, size = 0.15, color = "grey60"),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "Black"),
  axis.text.y = element_text(size = 5.5, color = "Black"),
  axis.title = element_blank(),
  axis.ticks = element_line(size = 0.1, color = "Black"),
  strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, 
                              size = 6, color = "Black"),
  strip.background = element_blank(),
  legend.title = element_text(size = 6),
  legend.position = "right",
  legend.box.spacing = unit(0.1, "lines"),
  legend.background = element_blank(),
  legend.title.align = 0,
  legend.text = element_text(size = 5),
  legend.key.size = unit(0.5, "lines"),
  legend.key = element_blank(),
  legend.spacing.y = unit(0.1, "lines"),
  panel.spacing.x = unit(0, "lines")
)

## dotplot function
plot_splitdot <- function(object, features, exclude.celltypes = c()) {
  
  # idents combined for celltype and covid
  Idents(object) <- object$celltype_covid
  
  # run seurat's dotplot
  p <- Seurat::DotPlot(
    object = object, 
    assay = "SCT", 
    features = features
    )
  
  
  # extract plotting data
  p.dat <- p$data
  p.dat$celltype <- gsub(p.dat$id, pattern = "(.+)_([a-z]{5})", replacement = "\\1")
  p.dat$covid <- gsub(p.dat$id, pattern = "(.+)_([a-z]{5})", replacement = "\\2")
  
  # subset for celltypes of interest
  p.sub <- p.dat[!(p.dat$celltype %in% exclude.celltypes), ]
  
  # redo scaling for selected celltype_covid idents
  p.wide <- pivot_wider(p.sub,
                        id_cols = features.plot, 
                        names_from = id, 
                        values_from = avg.exp) %>% 
    as.data.frame()
  
  rownames(p.wide) <- p.wide$features.plot
  p.wide$features.plot <- NULL
  p.scaled <- t(p.wide) %>% 
    scale(center = TRUE, scale = TRUE) %>% 
    t() %>% 
    as.data.frame()
  p.scaled$features.plot <- rownames(p.scaled)
  
  # melt rescaled avg.exp
  p.scaled.long <- pivot_longer(
    data = p.scaled, 
    cols = names(p.scaled)[names(p.scaled) != "features.plot"], 
    names_to = "id", 
    values_to = "scaled.new")
  
  # cap min and max scaled value; same as col.min and col.max in Seurat::DotPlot
  p.scaled.long$scaled.new[p.scaled.long$scaled.new >  2.5] <-  2.5
  p.scaled.long$scaled.new[p.scaled.long$scaled.new < -2.5] <- -2.5
  
  # join with p.dat to get pct.exp and scaled.new in the same df
  plot.dat <- dplyr::inner_join(x = p.dat, 
                                y = p.scaled.long, 
                                by = c("features.plot", "id"))
  
  # clustering for ordering
  cdat <- pivot_wider(data =  plot.dat, id_cols = "features.plot", 
                      names_from = id,
                      values_from = scaled.new) %>% 
    as.data.frame()

  rownames(cdat) <- cdat$features.plot
  cdat$features.plot <- NULL
  cor.matrix <- cor(t(cdat))
  
  # reorder correlation matrix based on clustering
  dd <- as.dist((1 - cor.matrix))
  hc <- hclust(dd, method = 'complete')
  cdat  <- cdat[hc$order, ]

  # reorder factor levels based on clustering
  plot.dat$features.plot <- factor(plot.dat$features.plot, 
                                   levels = rownames(cdat))
  
  # prepare annotation data (hpline for DE genes)
  de.ann <- de.genes
  de.ann <- do.call(rbind, de.ann)
  de.ann <- subset(de.ann, 
                   gene %in% plot.dat$features.plot &
                     celltype %in% plot.dat$celltype
  )
  de.ann$diff.exp <- de.ann$p_val_adj < 0.05
  de.ann$diff.exp <- tolower(as.character(de.ann$diff.exp))
  
  # plot
  q <- ggplot(data = plot.dat, 
              aes(x = covid, y = features.plot)) +
    geom_point(aes(fill = scaled.new, size = pct.exp), 
               shape = 21, stroke = 0, 
               color = "White") +
    scale_fill_gradient(low = "White", high = "Red", 
                        name = "avg.exp\n(scaled)") +
    scale_x_discrete(breaks = c("cntrl", "covid"), 
                     labels = c("Control", "Covid")) +
    scale_size_area(max_size = 3.5,
                    guide = guide_legend(override.aes = list(
                      color = "White", 
                      fill = "Black"))) +
    facet_grid(. ~ celltype) +
    ungeviz::geom_hpline(data = de.ann,
                         aes(y = gene, x = 1.5, color = diff.exp),
                         size = 0.30, width = 1, lineend = "butt") +
    scale_color_manual(breaks = "true",
                       labels = "< 0.05",
                       values = "Black",
                       name = "DE adj.p") +
    theme_dotplot
  
  return(q)
  
}

## exclude celltypes lists
exclude1 <- c("dec.Tcell_1", "dec.Tcell_2", "dec.Tcell_3", "dec.Bcells", "dec.Gran", "vil.Ery")
exclude2 <- c("dec.Tcell_3", "dec.Bcells", "dec.Gran", "vil.Ery")

## select genes to plot: top x number of genes from clusters of interest
de.genes <- lapply(de.genes,  # sort de.genes by logfc
                   function(x){x[order(x$avg_logFC, decreasing = TRUE), ]})

# function
select_genes <- function(n, exclude.celltypes) {
  de.to.plot <- do.call(
    rbind, 
    lapply(de.genes, function(x){x[1:min(n, nrow(x)), ]})
  )
  de.to.plot <- de.to.plot[which(!(de.to.plot$celltype %in% exclude.celltypes)), ]
  
  return(de.to.plot$gene %>% unique())
}

# plot
numgenes <- 5
features <- select_genes(n = numgenes, exclude.celltypes = exclude2)
splitdot.test <- plot_splitdot(
  object = seur, 
  features = features,
  exclude.celltypes = exclude2
) +
  labs(caption = paste0(
    length(features), " genes; top ", numgenes, " per cluster"))
cowplot::ggsave2(
  splitdot.test, 
  filename = paste0("results/99_paper-figures/fig4_single-cell/de-dotplot-options/de-dot_top", numgenes, "genes.pdf"), 
  width = 4.75, height = 6, units = "in"
  )

grid.test <-  plot_grid(umap.v, splitdot.noT, axis = "t", rel_widths = c(2, 4.5), labels = c("a", "b"))
cowplot::ggsave2(
  grid.test, 
  filename = "results/99_paper-figures/fig4_single-cell/grid-test.pdf", 
  width = 6.5, height = 6, units = "in"
)


grid.test2 <- umap.v + splitdot.noT + plot_layout(widths = c(2, 5), heights = c(1,0)) 

cowplot::ggsave2(
  grid.test2, 
  filename = "results/99_paper-figures/fig4_single-cell/grid-test2.pdf", 
  width = 6.5, height = 6, units = "in"
)

## Panel C: Interferome plot --------------------------------------------------
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
fig4$c <- ggplot(data = plot.dat, 
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
            size = 5/.pt, vjust = 0.5, hjust = 0.0, nudge_x = 0.5) +
  geom_text(data = plot.dat[, c("celltype", "pval")] %>% 
              unique(), # pvalues
            aes(y = celltype, x = 100, 
                label = formatC(pval, format = "e", digits = 0)),
            hjust = 0, size = 5/.pt, nudge_x = 1,
            color = unique(plot.dat[, c("celltype", "color")])$color) +
  annotate(geom = "richtext", x  = 101, y = 22, size = 5.5/.pt,
           label = "enrichment *p* value", hjust = 0,
           fill = NA, label.color = NA, # remove background and outline
           label.padding = grid::unit(0, "pt")) + # remove padding
  theme_classic() +
  theme(
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5.5),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 6, 
                               color = "black"),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks.y = element_line(size = 0.25),
    axis.ticks.x = element_line(size = 0.25),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 3))

ggsave(fig4$c, 
       filename = "results/99_paper-figures/fig4_single-cell/panelC_interferome.pdf", 
       width = 3.5, height = 3, units = "in")

## Panel D: GO plot -----------------------------------------------------------

### Arrange ===================================================================
