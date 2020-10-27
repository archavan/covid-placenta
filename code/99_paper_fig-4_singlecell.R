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
library(cowplot)
library(ggtext)
library(ungeviz) # for geom_hpline

### data ======================================================================
## Seurat
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")
celltypes <- seur@meta.data$annotation_merged %>% unique()

## ref data for annotation
ref <- read.csv("results/01_reference-atlas/vento-surya-pavli_joined.csv", stringsAsFactors = FALSE)

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

## CellphoneDB output
cpdb.cntrl <- read.table("results/06_cellphonedb/out/cntrl/count_network.txt", header = TRUE, sep = "\t", quote = "")
cpdb.covid <- read.table("results/06_cellphonedb/out/covid/count_network.txt", header = TRUE, sep = "\t", quote = "")

### Make individual plots =====================================================
## Panel A: UMAP --------------------------------------------------------------
## colours for clusters to match with stacked barplot from supp. fig
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

## set idents
Idents(seur) <- seur@meta.data$annotation_merged

# plotting fn
plot_umap <- function(object, title = NA, label = "none", 
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
    scale_color_manual(values = stack.colors) +
    annotate(geom = "text", x = Inf, y = -Inf, label = title,
             hjust = 1, vjust = -1, size = 6/.pt, fontface = 1) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_classic() +
    theme(
#      aspect.ratio = 1,
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(size = 5),
      axis.line = element_line(size = 0.25),
      axis.ticks = element_blank()
    )
  
  if(label == "none") {
    p <- p
  }
  
  if(label == "label") {
    p <- p + 
      geom_label_repel(data = ulab, 
                       aes(x = UMAP_1, y = UMAP_2, label = ident),
                       size = 5/.pt,
                       segment.size = 0.2, segment.color = "black",
                       force = 1,
                       fontface = 1, 
                       label.padding = unit(0.1, "lines"), 
                       label.r = 0,
                       label.size = 0, 
                       seed = 10, 
                       fill = alpha(c("white"), 0.65))
  }
  
  if(label == "text") {
    p <- p + 
      geom_text_repel(data = ulab, 
                      aes(x = UMAP_1, y = UMAP_2, label = ident),
                      size = 5/.pt,
                      segment.size = 0.2, segment.color = "black",
                      force = 1,
                      fontface = 1, 
                      seed = 10)
  }
  
  return(p)
  
}

## individual plots
umap.all <- plot_umap(object = seur, label = "label", title = "", )

## panel B: annotation by correlation: vento ----------------------------------
## average by cluster
Idents(seur) <- "seurat_clusters"
plac <- AverageExpression(seur) %>% .[["SCT"]]

## get data in shape 
plac$gene_name <- rownames(plac)
plac <- dplyr::inner_join(plac, 
                          dplyr::select(ref, "gene_name", contains("vento")), 
                          by = "gene_name")
plac <- dplyr::relocate(plac, gene_name)

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
  for(j in c("vento")){
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

## heatmap 
ann <- unique(seur@meta.data[, c("seurat_clusters", "annotation", "annotation_merged")])

plot.dat <- cormat %>% 
  filter(Var1_source %in% c("vento"), 
         Var2_source %in% c("clust"))

plot.dat$lab <- ann$annotation_merged[match(x = plot.dat$Var2, 
                                            table = ann$seurat_clusters)]
plot.dat$lab[!(plot.dat$Var1 %in% c("vento_HB"))] <- NA # this is a lazy hack to make geom_text add label only once. If this was filtered to Var1_source == vento, it would add 35 (number of vento clusters) labels on top of each other.

cor.plot <- ggplot(plot.dat, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = 'white', high = 'red',
                      name = "Spearman Correlation",
                      limits = c(0.40, 0.95),
                      breaks = c(0.40, 0.60, 0.80),
                      labels = c("0.4", "0.6", "0.8")) +
  geom_point(aes(Var1, Var2, alpha = top3),
             size = 1, shape = 19, stroke  = 0) +
  scale_alpha_manual(values = c(1, 0.5, 0.25), 
                     breaks = c(1, 2, 3),
                     name = "Top 3 match", na.value = 0) +
  guides(alpha = guide_legend(override.aes = list(size = 1.5),
                              title.position="top", 
                              title.hjust = 0),
         fill = guide_colorbar(title.position = "top", 
                               title.hjust = 1, 
                               frame.colour = "black", 
                               frame.linewidth = 0.25, 
                               ticks.colour = "black", 
                               ticks.linewidth = 0.25)) +
  coord_cartesian(clip = "off") +
  geom_text(aes(label = lab), x = 33, size = 5.25/.pt, 
            hjust = 0, color = "black") +
  labs(x = "Reference", y = "Query") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 5.25, 
                               hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 5.25, color = "black"),
    axis.ticks = element_line(size = 0.25),
    axis.ticks.length = unit(0.15, units = c('lines')),
    axis.title = element_text(size = 6),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.25, colour = "black"),
    legend.key.size = unit(0.5, units = c('lines')), 
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5),
    legend.position = "top",
    legend.box.spacing = unit(0, "lines"),
    legend.background = element_blank()
  )

## Panel C: DE dotplot --------------------------------------------------------
## create celltype_covid ident to use for seurat dotplot
seur$celltype_covid <- paste0(seur$annotation_merged, "_", seur$covid)

## theme for DE genes faceted dotplot
theme_dotplot <- theme(
  panel.background = element_blank(),
  text = element_text(color = "Black"),
  panel.grid = element_blank(),
  panel.border = element_rect(fill = NA, size = 0.15, color = "black"),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, color = "Black"),
  axis.text.y = element_text(size = 5, color = "Black"),
  axis.title = element_blank(),
  axis.ticks = element_line(size = 0.1, color = "Black"),
  strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, 
                              size = 5.25, color = "Black"),
  strip.background = element_blank(),
  legend.title = element_text(size = 5.25),
  # legend.position = "right",
  # legend.box.spacing = unit(0.1, "lines"),
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.box = "horizontal",
  legend.background = element_blank(),
  legend.box.spacing = unit(0, "lines"),
  legend.title.align = 0,
  legend.text = element_text(size = 5),
  legend.key.size = unit(0.5, "lines"),
  legend.key = element_blank(),
#  legend.spacing.y = unit(0, "lines"),
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
                        name = "avg. exp. (scaled)",
                        guide = guide_colorbar(
                          title.position = "top", 
                          title.hjust = 0.5
                        )
    ) +
    scale_x_discrete(breaks = c("cntrl", "covid"), 
                     labels = c("Control", "Covid")) +
    scale_size_area(max_size = 2.5,
                    name = "% cells with exp.",
                    guide = guide_legend(
                      override.aes = list(
                        color = "White", 
                        fill = "Grey30"
                      ),
                      title.position = "top", title.hjust = 0.5
                    )
    ) +
    facet_grid(. ~ celltype) +
    ungeviz::geom_hpline(data = de.ann,
                         aes(y = gene, x = 1.5, color = diff.exp),
                         size = 0.30, width = 1, lineend = "butt") +
    scale_color_manual(breaks = "true",
                       labels = "< 0.05",
                       values = "Black",
                       name = "DE adj.p",
                       guide = guide_legend(title.position = "top",
                                            title.hjust = 0.5)) +
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

## plot with T celss
numgenes <- 5
features <- select_genes(n = numgenes, exclude.celltypes = exclude2)
splitdot.top5 <- plot_splitdot(
  object = seur, 
  features = features,
  exclude.celltypes = exclude2
)

## plot without T cells
numgenes <- 5
features <- select_genes(n = numgenes, exclude.celltypes = exclude1)
splitdot.top5.noTcells <- plot_splitdot(
  object = seur, 
  features = features,
  exclude.celltypes = exclude1
)


## Panel D: Interferome plot --------------------------------------------------
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
p.ifome <- ggplot(data = plot.dat, 
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
  annotate(geom = "richtext", x  = 102, y = 22,
           size = 5/.pt,
           label = "enrichment<br>*p* value",
           hjust = 0, vjust = 0.25,
           fill = NA, label.color = NA, # remove background and outline
           label.padding = grid::unit(0, "pt")) + # remove padding
  theme_classic() +
  theme(
    plot.margin = margin(17, 23, 6, 6),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, 
                               size = 6, color = "black"),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks.y = element_line(size = 0.25),
    axis.ticks.x = element_line(size = 0.25),
    legend.direction = "horizontal",
    legend.position = c(0.5, 1),
    legend.justification = c(0.5, 0),
    legend.background = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5)
  )

## Panel E: CellPhoneDB results -----------------------------------------------
# create a dataset of differences
count.fc <- inner_join(x = cpdb.cntrl, y = cpdb.covid, 
                       by = c("SOURCE", "TARGET"),
                       suffix = c(".cntrl", ".covid"))
count.fc$fc <- count.fc$count.covid / count.fc$count.cntrl


# plot log fold change in number of interactions
cpdb <- count.fc %>% 
  filter(!(SOURCE %in% c("dec.Gran", "dec.Bcells", "vil.Ery", "dec.Tcell_3"))) %>% 
  filter(!(TARGET %in% c("dec.Gran", "dec.Bcells", "vil.Ery", "dec.Tcell_3"))) %>%
  select(SOURCE, TARGET, fc) %>%
  mutate(logfc = log(fc)) %>% 
  arrange(logfc) %>% 
  mutate(SOURCE = factor(SOURCE, levels = rev(.$SOURCE) %>% unique())) %>% 
  mutate(TARGET = factor(TARGET, levels = rev(.$TARGET) %>% unique())) %>% 
  
  ggplot(., aes(SOURCE, TARGET)) +
  geom_tile(aes(fill = logfc)) +
  scale_fill_gradient2(low = "#4393c3", mid = "white", high = "#d6604d", 
                       midpoint = 0,
                       name = "log(covid/control)") +
  theme_bw() +
  theme(
#    aspect.ratio = 1,
    panel.border = element_rect(size = 0.25, colour = "black"),
    axis.ticks = element_line(size = 0.25),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, color = "black"),
    axis.text.y = element_text(color = "black", size = 6),
    legend.key.size = unit(0.4, "lines"),
    legend.title = element_text(size = 5.5, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.position = "top",
    legend.box.spacing = unit(0, "lines"),
    legend.background = element_blank()
  )

### Arrange ===================================================================
## version with T cells (but not Tcell_3 -- too few cells)
ade.aligned <- align_plots(umap.all, p.ifome, cpdb, align = "v", axis = "lr")
col1 <- plot_grid(ade.aligned[[1]], 
                  ade.aligned[[2]],
                  ade.aligned[[3]],
                  ncol = 1,
                  rel_heights = c(2.75, 2.75, 3),
                  align = "none",
                  labels = c("A", "D", "E"),
                  label_size = 8,
                  label_fontface = "bold") 

bc.aligned <- align_plots(cor.plot, splitdot.top5, align = "v", axis = "lr")
col2 <- plot_grid(bc.aligned[[1]],
                  bc.aligned[[2]],
                  ncol = 1,
                  rel_heights = c(4, 5),
                  labels = c("B", "C"),
                  label_size = 8,
                  label_fontface = "bold")

composite <- plot_grid(col1, col2, nrow = 1, rel_widths = c(3, 4))

cowplot::ggsave2(
  filename="results/99_paper-figures/fig4_single-cell/fig04_composite.png",
  composite, 
  width = 7, height = 9, units = "in", type = "cairo", dpi = 600
)

 ## version without any T cells
bc.aligned.2 <- align_plots(cor.plot, splitdot.top5.noTcells, 
                            align = "v", axis = "lr")
col2.2 <- plot_grid(bc.aligned.2[[1]],
                    bc.aligned.2[[2]],
                    ncol = 1,
                    rel_heights = c(4, 5),
                    labels = c("B", "C"),
                    label_size = 8,
                    label_fontface = "bold")

composite.2 <- plot_grid(col1, col2.2, nrow = 1, rel_widths = c(3, 4))

cowplot::ggsave2(
  filename="results/99_paper-figures/fig4_single-cell/fig04_composite_noTcells.png",
  composite.2, 
  width = 7, height = 9, units = "in", type = "cairo", dpi = 600
)
  
### try different panel arrangements ==========================================
r1.aln <- align_plots(cor.plot + 
                        theme(plot.margin = margin(6, 50, 6, 6)), 
                      p.ifome, align = "h", axis = "tb")
r1 <- plot_grid(r1.aln[[1]], r1.aln[[2]],
                ncol = 2, 
                rel_widths = c(3.5, 3.25),
                labels = c("A", "D"),
                label_size = 8, 
                label_fontface = "bold")
c1.aln <- align_plots(umap.all, cpdb, align = "v", axis = "lr")
c1 <- plot_grid(c1.aln[[1]], c1.aln[[2]],
                ncol = 1, 
                rel_heights = c(0.9, 1.2),
                labels = c("B", "C"),
                label_size = 8, 
                label_fontface = "bold")
c2 <- plot_grid(splitdot.top5,
                ncol = 1,
                labels = c("E"),
                label_size = 8, 
                label_fontface = "bold")
r2 <- plot_grid(c1, c2, 
                ncol = 2, 
                rel_widths = c(2.6, 4.4))

comp.opt2 <- plot_grid(r1, r2, nrow = 2, rel_heights = c(3.8, 5.2))

cowplot::ggsave2(
  filename="results/99_paper-figures/fig4_single-cell/fig04_composite_option2.png",
  comp.opt2, 
  width = 7, height = 9, units = "in", type = "cairo", dpi = 600
)
### end =======================================================================
