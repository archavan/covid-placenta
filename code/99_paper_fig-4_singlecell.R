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

# go
library(GO.db)

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
## Panel A: UMAP --------------------------------------------------------------
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
    annotate(geom = "text", x = Inf, y = -Inf, label = title,
             hjust = 1, vjust = -1, size = 6/.pt, fontface = 1) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_classic() +
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 5),
          axis.line = element_line(size = 0.25),
          axis.ticks = element_blank())
  
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
                       label.size = 0, seed = 10)
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
umap.all <- plot_umap(object = seur, label = "text", title = "")
umap.cntrl <- plot_umap(object = seur,
                        split = TRUE, split.by = "covid", split.plot = "cntrl",
                        title = "Control")
umap.covid <- plot_umap(object = seur,
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
  # legend.position = "right",
  # legend.box.spacing = unit(0.1, "lines"),
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.box = "horizontal",
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
numgenes <- 7
features <- select_genes(n = numgenes, exclude.celltypes = exclude2)
splitdot.test <- plot_splitdot(
  object = seur, 
  features = features,
  exclude.celltypes = exclude2
)

cowplot::ggsave2(
  splitdot.test, 
  filename = paste0("results/99_paper-figures/fig4_single-cell/de-dotplot-options/de-dot_top", numgenes, "genes.pdf"), 
  width = 4.75, height = 6, units = "in"
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

ggsave(p.ifome, 
       filename = "results/99_paper-figures/fig4_single-cell/panelC_interferome.pdf", 
       width = 3.2, height = 3, units = "in")

## Panel D: GO plot -----------------------------------------------------------
## prepare data
# we  want only the following celltypes: dAPCs, monocytes, NK1, hoffbauers, cytotrophoblast
go.pdat <- de.go[c("dec.APC", "dec.Mono_1", "dec.Mono_2", "dec.NK_1", "vil.Hofb", "vil.VCT")]

go.pdat <- do.call(rbind, go.pdat)

go.pdat$pvalue <- NULL
go.pdat$geneID <- NULL
rownames(go.pdat) <- NULL

# calculate enrichment
go.pdat$b <- as.numeric(sapply(
  strsplit(go.pdat$GeneRatio, split = "/"), MARGIN = 1, FUN = '[['))

go.pdat$n <- as.numeric(sapply(
  strsplit(go.pdat$GeneRatio, split = "/"), MARGIN = 2, FUN = '[['))

go.pdat$B <- as.numeric(sapply(
  strsplit(go.pdat$BgRatio, split = "/"), MARGIN = 1, FUN = '[['))

go.pdat$N <- as.numeric(sapply(
  strsplit(go.pdat$BgRatio, split = "/"), MARGIN = 2, FUN = '[['))

go.pdat$enrichment <- (go.pdat$b / go.pdat$n) / (go.pdat$B / go.pdat$N)

head(go.pdat)

# function
tileplot_go <- function(go.pdat.sub) {
  
  go.pdat.sub <- unique(go.sub)
  
  # truncate long labels. do this first to retain factor levels
  go.pdat.sub$trunclab <- stringr::str_trunc(
    paste0(go.pdat.sub$Description, " ", go.pdat.sub$ID), 
    width = 65, side = "left")
  
  # clustering for ordering
  cdat <- pivot_wider(data = go.pdat.sub, id_cols = "trunclab", 
                      names_from = celltype, values_from = p.adjust)
  cdat[, 2:ncol(cdat)] <- apply(cdat[, 2:ncol(cdat)], 
                                MARGIN = 2, FUN = function(x){replace_na(x, 1)})
  
  cdat <- as.data.frame(cdat)
  rownames(cdat) <- cdat$trunclab
  cdat$trunclab <- NULL
  
  cor.matrix <- cor(t(cdat))
  
  dd <- as.dist((1 - cor.matrix)) 
  hc <- hclust(dd, method = 'complete')
  
  cdat <- cdat[hc$order, ]
  
  go.pdat.sub$trunclab <- factor(go.pdat.sub$trunclab, 
                                    levels = rownames(cdat))
  
  golab <- data.frame(
    x = 6,
    y = 1:length(levels(go.pdat.sub$trunclab)),
    lab = gsub("(.*)(GO:[0-9]*)", "\\2", levels(go.pdat.sub$trunclab))
  )
  
  # plot
  p <- ggplot(go.pdat.sub,
              aes(x = celltype, 
                  y = trunclab, 
                  fill = -log10(p.adjust))) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = colorRampPalette(c("White", "#4e79a7"))(10)[3],
                        high = colorRampPalette(c("White", "#4e79a7"))(10)[10],
                        name = "- log10 (adj. p)") +
    #coord_fixed() +
    theme_classic() +
    theme(
      plot.margin = margin(25, 6, 6, 6),
      axis.line = element_line(size = 0.25),
      axis.ticks = element_line(size = 0.25),
      panel.grid.major = element_line(size = 0.15),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                                 size = 5, color = "black"),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 5.25, color = "black"),
      legend.position = c(1.2, 1),
      legend.direction = "horizontal",
      legend.key.size = unit(0.5, "lines"),
      legend.justification = c(1, 0),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5),
      legend.title.align = 1
    )
  
  return(p)
}




##----

# write file and run it on revigo
write.table(
  go.pdat[, c("ID", "p.adjust")], 
  file = "results/99_paper-figures/fig4_single-cell/go-for-revigo.txt", 
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

revigo <- read.csv("results/99_paper-figures/fig4_single-cell/REVIGO.csv")

terminal <- function(terms, ontology=c("C", "P", "F"))
{
  FUN <- switch(match.arg(ontology),  
                C = GOCCPARENTS,
                P = GOBPPARENTS, 
                F = GOMFPARENTS)
  terminal <- terms
  seen <- c(terms, "all")
  while (length(terms)) {
    seen <- c(terms, seen)
    terms <- mappedRkeys(FUN[terms])
    terminal <- terminal[!terminal %in% terms]
    terms <- terms[!terms %in% seen]
  }
  terminal
}

leaf <- terminal(terms = as.character(revigo$term_ID[revigo$dispensability < 0.6]), ontology = "P")
length(leaf)

# plotting subset
to.keep <- revigo$term_ID[revigo$dispensability < 0.75]
# go.sub <- go.pdat[go.pdat$ID %in% leaf, ]
go.sub <- go.pdat[go.pdat$ID %in% to.keep, ]
# write to file to manually curate
go.sub$chars <- stringr::str_count(go.sub$Description)
write.csv(go.sub, "results/99_paper-figures/fig4_single-cell/go_revio_disp-below-0.75.csv", row.names = FALSE)



go.sub <- go.sub[order(go.sub$enrichment, decreasing = TRUE), ]
go.sub <- go.sub[1:40, ]

p.go <- tileplot_go(go.pdat.sub = go.sub)

go.sub$Description <- gsub("response", "resp.", go.sub$Description)
go.sub$Description <- gsub("endoplasmic reticulum", "ER", go.sub$Description)
go.sub$Description <- gsub("regulation", "reg", go.sub$Description)



ggsave(
  p, 
  filename = "results/99_paper-figures/fig4_single-cell/test_go.pdf",
  width = 3.5, height = 3, units = "in"
)

ggsave(
  p, 
  filename = "../results/04_de-genes-by-celltype/logfc_0.40/go/plots/go_dotplot_revigo.pdf",
  width = 5.75, height = 6, units = "in"
)





### Arrange ===================================================================
layout <- c(area(t = 1, l = 1, b = 3, r = 3),
            area(t = 1, l = 4, b = 8, r = 7),
            area(t = 4, l = 1, b = 6, r = 3),
            area(t = 7, l = 1, b = 9, r = 3))

umap.all + splitdot.test + p.ifome + p.go +
  plot_layout(guides = "keep", design = layout)

ac.aligned <- align_plots(umap.all, p.ifome, align = "v", axis = "lr")

col1 <- plot_grid(ac.aligned[[1]], 
                  ac.aligned[[2]], 
                  p.go,
                  ncol = 1,
                  rel_heights = c(2.5, 2.75, 3.5),
                  align = "none",
                  labels = c("A", "C", "D"),
                  label_size = 8,
                  label_fontface = "bold") 

col2 <- plot_grid(splitdot.test,
                  labels = c("B"),
                  label_size = 8,
                  label_fontface = "bold")

full.test <- plot_grid(col1, col2, nrow = 1, rel_widths = c(3, 4))

ggsave(full.test, filename="results/99_paper-figures/fig4_single-cell/grid-test_full.png", width = 7, height = 9, units = "in", type = "cairo", dpi = 600)
  
