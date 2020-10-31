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

### UMAP ======================================================================
## data -----------------------------------------------------------------------
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")
celltypes <- seur@meta.data$annotation_merged %>% unique()

## colours for clusters to match with stacked barplot from supp. fig ----------
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

## extract plotting data -------------------------------------------------------
Idents(seur) <- seur@meta.data$annotation_merged

umap.gg <- DimPlot(seur, label = TRUE, shuffle = TRUE) # save ggplot object, extract plotting data from it, and plot independently with ggplot. 
udat <- umap.gg$data # point coordinates
ulab <- umap.gg$layers[[2]]$data # label coordinates

ulab <- ulab[order(as.character(ulab$ident)), ]
ulab$newcolours <- stack.colors
ulab <- mutate(ulab, 
               "full_lab" = case_when(
                 ident == "dec.APC" ~ "dec.APC (Antigen presenting cell)",
                 ident == "dec.Bcells" ~ "dec.Bcells (B cell)",
                 ident == "dec.DSC" ~ "dec.DSC (Decidual stromal cell)",
                 ident == "dec.Endo" ~ "dec.Endo (Endothelial)",
                 ident == "dec.FB" ~ "dec.FB (Fibroblast)",
                 ident == "dec.Gran" ~ "dec.Gran (Granulocyte)",
                 ident == "dec.Mono_1" ~ "dec.Mono_1 (Monocyte)",
                 ident == "dec.Mono_2" ~ "dec.Mono_2",
                 ident == "dec.NK_1" ~ "dec.NK_1 (Natural killer)",
                 ident == "dec.NK_2" ~ "dec.NK_2",
                 ident == "dec.NK_3" ~ "dec.NK_3",
                 ident == "dec.SMC" ~ "dec.SMC (Smooth muscle cell)",
                 ident == "dec.Tcell_1" ~ "dec.Tcell_1 (T cell)",
                 ident == "dec.Tcell_2" ~ "dec.Tcell_2",
                 ident == "dec.Tcell_3" ~ "dec.Tcell_3",
                 ident == "vil.Ery" ~ "vil.Ery (Erythrocyte/blast)",
                 ident == "vil.EVT" ~  "vil.EVT (Extravillous trophoblast)",
                 ident == "vil.FB" ~ "vil.FB (Fibroblast)",
                 ident == "vil.Hofb" ~ "vil.Hofb (Hofbauer)",
                 ident == "vil.SCT" ~ "vil.SCT (Syncytial trophoblast)",
                 ident == "vil.VCT" ~ "vil.VCT (Villous cytotrophoblast)"
               ))

udat$ident <- factor(udat$ident, levels = ulab$ident)

## plot -----------------------------------------------------------------------
umap <- ggplot() +
  geom_point(data = udat, 
             aes(x = UMAP_1, y = UMAP_2, color = ident),
             size = 0.1) +
  scale_color_manual(values = alpha(ulab$newcolours, 1), 
                     breaks = ulab$ident, 
                     labels = ulab$full_lab) +
  guides(color = guide_legend(override.aes = list(size = 1.5),
                              keyheight = unit(0.5, "lines"), 
                              keywidth = unit(0, "lines"),
                              ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
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
  annotate(geom = "text", x = Inf, y = -14.25, label = "dec = decidual",
           hjust = 1, color = "black", size = 5/.pt) +
  annotate(geom = "text", x = Inf, y = -15.25, label = "vil = villous",
           hjust = 1, color = "black", size = 5/.pt) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 5),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 5),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.box.spacing = unit(0, "lines")
  )

cowplot::ggsave2(
  umap,
  filename = "results/99_paper-figures/fig4_single-cell/04a_umap.pdf",
  width = 4, height = 2.75, units = "in"
)

### DE dotplot ================================================================
## data -----------------------------------------------------------------------
de.genes <- list()
for(i in celltypes){
  de.genes[[i]] <- read.csv(
    paste0("results/04_de-genes-by-celltype/logfc_0.40/de/files/de-genes_", i, ".csv")
  )
  de.genes[[i]] <- dplyr::rename(de.genes[[i]], "gene" = "X")
  de.genes[[i]]$celltype <- i
}

## create celltype_covid ident to use for seurat dotplot ----------------------
seur$celltype_covid <- paste0(seur$annotation_merged, "_", seur$covid)

## theme for DE genes faceted dotplot -----------------------------------------
theme_dotplot <- theme(
  panel.background = element_blank(),
  text = element_text(color = "Black"),
  panel.grid = element_blank(),
  panel.border = element_rect(fill = NA, size = 0.25, color = "grey"),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, color = "Black"),
  axis.text.y = element_text(size = 5, color = "Black"),
  axis.title = element_blank(),
  axis.ticks = element_line(size = 0.25, color = "grey"),
  strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, 
                              size = 5.25, color = "Black"),
  strip.background = element_blank(),
  legend.title = element_text(size = 5.25),
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.box = "horizontal",
  legend.background = element_blank(),
  legend.box.spacing = unit(0, "lines"),
  legend.title.align = 0,
  legend.text = element_text(size = 5),
  legend.key.size = unit(0.5, "lines"),
  legend.key = element_blank(),
  panel.spacing.x = unit(0, "lines")
)

## dotplot function -----------------------------------------------------------
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
    scale_fill_gradient(low = "White", high = "#d94801", 
                        name = "avg. exp. (scaled)",
                        guide = guide_colorbar(
                          title.position = "top", 
                          title.hjust = 0.5, 
                          ticks.colour = "black"
                        )
    ) +
    scale_x_discrete(breaks = c("cntrl", "covid"), 
                     labels = c("ctrl", "COVID")) +
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

## exclude celltypes lists ----------------------------------------------------
exclude <- c("dec.Tcell_3", "dec.Bcells", "dec.Gran", "vil.Ery")

## select genes to plot: top x number of genes from clusters of interest ------
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

## plot -----------------------------------------------------------------------
numgenes <- 5
features <- select_genes(n = numgenes, exclude.celltypes = exclude)
splitdot.top5 <- plot_splitdot(
  object = seur, 
  features = features,
  exclude.celltypes = exclude
)

cowplot::ggsave2(
  splitdot.top5,
  filename = "results/99_paper-figures/fig4_single-cell/04b_splitdot-top5genes.pdf",
  width = 4, height = 5.5, units = "in"
)

### Interferome plot ==========================================================
## data -----------------------------------------------------------------------
ifome.n <- read.csv("results/04_de-genes-by-celltype/logfc_0.40/interferome/files/number-of-DE-genes-in-interferome.csv")

## color by pval --------------------------------------------------------------
ifome.n$color <- ifelse(test = ifome.n$pval < 0.05, yes = "Red", no = "Grey30")

## tidy -----------------------------------------------------------------------
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

## plot -----------------------------------------------------------------------
plot.dat <- ifome.n.long
p.ifome <- ggplot(data = plot.dat, 
                  aes(y = celltype, x = pct)) +
  geom_bar(aes(fill = interferome),
           position = "stack", 
           stat = "identity", 
           width = 0.8) +
  scale_fill_manual(values = c("grey90", "grey40"), # #4e79a7 or #74add1
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
            size = 4.5/.pt, vjust = 0.5, hjust = 0.0, nudge_x = 0.5) +
  geom_text(data = plot.dat[, c("celltype", "pval")] %>% 
              unique(), # pvalues
            aes(y = celltype, x = 100, 
                label = formatC(pval, format = "e", digits = 0)),
            hjust = 0, 
            size = 5/.pt, 
            nudge_x = 1,
            color = unique(plot.dat[, c("celltype", "color")])$color) +
  annotate(geom = "richtext", x  = 102, y = 22,
           size = 5/.pt,
           label = "enrichment<br>*p* value",
           hjust = 0, vjust = 0.25,
           fill = NA, label.color = NA, # remove background and outline
           label.padding = grid::unit(0, "pt")) + # remove padding
  theme_classic() +
  theme(
    plot.margin = margin(6, 23, 6, 6),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, 
                               size = 5, color = "black"),
    axis.text.x = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5.25),
    axis.title.y = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks.y = element_line(size = 0.25),
    axis.ticks.x = element_line(size = 0.25),
    legend.direction = "horizontal",
    legend.position = "top",
    legend.box.spacing = unit(0, "lines"),
    legend.background = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 5.25),
    legend.text = element_text(size = 5)
  )

cowplot::ggsave2(
  p.ifome,
  filename = "results/99_paper-figures/fig4_single-cell/04c_interferome.pdf",
  width = 3, height = 2.75, units = "in"
)

### Metascape results ==========================================================
## data -----------------------------------------------------------------------
meta <- read.csv("results/07_metascape/metascape_fewer-celltypes/Enrichment_heatmap/HeatmapSelectedGO.csv")

## prepare for plotting -------------------------------------------------------
meta$Description_new <- gsub(pattern = "positive", 
                             replacement = "+ve", 
                             x = meta$Description)
meta$Description_new[meta$Description_new == "Hemostasis"] <- "hemostasis"
meta$Description_new[meta$Description_new == "Adaptive Immune System"] <- "adaptive immune system"
meta$Description_new[meta$Description_new == "Selenocysteine synthesis"] <- "selenocysteine synthesis"

meta.long <- pivot_longer(data = meta, 
                          cols = 3:19, 
                          names_to = "celltype", 
                          names_prefix = "X_LogP_", 
                          values_to = "logP")

meta.long$Description_new <- factor(
  meta.long$Description_new,
  levels = rev(c( 
    "regulation of MAPK cascade",                
    "hemostasis",                  
    "regulated exocytosis",                         
    "+ve regulation of cellular component movement",
    "activation of immune response",               
    "+ve regulation of cell death",                 
    "response to interferon-gamma",                
    "response to toxic substance",                 
    "viral life cycle",                           
    "regulation of multi-organism process",       
    "response to wounding",        
    "+ve regulation of cytokine production",        
    "regulation of cell adhesion",       
    "apoptotic signaling pathway",                 
    "Signaling by Interleukins",                 
    "regulation of cellular response to stress", 
    "supramolecular fiber organization",         
    "adaptive immune system",                      
    "leukocyte migration",                        
    "selenocysteine synthesis"                    
  ))
) # same order as in metascape output, which uses clustering

meta.long$celltype <- factor(
  meta.long$celltype,
  levels = c(
    "vil.EVT", "vil.SCT", "vil.Hofb", "dec.DSC", "dec.Endo", "dec.SMC", "dec.NK_3", "dec.NK_2", "dec.NK_1", "dec.Tcell_1", "dec.APC", "dec.FB", "vil.FB", "dec.Mono_1", "vil.VCT", "dec.Mono_2", "dec.Tcell_2"
  )
)

## plot -----------------------------------------------------------------------
p.meta <- ggplot(meta.long, aes(x = celltype, y = Description_new)) +
  geom_tile(aes(fill = -logP)) +
  scale_fill_gradient(
    name = "- log(p)",
    low = "white", high = "#d94801",
    limits = c(0, 20), 
    oob = scales::squish,
    breaks = c(0, 10, 20),
    labels = c("0", "10", ">20"),
    guide = guide_colorbar(ticks.colour = "black")
    ) +
  coord_fixed() +
  scale_y_discrete(position = "right") +
  theme(
    panel.border = element_rect(size = 0.25, colour = "black", fill = NA),
    axis.ticks = element_line(size = 0.25),
    axis.ticks.length = unit(0.1, "lines"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, color = "black"),
    axis.text.y = element_text(color = "black", size = 5),
    legend.key.width = unit(0.4, "lines"),
    legend.key.height = unit(0.25, "lines"),
    legend.title = element_text(size = 5.5, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.position = "top",
    legend.box.spacing = unit(0, "lines"),
    legend.background = element_blank()
  )

cowplot::ggsave2(
  p.meta,
  filename = "results/99_paper-figures/fig4_single-cell/04d_metascape.pdf",
  width = 3, height = 2.5, units = "in"
)

### CellPhoneDB results =======================================================
## data -----------------------------------------------------------------------
cpdb.cntrl <- read.table("results/06_cellphonedb/out/cntrl/count_network.txt", header = TRUE, sep = "\t", quote = "")
cpdb.covid <- read.table("results/06_cellphonedb/out/covid/count_network.txt", header = TRUE, sep = "\t", quote = "")

## create a dataset of differences --------------------------------------------
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
  scale_fill_gradient2(low = "#4393c3", mid = "white", high = "#d94801", 
                       midpoint = 0,
                       name = "log(covid/control)",
                       guide = guide_colorbar(ticks.colour = "black")) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(size = 0.25, colour = "black"),
    axis.ticks = element_line(size = 0.25),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, color = "black"),
    axis.text.y = element_text(color = "black", size = 5),
    legend.key.width = unit(0.4, "lines"),
    legend.key.height = unit(0.25, "lines"),
    legend.title = element_text(size = 5.25, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.position = "top",
    legend.box.spacing = unit(0, "lines"),
    legend.background = element_blank()
  )

cowplot::ggsave2(
  cpdb,
  filename = "results/99_paper-figures/fig4_single-cell/04e_cellphonedb.pdf",
  width = 3, height = 3, units = "in"
)

### Arrange ===================================================================
## version with T cells (but not Tcell_3 -- too few cells)
ac.aligned <- align_plots(umap, p.ifome, align = "h", axis = "tb")
r1 <- plot_grid(ac.aligned[[1]], 
                ac.aligned[[2]],
                ncol = 2,
                rel_widths = c(4, 3),
                align = "none",
                labels = c("A", "C"),
                label_size = 8,
                label_fontface = "bold") 

c1 <- plot_grid(splitdot.top5,
                ncol = 1,
                align = "none",
                labels = c("B"),
                label_size = 8,
                label_fontface = "bold")

c2 <- plot_grid(p.meta, cpdb,
                ncol = 1,
                rel_heights = c(2.5, 3),
                align = "v",
                axis = "l",
                labels = c("D", "E"),
                label_size = 8,
                label_fontface = "bold")

r2 <- plot_grid(c1,
                c2,
                ncol = 2,
                rel_heights = c(4, 3),
                align = "none")

composite <- plot_grid(r1, r2, nrow = 2, rel_heights = c(2.75, 5.5))

cowplot::ggsave2(
  filename="results/99_paper-figures/fig4_single-cell/04_composite_v1.png",
  composite, 
  width = 7, height = 8.25, units = "in", type = "cairo", dpi = 600
)

cowplot::ggsave2(
  filename="results/99_paper-figures/fig4_single-cell/04_composite_v1.pdf",
  composite, 
  width = 7, height = 8.25, units = "in"
)

### end =======================================================================
