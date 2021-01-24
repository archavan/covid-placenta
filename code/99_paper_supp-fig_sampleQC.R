# Make supplementary figures related to sample QC

### packages ==================================================================
library(tidyverse)

library(Seurat) # single cell

library(ggrepel) # plotting
library(patchwork)
library(ggthemes)
library(cowplot)
library(grid)

### data ======================================================================
seur <- readRDS("data/seurat-object_annotated.rds")

### uUMI per sample ===========================================================
n.umi <- ggplot(data = seur@meta.data, 
                aes(x = sample_id, y = nCount_RNA)) +
  geom_violin(fill = "grey", size = 0.35) +
  facet_grid(. ~ covid + tissue, 
             space = "free", 
             scales = "free",
             labeller = as_labeller(c(cntrl = "", # will be added later with textGrob
                                      covid = "", 
                                      decidua = "Decidua", 
                                      villi = "Villi"))) +
  labs(tag = "A", x = "Sample", y = "nUMI") +
  theme_classic(base_line_size = 0.25) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 5, color = "black", 
                               angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5.25),
    strip.text = element_text(size = 5.25),
    strip.background = element_blank(),
    plot.tag = element_text(size = 8, color = "black", face = 2)
  )

n.umi

### nGene per sample ==========================================================
n.gene <- ggplot(data = seur@meta.data, 
                aes(x = sample_id, y = nFeature_RNA)) +
  geom_violin(fill = "grey", size = 0.35) +
  facet_grid(. ~ covid + tissue, 
             space = "free", 
             scales = "free",
             labeller = as_labeller(c(cntrl = "", # will be added later with textGrob
                                      covid = "", 
                                      decidua = "Decidua", 
                                      villi = "Villi"))) +
  labs(tag = "B", x = "Sample", y = "nGene") +
  theme_classic(base_line_size = 0.25) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 5, color = "black", 
                               angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5.25),
    strip.text = element_text(size = 5.25),
    strip.background = element_blank(),
    plot.tag = element_text(size = 8, color = "black", face = 2)
  )

n.gene

## stacked bar: cells/celltype/sample -----------------------------------------
cell.sample <- seur@meta.data %>% 
  group_by(sample_id, annotation_merged, covid, tissue) %>% 
  summarize(n_cells = n()) %>% 
  ungroup() %>% 
  group_by(sample_id, covid, tissue) %>% 
  mutate(frac = n_cells / sum(n_cells)) %>% 
  ungroup()

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

clust.order <- c("dec.DSC", "dec.Endo", "dec.SMC", "dec.FB", 
                 "vil.FB", "vil.EVT", "vil.SCT", "vil.VCT", "vil.Ery", "vil.Hofb", 
                 "APC", "Bcell", "Gran", "Mono_1", "Mono_2", "NK_1", "NK_2", "NK_3", 
                 "Tcell_1", "Tcell_2", "Tcell_3")

cell.sample$annotation_merged <- factor(cell.sample$annotation_merged, levels = clust.order)

bar.sample <- ggplot(data = cell.sample %>% 
                       filter(frac > 0)) +
  geom_bar(aes(x = sample_id, y = frac, fill = annotation_merged),
           stat = "identity", width = 0.8) +
  scale_fill_manual(name = "Cell type",
                    values = alpha(stack.colors, 0.8)) +
  facet_grid(. ~ covid + tissue, 
             space = "free_x", 
             scales = "free_x",
             labeller = as_labeller(c(cntrl = "",
                                      covid = "",
                                      decidua = "Decidua",
                                      villi = "Villi"))
             ) +
  scale_y_continuous(expand = c(0, 0.05)) +
  labs(x = "Sample", y = "Fraction of cells", tag = "C") +
  theme_classic(base_line_size = 0.25) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 5, color = "black", 
                               angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5.25),
    strip.text = element_text(size = 5.25),
    strip.background = element_blank(),
    plot.tag = element_text(size = 8, color = "black", face = 2),
    legend.key.size = unit(0.7, "lines"),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )

### umap colored by sample ====================================================
theme_umap <- theme_classic(base_line_size = 0.25) +
  theme(aspect.ratio = 1,
        plot.margin = margin(6, 6, 6, 6),
        plot.tag = element_text(size = 8, face = 2),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.5, "lines"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

Idents(seur) <- "sample_id"
umap.sample.gg <- DimPlot(seur, label = FALSE, shuffle = TRUE) 

umap.sample <- ggplot() +
  geom_point(data = umap.sample.gg$data, 
             aes(x = UMAP_1, y = UMAP_2, color = ident),
             size = 0.1) +
  labs(x = "UMAP 1", y = "UMAP 2", tag = "E") +
  scale_color_manual(name = "Sample",
                     values = c("#ca5e4a",
                                "#5ba965",
                                "#c55a9f",
                                "#ad963d",
                                "#777acd"),
                     guide = guide_legend(override.aes = list(size = 1))) +
  theme_umap

### UMAP split by covid and control ===========================================
Idents(seur) <- "covid"
umap.split.gg <- DimPlot(seur, shuffle = TRUE, label = FALSE)

umap.cntrl <- ggplot(
  data = umap.split.gg$data[umap.split.gg$data$ident == "cntrl", ],
  aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, color = "grey40") +
  annotate(geom = "text", label = "CNTRL", 
           x = Inf, y = -Inf, hjust = 1.25, vjust = -0.5, size = 6/.pt) +
  labs(tag = "D") +
  theme_umap

umap.covid <- ggplot(
  data = umap.split.gg$data[umap.split.gg$data$ident == "covid", ],
  aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, color = "grey40") +
  annotate(geom = "text", label = "COVID", 
           x = Inf, y = -Inf, hjust = 1.25, vjust = -0.5, size = 6/.pt) +
  theme_umap

### arrange ===================================================================
row1 <- n.umi + n.gene
row2 <- bar.sample 
row3 <- umap.cntrl + umap.covid + umap.sample + plot_layout(widths = c(1, 1, 1))

composite <- row1 / row2 / row3 +
  plot_layout(heights = c(1.75, 2.5, 2)) &
  theme(legend.justification = "left")

# It's not easy to label nested facet labels appropriately using ggplot. The higher level is printed twice (ctrl is printed for both decidua and villi; same for COVID). We printed an empty string for these instead above when making the plot. Now using textGrob and linesGrobs, and cooridinates in inches (for a 6x6 inch final figure dimensions), we will add the lines and text. (I first saved the plot as pdf, opened in illustrator, grabbed all inch coordniates that are used below.)

l1 <- linesGrob(
  x = unit(c(  0.6040,   0.6040,   1.5037,   1.5037), units = "in"), 
  y = unit(c(6-0.5069, 6-0.4433, 6-0.4433, 6-0.5069), units = "in"),
  gp = gpar(colour = "black", lwd = 0.5) 
)
t1 <- textGrob(label = "CNTRL", 
               x = unit(mean(c(0.6040, 1.5037)), "in"), 
               y = unit(6-0.3749, "in"), 
               just = c(0.5, 0), 
               gp = gpar(fontsize = 5, col = "black"))

l2 <- linesGrob(
  x = unit(c(  1.5799,   1.5799,   2.3269,   2.3269), units = "in"), 
  y = unit(c(6-0.5069, 6-0.4433, 6-0.4433, 6-0.5069), units = "in"),
  gp = gpar(colour = "black", lwd = 0.5) 
)
t2 <- textGrob(label = "COVID", 
               x = unit(mean(c(1.5799, 2.3269)), "in"), 
               y = unit(6-0.3749, "in"), 
               just = c(0.5, 0), 
               gp = gpar(fontsize = 5, col = "black"))

l3 <- linesGrob(
  x = unit(c(  2.8856,   2.8856,   3.7851,   3.7851), units = "in"), 
  y = unit(c(6-0.5069, 6-0.4433, 6-0.4433, 6-0.5069), units = "in"),
  gp = gpar(colour = "black", lwd = 0.5) 
)
t3 <- textGrob(label = "CNTRL", 
               x = unit(mean(c(2.8856, 3.7851)), "in"), 
               y = unit(6-0.3749, "in"), 
               just = c(0.5, 0), 
               gp = gpar(fontsize = 5, col = "black"))

l4 <- linesGrob(
  x = unit(c(  3.8613,   3.8613,   4.6083,   4.6083), units = "in"), 
  y = unit(c(6-0.5069, 6-0.4433, 6-0.4433, 6-0.5069), units = "in"),
  gp = gpar(colour = "black", lwd = 0.5) 
)
t4 <- textGrob(label = "COVID", 
               x = unit(mean(c(3.8613, 4.6083)), "in"), 
               y = unit(6-0.3749, "in"), 
               just = c(0.5, 0), 
               gp = gpar(fontsize = 5, col = "black"))

l5 <- linesGrob(
  x = unit(c(  0.6040,   0.6040,   2.7608,   2.7608), units = "in"), 
  y = unit(c(6-2.5694, 6-2.5058, 6-2.5058, 6-2.5694), units = "in"),
  gp = gpar(colour = "black", lwd = 0.5) 
)
t5 <- textGrob(label = "CNTRL", 
               x = unit(mean(c(0.6040, 2.7608)), "in"), 
               y = unit(6-2.4375, "in"), 
               just = c(0.5, 0), 
               gp = gpar(fontsize = 5, col = "black"))

l6 <- linesGrob(
  x = unit(c(  2.8369,   2.8369,   4.6083,   4.6083), units = "in"), 
  y = unit(c(6-2.5694, 6-2.5058, 6-2.5058, 6-2.5694), units = "in"),
  gp = gpar(colour = "black", lwd = 0.5) 
)
t6 <- textGrob(label = "COVID", 
               x = unit(mean(c(2.8369, 4.6083)), "in"), 
               y = unit(6-2.4375, "in"), 
               just = c(0.5, 0), 
               gp = gpar(fontsize = 5, col = "black"))

composite <- ggdraw(composite) + 
  draw_grob(l1) + draw_grob(t1) + 
  draw_grob(l2) + draw_grob(t2) +
  draw_grob(l3) + draw_grob(t3) +
  draw_grob(l4) + draw_grob(t4) +
  draw_grob(l5) + draw_grob(t5) +
  draw_grob(l6) + draw_grob(t6)

cowplot::ggsave2(
  composite, 
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC_composite_v3.png",
  width = 6, height = 6, units = "in", type = "cairo", dpi = 600
)

cowplot::ggsave2(
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC_composite_v3.pdf",
  composite, 
  width = 6, height = 6, units = "in"
)

### end  ======================================================================
