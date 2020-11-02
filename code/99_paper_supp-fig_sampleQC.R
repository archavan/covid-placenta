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

### uUMI per sample ===========================================================
n.umi <- ggplot(data = seur@meta.data, 
                aes(x = orig.ident_renamed, y = nCount_RNA)) +
  geom_violin(fill = "grey") +
  facet_grid(. ~ covid + tissue, 
             space = "free", 
             scales = "free",
             labeller = as_labeller(c(cntrl = "ctrl", 
                                      covid = "COVID", 
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
                aes(x = orig.ident_renamed, y = nFeature_RNA)) +
  geom_violin(fill = "grey") +
  facet_grid(. ~ covid + tissue, 
             space = "free", 
             scales = "free",
             labeller = as_labeller(c(cntrl = "ctrl", 
                                      covid = "COVID", 
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
cell.sample <- table(seur@meta.data$orig.ident_renamed,
                     seur@meta.data$annotation_merged,
                     seur@meta.data$covid,
                     seur@meta.data$tissue) %>% 
  as.data.frame()

cell.sample <- rename(cell.sample,
                      "sample" = "Var1",
                      "celltype" = "Var2",
                      "covid" = "Var3",
                      "tissue" = "Var4")

cell.sample$frac <- NA
for(i in 1:nrow(cell.sample)) {
  cell.sample$frac[i] <- cell.sample$Freq[i] / 
    sum(cell.sample$Freq[cell.sample$sample == cell.sample$sample[i]])
}

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


bar.sample <- ggplot(data = cell.sample %>% 
                       filter(frac > 0)) +
  geom_bar(aes(x = sample, y = frac, fill = celltype),
           stat = "identity", width = 0.8) +
  scale_fill_manual(name = "Cell type",
                    values = alpha(stack.colors, 0.8)) +
  facet_grid(. ~ covid + tissue, 
             space = "free_x", 
             scales = "free_x",
             labeller = as_labeller(c(cntrl = "ctrl",
                                      covid = "COVID",
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

bar.sample

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

Idents(seur) <- "orig.ident_renamed"
umap.sample.gg <- DimPlot(seur, label = FALSE, shuffle = TRUE) 

umap.sample <- ggplot() +
  geom_point(data = umap.sample.gg$data, 
             aes(x = UMAP_1, y = UMAP_2, color = ident),
             size = 0.1) +
  labs(x = "UMAP 1", y = "UMAP 2", tag = "E") +
  scale_color_manual(name = "Sample",
                     values = c("#c65c8a",
                                "#7cb843",
                                "#a361c7",
                                "#50a166",
                                "#ca5842",
                                "#48bcc1",
                                "#d19a44",
                                "#6683cb",
                                "#848039"),
                     guide = guide_legend(override.aes = list(size = 1))) +
  theme_umap

### UMAP split by covid and control ===========================================
Idents(seur) <- "covid"
umap.split.gg <- DimPlot(seur, shuffle = TRUE, label = FALSE)

umap.cntrl <- ggplot(
  data = umap.split.gg$data[umap.split.gg$data$ident == "cntrl", ],
  aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, color = "grey40") +
  annotate(geom = "text", label = "Control", 
           x = Inf, y = -Inf, hjust = 1.25, vjust = -0.5, size = 6/.pt) +
  labs(tag = "D") +
  theme_umap

umap.covid <- ggplot(
  data = umap.split.gg$data[umap.split.gg$data$ident == "covid", ],
  aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, color = "grey40") +
  annotate(geom = "text", label = "Covid", 
           x = Inf, y = -Inf, hjust = 1.25, vjust = -0.5, size = 6/.pt) +
  theme_umap

### arrange ===================================================================
row1 <- n.umi + n.gene
row2 <- bar.sample 
row3 <- umap.cntrl + umap.covid + umap.sample + plot_layout(widths = c(1, 1, 1))

composite <- row1 / row2 / row3 +
  plot_layout(heights = c(1.75, 2.5, 1.85)) &
  theme(legend.justification = "left")

cowplot::ggsave2(
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC_composite_v1.pdf",
  composite, 
  width = 6, height = 6, units = "in"
)

cowplot::ggsave2(
  composite, 
  filename = "results/99_paper-figures/supp-fig_sampleQC/supp-fig_sampleQC_composite_v1.png",
  width = 6, height = 6, units = "in", type = "cairo", dpi = 600
)

### end  ======================================================================
