# CellPhoneDB dotplot:
# 1. Only plot interactions among immune cells
# 2. Only plot interactions present in covid but absent in control.

### packages ==================================================================
library(tidyverse)
library(ggtext)
library(glue)

### data ======================================================================
cntrl <- list()
cntrl$means <- read.table("results/06_cellphonedb/out/cntrl/means.txt", header = TRUE, sep = "\t", quote = "")
cntrl$pvals <- read.table("results/06_cellphonedb/out/cntrl/pvalues.txt", header = TRUE, sep = "\t", quote = "")

covid <- list()
covid$means <- read.table("results/06_cellphonedb/out/covid/means.txt", header = TRUE, sep = "\t", quote = "")
covid$pvals <- read.table("results/06_cellphonedb/out/covid/pvalues.txt", header = TRUE, sep = "\t", quote = "")

### melt and join means and pval data =========================================
cntrl$means.long <- pivot_longer(data = cntrl$means, 
                                 cols = 12:ncol(cntrl$means), 
                                 names_to = "celltype.pair", 
                                 values_to = "mean")
cntrl$pvals.long <- pivot_longer(data = cntrl$pvals, 
                                 cols = 12:ncol(cntrl$pvals), 
                                 names_to = "celltype.pair", 
                                 values_to = "pval")
cntrl$means.pval <- inner_join(cntrl$means.long, cntrl$pvals.long) %>% 
  select("id_cp_interaction", "interacting_pair", "celltype.pair", "mean", "pval")


covid$means.long <- pivot_longer(data = covid$means, 
                                 cols = 12:ncol(covid$means), 
                                 names_to = "celltype.pair", 
                                 values_to = "mean")
covid$pvals.long <- pivot_longer(data = covid$pvals, 
                                 cols = 12:ncol(covid$pvals), 
                                 names_to = "celltype.pair", 
                                 values_to = "pval")
covid$means.pval <- inner_join(covid$means.long, covid$pvals.long) %>% 
  select("id_cp_interaction", "interacting_pair", "celltype.pair", "mean", "pval")

### exclude non-significant interactions ======================================
cntrl$means.pval <- filter(cntrl$means.pval, pval < 0.01)
covid$means.pval <- filter(covid$means.pval, pval < 0.01)

### interactions specific to covid ============================================
covid$means.pval <- covid$means.pval %>%  
  mutate("interaction_celltype.pair" = paste0(id_cp_interaction, "_", celltype.pair)) 
cntrl$means.pval <- cntrl$means.pval %>% 
  mutate("interaction_celltype.pair" = paste0(id_cp_interaction, "_", celltype.pair))

covid$unique <- covid$means.pval %>% 
  filter(!(interaction_celltype.pair %in% cntrl$means.pval$interaction_celltype.pair)) %>% 
  select(-c("interaction_celltype.pair"))

### filter to keep only immune cells ==========================================
celltypes.keep <- c("dec.NK_1", "dec.NK_2", "dec.Tcell_1", "dec.Tcell_2", "dec.APC", "dec.Mono_1", "dec.Mono_2")

pairs.keep <- c(
  paste0(celltypes.keep[1], ".", celltypes.keep),
  paste0(celltypes.keep[2], ".", celltypes.keep),
  paste0(celltypes.keep[3], ".", celltypes.keep),
  paste0(celltypes.keep[4], ".", celltypes.keep),
  paste0(celltypes.keep[5], ".", celltypes.keep),
  paste0(celltypes.keep[6], ".", celltypes.keep),
  paste0(celltypes.keep[7], ".", celltypes.keep)
)

covid$unique.immune <- covid$unique %>% 
  filter(celltype.pair %in% pairs.keep)

### get celltype and gene pair names in order =================================
# celltype pair and gene pair labels when plotted as they are, are hard to read. We want to separate the labels so that we can plot them in different colors or add a prominent spacer that can help with reading the labels. 

## separate celltype pair into individual celltypes ----------------------------
covid$unique.immune$celltype1 <- gsub("(dec\\..*)(\\.)(dec.*)", 
                                      "\\1", 
                                      covid$unique.immune$celltype.pair)
covid$unique.immune$celltype2 <- gsub("(dec\\..*)(\\.)(dec.*)", 
                                      "\\3", 
                                      covid$unique.immune$celltype.pair)

## separate interacting pair into individual molecules ------------------------
covid$unique.immune$gene1 <- gsub("(.*)(_)(.*)", 
                                  "\\1", 
                                  covid$unique.immune$interacting_pair)
covid$unique.immune$gene2 <- gsub("(.*)(_)(.*)", 
                                  "\\3", 
                                  covid$unique.immune$interacting_pair)

### prepare data for coloring pair items differently ==========================
# using ggtext: https://github.com/wilkelab/ggtext

# add colour to genes
covid$unique.immune$interacting_pair <- glue("<span style='color:#a50026
'>{covid$unique.immune$gene1}</span> <span style='color:#bdbdbd'>&</span> <span style='color:#313695'>{covid$unique.immune$gene2}</span>")

# add colour to celltypes
covid$unique.immune$celltype.pair <- glue("<span style='color:#a50026
'>{covid$unique.immune$celltype1}</span> <span style='color:#bdbdbd'>&</span> <span style='color:#313695'>{covid$unique.immune$celltype2}</span>")

### plot ======================================================================
cpdb <- ggplot(covid$unique.immune, aes(celltype.pair, interacting_pair)) +
  geom_point(aes(color = log(mean), size = -log10(pval+0.001)), shape = 19) +
  scale_color_gradient(
    low = "#fdd0a2", high = "#d94801",
    name = "log(mean expr.)",
    guide = guide_colorbar(ticks = TRUE, 
                           ticks.colour = "black",
                           frame.colour = "black", 
                           frame.linewidth = 0.25)
  ) +
  scale_size(range = c(0, 1.5), name = "-log<sub>10</sub>(*p* val.)") +
  theme_bw() +
  theme(
    panel.grid = element_line(size = 0.2),
    panel.border = element_rect(size = 0.25, color = "black"),
    axis.ticks = element_line(size = 0.25, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_markdown(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_markdown(size = 5),
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_markdown(size = 5),
    legend.text = element_text(size = 5)
  )

cowplot::ggsave2(
  cpdb,
  filename = "results/99_paper-figures/supp-fig_cellphonedb/supp-fig_cellphonedb_v1.pdf",
  width = 5.75 , height = 10, units = "in"
)

### end =======================================================================
