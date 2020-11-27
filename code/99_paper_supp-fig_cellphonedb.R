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

### update sep characters =====================================================
# CellPhoneDB, when creating names for columns in means and pvals output files, join cell type names with a period, e.g. celltype1.celltype2. But our celltype names have periods (as well as underscores) in them. We’ll replace the periods separating cell types with a double underscore, "__".

# function
replace_sep.char <- function(x) {
  gsub("([a-z]{3}\\.[A-z]+_*\\d*)(\\.)([a-z]{3}\\.[A-z]+_*\\d*)", 
       "\\1__\\3", 
       x)
}

# replace
colnames(cntrl$means) <- replace_sep.char(names(cntrl$means))
colnames(cntrl$pvals) <- replace_sep.char(names(cntrl$pvals))

colnames(covid$means) <- replace_sep.char(names(covid$means))
colnames(covid$pvals) <- replace_sep.char(names(covid$pvals))

### Update cluster annotations ================================================
# Since we ran CellPhoneDB, we decided to update immune cell cluster labels to by removing their “dec.” prefix to convey the uncertaintly about their tissue of origin. We need to replace the old labels with the new labels.

# function to replace labels
replace_labels <- function(x) {
  
  x <- gsub("dec.APC", "APC", x)
  x <- gsub("dec.Bcells", "Bcell", x)
  x <- gsub("dec.Tcell_1", "Tcell_1", x)
  x <- gsub("dec.Tcell_2", "Tcell_2", x)
  x <- gsub("dec.Tcell_3", "Tcell_3", x)
  x <- gsub("dec.NK_1", "NK_1", x)
  x <- gsub("dec.NK_2", "NK_2", x)
  x <- gsub("dec.NK_3", "NK_3", x)
  x <- gsub("dec.Mono_1", "Mono_1", x)
  x <- gsub("dec.Mono_2", "Mono_2", x)
  x <- gsub("dec.Gran", "Gran", x)
  
  return(x)
  
}

# replace
colnames(cntrl$means) <- replace_labels(names(cntrl$means))
colnames(cntrl$pvals) <- replace_labels(names(cntrl$pvals))

colnames(covid$means) <- replace_labels(names(covid$means))
colnames(covid$pvals) <- replace_labels(names(covid$pvals))

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
celltypes.keep <- c("NK_1", "NK_2", "Tcell_1", "Tcell_2", "APC", "Mono_1", "Mono_2")

pairs.keep <- c(
  paste0(celltypes.keep[1], "__", celltypes.keep),
  paste0(celltypes.keep[2], "__", celltypes.keep),
  paste0(celltypes.keep[3], "__", celltypes.keep),
  paste0(celltypes.keep[4], "__", celltypes.keep),
  paste0(celltypes.keep[5], "__", celltypes.keep),
  paste0(celltypes.keep[6], "__", celltypes.keep),
  paste0(celltypes.keep[7], "__", celltypes.keep)
)

covid$unique.immune <- covid$unique %>% 
  filter(celltype.pair %in% pairs.keep)

### get celltype and gene pair names in order =================================
# celltype pair and gene pair labels when plotted as they are, are hard to read. We want to separate the labels so that we can plot them in different colors or add a prominent spacer that can help with reading the labels. 

## separate celltype pair into individual celltypes ----------------------------
covid$unique.immune$celltype1 <- gsub("([A-z]+_*\\d*)(__)([A-z]+_*\\d*)", 
                                      "\\1", 
                                      covid$unique.immune$celltype.pair)

covid$unique.immune$celltype2 <- gsub("([A-z]+_*\\d*)(__)([A-z]+_*\\d*)", 
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
covid$unique.immune$interacting_pair <- glue("<span style='color:#a50026'>{covid$unique.immune$gene1}</span> <span style='color:#bdbdbd'>&</span> <span style='color:#313695'>{covid$unique.immune$gene2}</span>")

# add colour to celltypes
covid$unique.immune$celltype.pair <- glue("<span style='color:#a50026'>{covid$unique.immune$celltype1}</span> <span style='color:#bdbdbd'>&</span> <span style='color:#313695'>{covid$unique.immune$celltype2}</span>")

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
  filename = "results/99_paper-figures/supp-fig_cellphonedb/supp-fig_cellphonedb_v2.pdf",
  width = 5.75 , height = 10, units = "in"
)

cowplot::ggsave2(
  cpdb,
  filename = "results/99_paper-figures/supp-fig_cellphonedb/supp-fig_cellphonedb_v2.png",
  width = 5.75 , height = 10, units = "in", type = "cairo", dpi = 600
)


### end =======================================================================
