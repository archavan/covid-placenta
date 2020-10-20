# Following the instructions from here: https://www.cellphonedb.org/faq-and-troubleshooting

setwd("/home/arc78/scratch60/covid-placenta")

### packages ==================================================================
library(Seurat)

### data ======================================================================
seur <- readRDS("results/02_annotation/seurat-object_annotated.rds")

### subset by covid status ====================================================
seur.cntrl <- subset(x = seur, subset = covid == "cntrl")
seur.covid <- subset(x = seur, subset = covid == "covid")

### extract raw data, normalize it, and write it ==============================
# control
count.raw_cntrl <- seur.cntrl@assays$RNA@counts
count.norm_cntrl <- apply(count.raw_cntrl, 2, function(x){(x/sum(x)) * 10000})
write.table(count.norm_cntrl, 
            "results/06_cellphonedb/01_input/cntrl_count.txt", 
            sep = "\t", quote = FALSE)

# covid
count.raw_covid <- seur.covid@assays$RNA@counts
count.norm_covid <- apply(count.raw_covid, 2, function(x){(x/sum(x)) * 10000})
write.table(count.norm_covid, 
            "results/06_cellphonedb/01_input/covid_count.txt", 
            sep = "\t", quote = FALSE)

### write metadata ============================================================
# control
meta_cntrl <- cbind(rownames(seur.cntrl@meta.data), 
                    seur.cntrl@meta.data[, "annotation_merged", drop = FALSE])
write.table(meta_cntrl,
            "results/06_cellphonedb/01_input/cntrl_meta.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# covid
meta_covid <- cbind(rownames(seur.covid@meta.data), 
                    seur.covid@meta.data[, "annotation_merged", drop = FALSE])
write.table(meta_covid,
            "results/06_cellphonedb/01_input/covid_meta.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

### end =======================================================================
