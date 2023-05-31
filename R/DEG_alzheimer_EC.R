library(scRNAseq)
library(scater)
library(scran)
library(Glimma)
library(edgeR)
library(Seurat)
library(dplyr)

so.renamed <- readRDS("EC_neurons_annoted_from_SFG.rds")

table(Idents(so.renamed), so.renamed$group_id)

new_labels <- paste0(so.renamed$group_id, "_", Idents(so.renamed))

Idents(so.renamed) <- new_labels

f.markers_RORB <- FindMarkers(so.renamed,
                              ident.1 = "AD_RORB",
                              ident.2 = "Control_RORB",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_RORB$gene <- rownames(f.markers_RORB)
readr::write_csv(f.markers_RORB, "EC_DEG/gene_markers_per_markers_RORB.csv")


f.markers_Ex1 <- FindMarkers(so.renamed,
                             ident.1 = "AD_Ex_1",
                             ident.2 = "Control_Ex_1",
                             min.cells.group = 1,
                             min.cells.feature = 1,
                             min.pct = 0,
                             logfc.threshold = 0,
                             only.pos = FALSE)


f.markers_Ex1$gene <- rownames(f.markers_Ex1)
readr::write_csv(f.markers_Ex1, "EC_DEG/gene_markers_per_markers_Ex1.csv")

f.markers_Vip <- FindMarkers(so.renamed,
                             ident.1 = "AD_Vip",
                             ident.2 = "Control_Vip",
                             min.cells.group = 1,
                             min.cells.feature = 1,
                             min.pct = 0,
                             logfc.threshold = 0,
                             only.pos = FALSE)


f.markers_Vip$gene <- rownames(f.markers_Vip)
readr::write_csv(f.markers_Vip, "EC_DEG/gene_markers_per_markers_Vip.csv")
# high <- f.markers_Vip %>% filter(avg_log2FC > 0.5) %>%
#   filter(p_val_adj <= 0.001)

f.markers_Pv <- FindMarkers(so.renamed,
                            ident.1 = "AD_Pv",
                            ident.2 = "Control_Pv",
                            min.cells.group = 1,
                            min.cells.feature = 1,
                            min.pct = 0,
                            logfc.threshold = 0,
                            only.pos = FALSE)


f.markers_Pv$gene <- rownames(f.markers_Pv)
readr::write_csv(f.markers_Pv, "EC_DEG/gene_markers_per_markers_Pv.csv")

f.markers_Ex2 <- FindMarkers(so.renamed,
                             ident.1 = "AD_Ex_2",
                             ident.2 = "Control_Ex_2",
                             min.cells.group = 1,
                             min.cells.feature = 1,
                             min.pct = 0,
                             logfc.threshold = 0,
                             only.pos = FALSE)


f.markers_Ex2$gene <- rownames(f.markers_Ex2)
readr::write_csv(f.markers_Ex2, "EC_DEG/gene_markers_per_markers_Ex2.csv")

f.markers_Non_Vip <- FindMarkers(so.renamed,
                                 ident.1 = "AD_Non-Vip",
                                 ident.2 = "Control_Non-Vip",
                                 min.cells.group = 1,
                                 min.cells.feature = 1,
                                 min.pct = 0,
                                 logfc.threshold = 0,
                                 only.pos = FALSE)


f.markers_Non_Vip$gene <- rownames(f.markers_Non_Vip)
readr::write_csv(f.markers_Non_Vip, "EC_DEG/gene_markers_per_markers_Non_Vip.csv")

## Sst
f.markers_Non_Vip <- FindMarkers(so.renamed,
                                 ident.1 = "AD_Sst",
                                 ident.2 = "Control_Sst",
                                 min.cells.group = 1,
                                 min.cells.feature = 1,
                                 min.pct = 0,
                                 logfc.threshold = 0,
                                 only.pos = FALSE)


f.markers_Non_Vip$gene <- rownames(f.markers_Non_Vip)
readr::write_csv(f.markers_Non_Vip, "EC_DEG/gene_markers_per_markers_Sst.csv")


## Ex_3
f.markers_Non_Vip <- FindMarkers(so.renamed,
                                 ident.1 = "AD_Ex_3",
                                 ident.2 = "Control_Ex_3",
                                 min.cells.group = 1,
                                 min.cells.feature = 1,
                                 min.pct = 0,
                                 logfc.threshold = 0,
                                 only.pos = FALSE)


f.markers_Non_Vip$gene <- rownames(f.markers_Non_Vip)
readr::write_csv(f.markers_Non_Vip, "EC_DEG/gene_markers_per_markers_Ex_3.csv")


## Ex_4
f.markers_Non_Vip <- FindMarkers(so.renamed,
                                 ident.1 = "AD_Ex_4",
                                 ident.2 = "Control_Ex_4",
                                 min.cells.group = 1,
                                 min.cells.feature = 1,
                                 min.pct = 0,
                                 logfc.threshold = 0,
                                 only.pos = FALSE)


f.markers_Non_Vip$gene <- rownames(f.markers_Non_Vip)
readr::write_csv(f.markers_Non_Vip, "EC_DEG/gene_markers_per_markers_Ex_4.csv")

## Ex_5
f.markers_Non_Vip <- FindMarkers(so.renamed,
                                 ident.1 = "AD_Ex_5",
                                 ident.2 = "Control_Ex_5",
                                 min.cells.group = 1,
                                 min.cells.feature = 1,
                                 min.pct = 0,
                                 logfc.threshold = 0,
                                 only.pos = FALSE)


f.markers_Non_Vip$gene <- rownames(f.markers_Non_Vip)
readr::write_csv(f.markers_Non_Vip, "EC_DEG/gene_markers_per_markers_Ex_5.csv")



