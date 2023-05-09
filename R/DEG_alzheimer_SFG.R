library(scRNAseq)
library(scater)
library(scran)
library(Glimma)
library(edgeR)
library(Seurat)
library(dplyr)

so.renamed <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")

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
readr::write_csv(f.markers_RORB, "gene_markers_per_markers_RORB.csv")


f.markers_Ex1 <- FindMarkers(so.renamed,
                              ident.1 = "AD_Ex_1",
                              ident.2 = "Control_Ex_1",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_Ex1$gene <- rownames(f.markers_Ex1)
readr::write_csv(f.markers_Ex1, "gene_markers_per_markers_Ex1.csv")

f.markers_Vip <- FindMarkers(so.renamed,
                             ident.1 = "AD_Vip",
                             ident.2 = "Control_Vip",
                             min.cells.group = 1,
                             min.cells.feature = 1,
                             min.pct = 0,
                             logfc.threshold = 0,
                             only.pos = FALSE)


f.markers_Vip$gene <- rownames(f.markers_Vip)
readr::write_csv(f.markers_Vip, "gene_markers_per_markers_Vip.csv")
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
readr::write_csv(f.markers_Pv, "gene_markers_per_markers_Pv.csv")

f.markers_Ex2 <- FindMarkers(so.renamed,
                            ident.1 = "AD_Ex_2",
                            ident.2 = "Control_Ex_2",
                            min.cells.group = 1,
                            min.cells.feature = 1,
                            min.pct = 0,
                            logfc.threshold = 0,
                            only.pos = FALSE)


f.markers_Ex2$gene <- rownames(f.markers_Ex2)
readr::write_csv(f.markers_Ex2, "gene_markers_per_markers_Ex2.csv")

f.markers_Non_Vip <- FindMarkers(so.renamed,
                             ident.1 = "AD_Non-Vip",
                             ident.2 = "Control_Non-Vip",
                             min.cells.group = 1,
                             min.cells.feature = 1,
                             min.pct = 0,
                             logfc.threshold = 0,
                             only.pos = FALSE)


f.markers_Non_Vip$gene <- rownames(f.markers_Non_Vip)
readr::write_csv(f.markers_Non_Vip, "gene_markers_per_markers_Non_Vip.csv")

rorb_grpah <- readr::read_csv("../GRN_atac/Nets/Graph_rorb_PAPER_uniques.csv")

markers <- FindMarkers(so.renamed, ident.1 = "g1", group.by = 'group_id', subset.ident = "2")



## Grafica DEG para las células RORB+

f.markers_RORB <- readr::read_csv("gene_markers_per_markers_RORB.csv")
library(dplyr)
library(ggplot2)

f.markers_RORB <- f.markers_RORB[!stringr::str_detect(f.markers_RORB$gene, "MT"),]


high <- f.markers_RORB %>% filter(avg_log2FC > 0.5) %>%
  filter(p_val_adj <= 0.001)

high <- high[!stringr::str_detect(high$gene, "MT"),] %>% head(15)

low <- f.markers_RORB %>% filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj <= 0.001)

low <- low[!stringr::str_detect(low$gene, "MT"),] %>% head(15)

# rorb_tg <- rorb_grpah[rorb_grpah$TG %in% low$gene,]
# rorb_tf <- rorb_grpah[rorb_grpah$TF %in% low$gene,]


# saveRDS(prueba, "resultados/DE_Ex_SFG_ctrl_AD.rds")

jpeg("images/DEG_RORB_AD_ct.jpeg", units="in", width=15, height=10, res=300)

ggplot(f.markers_RORB) +
  geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj)),
             color = "gray", cex = 2, alpha = 0.5) +
  theme_classic() +
  geom_hline(yintercept = 3,
             col = "red",
             linetype = "dotted",
             size = 1) +
  geom_vline(xintercept = c(0.5, -0.5),
             col = "red",
             linetype = "dashed",
             size = 0.5) +
  geom_point(data = high, aes(x = avg_log2FC, y = -log10(p_val_adj)),
             color = "dark red", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = high, aes(x = avg_log2FC, y = -log10(p_val_adj)), label = high$gene,
                           color = "black", size = 5, max.overlaps = 20) +
  geom_point(data = low, aes(x = avg_log2FC, y = -log10(p_val_adj)),
             color = "dark blue", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = low, aes(x = avg_log2FC, y = -log10(p_val_adj)), label = low$gene,
                           color = "black", size = 5,  max.overlaps = 158) +
  theme(text = element_text(size = 22)) +
  ylab("-Log10(FDR)") + xlab("LogFC")

dev.off()


# Enrich plots

input_pathfinder <- data.frame(Gene.symbol = f.markers_RORB$gene, logFC = f.markers_RORB$avg_log2FC,
                               adj.P.Val = f.markers_RORB$p_val_adj)

library(pathfindR)
output_rorb <- run_pathfindR(input_pathfinder, output_dir = "results_pathfinder/", gene_sets = "GO-All")


#
#
#
# contr <- makeContrasts("Inh_Ct - Inh_AD", levels = design)
#
# pb_Qfit <- glmQLFit(pb_dge, design)
# pb_Qlrt <- glmQLFTest(pb_Qfit, contrast = contr)
#
# prueba <- pb_Qlrt$table
# fdr <- p.adjust(prueba$PValue, "fdr")
# prueba$FDR <- fdr
# prueba$genes <- row.names(prueba)
#
# high <- prueba %>% filter(logFC > 2) %>%
#   filter(PValue <= 0.005)
#
# low <- prueba %>% filter(logFC < -2) %>%
#   filter(PValue <= 0.005)
#
# saveRDS(prueba, "resultados/DE_Inh_SFG_ctrl_AD.rds")
#
# jpeg("images/Inh_neurons_SFG_DE_ct_AD.jpeg", units="in", width=15, height=10, res=300)
#
# ggplot(prueba) +
#   geom_point(aes(x = logFC, y = -log10(PValue)),
#              color = "gray", cex = 2, alpha = 0.5) +
#   theme_classic() +
#   geom_hline(yintercept = 2.30103,
#              col = "red",
#              linetype = "dotted",
#              size = 1) +
#   geom_vline(xintercept = c(2, -2),
#              col = "red",
#              linetype = "dashed",
#              size = 0.5) +
#   geom_point(data = high, aes(x = logFC, y = -log10(PValue)),
#              color = "dark red", cex = 1.5, alpha = 0.5) +
#   ggrepel::geom_text_repel(data = high, aes(x = logFC, y = -log10(PValue)), label = high$genes,
#                            color = "black", size = 5, max.overlaps = 20) +
#   geom_point(data = low, aes(x = logFC, y = -log10(PValue)),
#              color = "dark blue", cex = 1.5, alpha = 0.5) +
#   ggrepel::geom_text_repel(data = low, aes(x = logFC, y = -log10(PValue)), label = low$genes,
#                            color = "black", size = 5,  max.overlaps = 50) +
#   theme(text = element_text(size = 22)) +
#   ylab("-Log10(PValue)") + xlab("LogFC") +
#   coord_flip()
#
# dev.off()
#
#
# contr <- makeContrasts("Pvalb_Ct - Pvalb_AD", levels = design)
#
# pb_Qfit <- glmQLFit(pb_dge, design)
# pb_Qlrt <- glmQLFTest(pb_Qfit, contrast = contr)
#
# prueba <- pb_Qlrt$table
# fdr <- p.adjust(prueba$PValue, "fdr")
# prueba$FDR <- fdr
# prueba$genes <- row.names(prueba)
#
# high <- prueba %>% filter(logFC > 2) %>%
#   filter(PValue <= 0.005)
#
# low <- prueba %>% filter(logFC < -2) %>%
#   filter(PValue <= 0.005)
#
# saveRDS(prueba, "resultados/DE_Pvalb_SFG_ctrl_AD.rds")
#
# jpeg("images/Pvalb_neurons_SFG_DE_ct_AD.jpeg", units="in", width=15, height=10, res=300)
#
# ggplot(prueba) +
#   geom_point(aes(x = logFC, y = -log10(PValue)),
#              color = "gray", cex = 2, alpha = 0.5) +
#   theme_classic() +
#   geom_hline(yintercept = 2.30103,
#              col = "red",
#              linetype = "dotted",
#              size = 1) +
#   geom_vline(xintercept = c(2, -2),
#              col = "red",
#              linetype = "dashed",
#              size = 0.5) +
#   geom_point(data = high, aes(x = logFC, y = -log10(PValue)),
#              color = "dark red", cex = 1.5, alpha = 0.5) +
#   ggrepel::geom_text_repel(data = high, aes(x = logFC, y = -log10(PValue)), label = high$genes,
#                            color = "black", size = 5, max.overlaps = 20) +
#   geom_point(data = low, aes(x = logFC, y = -log10(PValue)),
#              color = "dark blue", cex = 1.5, alpha = 0.5) +
#   ggrepel::geom_text_repel(data = low, aes(x = logFC, y = -log10(PValue)), label = low$genes,
#                            color = "black", size = 5,  max.overlaps = 50) +
#   theme(text = element_text(size = 22)) +
#   ylab("-Log10(PValue)") + xlab("LogFC") +
#   coord_flip()
#
# dev.off()
#
#
# contr <- makeContrasts("Vip_Ct - Vip_AD", levels = design)
#
# pb_Qfit <- glmQLFit(pb_dge, design)
# pb_Qlrt <- glmQLFTest(pb_Qfit, contrast = contr)
#
# prueba <- pb_Qlrt$table
# fdr <- p.adjust(prueba$PValue, "fdr")
# prueba$FDR <- fdr
# prueba$genes <- row.names(prueba)
#
# high <- prueba %>% filter(logFC > 2) %>%
#   filter(PValue <= 0.005)
#
# low <- prueba %>% filter(logFC < -2) %>%
#   filter(PValue <= 0.005)
#
# saveRDS(prueba, "resultados/DE_Vip_SFG_ctrl_AD.rds")
#
# jpeg("images/Vip_neurons_SFG_DE_ct_AD.jpeg", units="in", width=15, height=10, res=300)
#
# ggplot(prueba) +
#   geom_point(aes(x = logFC, y = -log10(PValue)),
#              color = "gray", cex = 2, alpha = 0.5) +
#   theme_classic() +
#   geom_hline(yintercept = 2.30103,
#              col = "red",
#              linetype = "dotted",
#              size = 1) +
#   geom_vline(xintercept = c(2, -2),
#              col = "red",
#              linetype = "dashed",
#              size = 0.5) +
#   geom_point(data = high, aes(x = logFC, y = -log10(PValue)),
#              color = "dark red", cex = 1.5, alpha = 0.5) +
#   ggrepel::geom_text_repel(data = high, aes(x = logFC, y = -log10(PValue)), label = high$genes,
#                            color = "black", size = 5, max.overlaps = 20) +
#   geom_point(data = low, aes(x = logFC, y = -log10(PValue)),
#              color = "dark blue", cex = 1.5, alpha = 0.5) +
#   ggrepel::geom_text_repel(data = low, aes(x = logFC, y = -log10(PValue)), label = low$genes,
#                            color = "black", size = 5,  max.overlaps = 50) +
#   theme(text = element_text(size = 22)) +
#   ylab("-Log10(PValue)") + xlab("LogFC") +
#   coord_flip()
#
# dev.off()
#
