library(scRNAseq)
library(scater)
library(scran)
library(Glimma)
library(edgeR)
library(Seurat)
library(dplyr)

so.renamed <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")
so.sfg <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")


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
readr::write_csv(f.markers_Non_Vip, "gene_markers_per_markers_Sst.csv")


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
readr::write_csv(f.markers_Non_Vip, "gene_markers_per_markers_Ex_3.csv")


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
readr::write_csv(f.markers_Non_Vip, "gene_markers_per_markers_Ex_4.csv")

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
readr::write_csv(f.markers_Non_Vip, "gene_markers_per_markers_Ex_5.csv")

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

saveRDS(output_rorb, "results_pathfinder/output_rorb.rds")

jpeg("images/enrichChart_RORB_ADvsCt.jpeg", units="in", width=15, height=10, res=300)
enrichment_chart(output_rorb, top_terms = 20) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

# term_gene_heatmap(output_rorb, use_description = TRUE, num_terms = 5) +
#   theme(text = element_text(size = 20))

jpeg("images/UpSet_RORB_ADvsCt.jpeg", units="in", width=15, height=10, res=300)
UpSet_plot(output_rorb, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

visualize_terms(
  result_df = example_pathfindR_output,
  input_processed = input_processed,
  hsa_KEGG = TRUE
)


get_pathfinder <- function(path="KEGG"){
  suppressMessages(library(pathfindR))
  suppressMessages(library(ggplot2))

  output_1 <- list()

  for(i in list.files("to_analyse/")){
    cell_type <- i

    f.markers <- readr::read_csv(paste0("to_analyse/",i))
    markers_ex1 <- data.frame(Gene.symbol = f.markers$gene,
                              logFC = f.markers$avg_log2FC,
                              adj.P.Val = f.markers$p_val_adj)

    dir_path <- paste0("results_pathfinder","_", cell_type, "_",path)
    output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = path)

    saveRDS(output_, paste0(cell_type, "_", path, ".rds"))

    name_enrichment_chart <- paste0("images/", cell_type,"_",path, "_","enrichment_chart.jpeg")
    jpeg(name_enrichment_chart, units="in", width=25, height=17, res=300)
    print(enrichment_chart(output_, top_terms = 15) +
            theme(text = element_text(size = 25),
                  axis.text.x = element_text(size = 20),
                  axis.text.y = element_text(size = 20)))
    dev.off()

    name_term_gene_heatmap <- paste0("images/", cell_type,"_",path, "_",path, "_","gene_heatmap.jpeg")
    jpeg(name_term_gene_heatmap, units="in", width=25, height=17, res=300)
    print(term_gene_heatmap(output_, use_description = TRUE, num_terms = 15) +
            theme(text = element_text(size = 20)))
    dev.off()

    name_term_gene_graph <- paste0("images/", cell_type,"_",path, "_","gene_graph.jpeg")
    jpeg(name_term_gene_graph, units="in", width=25, height=17, res=300)
    print(term_gene_graph(output_, use_description = TRUE, num_terms = 15) +
            theme(text = element_text(size = 20)))
    dev.off()

    name_UpSet_plot <- paste0("images/", cell_type,"_",path, "_","UpSet_plot.jpeg")
    jpeg(name_UpSet_plot, units="in", width=25, height=17, res=300)
    print(UpSet_plot(output_, num_terms = 10, use_description = TRUE,
                     method = "barplot") +
            theme(text = element_text(size = 20),
                  axis.text.x = element_text(size = 20),
                  axis.text.y = element_text(size = 20)))
    dev.off()
  }



  return(output_)
}


get_pathfinder()
get_pathfinder(path = "Reactome")
get_pathfinder(path = "GO-All")

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


##################### DEG per celltype

ex1 <- readr::read_csv("to_analyse/gene_markers_per_markers_Ex1.csv")
ex1$cluster <- "ex1"
ex2 <- readr::read_csv("to_analyse/gene_markers_per_markers_Ex2.csv")
ex2$cluster <- "ex2"
ex3 <- readr::read_csv("to_analyse/gene_markers_per_markers_Ex_3.csv")
ex3$cluster <- "ex3"

ex4 <- readr::read_csv("to_analyse/gene_markers_per_markers_Ex_4.csv")
ex4$cluster <- "ex4"

ex5 <- readr::read_csv("to_analyse/gene_markers_per_markers_Ex_5.csv")
ex5$cluster <- "ex5"

vip <- readr::read_csv("to_analyse/gene_markers_per_markers_Vip.csv")
vip$cluster <- "vip"

nonvip <- readr::read_csv("to_analyse/gene_markers_per_markers_Non_Vip.csv")
nonvip$cluster <- "nonvip"

sst <- readr::read_csv("to_analyse/gene_markers_per_markers_Sst.csv")
sst$cluster <- "sst"

pv <- readr::read_csv("to_analyse/gene_markers_per_markers_Pv.csv")
pv$cluster <- "pv"

rorb <- readr::read_csv("to_analyse/gene_markers_per_markers_RORB.csv")
rorb$cluster <- "rorb"


all_up <- rbind(ex1, ex2, ex3, ex4, ex5, vip, nonvip, sst, pv, rorb) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC >= 0.5)

all_down <- rbind(ex1, ex2, ex3, ex4, ex5, vip, nonvip, sst, pv, rorb) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC <= -0.5)

table(all_down$cluster)
table(all_up$cluster)

readr::write_csv(all_up, "DEG_up_ADvsCt_perCelltype.csv")

readr::write_csv(all_down, "DEG_down_ADvsCt_perCelltype.csv")


suppressMessages(library(pathfindR))
suppressMessages(library(ggplot2))


select_genes_rorb <- readr::read_csv("deg_unicos_RORB.csv")

rorb <- rorb[rorb$gene %in% select_genes_rorb$gene,]

get_pathfinder <- function(f.markers){

  path <- c("KEGG", "Reactome", "GO-All")
  suppressMessages(library(pathfindR))
  suppressMessages(library(ggplot2))

  output_1 <- list()

  for(i in path){
    cell_type <- "RORB"
    markers_ex1 <- data.frame(Gene.symbol = f.markers$gene,
                              logFC = f.markers$avg_log2FC,
                              adj.P.Val = f.markers$p_val)

    dir_path <- paste0("results_pathfinder","_", cell_type, "_",i, "_unicos")
    output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = i)

    saveRDS(output_, paste0(cell_type, "_", i, ".rds"))

    name_enrichment_chart <- paste0("images/", cell_type,"_",i, "_","enrichment_chart.jpeg")
    jpeg(name_enrichment_chart, units="in", width=25, height=17, res=300)
    print(enrichment_chart(output_, top_terms = 15) +
            theme(text = element_text(size = 25),
                  axis.text.x = element_text(size = 25),
                  axis.text.y = element_text(size = 25)))
    dev.off()

    name_term_gene_heatmap <- paste0("images/", cell_type,"_",i, "_","gene_heatmap.jpeg")
    jpeg(name_term_gene_heatmap, units="in", width=25, height=17, res=300)
    print(term_gene_heatmap(output_, use_description = TRUE, num_terms = 15) +
            theme(text = element_text(size = 25)))
    dev.off()

    name_term_gene_graph <- paste0("images/", cell_type,"_",i, "_","gene_graph.jpeg")
    jpeg(name_term_gene_graph, units="in", width=25, height=17, res=300)
    print(term_gene_graph(output_, use_description = TRUE, num_terms = 15) +
            theme(text = element_text(size = 25)))
    dev.off()

    name_UpSet_plot <- paste0("images/", cell_type,"_",i, "_","UpSet_plot.jpeg")
    jpeg(name_UpSet_plot, units="in", width=25, height=17, res=300)
    print(UpSet_plot(output_, num_terms = 10, use_description = TRUE,
                     method = "barplot") +
            theme(text = element_text(size = 25),
                  axis.text.x = element_text(size = 25),
                  axis.text.y = element_text(size = 25)))
    dev.off()
  }



  return(output_)
}


get_pathfinder(rorb)


library(ggplot2)

so.renamed <- ScaleData(so.renamed)
# Single cell heatmap of feature expression
DoHeatmap(subset(so.renamed), features = c("ZBTB7A", "JUND", "SMARCC2", "NFIB", "THRB", "NFE2L1", "NFIC", "EEF1A2", "NEUROD2", "SYN1"),
          size = 3,group.by = "ident")

# Plot single genes

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(so.renamed, features = c("SYN1", "STX1A", "STX1B", "CPLX1", "DLG4", "PPFIA3"),
          idents = c("RORB"), group.by = "dataset")

RidgePlot(so.renamed, features = c("ZBTB7A", "JUND", "SMARCC2", "NFIB", "THRB", "NFE2L1", "NFIC", "EEF1A2", "NEUROD2", "SYN1"),
          idents = c("RORB"), group.by = "group_id")

RidgePlot(so.renamed, features = c("ZBTB7A", "JUND", "SMARCC2", "NFIB", "THRB", "NFE2L1", "NFIC", "EEF1A2", "NEUROD2", "SYN1"),
          idents = c("Sst"), group.by = "braak")

# Nuclear events
library("wesanderson")
RidgePlot(so.renamed, features = c("EGR1","PEG3", "RORB", "JUND", "PRNP","CUX1", "YWHAZ"),
          idents = c("RORB"), group.by = "braak") +
  theme(text = element_text(size = 25),
        title = element_text(size = 25)) +
  labs(x = NULL, y = NULL, color = NULL, alt_insight = NULL) +
  scale_color_manual(values=wes_palette("Rushmore1"))

RidgePlot(so.renamed, features = c("EGR1"),
          idents = c("RORB"), group.by = "braak") +
  theme(text = element_text(size = 25),
        title = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25)) +
  labs(y = NULL, color = NULL, alt_insight = NULL)+
  scale_color_manual(values=wes_palette("Rushmore1"))


# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(so.renamed, features = "SYN1")

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(so.renamed, features = c("PEG3", "RORB", "JUND", "PRNP","CUX1", "YWHAZ"),
        idents = c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5", "RORB")) + RotatedAxis()

DotPlot(so.renamed, features = c("PEG3", "RORB", "JUND", "PRNP","CUX1", "YWHAZ"),
        idents = c("RORB"),
        group.by = "braak",
        cols = c("lightgrey", "red"))+ RotatedAxis() +
  theme(text = element_text(size = 20))

RidgePlot(so.renamed, features = c("PRNP"),
          idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(y = NULL, color = NULL, alt_insight = NULL) +
  scale_color_manual(values=wes_palette("Rushmore1"))

VlnPlot(so.renamed, features = c("PRNP"),
          idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(color = NULL, alt_insight = NULL) +
  scale_color_manual(values=wes_palette("Rushmore1"))


RidgePlot(so.sfg, features = c("FOXP1"),
          idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(y = NULL, color = NULL, alt_insight = NULL)

VlnPlot(so.sfg, features = c("FOXP1"),
        idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(color = NULL, alt_insight = NULL)


RidgePlot(so.sfg, features = c("PKM"),
          idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(y = NULL, color = NULL, alt_insight = NULL)

VlnPlot(so.sfg, features = c("PKM"),
        idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(color = NULL, alt_insight = NULL)

