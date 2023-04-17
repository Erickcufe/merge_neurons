library(Seurat)
library(dplyr)
library(ggplot2)

so <- readRDS("anotation_SFG.rds")
Idents(so)

f.markers <- FindAllMarkers(so, min.pct = 0.25, logfc.threshold = 0.25)

f.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

library(viridis)

jpeg("images/markers_CellTypes.jpeg", units="in", width=15, height=10, res=300)
DoHeatmap(so, features = top10$gene) +
  theme(text = element_text(size = 18)) +
  scico::scale_fill_scico(palette = "imola")
dev.off()


so_neuron_merge <- ScaleData(so_neuron_merge)
f.markers <- FindAllMarkers(so_neuron_merge_copy, min.pct = 0.25, logfc.threshold = 0.25)
# clusters_markers <- read.csv("../Datos_scRNA/neurons_integrated/SFG/so_neuron_merge_all-ref_MGZS_20PC_res0.4_genes-RNA.csv")
f.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

jpeg("images/SFG_markers_Clusters.jpeg", units="in", width=25, height=17, res=300)
DoHeatmap(so_neuron_merge, features = top10$gene) +
  theme(text = element_text(size = 18)) +
  scico::scale_fill_scico(palette = "imola")
dev.off()


library(pathfindR)

rorb.markers <- f.markers[f.markers$cluster=="RORB+",]
rorb.markers <- data.frame(Gene.symbol = rownames(rorb.markers),
                           logFC = rorb.markers$avg_log2FC,
                           adj.P.Val = rorb.markers$p_val_adj)

genes <- sapply(rorb.markers$Gene.symbol, FUN = function(x){
  a <- stringr::str_split(string = x, pattern = "[.]")
  b <- a[[1]][[1]]
})
rorb.markers$Gene.symbol <- genes

output_rorb <- run_pathfindR(rorb.markers, output_dir = "results_pathfinder/", gene_sets = "GO-All")

enrichment_chart(output_rorb)
term_gene_heatmap(output_rorb, use_description = TRUE, num_terms = 20) +
  theme(text = element_text(size = 20)) + NoLegend()

term_gene_graph(output_rorb, use_description = TRUE, num_terms = 20)

UpSet_plot(output_rorb, num_terms = 15, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20)) + NoLegend()

library(clusterProfiler)
library(enrichplot)

library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db)

tum_down  <- subset(rorb.markers,
                    rorb.markers$logFC < -0.5
                    & rorb.markers$adj.P.Val < 0.05)
# rorb_down <- tum_down[tum_down$cluster=="RORB+",]
rorb_down_genes <- tum_down$Gene.symbol

tum_up <- subset(rorb.markers,
                 rorb.markers$logFC > 0.5
                 & rorb.markers$adj.P.Val < 0.05)
# rorb_up <- tum_up[tum_up$cluster=="RORB+",]
rorb_up_genes <- tum_up$Gene.symbol

# DOWN
rorb_vs_norm_go <- clusterProfiler::enrichGO(rorb_down_genes,
                                            "org.Hs.eg.db",
                                            keyType = "SYMBOL",
                                            ont = "BP",
                                            minGSSize = 50)

enr_go <- clusterProfiler::simplify(rorb_vs_norm_go)
View(enr_go@result)

jpeg("images/GO_down_RORB.jpeg", units="in", width=15, height=10, res=300)
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
                     showCategory = 30, cex_label_category = 0.7)
dev.off()



# UP
rorb_vs_norm_go_UP <- clusterProfiler::enrichGO(rorb_up_genes,
                                             "org.Hs.eg.db",
                                             keyType = "SYMBOL",
                                             ont = "BP",
                                             minGSSize = 50)

enr_go_up <- clusterProfiler::simplify(rorb_vs_norm_go_UP)

jpeg("images/GO_up_RORB.jpeg", units="in", width=15, height=10, res=300)
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go_up),
                     showCategory = 30, cex_label_category = 0.7)
dev.off()

