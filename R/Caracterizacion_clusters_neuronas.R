library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

# Caracterizacion y eleccion de clusters para paper

so <- readRDS("anotacion_Parcial_neuronas.rds")
f.markers <- FindAllMarkers(so, min.pct = 0.25, logfc.threshold = 0.25)
readr::write_csv(f.markers, "gene_markers_per_cluster.csv")
f.markers <- readr::read_csv("gene_markers_per_cluster.csv")

library(pathfindR)

get_pathfinder <- function(cell_type, markers = f.markers, path="KEGG"){
  suppressMessages(library(pathfindR))
  markers_ex1 <- f.markers[f.markers$cluster==cell_type,]
  markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                            logFC = markers_ex1$avg_log2FC,
                            adj.P.Val = markers_ex1$p_val_adj)

  dir_path <- paste0("results_pathfinder","_", cell_type, "_",path)
  output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = path)

  name_enrichment_chart <- paste0("images/", cell_type, "_","enrichment_chart.jpeg")
  jpeg(name_enrichment_chart, units="in", width=25, height=17, res=300)
  enrichment_chart(output_) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
  dev.off()

  name_term_gene_heatmap <- paste0("images/", cell_type, "_","gene_heatmap.jpeg")
  jpeg(name_term_gene_heatmap, units="in", width=25, height=17, res=300)
  term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
    theme(text = element_text(size = 20))
  dev.off()

  name_term_gene_graph <- paste0("images/", cell_type, "_","gene_graph.jpeg")
  jpeg(name_term_gene_graph, units="in", width=25, height=17, res=300)
  term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
    theme(text = element_text(size = 20))
  dev.off()

  name_UpSet_plot <- paste0("images/", cell_type, "_","UpSet_plot.jpeg")
  jpeg(name_UpSet_plot, units="in", width=25, height=17, res=300)
  UpSet_plot(output_, num_terms = 10, use_description = TRUE,
             method = "barplot") +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
  dev.off()

  return(output_)
}


results <- purrr::map(unique(f.markers$cluster), get_pathfinder)
names(results) <- unique(f.markers$cluster)
save(results, file = "resultados_PathFindR.rda")


## Ex_1
markers_ex1 <- f.markers[f.markers$cluster=="Ex_1",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Ex_1")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Ex_1_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Ex_1_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_1_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_1_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
          method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()


#### Ex_2

markers_ex1 <- f.markers[f.markers$cluster=="Ex_2",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Ex_2")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Ex_2_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Ex_2_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_2_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_2_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()


#### Ex_5

markers_ex1 <- f.markers[f.markers$cluster=="Ex_5",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Ex_5")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Ex_5_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Ex_5_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_5_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_5_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()


#### Ex_6

markers_ex1 <- f.markers[f.markers$cluster=="Ex_6",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Ex_6")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Ex_6_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Ex_6_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_6_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_6_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()


#### Ex_8

markers_ex1 <- f.markers[f.markers$cluster=="Ex_8",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Ex_8")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Ex_8_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Ex_8_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_8_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Ex_8_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()



#### Inh_5

markers_ex1 <- f.markers[f.markers$cluster=="Inh_5",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Inh_5")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Inh_5_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Inh_5_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_5_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_5_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()




#### Inh_1

markers_ex1 <- f.markers[f.markers$cluster=="Inh_1",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Inh_1")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Inh_1_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Inh_1_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_1_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_1_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()




#### Inh_2

markers_ex1 <- f.markers[f.markers$cluster=="Inh_2",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Inh_2")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Inh_2_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Inh_2_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_2_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_2_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()


#### Inh_3

markers_ex1 <- f.markers[f.markers$cluster=="Inh_3",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Inh_3")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Inh_3_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Inh_3_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_3_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_3_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

#### Inh_4

markers_ex1 <- f.markers[f.markers$cluster=="Inh_4",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Inh_4")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Inh_4_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Inh_4_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_4_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_4_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()



#### Inh_6

markers_ex1 <- f.markers[f.markers$cluster=="Inh_6",]
markers_ex1 <- data.frame(Gene.symbol = markers_ex1$gene,
                          logFC = markers_ex1$avg_log2FC,
                          adj.P.Val = markers_ex1$p_val_adj)

dir_path <- paste0("results_pathfinder", "Inh_6")
output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = "KEGG")

jpeg("images/Inh_6_enrichment_chart.jpeg", units="in", width=25, height=17, res=300)
enrichment_chart(output_) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()

jpeg("images/Inh_6_term_gene_heatmap.jpeg", units="in", width=25, height=17, res=300)
term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_6_term_gene_graph.jpeg", units="in", width=25, height=17, res=300)
term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
  theme(text = element_text(size = 20))
dev.off()

jpeg("images/Inh_6_UpSet_plot.jpeg", units="in", width=25, height=17, res=300)
UpSet_plot(output_, num_terms = 10, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
dev.off()





results_paths <- sapply(c("Ex_1", "Ex_2", "Ex_5", "Ex_6"), get_pathfinder)

enrichment_chart(results_paths[[1]])
term_gene_heatmap(output_rorb, use_description = TRUE, num_terms = 20) +
  theme(text = element_text(size = 20)) + NoLegend()

term_gene_graph(output_rorb, use_description = TRUE, num_terms = 10)

UpSet_plot(output_rorb, num_terms = 15, use_description = TRUE,
           method = "barplot") +
  theme(text = element_text(size = 20)) + NoLegend()
