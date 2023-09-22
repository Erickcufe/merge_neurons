make_DEG_and_PathFindR <- function(so, cell_type, braak, directory = "SFG_DEG", braak_1= 0){

  suppressMessages(library(Seurat))

  new_labels <- paste0(so$braak, "_", Idents(so))
  Idents(so) <- new_labels

  ident.2 <- paste0(braak_1,"_",cell_type)
  ident.1 <- paste0(braak,"_",cell_type)

  f.markers <- FindMarkers(so,
                          ident.1 = ident.1,
                          ident.2 = ident.2,
                          min.cells.group = 3,
                          min.cells.feature = 3,
                          min.pct = 0.1,
                          logfc.threshold = 0.25,
                          only.pos = FALSE)

  cell_type <- paste0(braak, "_", 0, "_", cell_type,"_", directory)

  f.markers$gene <- rownames(f.markers)

  path <- paste0(directory, "/", cell_type, ".csv")

  readr::write_csv(f.markers, file = path)


    paths <- c("KEGG", "Reactome", "GO-All")
    suppressMessages(library(pathfindR))
    suppressMessages(library(ggplot2))

    output_1 <- list()

    for(i in paths){
      # cell_type <- "RORB"
      markers_ex1 <- data.frame(Gene.symbol = f.markers$gene,
                                logFC = f.markers$avg_log2FC,
                                adj.P.Val = f.markers$p_val_adj)

      dir_path <- paste0("results_pathfinder","_", cell_type, "_",i)
      output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = i)

      saveRDS(output_, paste0(cell_type, "_", i, ".rds"))

      name_enrichment_chart <- paste0("images/", cell_type,"_",i, "_","enrichment_chart.jpeg")
      jpeg(name_enrichment_chart, units="in", width=25, height=17, res=300)
      print(enrichment_chart(output_, top_terms = 15) +
              theme(text = element_text(size = 25),
                    axis.text.x = element_text(size = 25),
                    axis.text.y = element_text(size = 25)) +
              scale_color_gradient(low = "#f7fcb9", high = "#31a354")
            )
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
}

so_sfg <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")

# Braak 2 vs 0
purrr::map(c("Ex_4", "Ex_5",
             "RORB", "Pv", "Sst", "Vip", "Non-Vip"), purrr::safely(.f = make_DEG_and_PathFindR), braak = 2,
           directory = "SFG_DEG", so = so_sfg,
           .progress = TRUE)

# Braak 6 vs 0
purrr::map(c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
             "RORB", "Pv", "Sst", "Vip", "Non-Vip"), purrr::safely(.f = make_DEG_and_PathFindR), braak = 6,
           braak_1= 0,
           directory = "SFG_DEG", so = so_sfg,
           .progress = TRUE)

# Braak 6 vs 2

purrr::map(c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
             "RORB", "Pv", "Sst", "Vip", "Non-Vip"), purrr::safely(.f = make_DEG_and_PathFindR), braak = 6,
           braak_1= 2,
           directory = "SFG_DEG", so = so_sfg,
           .progress = TRUE)

so_ec <- readRDS("EC_neurons_annoted_from_SFG.rds")

# Braak 2 vs 0
purrr::map(c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
             "RORB", "Pv", "Sst", "Vip", "Non-Vip"), purrr::safely(.f = make_DEG_and_PathFindR), braak = 2,
           directory = "EC_DEG", so = so_ec,
           .progress = TRUE)

# Braak 6 vs 0
purrr::map(c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
             "RORB", "Pv", "Sst", "Vip", "Non-Vip"), purrr::safely(.f = make_DEG_and_PathFindR), braak = 6,
           directory = "EC_DEG", so = so_ec,
           .progress = TRUE)
